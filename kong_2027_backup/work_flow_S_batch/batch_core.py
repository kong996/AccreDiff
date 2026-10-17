"""Execution engine used by run_batch.ipynb for the sulfur A–E workflow."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from contextlib import contextmanager
import ast
import copy
import pprint

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_DIR = Path(__file__).resolve().parent.parent / "work_flow_S"
STAGES = (
    ("A_01", "A_01_meteorite_bulk_composition_s.ipynb", ("Meteorites_composition.csv", "Meteorites_inventory_per_100g.csv", "Meteorites_mass_coverage.csv")),
    ("B_01", "B_01_embryo_bulk_composition_s.ipynb", ("03_embryo_bulk_composition.csv",)),
    ("B_02", "B_02_accretion_history_s.ipynb", ("04_accretion_history_S.xlsx", "04_planet_order_S.csv")),
    ("C_01", "C_01_small_impact_s.ipynb", ("05_small_global_impact_events_S.xlsx", "05_impact_parameters_S.json")),
    ("C_02", "C_02_giant_impact_s.ipynb", ("05_impact_events_S.xlsx",)),
    ("D_01", "D_01_differentiation_s.ipynb", ("06_differentiation_S.xlsx",)),
    ("E_01", "E_01_bulk_comparison_s.ipynb", ("07_bulk_comparison_S.csv",)),
)
STAGE_BY_NAME = {stage: (filename, outputs) for stage, filename, outputs in STAGES}
INPUTS = {
    'A_01': (),
    'B_01': ('Meteorites_composition.csv', 'Meteorites_inventory_per_100g.csv', 'Meteorites_mass_coverage.csv'),
    'B_02': ('03_embryo_bulk_composition.csv',),
    'C_01': ('04_accretion_history_S.xlsx', '04_planet_order_S.csv'),
    'C_02': ('05_small_global_impact_events_S.xlsx', '05_impact_parameters_S.json'),
    'D_01': ('Meteorites_composition.csv', '03_embryo_bulk_composition.csv', '05_impact_events_S.xlsx'),
    'E_01': ('04_planet_order_S.csv', '06_differentiation_S.xlsx'),
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    os.replace(temporary, path)


def resolve_path(value: str) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else PROJECT_ROOT / path).resolve()


def check_id(value: str, label: str) -> None:
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value) or value in {".", ".."}:
        raise ValueError(f"Invalid {label}: {value!r}")


def parameter_names(stage: str) -> set[str]:
    """The source notebook's first code cell is the single parameter contract."""
    notebook = json.loads((NOTEBOOK_DIR / STAGE_BY_NAME[stage][0]).read_text())
    cells = [cell for cell in notebook['cells'] if cell['cell_type'] == 'code']
    if not cells or 'parameters' not in cells[0].get('metadata', {}).get('tags', []):
        raise ValueError(f'{stage}: first code cell must have parameters tag')
    names = set()
    for node in ast.parse(''.join(cells[0]['source'])).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
            raise ValueError(f'{stage}: parameter cell must contain assignments only')
        names.add(node.targets[0].id)
    return names


PATH_KEYS = {'PROJECT_ROOT', 'WORKFLOW_DIR', 'MODEL_ID', 'MODEL_DIR', 'METEORITE_FILE',
             'METEORITE_INVENTORY_FILE', 'METEORITE_COVERAGE_FILE',
             'OUTPUT_FILE', 'OUTPUT_HISTORY_FILE', 'OUTPUT_ORDER_FILE', 'ACCRETION_FILE',
             'OUTPUT_INVENTORY_FILE', 'OUTPUT_COVERAGE_FILE',
             'OUTPUT_EVENTS_FILE', 'OUTPUT_PARAMETERS_FILE', 'INTERMEDIATE_FILE',
             'PARAMETER_FILE', 'COMPOSITION_FILE', 'IMPACT_FILE', 'DIFFERENTIATION_FILE',
             'ORDER_FILE'}
AUTO_KEYS = PATH_KEYS  # resolved and printed by this driver, never scientific defaults


def validate_plan(config: dict) -> dict:
    required = {'data_root', 'output_root', 'models', 'scenario_id', 'execution_plan',
                'timeouts'}
    if set(config) != required:
        raise ValueError(f'Configuration keys differ: {sorted(set(config) ^ required)}')
    result = copy.deepcopy(config)
    result['data_root'] = resolve_path(str(config['data_root']))
    result['output_root'] = resolve_path(str(config['output_root']))
    if not result['data_root'].is_dir():
        raise NotADirectoryError(result['data_root'])
    if result['output_root'] == result['data_root'] or result['output_root'].is_relative_to(result['data_root']):
        raise ValueError('Output root must be outside N-body data')
    check_id(result['scenario_id'], 'scenario ID')
    if not isinstance(result['models'], list) or not result['models'] or len(set(result['models'])) != len(result['models']):
        raise ValueError('models must be a nonempty unique list')
    for model in result['models']:
        check_id(model, 'model ID')
        if not (result['data_root'] / model).is_dir():
            raise NotADirectoryError(result['data_root'] / model)
    if set(result['timeouts']) - {'default', *STAGE_BY_NAME} or not all(isinstance(v, int) and v > 0 for v in result['timeouts'].values()) or 'default' not in result['timeouts']:
        raise ValueError('Invalid stage timeout')
    plan = result['execution_plan']
    if not isinstance(plan, list) or [item.get('stage') for item in plan] != [x[0] for x in STAGES]:
        raise ValueError('Execution plan must list A_01 through E_01 in order')
    for item in plan:
        stage = item['stage']
        if set(item) != {'stage', 'notebook', 'enabled', 'parameters'} or item['notebook'] != STAGE_BY_NAME[stage][0] or not isinstance(item['enabled'], bool):
            raise ValueError(f'{stage}: malformed plan entry')
        keys = set(item['parameters'])
        expected = parameter_names(stage) - AUTO_KEYS
        if keys != expected:
            raise ValueError(f'{stage}: missing/unknown parameters {sorted(keys ^ expected)}')
    active = [item['stage'] for item in plan if item['enabled']]
    if not active or active != [x[0] for x in STAGES[:len(active)]]:
        raise ValueError('Enabled stages must be a contiguous prefix from A_01')
    values = {item['stage']: item['parameters'] for item in plan}
    if values['B_01']['GAS_TIME'] != values['B_02']['GAS_TIME']:
        raise ValueError('B_01/B_02 GAS_TIME differ')
    for key in ('GLOBAL_MASS_THRESHOLD', 'IMPACT_RATIO_RANGE'):
        if values['C_01'][key] != values['C_02'][key]:
            raise ValueError(f'C_01/C_02 {key} differ')
    sulfur = values['A_01']['SULFUR_WT']
    if set(sulfur) != {'EF', 'EC', 'OC', 'CI'} or any(not isinstance(v,(int,float)) or not np.isfinite(v) or v < 0 for v in sulfur.values()):
        raise ValueError('Invalid A_01 SULFUR_WT')
    if any(not isinstance(values['A_01'][key],(int,float)) or not np.isfinite(values['A_01'][key])
           or values['A_01'][key] <= 0 for key in ('EF_ENRICHMENT_FACTOR','EF_DEPLETION_FACTOR')):
        raise ValueError('A_01 EF enrichment/depletion factors must be positive')
    if not re.fullmatch(r'\d+(?:\.\d+)? Myr',values['B_01']['GAS_TIME']):
        raise ValueError('GAS_TIME must be an available snapshot label')
    if not re.fullmatch(r'\d+(?:\.\d+)? Myr',values['B_01']['INITIAL_TIME']):
        raise ValueError('B_01 INITIAL_TIME must be a snapshot label')
    b=values['B_01']
    if not 0 < b['EF_MAX_AU'] < b['EC_MAX_AU'] < b['OC_MAX_AU'] or not b['EC_MAX_AU'] < b['CI_MIN_AU'] <= b['OC_MAX_AU']:
        raise ValueError('B_01 initial radial composition boundaries are inconsistent')
    if values['B_02']['PLANET_MASS_THRESHOLD'] <= 0 or values['B_01']['EMBRYO_DISPLAY_THRESHOLD'] < 0:
        raise ValueError('Invalid mass threshold')
    for stage,key in [('C_01','P_RATIO'),('C_02','PARTIAL_PRESSURE_RATIO')]:
        if not 0 < values[stage][key] <= 1:
            raise ValueError(f'{stage} {key} outside (0,1]')
    if values['C_01']['GLOBAL_MASS_THRESHOLD'] <= 0 or values['C_01']['COOLING_DURATION_MYR'] <= 0:
        raise ValueError('Invalid C_01 mass/cooling parameters')
    ratios = values['C_01']['IMPACT_RATIO_RANGE']
    if len(ratios) != 2 or not 0 <= ratios[0] < ratios[1] <= 1:
        raise ValueError('Invalid impact ratio interval')
    if values['D_01']['ROOT_SEARCH'] not in {'fast','thorough'}:
        raise ValueError('D_01 ROOT_SEARCH must be fast or thorough')
    for key,value in values['D_01'].items():
        if key == 'ROOT_SEARCH':
            continue
        if not isinstance(value,(int,float)) or not np.isfinite(value):
            raise ValueError(f'D_01 {key} must be finite')
        if key.endswith(('FRACTION','FACTOR')) and not 0 <= value <= 1:
            raise ValueError(f'D_01 {key} outside [0,1]')
        if not key.endswith(('FRACTION','FACTOR')) and value <= 0:
            raise ValueError(f'D_01 {key} must be positive')
    if not isinstance(values['E_01']['SAVE_FIGURES'],bool):
        raise ValueError('E_01 SAVE_FIGURES must be Boolean')
    figure_size=values['E_01']['FIGURE_SIZE_MM']
    if len(figure_size)!=2 or any(not isinstance(x,(int,float)) or x<=0 for x in figure_size):
        raise ValueError('E_01 FIGURE_SIZE_MM must contain positive dimensions')
    limits=values['E_01']['PLOT_DEFAULT_YLIMS']
    if set(limits)!={'Mass','Mg#','CMF'} or any(len(x)!=2 or x[0]>=x[1] for x in limits.values()):
        raise ValueError('E_01 plot limits must be increasing pairs for Mass, Mg# and CMF')
    return result


@contextmanager
def exclusive_output_lock(output_root: Path):
    output_root.mkdir(parents=True, exist_ok=True)
    with (output_root / '.batch.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another batch run is writing to {output_root}') from exc
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns - set(frame.columns)
    if frame.empty or missing:
        raise ValueError(f"{label}: empty table or missing columns {sorted(missing)}")


def require_workbook(path: Path, columns: set[str]) -> dict[str, pd.DataFrame]:
    sheets = pd.read_excel(path, sheet_name=None)
    if not sheets:
        raise ValueError(f"No planet sheets: {path}")
    for name, frame in sheets.items():
        require_columns(frame, columns, f"{path.name}/{name}")
        if "Time" in frame and not frame["Time"].is_monotonic_increasing:
            raise ValueError(f"{path.name}/{name}: events out of time order")
    return sheets


def validate_output(stage: str, tables: Path) -> dict:
    outputs = STAGE_BY_NAME[stage][1]
    for filename in outputs:
        path = tables / filename
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty {stage} output: {path}")
    info = {}
    if stage == "A_01":
        frame = pd.read_csv(tables / outputs[0], index_col=0)
        if not {'EF', 'EC', 'OC', 'CI'} <= set(frame.columns) or 'S' not in frame.index:
            raise ValueError("A_01 meteorite table lacks source or sulfur entries")
        inventory = pd.read_csv(tables / outputs[1], index_col=0)
        coverage = pd.read_csv(tables / outputs[2], index_col=0)
        if not frame.index.equals(inventory.index) or not frame.columns.equals(inventory.columns):
            raise ValueError("A_01 display and inventory tables have different components")
        if not {'tracked_mass_fraction', 'untracked_mass_fraction'} <= set(coverage.columns):
            raise ValueError("A_01 lacks mass coverage columns")
        if set(coverage.index) != set(inventory.columns):
            raise ValueError("A_01 coverage and inventory source labels differ")
        if not np.isfinite(inventory.to_numpy(dtype=float)).all() or (inventory < 0).any().any():
            raise ValueError("A_01 inventory must be finite and nonnegative")
        if not np.allclose(coverage['tracked_mass_fraction'] + coverage['untracked_mass_fraction'], 1, atol=1e-10):
            raise ValueError("A_01 source mass does not close")
    elif stage == "B_01":
        frame = pd.read_csv(tables / outputs[0])
        require_columns(frame, {'i', 'm_e', 'S', 'S_sil', 'S_met',
            'tracked_mass_fraction', 'untracked_mass_fraction', 'source_mass_relative_error'}, stage)
        if frame['i'].duplicated().any() or not np.allclose(frame['S'], frame['S_sil'] + frame['S_met'], atol=1e-10):
            raise ValueError("B_01 particle IDs or sulfur inventory are inconsistent")
        if not np.allclose(frame['tracked_mass_fraction'] + frame['untracked_mass_fraction'], 1, atol=1e-7):
            raise ValueError("B_01 modeled and unmodeled mass do not close")
    elif stage == "B_02":
        sheets = require_workbook(tables / outputs[0], {'Time', 'indexi', 'indexj', 'Mass'})
        order = pd.read_csv(tables / outputs[1])
        require_columns(order, {'planet', 'final_particle_id', 'a_final_au', 'm_e_final'}, stage)
        if order['planet'].tolist() != list(sheets) or not order['a_final_au'].is_monotonic_increasing:
            raise ValueError("B_02 planet sheets and orbital ordering differ")
    elif stage in {'C_01', 'C_02'}:
        sheets = require_workbook(tables / outputs[0], {'Time', 'target_id', 'impactor_id', 'events', 'P_equil', 'T_equil', 'k_mantle'})
        allowed = {'small', 'global'} if stage == 'C_01' else {'small', 'global', 'partial'}
        for name, frame in sheets.items():
            if frame['events'].isna().any() or not set(frame['events']) <= allowed:
                raise ValueError(f"{stage}/{name}: unexpected event type")
            pt = frame[['P_equil', 'T_equil']].to_numpy(dtype=float)
            if not np.isfinite(pt).all() or (pt <= 0).any():
                raise ValueError(f"{stage}/{name}: invalid pressure or temperature")
        if stage == 'C_01':
            params = json.loads((tables / outputs[1]).read_text())
            if not {'P_RATIO', 'GLOBAL_MASS_THRESHOLD', 'IMPACT_RATIO_RANGE', 'COOLING_DURATION_MYR'} <= set(params):
                raise ValueError("C_01 parameter record is incomplete")
    elif stage == "D_01":
        sheets = require_workbook(tables / outputs[0], {'Time', 'm_e', 'accepted', 'accepted_approximate',
            'oxygen_converged', 'solution_strategy', 'O_KD_log_residual', 'S', 'S_sil', 'S_met', 'max_conservation_error',
            'tracked_mass_fraction', 'untracked_mass_fraction'})
        event_rows = []
        for planet, frame in sheets.items():
            if frame[['accepted', 'accepted_approximate']].isna().any().any() or not frame['accepted'].astype(bool).all():
                raise ValueError(f"D_01/{planet}: unaccepted differentiation event")
            if not np.allclose(frame['S'], frame['S_sil'] + frame['S_met'], atol=1e-10):
                raise ValueError(f"D_01/{planet}: sulfur phase closure failed")
            if not np.allclose(
                frame['tracked_mass_fraction'] + frame['untracked_mass_fraction'],
                1, rtol=0, atol=1e-7,
            ):
                raise ValueError(f"D_01/{planet}: complete-body mass does not close")
            if not np.isfinite(frame[['O_KD_log_residual', 'max_conservation_error']].to_numpy(dtype=float)).all():
                raise ValueError(f"D_01/{planet}: nonfinite residual")
            if (frame['m_e'] <= 0).any() or frame['solution_strategy'].isna().any():
                raise ValueError(f"D_01/{planet}: invalid mass or missing solution strategy")
            event_rows.append(frame.assign(planet=planet))
        events = pd.concat(event_rows, ignore_index=True)
        info = {"events": len(events), "accepted_approximate": int(events['accepted_approximate'].astype(bool).sum()),
                "max_abs_O_KD_log_residual": float(events['O_KD_log_residual'].abs().max()),
                "max_conservation_error": float(events['max_conservation_error'].abs().max())}
    elif stage == "E_01":
        frame = pd.read_csv(tables / outputs[0])
        require_columns(frame, {'planet', 'a_final_au', 'Mass', 'Mg_number', 'CMF', 'S', 'S_sil', 'S_met',
            'tracked_mass_fraction', 'untracked_mass_fraction', 'tracked_mass_fraction_reconstructed',
            'CMF_full_lower_bound', 'CMF_full_upper_bound'}, stage)
        values = frame[['a_final_au', 'Mass', 'Mg_number', 'CMF', 'S', 'S_sil', 'S_met']].to_numpy(dtype=float)
        if not np.isfinite(values).all() or (frame[['a_final_au', 'Mass']] <= 0).any().any():
            raise ValueError('E_01 contains nonfinite or nonpositive orbit/mass values')
        if frame['planet'].duplicated().any() or not frame['a_final_au'].is_monotonic_increasing:
            raise ValueError('E_01 planet identity or orbital order is invalid')
        if not np.allclose(frame['S'], frame['S_sil'] + frame['S_met'], atol=1e-10):
            raise ValueError("E_01 sulfur phase closure failed")
        if not np.allclose(
            frame['tracked_mass_fraction'] + frame['untracked_mass_fraction'],
            1, rtol=0, atol=1e-7,
        ):
            raise ValueError("E_01 complete-body mass does not close")
        if not frame['Mg_number'].between(0, 1).all() or not frame['CMF'].between(0, 1).all():
            raise ValueError("E_01 Mg# or CMF outside physical range")
        if not np.allclose(frame['tracked_mass_fraction_reconstructed'],
                           frame['tracked_mass_fraction'], rtol=0, atol=1e-7):
            raise ValueError("E_01 reconstructed mass coverage differs from differentiation history")
        if not ((frame['CMF_full_lower_bound'] <= frame['CMF'])
                & (frame['CMF'] <= frame['CMF_full_upper_bound'])).all():
            raise ValueError("E_01 full-body CMF bounds do not contain modeled CMF")
    return info


def notebook_text_output(path: Path) -> str:
    if not path.is_file():
        return ''
    notebook = json.loads(path.read_text())
    chunks = []
    for index, cell in enumerate(notebook.get('cells', [])):
        if cell.get('cell_type') != 'code':
            continue
        for output in cell.get('outputs', []):
            payload = output.get('text') or output.get('data', {}).get('text/plain')
            if payload:
                chunks.append(f'\n--- cell {index} ---\n' + ''.join(payload))
            if output.get('traceback'):
                chunks.append('\n'.join(output['traceback']))
    return ''.join(chunks)


OUTPUT_PARAMETER_KEYS = {
    'A_01': {'OUTPUT_FILE':'Meteorites_composition.csv',
             'OUTPUT_INVENTORY_FILE':'Meteorites_inventory_per_100g.csv',
             'OUTPUT_COVERAGE_FILE':'Meteorites_mass_coverage.csv'},
    'B_01': {'OUTPUT_FILE':'03_embryo_bulk_composition.csv'},
    'B_02': {'OUTPUT_HISTORY_FILE':'04_accretion_history_S.xlsx','OUTPUT_ORDER_FILE':'04_planet_order_S.csv'},
    'C_01': {'OUTPUT_EVENTS_FILE':'05_small_global_impact_events_S.xlsx','OUTPUT_PARAMETERS_FILE':'05_impact_parameters_S.json'},
    'C_02': {'OUTPUT_FILE':'05_impact_events_S.xlsx'},
    'D_01': {'OUTPUT_FILE':'06_differentiation_S.xlsx'},
    'E_01': {'OUTPUT_FILE':'07_bulk_comparison_S.csv'},
}
INPUT_PARAMETER_KEYS = {
    'B_01': {'METEORITE_FILE':'Meteorites_composition.csv',
             'METEORITE_INVENTORY_FILE':'Meteorites_inventory_per_100g.csv',
             'METEORITE_COVERAGE_FILE':'Meteorites_mass_coverage.csv'},
    'C_01': {'ACCRETION_FILE':'04_accretion_history_S.xlsx'},
    'C_02': {'INTERMEDIATE_FILE':'05_small_global_impact_events_S.xlsx','PARAMETER_FILE':'05_impact_parameters_S.json'},
    'D_01': {'COMPOSITION_FILE':'03_embryo_bulk_composition.csv','IMPACT_FILE':'05_impact_events_S.xlsx','METEORITE_FILE':'Meteorites_composition.csv'},
    'E_01': {'DIFFERENTIATION_FILE':'06_differentiation_S.xlsx','ORDER_FILE':'04_planet_order_S.csv'},
}


def resolved_parameters(stage: str, science: dict, workspace: Path, model_id: str, model_dir: Path) -> dict:
    tables = workspace / 'Tables'
    params = {'PROJECT_ROOT':str(PROJECT_ROOT),'WORKFLOW_DIR':str(workspace)}
    if stage in {'B_01','B_02'}:
        params.update(MODEL_ID=model_id,MODEL_DIR=str(model_dir))
    params.update(science)
    params.update({key:str(tables/name) for key,name in INPUT_PARAMETER_KEYS.get(stage,{}).items()})
    params.update({key:str(tables/name) for key,name in OUTPUT_PARAMETER_KEYS[stage].items()})
    if set(params) != parameter_names(stage):
        raise ValueError(f'{stage}: resolved parameter mismatch {set(params)^parameter_names(stage)}')
    return params


def parameterize(source: Path, destination: Path, params: dict) -> None:
    nb = json.loads(source.read_text())
    code = next(c for c in nb['cells'] if c['cell_type']=='code')
    if 'parameters' not in code.get('metadata',{}).get('tags',[]):
        raise ValueError('First code cell not tagged parameters')
    lines = ['# Actual parameters used for this saved run.\n']
    for key,value in params.items():
        lines.append(f'{key} = {pprint.pformat(value, width=90, sort_dicts=False)}\n')
    code['source']=lines
    for c in nb['cells']:
        if c['cell_type']=='code':
            c['execution_count']=None;c['outputs']=[]
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def preview_plan(config: dict) -> None:
    config = validate_plan(config)
    print(f"Scenario {config['scenario_id']}; models {config['models']}")
    print(f"Data: {config['data_root']}\nOutput: {config['output_root']}")
    print('每次从 A_01 重新计算；同模型成功后覆盖旧结果，失败时保留旧结果。')
    shared_workspace = (config['output_root'] / '.staging' / '<attempt_id>' /
                        'shared' / 'stage_workspaces' / 'A_01')
    shared_params = resolved_parameters('A_01', config['execution_plan'][0]['parameters'],
                                        shared_workspace, '_shared', config['data_root'])
    print(f'\nA_01 {config["execution_plan"][0]["notebook"]} (本批次共用)')
    for key, value in shared_params.items():
        print(f'  {key} = {value!r}')
    for model_id in config['models']:
        model_dir = config['data_root'] / model_id
        final_dir = config['output_root'] / 'models' / model_id
        print(f'\nMODEL {model_id}: input={model_dir}; result={final_dir}')
        for item in config['execution_plan'][1:]:
            if not item['enabled']:
                break
            stage = item['stage']
            workspace = (config['output_root'] / '.staging' / '<attempt_id>' /
                         'models' / model_id / 'stage_workspaces' / stage)
            params = resolved_parameters(stage, item['parameters'], workspace, model_id, model_dir)
            print(f'\n{stage} {item["notebook"]}')
            for key, value in params.items():
                print(f'  {key} = {value!r}')


def record_stage(destination: Path, stage: str, result: dict) -> None:
    path=destination/'manifest.json'
    manifest=json.loads(path.read_text()) if path.exists() else {'stages':{},'history':[]}
    manifest['stages'][stage]=result
    manifest['history'].append({'stage':stage,**result})
    atomic_json(path,manifest)


def run_stage(stage: str, destination: Path, science: dict, *, model_id: str,
              model_dir: Path, timeout: int) -> dict:
    filename,outputs=STAGE_BY_NAME[stage]
    source=NOTEBOOK_DIR/filename
    tables=destination/'tables';tables.mkdir(parents=True,exist_ok=True)
    workspace=destination/'stage_workspaces'/stage
    stage_tables=workspace/'Tables';stage_tables.mkdir(parents=True,exist_ok=True)
    for name in INPUTS[stage]:
        src=tables/name
        if not src.exists():raise FileNotFoundError(f'{stage} missing upstream input: {src}')
        shutil.copy2(src,stage_tables/name)
    params=resolved_parameters(stage,science,workspace,model_id,model_dir)
    param_nb=destination/'parameterized_notebooks'/filename
    exec_nb=destination/'executed_notebooks'/filename
    parameterize(source,param_nb,params)
    started=datetime.now(timezone.utc).isoformat()
    log_path=destination/'logs'/f'{stage}.log';log_path.parent.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['MPLCONFIGDIR']=str(workspace/'matplotlib')
    exec_nb.parent.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,'-m','jupyter','nbconvert','--to','notebook','--execute',
         f'--ExecutePreprocessor.timeout={timeout}','--output',str(exec_nb),str(param_nb)]
    start=time.monotonic()
    try:
        process=subprocess.run(cmd,cwd=PROJECT_ROOT,env=env,text=True,stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,timeout=timeout+90)
        log_path.write_text(process.stdout+notebook_text_output(exec_nb))
        if process.returncode:raise RuntimeError(f'{stage} execution failed; see {log_path}')
        info=validate_output(stage,stage_tables)
        for name in outputs:shutil.copy2(stage_tables/name,tables/name)
        checksums={name:sha256_file(tables/name) for name in outputs}
        result={'status':'passed','started_at':started,
                'finished_at':datetime.now(timezone.utc).isoformat(),'seconds':round(time.monotonic()-start,2),
                'log':str(log_path),'parameterized_notebook':str(param_nb),'executed_notebook':str(exec_nb),
                'outputs':checksums,**info}
    except Exception as exc:
        if isinstance(exc,subprocess.TimeoutExpired):log_path.write_text(str(exc))
        result={'status':'failed','started_at':started,
                'finished_at':datetime.now(timezone.utc).isoformat(),'error':str(exc),'log':str(log_path),
                'parameterized_notebook':str(param_nb),'executed_notebook':str(exec_nb) if exec_nb.exists() else None}
    record_stage(destination,stage,result)
    return result


def collect_summary(latest_dir: Path, outcomes: list[dict], attempt_id: str) -> None:
    status=pd.DataFrame(outcomes)
    status.to_csv(latest_dir/'batch_status.csv',index=False)
    planets=[];events=[];impacts=[]
    for row in outcomes:
        if row['status']!='passed':continue
        tables=Path(row['run_dir'])/'tables'
        if row['last_stage']=='E_01':
            frame=pd.read_csv(tables/'07_bulk_comparison_S.csv')
            frame.insert(0,'attempt_id',attempt_id);frame.insert(1,'scenario_id',row['scenario_id']);frame.insert(2,'model_id',row['model_id'])
            planets.append(frame)
        if (tables/'06_differentiation_S.xlsx').exists():
            for planet,frame in pd.read_excel(tables/'06_differentiation_S.xlsx',sheet_name=None).items():
                frame.insert(0,'attempt_id',attempt_id);frame.insert(1,'scenario_id',row['scenario_id']);frame.insert(2,'model_id',row['model_id']);frame.insert(3,'planet',planet)
                frame.insert(4,'event_number',range(1,len(frame)+1));events.append(frame)
        if (tables/'05_impact_events_S.xlsx').exists():
            for planet,frame in pd.read_excel(tables/'05_impact_events_S.xlsx',sheet_name=None).items():
                counts=frame['events'].value_counts().to_dict()
                impacts.append({'attempt_id':attempt_id,'scenario_id':row['scenario_id'],'model_id':row['model_id'],'planet':planet,**counts})
    if planets:pd.concat(planets,ignore_index=True).to_csv(latest_dir/'planets.csv',index=False)
    if events:
        df=pd.concat(events,ignore_index=True);df.to_csv(latest_dir/'solver_events.csv',index=False)
        df.groupby(['attempt_id','scenario_id','model_id','planet']).agg(
            total_events=('accepted','size'),accepted_approximate=('accepted_approximate','sum'),
            max_abs_oxygen_kd_residual=('O_KD_log_residual',lambda x:x.abs().max())
        ).reset_index().to_csv(latest_dir/'solver_summary.csv',index=False)
    if impacts:pd.DataFrame(impacts).fillna(0).to_csv(latest_dir/'impact_counts.csv',index=False)


def recover_staging(output_root: Path) -> None:
    """Restore a previous successful model after an interrupted promotion."""
    staging_root = output_root / '.staging'
    if not staging_root.exists():
        return
    models_root = output_root / 'models'
    models_root.mkdir(parents=True, exist_ok=True)
    for attempt in staging_root.iterdir():
        if not attempt.is_dir():
            raise RuntimeError(f'Unexpected staging entry: {attempt}')
        backups = attempt / 'backup'
        if backups.exists():
            for backup in backups.iterdir():
                live = models_root / backup.name
                if live.exists():
                    shutil.rmtree(backup)
                else:
                    os.replace(backup, live)
        shutil.rmtree(attempt)


def rebase_notebooks(staged_model: Path, final_model: Path, shared: Path) -> None:
    """Keep saved notebook paths valid after a successful directory promotion."""
    old_prefixes = (str(staged_model), str(shared))
    for folder in ('parameterized_notebooks', 'executed_notebooks'):
        for path in (staged_model / folder).glob('*.ipynb'):
            notebook = json.loads(path.read_text())
            cell = next(c for c in notebook['cells'] if c['cell_type'] == 'code')
            original_source = ''.join(cell['source'])
            rebased_source = original_source.replace(old_prefixes[0], str(final_model))
            rebased_source = rebased_source.replace(old_prefixes[1], str(final_model))
            if original_source == rebased_source:
                raise ValueError(f'Saved notebook has no staging path to rebase: {path}')
            cell['source'] = rebased_source.splitlines(keepends=True)
            notebook.setdefault('metadata', {})['batch_execution_paths'] = {
                'model_workspace': old_prefixes[0], 'shared_workspace': old_prefixes[1],
                'saved_model_directory': str(final_model),
            }
            # Check every final stage input/output path against its staged counterpart.
            for node in ast.parse(rebased_source).body:
                if not isinstance(node, ast.Assign):
                    continue
                name = node.targets[0].id
                if name not in PATH_KEYS or name in {'PROJECT_ROOT', 'MODEL_ID', 'MODEL_DIR'}:
                    continue
                value = ast.literal_eval(node.value)
                if isinstance(value, str) and value.startswith(str(final_model)):
                    staged = value.replace(str(final_model), old_prefixes[0], 1)
                    if not Path(staged).exists():
                        staged = value.replace(str(final_model), old_prefixes[1], 1)
                    if not Path(staged).exists():
                        raise FileNotFoundError(f'{path.name}: missing retained input/output {value}')
            path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + '\n')
    manifest_path = staged_model / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    for entry in list(manifest.get('stages', {}).values()) + manifest.get('history', []):
        for key in ('parameterized_notebook', 'executed_notebook'):
            if isinstance(entry.get(key), str):
                entry[key] = entry[key].replace(old_prefixes[0], str(final_model)).replace(
                    old_prefixes[1], str(final_model))
        if 'shared_source' in entry:
            entry['shared_source'] = 'computed once for this batch; retained under stage_workspaces/A_01'
        # Latest-attempt logs are overwritten by future batches and cannot be
        # cited as permanent records of this committed model result.
        entry.pop('log', None)
    manifest['saved_model_directory'] = str(final_model)
    manifest['logs_retained_with_model'] = False
    atomic_json(manifest_path, manifest)


def commit_model(staged_model: Path, final_model: Path, shared: Path) -> None:
    """Replace one model only after its complete staged result is validated."""
    rebase_notebooks(staged_model, final_model, shared)
    backup = staged_model.parent.parent / 'backup' / final_model.name
    backup.parent.mkdir(parents=True, exist_ok=True)
    final_model.parent.mkdir(parents=True, exist_ok=True)
    had_old = final_model.exists()
    if had_old:
        os.replace(final_model, backup)
    try:
        os.replace(staged_model, final_model)
        for folder in ('parameterized_notebooks', 'executed_notebooks'):
            for path in (final_model / folder).glob('*.ipynb'):
                notebook = json.loads(path.read_text())
                source = ''.join(next(c for c in notebook['cells'] if c['cell_type'] == 'code')['source'])
                if str(staged_model) in source or str(shared) in source:
                    raise ValueError(f'Unrebased staging path in {path}')
                for node in ast.parse(source).body:
                    if not isinstance(node, ast.Assign):
                        continue
                    key = node.targets[0].id
                    if key not in PATH_KEYS or key in {'PROJECT_ROOT', 'MODEL_ID', 'MODEL_DIR'}:
                        continue
                    value = ast.literal_eval(node.value)
                    if isinstance(value, str) and value.startswith(str(final_model)):
                        if not Path(value).exists():
                            raise FileNotFoundError(f'Committed notebook path is missing: {value}')
    except Exception:
        if final_model.exists():
            shutil.rmtree(final_model)
        if had_old and backup.exists():
            os.replace(backup, final_model)
        raise
    if had_old:
        shutil.rmtree(backup)


def copy_attempt_logs(source: Path, latest_dir: Path, model_id: str) -> None:
    source_logs = source / 'logs'
    if not source_logs.is_dir():
        return
    target = latest_dir / 'logs' / model_id
    target.mkdir(parents=True, exist_ok=True)
    for path in source_logs.glob('*.log'):
        shutil.copy2(path, target / path.name)


def execute(raw_config: dict) -> Path:
    config = validate_plan(raw_config)
    output_root = config['output_root']
    with exclusive_output_lock(output_root):
        recover_staging(output_root)
        attempt_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        attempt_root = output_root / '.staging' / attempt_id
        attempt_root.mkdir(parents=True)
        latest_dir = output_root / 'latest_run'
        if latest_dir.exists():
            shutil.rmtree(latest_dir)
        latest_dir.mkdir(parents=True)
        active = [item for item in config['execution_plan'] if item['enabled']]
        science = {item['stage']: item['parameters'] for item in active}
        attempt_record = {
            'attempt_id': attempt_id, 'started_at': datetime.now(timezone.utc).isoformat(),
            'status': 'running', 'models': config['models'], 'scenario_id': config['scenario_id'],
            'data_root': str(config['data_root']), 'output_root': str(output_root),
            'last_stage': active[-1]['stage'], 'execution_plan': config['execution_plan'],
        }
        atomic_json(latest_dir / 'attempt.json', attempt_record)
        outcomes = []
        shared = attempt_root / 'shared'
        try:
            try:
                a = run_stage('A_01', shared, science['A_01'], model_id='_shared',
                              model_dir=config['data_root'],
                              timeout=config['timeouts'].get('A_01', config['timeouts']['default']))
            except Exception as exc:
                a = {'status': 'failed', 'error': str(exc)}
                shared.mkdir(parents=True, exist_ok=True)
                record_stage(shared, 'A_01', a)
                (shared / 'logs').mkdir(exist_ok=True)
                (shared / 'logs' / 'A_01_driver.log').write_text(repr(exc) + '\n')
            copy_attempt_logs(shared, latest_dir, '_shared')
            print(f'A_01: {a["status"]}', flush=True)
            for model_id in config['models']:
                staged = attempt_root / 'models' / model_id
                staged.mkdir(parents=True)
                final = output_root / 'models' / model_id
                model_dir = config['data_root'] / model_id
                atomic_json(staged / 'parameters.json', {
                    'attempt_id': attempt_id, 'model_id': model_id,
                    'scenario_id': config['scenario_id'], 'data_root': str(config['data_root']),
                    'model_dir': str(model_dir), 'last_stage': active[-1]['stage'],
                    'composition_basis': 'mol_per_100g_whole_body',
                    'execution_plan': config['execution_plan'],
                })
                outcome = {'attempt_id': attempt_id, 'model_id': model_id,
                           'scenario_id': config['scenario_id'], 'run_dir': str(final),
                           'last_stage': active[-1]['stage'], 'status': 'failed',
                           'failed_stage': '', 'error': '', 'old_result_retained': False}
                if a['status'] == 'failed':
                    outcome.update(failed_stage='A_01', error=a['error'],
                                   old_result_retained=final.exists())
                    outcomes.append(outcome)
                    continue
                tables = staged / 'tables'
                tables.mkdir(exist_ok=True)
                for source_name in STAGE_BY_NAME['A_01'][1]:
                    shutil.copy2(shared / 'tables' / source_name, tables / source_name)
                shutil.copytree(shared / 'stage_workspaces' / 'A_01',
                                staged / 'stage_workspaces' / 'A_01')
                for folder in ('parameterized_notebooks', 'executed_notebooks'):
                    source = shared / folder / STAGE_BY_NAME['A_01'][0]
                    target = staged / folder / source.name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
                record_stage(staged, 'A_01', {'status': a['status'],
                             'shared_source': str(shared), 'outputs': a.get('outputs', {})})
                failed = None
                for item in active[1:]:
                    stage = item['stage']
                    try:
                        result = run_stage(stage, staged, science[stage], model_id=model_id,
                                           model_dir=model_dir,
                                           timeout=config['timeouts'].get(stage, config['timeouts']['default']))
                    except Exception as exc:
                        result = {'status': 'failed', 'error': str(exc)}
                        record_stage(staged, stage, result)
                        (staged / 'logs').mkdir(exist_ok=True)
                        (staged / 'logs' / f'{stage}_driver.log').write_text(repr(exc) + '\n')
                    print(f'{model_id} / {stage}: {result["status"]}', flush=True)
                    if result['status'] == 'failed':
                        failed = (stage, result['error'])
                        break
                copy_attempt_logs(staged, latest_dir, model_id)
                if (staged / 'logs').exists():
                    shutil.rmtree(staged / 'logs')
                if failed:
                    outcome.update(failed_stage=failed[0], error=failed[1],
                                   old_result_retained=final.exists())
                else:
                    manifest_path = staged / 'manifest.json'
                    manifest = json.loads(manifest_path.read_text())
                    manifest.update(attempt_id=attempt_id, scenario_id=config['scenario_id'],
                                    last_stage=active[-1]['stage'],
                                    completed_at=datetime.now(timezone.utc).isoformat())
                    atomic_json(manifest_path, manifest)
                    try:
                        commit_model(staged, final, shared)
                        outcome.update(status='passed', old_result_retained=False)
                    except Exception as exc:
                        outcome.update(failed_stage='commit', error=str(exc),
                                       old_result_retained=final.exists())
                outcomes.append(outcome)
            collect_summary(latest_dir, outcomes, attempt_id)
            attempt_record.update(status='complete',
                                  finished_at=datetime.now(timezone.utc).isoformat(),
                                  passed_models=sum(row['status']=='passed' for row in outcomes),
                                  failed_models=sum(row['status']=='failed' for row in outcomes))
            atomic_json(latest_dir / 'attempt.json', attempt_record)
            print(f'Batch status: {latest_dir / "batch_status.csv"}')
            return latest_dir
        except Exception as exc:
            attempt_record.update(status='failed', error=str(exc),
                                  finished_at=datetime.now(timezone.utc).isoformat())
            atomic_json(latest_dir / 'attempt.json', attempt_record)
            raise
        finally:
            backups = attempt_root / 'backup'
            if attempt_root.exists() and not (backups.exists() and any(backups.iterdir())):
                shutil.rmtree(attempt_root)
