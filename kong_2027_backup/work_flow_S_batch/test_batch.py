"""Contracts for visible parameters, path expansion, and validated outputs."""
import importlib.util
import json
import ast
from pathlib import Path

import pandas as pd
import pytest
import numpy as np

runner = Path(__file__).with_name('batch_core.py')
spec = importlib.util.spec_from_file_location('sulfur_batch_runner', runner)
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)


def example_plan():
    notebook = json.loads(Path(__file__).with_name('run_batch.ipynb').read_text())
    namespace = {'batch': batch}
    exec(''.join(notebook['cells'][3]['source']), namespace)
    return namespace['CONFIG']


@pytest.fixture(autouse=True)
def isolated_nbody_directories(tmp_path, monkeypatch):
    """Contract tests need directory identities, not private N-body datasets."""
    load_plan = example_plan
    data_root = tmp_path / 'nbody_inputs'
    for model_id in load_plan()['models']:
        (data_root / model_id).mkdir(parents=True)

    def isolated_plan():
        config = load_plan()
        config['data_root'] = str(data_root)
        return config

    monkeypatch.setitem(globals(), 'example_plan', isolated_plan)


def test_backup_paths_and_stage_workspace_overrides(tmp_path, monkeypatch):
    backup = runner.resolve().parent.parent
    assert batch.NOTEBOOK_DIR == backup / 'work_flow_S'
    assert batch.validate_plan(example_plan())['output_root'].is_relative_to(backup)
    for cwd in (batch.PROJECT_ROOT, backup / 'work_flow_S', runner.parent.resolve()):
        monkeypatch.chdir(cwd)
        for filename in batch.NOTEBOOK_DIR.glob('*.ipynb'):
            cells = json.loads(filename.read_text())['cells']
            source = ''.join(cells[2]['source'])
            # Execute just path resolution, before imports or data access.
            bootstrap = source[source.index('working_dir ='):source.index('sys.path.insert')]
            for override in (None, tmp_path / filename.stem):
                env = {'Path': Path, 'PROJECT_ROOT': None, 'WORKFLOW_DIR': override}
                exec(bootstrap, env)
                assert env['project_root'] == batch.PROJECT_ROOT
                assert env['workflow_dir'] == (override or backup / 'work_flow_S')


def test_mass_based_endmember_handoff(tmp_path, monkeypatch):
    # Run actual A table-building cells and B mixing cell on one synthetic
    # particle, using real S8 equilibrium and temporary outputs only.
    monkeypatch.chdir(batch.PROJECT_ROOT)
    a = json.loads((batch.NOTEBOOK_DIR / batch.STAGE_BY_NAME['A_01'][0]).read_text())
    env = {'display': lambda *args: None}
    exec(''.join(a['cells'][1]['source']), env)
    env['WORKFLOW_DIR'] = tmp_path
    for i in (2, 4, 6, 8, 12, 13):
        exec(''.join(a['cells'][i]['source']), env)
    batch.validate_output('A_01', tmp_path / 'Tables')
    choices = ['EF', 'EC', 'OC', 'CI']
    weights = np.array([0.0, 0.5, 0.5, 0.0])
    env.update(choices=choices, METEORITE_FILE=None,
               METEORITE_INVENTORY_FILE=None, METEORITE_COVERAGE_FILE=None,
               df_gas_all=pd.DataFrame([dict(i=1, P=0.1, T=2000.0,
                   **dict(zip([f'{c}_ratio' for c in choices], weights)))]))
    b = json.loads((batch.NOTEBOOK_DIR / batch.STAGE_BY_NAME['B_01'][0]).read_text())
    exec(''.join(b['cells'][19]['source']), env)
    row = env['df_gas_all'].iloc[0]
    inventory = env['df_inventory']
    expected_s = inventory.loc['S', choices].to_numpy() @ weights
    assert row['S'] == pytest.approx(expected_s)
    assert row['Al2O3'] == pytest.approx(inventory.loc['Al2O3', choices].to_numpy() @ weights)
    assert row['S_sil'] + row['S_met'] == pytest.approx(expected_s)
    assert row['tracked_mass_fraction'] + row['untracked_mass_fraction'] == pytest.approx(1.0)
    # Prevent accidental reintroduction of unit-mole normalization.
    wrong_s = (inventory.div(inventory.sum(axis=0), axis=1)
               .loc['S', choices].to_numpy() @ weights)
    assert not np.isclose(row['S'], wrong_s)

    # Exercise D/E mass bookkeeping without rerunning an N-body history or
    # altering its chemical-solver strategy. Identity chemistry isolates the
    # partition/merge bookkeeping for all three impact types.
    d = json.loads((batch.NOTEBOOK_DIR / batch.STAGE_BY_NAME['D_01'][0]).read_text())
    denv = dict(env)
    exec(''.join(d['cells'][1]['source']), denv)
    exec(''.join(d['cells'][6]['source']), denv)
    definitions = [node for node in ast.parse(''.join(d['cells'][8]['source'])).body
                   if isinstance(node, ast.FunctionDef)]
    exec(compile(ast.Module(body=definitions, type_ignores=[]), 'D helpers', 'exec'), denv)
    denv['solve_sulfur_equilibrium'] = lambda incoming, **kwargs: (dict(incoming), {})
    initial = pd.DataFrame([row.to_dict(), row.to_dict()], index=[1, 2])
    initial['m_e'] = [0.01, 0.02]
    e = json.loads((batch.NOTEBOOK_DIR / batch.STAGE_BY_NAME['E_01'][0]).read_text())
    for event in ('small', 'global', 'partial'):
        events = pd.DataFrame([dict(Time=1.0, target_id=1, impactor_id=2,
            m_target=0.01, m_impactor=0.02, Mass=0.03, events=event,
            k_mantle=0.5, P_equil=5.0, T_equil=2500.0)])
        history, _ = denv['process_planet_history'](events, initial, planet_key='probe')
        assert history['untracked_mass_fraction'].iloc[0] == pytest.approx(row['untracked_mass_fraction'])
        eenv = dict(env)
        exec(''.join(e['cells'][1]['source']), eenv)
        eenv['df_planets'] = history.assign(planet='probe')
        exec(''.join(e['cells'][16]['source']), eenv)
        metrics = eenv['metrics'].iloc[0]
        assert metrics['tracked_mass_fraction_reconstructed'] == pytest.approx(row['tracked_mass_fraction'])
        assert metrics['CMF_full_lower_bound'] <= metrics['CMF'] <= metrics['CMF_full_upper_bound']


def test_complete_parameter_contract_and_shared_values():
    config = batch.validate_plan(example_plan())
    assert len(config['execution_plan']) == 7
    for entry in config['execution_plan']:
        names = batch.parameter_names(entry['stage'])
        assert set(entry['parameters']) == names - batch.AUTO_KEYS
    changed = example_plan()
    changed['execution_plan'][4]['parameters']['GLOBAL_MASS_THRESHOLD'] = 0.2
    with pytest.raises(ValueError, match='differ'):
        batch.validate_plan(changed)
    changed = example_plan()
    changed['execution_plan'][5]['parameters']['TYPO'] = 3
    with pytest.raises(ValueError, match='missing/unknown'):
        batch.validate_plan(changed)


def test_parameterized_copy_shows_resolved_paths(tmp_path):
    source = batch.NOTEBOOK_DIR / batch.STAGE_BY_NAME['C_02'][0]
    params = batch.resolved_parameters('C_02', example_plan()['execution_plan'][4]['parameters'],
                                       tmp_path / 'workspace', 'model_01', tmp_path / 'model_01')
    target = tmp_path / 'C_02.ipynb'
    batch.parameterize(source, target, params)
    original = json.loads(source.read_text())
    copied = json.loads(target.read_text())
    first = ''.join(copied['cells'][1]['source'])
    assert "PARTIAL_PRESSURE_RATIO = 0.6" in first
    assert str(tmp_path / 'workspace' / 'Tables') in first
    assert 'parameters' in copied['cells'][1]['metadata']['tags']
    assert ''.join(original['cells'][1]['source']) != first


def test_missing_upstream_output_stops_before_notebook_execution(tmp_path):
    with pytest.raises(FileNotFoundError, match='Meteorites_composition.csv'):
        batch.run_stage('B_01', tmp_path, example_plan()['execution_plan'][1]['parameters'],
                        model_id='model_01', model_dir=tmp_path / 'model_01', timeout=1)
    assert not (tmp_path / 'executed_notebooks').exists()


def test_final_output_validation_checks_sulfur_closure_and_orbital_order(tmp_path):
    path = tmp_path / '07_bulk_comparison_S.csv'
    frame = pd.DataFrame({'planet': ['planets_1', 'planets_2'], 'a_final_au': [0.8, 1.5],
                          'Mass': [0.7, 0.2], 'Mg_number': [0.85, 0.8], 'CMF': [0.3, 0.2],
                          'S': [0.10, 0.20], 'S_sil': [0.03, 0.05], 'S_met': [0.07, 0.15],
                          'tracked_mass_fraction': [0.9, 0.8],
                          'untracked_mass_fraction': [0.1, 0.2],
                          'tracked_mass_fraction_reconstructed': [0.9, 0.8],
                          'CMF_full_lower_bound': [0.27, 0.16],
                          'CMF_full_upper_bound': [0.37, 0.36]})
    frame.to_csv(path, index=False)
    assert batch.validate_output('E_01', tmp_path) == {}
    frame.loc[1, 'a_final_au'] = 0.5
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match='orbital order'):
        batch.validate_output('E_01', tmp_path)


def test_model_commit_replaces_only_its_own_result_and_rebases_notebooks(tmp_path):
    attempt = tmp_path / '.staging' / 'attempt_1'
    shared = attempt / 'shared'
    staged = attempt / 'models' / 'model_01'
    final = tmp_path / 'models' / 'model_01'
    other = tmp_path / 'models' / 'model_02'
    (final / 'tables').mkdir(parents=True)
    (final / 'tables' / 'old.csv').write_text('old')
    other.mkdir(parents=True)
    (other / 'marker.txt').write_text('keep')
    stage_tables = staged / 'stage_workspaces' / 'A_01' / 'Tables'
    stage_tables.mkdir(parents=True)
    for source_name in batch.STAGE_BY_NAME['A_01'][1]:
        (stage_tables / source_name).write_text('new')
    (staged / 'tables').mkdir()
    for source_name in batch.STAGE_BY_NAME['A_01'][1]:
        (staged / 'tables' / source_name).write_text('new')
    stage = 'A_01'
    filename = batch.STAGE_BY_NAME[stage][0]
    params = batch.resolved_parameters(stage, example_plan()['execution_plan'][0]['parameters'],
                                       staged / 'stage_workspaces' / stage, '_shared', tmp_path)
    for folder in ('parameterized_notebooks', 'executed_notebooks'):
        batch.parameterize(batch.NOTEBOOK_DIR / filename, staged / folder / filename, params)
    batch.atomic_json(staged / 'manifest.json', {
        'stages': {stage: {'status': 'passed', 'log': str(staged / 'logs' / 'A_01.log')}},
        'history': [],
    })
    batch.commit_model(staged, final, shared)
    assert not (final / 'tables' / 'old.csv').exists()
    assert (final / 'tables' / 'Meteorites_composition.csv').read_text() == 'new'
    assert (other / 'marker.txt').read_text() == 'keep'
    saved = json.loads((final / 'executed_notebooks' / filename).read_text())
    first = ''.join(saved['cells'][1]['source'])
    assert str(final / 'stage_workspaces' / stage) in first
    assert str(staged) not in first
    assert saved['metadata']['batch_execution_paths']['model_workspace'] == str(staged)
    manifest = json.loads((final / 'manifest.json').read_text())
    assert 'log' not in manifest['stages'][stage]
    assert manifest['logs_retained_with_model'] is False


def test_interrupted_promotion_restores_old_model(tmp_path):
    attempt = tmp_path / '.staging' / 'interrupted'
    backup = attempt / 'backup' / 'model_01'
    backup.mkdir(parents=True)
    (backup / 'old.txt').write_text('prior success')
    staged = attempt / 'models' / 'model_01'
    staged.mkdir(parents=True)
    (staged / 'incomplete.txt').write_text('discard')
    batch.recover_staging(tmp_path)
    assert (tmp_path / 'models' / 'model_01' / 'old.txt').read_text() == 'prior success'
    assert not attempt.exists()


def test_current_batch_summary_counts_impacts_without_old_run_id(tmp_path):
    model = tmp_path / 'models' / 'model_01'
    tables = model / 'tables'
    tables.mkdir(parents=True)
    with pd.ExcelWriter(tables / '05_impact_events_S.xlsx') as writer:
        pd.DataFrame({'events': ['small', 'global', 'small']}).to_excel(
            writer, sheet_name='planets_1', index=False)
    latest = tmp_path / 'latest_run'
    latest.mkdir()
    outcome = {'attempt_id': 'attempt_1', 'model_id': 'model_01',
               'scenario_id': 'baseline', 'run_dir': str(model),
               'last_stage': 'C_02', 'status': 'passed',
               'failed_stage': '', 'error': '', 'old_result_retained': False}
    batch.collect_summary(latest, [outcome], 'attempt_1')
    counts = pd.read_csv(latest / 'impact_counts.csv')
    assert counts.loc[0, 'attempt_id'] == 'attempt_1'
    assert counts.loc[0, 'small'] == 2
    assert counts.loc[0, 'global'] == 1
