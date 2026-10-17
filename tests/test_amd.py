"""Analytic AMD values and selection semantics for the public diagnostics."""
import math

import pytest

from accrediff import AMDProcessor, PlanetarySystem


def write_orbits(path, bodies):
    # Multiple spaces and tabs also occur in whitespace-delimited tables.
    path.write_text('\n'.join(
        f'{time}  {index}\t{a} {e} {inc} 0 0 0 0 0 {mass} 0'
        for index, (time, a, e, inc, mass) in enumerate(bodies)
    ) + '\n')
    return path


@pytest.mark.parametrize('normalized', [False, True])
def test_constructor_row_selection_is_applied_once(tmp_path, normalized):
    path = write_orbits(tmp_path / 'orbits.txt', [
        (0, 5, 0.5, 10, 1), (0, 9, 0.5, 10, 1),
        (0, 1, 0, 0, 0.001), (0, 4, 0.6, 60, 0.002),
        (0, 16, 0.9, 80, 1),
    ])
    expected = 0.002 * math.sqrt(4 * 1.002) * (1 - 0.8 * 0.5)
    if normalized:
        expected /= 0.003
    system = PlanetarySystem(path, skip_rows=2, limit_rows=2)
    assert len(system.planet_data) == 2
    assert system.calculate_amd(normalized=normalized) == pytest.approx(expected)
    value, time = AMDProcessor(skip_rows=2, limit_rows=2).process_file(
        path, normalized=normalized,
    )
    assert value == pytest.approx(expected)
    assert time == 0.0


def test_per_call_selection_overrides_original_rows(tmp_path):
    path = write_orbits(tmp_path / 'orbits.txt', [
        (0, 1, 0, 0, 0.001), (0, 4, 0.6, 0, 0.002),
        (0, 9, 0, 0, 0.003),
    ])
    system = PlanetarySystem(path, skip_rows=2, limit_rows=1)
    assert system.calculate_amd() == 0.0
    expected = 0.002 * math.sqrt(4 * 1.002) * 0.2
    assert system.calculate_amd(skip_rows=1, limit_rows=1) == pytest.approx(expected)
    value, _ = AMDProcessor(skip_rows=2, limit_rows=1).process_file(
        path, skip_rows=1, limit_rows=1,
    )
    assert value == pytest.approx(expected)


def test_time_selection_and_processor_timestamp_agree(tmp_path):
    path = write_orbits(tmp_path / 'orbits.txt', [
        (1, 1, 0, 0, 0.001), (2, 4, 0.6, 0, 0.002),
        (3, 9, 0, 0, 0.003),
    ])
    expected = 0.002 * math.sqrt(4 * 1.002) * 0.2
    system = PlanetarySystem(path, skip_rows=1)
    assert system.calculate_amd() == pytest.approx(expected)
    assert system.calculate_amd(time_point=1) == 0.0
    assert system.calculate_amd(time_point=3) == 0.0
    value, time = AMDProcessor(skip_rows=1).process_file(path)
    assert value == pytest.approx(expected)
    assert time == 2.0


@pytest.mark.parametrize('selection', [{'skip_rows': 5}, {'limit_rows': 0}])
def test_empty_selection_returns_zero(tmp_path, selection):
    path = write_orbits(tmp_path / 'orbits.txt', [(1, 1, 0.6, 0, 0.001)])
    assert PlanetarySystem(path, **selection).calculate_amd() == 0.0
    assert AMDProcessor(**selection).process_file(path) == (0.0, 0.0)


def test_negative_row_limit_is_rejected(tmp_path):
    path = write_orbits(tmp_path / 'orbits.txt', [(1, 1, 0, 0, 0.001)])
    with pytest.raises(ValueError, match='limit_rows'):
        PlanetarySystem(path, limit_rows=-1)
