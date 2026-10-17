"""Angular momentum deficit diagnostics for GENGA-style orbital tables."""

import pandas as pd # type: ignore
import numpy as np


class PlanetarySystem:
    """Read orbital rows; inclinations are degrees and units must match G.

    Row selection is relative to the original table. Per-call selections
    override constructor defaults rather than filtering an already selected
    table. ``planet_data`` exposes the constructor's selected rows.
    """

    def __init__(self, data_file, G=1, star_mass=1, skip_rows=0, limit_rows=None):
        self.G = G
        self.star_mass = star_mass
        self.skip_rows = skip_rows
        self.limit_rows = limit_rows
        self._raw_planet_data = pd.read_csv(
            data_file, sep=r'\s+', header=None,
            names=['t', 'i', 'a', 'e', 'inc', 'Omega', 'w', 'T', 'E', 'M', 'm', 'r']
        )
        self.planet_data = self._filter_rows(self._raw_planet_data, self.skip_rows, self.limit_rows)

    @staticmethod
    def _filter_rows(data, skip_rows=0, limit_rows=None):
        if skip_rows is None:
            skip_rows = 0
        if limit_rows is not None and limit_rows < 0:
            raise ValueError("limit_rows must be non-negative")

        filtered = data.iloc[skip_rows:].reset_index(drop=True) if skip_rows > 0 else data.copy()
        if limit_rows is not None:
            filtered = filtered.iloc[:limit_rows].reset_index(drop=True)
        return filtered

    def calculate_amd(self, time_point=None, normalized=False, skip_rows=None, limit_rows=None):
        """Return AMD for one time, defaulting to the first selected time.

        ``normalized=True`` divides AMD by the selected bodies' total mass;
        it does not normalize by their circular angular momentum.
        """
        if skip_rows is None:
            skip_rows = self.skip_rows
        if limit_rows is None:
            limit_rows = self.limit_rows

        data = self._filter_rows(self._raw_planet_data, skip_rows=skip_rows, limit_rows=limit_rows)

        if data.empty:
            return 0.0

        if time_point is None:
            time_point = data['t'].iloc[0]

        selected_data = data[data['t'] == time_point].copy()
        if selected_data.empty:
            return 0.0

        i_rad = np.radians(selected_data['inc'])

        h = np.sqrt(selected_data['a'] * (1 - selected_data['e']**2) * self.G * (self.star_mass + selected_data['m']))
        h_circular = np.sqrt(selected_data['a'] * self.G * (self.star_mass + selected_data['m']))

        amd = np.sum(selected_data['m'] * (h_circular - h * np.cos(i_rad)))

        total_mass = np.sum(selected_data['m'])
        if normalized:
            return amd / total_mass if total_mass != 0 else 0.0
        return amd


class AMDProcessor:
    def __init__(self, G=1, star_mass=1, skip_rows=0, limit_rows=None):
        self.G = G
        self.star_mass = star_mass
        self.skip_rows = skip_rows
        self.limit_rows = limit_rows

    def process_file(self, file_path, normalized=False, skip_rows=None, limit_rows=None):
        if skip_rows is None:
            skip_rows = self.skip_rows
        if limit_rows is None:
            limit_rows = self.limit_rows

        system = PlanetarySystem(
            file_path,
            G=self.G,
            star_mass=self.star_mass,
            skip_rows=skip_rows,
            limit_rows=limit_rows,
        )
        amd = system.calculate_amd(normalized=normalized)
        time_stamp = float(system.planet_data['t'].iloc[0]) if not system.planet_data.empty else 0.0
        return amd, time_stamp
