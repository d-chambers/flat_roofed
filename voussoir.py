"""
Voussoir beam analytical solution.

See Diederichs & Kaiser (1999). 
"""
from dataclasses import dataclass
from functools import cache

import numpy as np
import pandas as pd

import local


# minimum properties (chapter 2.9 of his thesis)
ABOUSELEIMAN_MIN = dict(
    span=12.5,  # (m)
    thickness=3,  # (m)
    stiffness=13.8e9,  # (Pa)
    rockmass_ucs=30.7e6 * 0.5,  # (Pa)
    density=24000/9.81,  # (kg/m**3)
    joint_stiffness=1500e9,  # (Pa/m)
    joint_spacing=0.8,  # (m)
    joint_friction_angle=28,  # (degrees)
    bolt_length=3.9,
)

ABOUSELEIMAN_MAX = dict(
    span=12.5,  # (m)
    thickness=4.5,  # (m)
    stiffness=13.8e9,  # (Pa)
    rockmass_ucs=30.7e6,  # (Pa)
    density=24000/9.81,  # (kg/m**3)
    joint_stiffness=4000e9,  # (Pa/m)
    joint_spacing=0.8,  # (m)
    joint_friction_angle=45,  # (degrees)
    bolt_length=3.9,
)


@dataclass
class DiederichBeam:
    """
    A class for analyzing Voussoir beams using Diederichs and Kaiser (1999).

    Parameters
    ----------
    span
        The excavation span in m.
    thickness
        The thickness of the beam in m.
    rockmass_ucs
        The unconfined compressive strength of the beam material (in Pa).
        Note: The authors recommend multiplying laboratory UCS by a number
        between 0.5 and 0.3 to account for rockmass strength reduction.
    density
        The density of the beam material (in kg/m**3).
    joint_stiffness
        The horizontal stiffness of the joints in Pa/m.
    joint_spacing
        The spacing of the vertical joints.
    joint_friction_angle
        Internal friction angle for beam material (degrees)
    rockmass_stiffness
        The effective stiffness of the rockmass. If None, this is calculated
        according to stiffness and joint stiffness.
    pressure
        The total pressure applied to the beam and overburden (per 1 m slice).
        If the density is None, the beam weight is taken out. If not, this
        is the weight above the beam.
    """
    span: float = 15.0  # (m)
    thickness: float = 1.  # (m)
    stiffness: float = 1.0e10  # (Pa)
    rockmass_ucs: float = 30_000_000  # (Pa)
    density: float | None = 3000  # (kg/m**3)
    joint_stiffness: float = 5.0e9  # (Pa/m)
    joint_spacing: float = 1  # (m)
    joint_friction_angle: float = 41  # (degrees)
    rockmass_stiffness: float | None = None
    pressure: float | None = None

    _gravity = local.GRAVITATIONAL_CONSTANT
    _maximum_buckling_limit = 35  # this is when snap-through is assumed to occur
    _test_columns = (
        "iteration_count",
        "N",
        "buckling_limit",
        "midspan_displacement",
        "Fm",
        "span",
    )

    def get_rockmass_stiffness(self):
        """Get the rockmass stiffness from material and joint stiffness"""
        inv_stiff = 1. / self.stiffness
        inv_joint = 1. / (self.joint_stiffness * self.joint_spacing)
        return 1. / (inv_stiff + inv_joint)

    def get_effective_weight(self):
        """Get the effective weight."""
        pressure = self.pressure if self.pressure is not None else 0
        return self._gravity * self.density + np.sum(pressure) / self.thickness

    def solve(self):
        """
        Return a dataframe of:

        [iteration_count, Np, buckling_limit, midspan_displacement, Fm, span]

        After iteratively solving for beam parameters.
        """
        stiffness = self.rockmass_stiffness or self.get_rockmass_stiffness()
        weight = self.get_effective_weight()
        N_min, N_max = 1, 0
        N = np.arange(0.01, 1.01, 0.01).T
        Fmp = 1e100
        Fm, mm, Np, disp = 0, 0, 0, 0
        for count, n in enumerate(N):
            Zo = self.thickness * (1 - (2 / 3) * n)
            L = self.span + (8 * Zo * Zo / (3 * self.span))
            del_l = 0
            del_l_prev = 100
            k = 0
            flag = 0
            while abs(del_l - del_l_prev) > 0.000001:
                zchk = ((8 * Zo * Zo / (3 * self.span)) - del_l)
                if zchk >= 0:
                    k += 1
                    Z = np.sqrt(3 * self.span * zchk / 8)
                    Fm = weight * (self.span ** 2) / (4 * n * Z)
                    Fav = Fm * ((2 / 3) + n) / 3
                    del_l_prev = del_l
                    del_l = (Fav / stiffness) * L
                    mm += 1
                else:
                    del_l = del_l_prev
                    flag = 1

            if flag == 0 and k > 1:
                N_min = min([n, N_min])
                N_max = max([n, N_max])

            if Fmp > Fm and k > 1 and flag == 0:
                Fmp = Fm
                Np = n
                Zp = Z
                Zop = Zo
                disp = Zop - Zp

        friction_coef = np.tan(np.deg2rad(self.joint_friction_angle))
        buckling_limit = np.min([100 * (1 - (N_max - N_min)), 100])
        results = dict(
            iteration_count=(count + 1) / 100,
            N=Np,
            buckling_limit=buckling_limit,
            buckling_fs=self._maximum_buckling_limit / buckling_limit,
            crushing_fs=self.rockmass_ucs / Fmp,
            sliding_fs=(Fmp * Np) / (weight * self.span) * friction_coef,
            midspan_displacement=-disp,
            Fm=Fmp,
            span=self.span,
        )
        results.update(self.diederich_inputs)
        return pd.Series(results)

    @property
    def diederich_inputs(self) -> dict:
        """Return a dict of input values"""
        params = [x for x in list(DiederichBeam.__annotations__) if not x.startswith('_')]
        out = {x: getattr(self, x) for x in params}
        return out

    @classmethod
    def _get_test_output(cls):
        """Get the output for testing against matlab defaults."""
        out = []
        for span in np.arange(1, 21):
            vb = cls(
                span=span,  # (m)
                thickness=1.,  # (m)
                stiffness=1.0e10,  # (Pa)
                rockmass_ucs=30_000_000,  # (Pa)
                density=3000,  # (kg/m**3)
                joint_stiffness=5.0e9,  # (Pa)
                joint_spacing=1,  # (m)
            )
            solution = vb.solve()[list(cls._test_columns)]
            out.append(solution)
        return pd.DataFrame(out)

    def __hash__(self):
        return hash(id(self))


@dataclass
class AbousleimanBeam(DiederichBeam):
    """
    A class for analyzing Voussoir beams using Abousleiman (2021).

    This allows the analytical solution to account for bolting and horizontal
    bedding. Inputs are the same as DiederichBeam with the following additional
    parameters:

    Parameters
    ---------
    bolt_length
        The bolt length, in m. If None, assume bolt_length is the same as
        thickness.
    horizontal_joint_spacing
        The spacing of bedding planes, in m. If None, assume same as
        thickness (no bedding)
    brittleness_factor
        A factor describing the post-peak behavior of the material. Values range
        from 0.6 (very brittle) to 1.25 (perfectly plastic). This is B in eq 10.
    """
    bolt_length: float | None = None
    horizontal_joint_spacing: float | None = None
    brittleness_factor: float = 1.0

    def __post_init__(self):
        """Set default values if left None, check brittleness_factor"""
        self.bolt_length = self.bolt_length or self.thickness
        self.horizontal_joint_spacing = self.horizontal_joint_spacing or self.thickness
        if not (0.6 <= self.brittleness_factor <= 1.25):
            msg = 'brittleness_factor must be between 0.6 and 1.25'
            raise ValueError(msg)

    @cache
    def get_modified_beam_stiffness(self):
        """Return the modified beam stiffness (eq 7)"""
        bolt_joint_ratio = self.bolt_length / self.horizontal_joint_spacing
        num_layers = np.ceil(bolt_joint_ratio)
        beam_stiffness = self.get_rockmass_stiffness()
        modified_beam_stiffness = (1.1 / (num_layers ** 2)) * beam_stiffness
        return modified_beam_stiffness

    @cache
    def get_modified_thickness(self):
        """Return the modified thickness (EQ 9)"""
        mod_stiff = self.get_modified_beam_stiffness()
        # need to convert mod_stiff to GPa for eq to work.
        mod_stiff = mod_stiff / 1e9
        thick = self.horizontal_joint_spacing
        out = (0.19 * 1 / np.sqrt(mod_stiff) + 1.05) * thick
        # print(out, mod_stiff)
        return out

    def get_inputs_1(self):
        """Get the first inputs for solving for buckling limit/displacement."""
        vals = dict(self.diederich_inputs)
        vals['rockmass_stiffness'] = self.get_modified_beam_stiffness()
        vals['thickness'] = self.bolt_length
        return vals

    def get_inputs_2(self):
        """Get the second inputs for solving for max stress in beam."""
        vals = dict(self.diederich_inputs)
        vals['rockmass_stiffness'] = self.get_rockmass_stiffness()
        vals['thickness'] = self.get_modified_thickness()
        return vals

    def solve(self):
        """Solve getting new factors of safety."""
        # populate run inputs
        out = dict(self.diederich_inputs)
        out['bolt_length'] = self.bolt_length
        out['horizontal_joint_spacing'] = self.horizontal_joint_spacing
        out['brittleness_factor'] = self.brittleness_factor
        # next get buckling limit and displacement
        valid_outputs1 = (
            'buckling_limit', 'buckling_fs', 'midspan_displacement',
        )
        results_1 = DiederichBeam(**self.get_inputs_1()).solve()
        out.update({x: results_1[x] for x in valid_outputs1})
        # now solve with modified thickness to get max stress and crushing limit
        results_2 = DiederichBeam(**self.get_inputs_2()).solve()
        out['Fm'] = results_2['Fm']
        # print(out['Fm'])
        out['sliding_fs'] = results_2['sliding_fs']
        out['crushing_fs'] = self.get_crushing_fs(out)
        return pd.Series(out)

    def __hash__(self):
        return hash(id(self))

    def get_crushing_fs(self, out):
        """Get the crushing factor of safety."""
        return (self.rockmass_ucs / out['Fm']) * self.brittleness_factor


def test_with_matlab_results():
    """Ensure matlab results are consistent."""
    results = DiederichBeam()._get_test_output()
    # test against matlab code output
    expected = pd.read_csv("test/matlab_output.csv", names=results.columns)
    assert np.allclose(results.values, expected.values)


def test_abouslieman_case_1():
    """tests for Abousleiman cases from chapter 2.9 of his thesis."""
    # no load case
    params = dict(ABOUSELEIMAN_MIN)
    params['bolt_length'] = 3.0
    params['brittleness_factor'] = 0.6
    params['horizontal_joint_spacing'] = params['thickness'] / 3
    params['pressure'] = 0
    solution_1 = AbousleimanBeam(**params).solve()

    # loaded case
    params = dict(params)
    params['pressure'] = 120_000 * (2/3)
    breakpoint()
    solution_2 = AbousleimanBeam(**params).solve()
    breakpoint()



    breakpoint()
    # solution = beam1.solve()


if __name__ == "__main__":
    test_with_matlab_results()
    test_abouslieman_case_1()
