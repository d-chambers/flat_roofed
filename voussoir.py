"""
Voussoir beam analytical solution.

See Diederichs & Kaiser (1999). 
"""
import random
from dataclasses import dataclass
from functools import cache

import pandas as pd
import numpy as np

import local

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
        The external pressure applied down on the beam. See eq 9-12 for applying
        different types of pressure geometries. A list can be provided for
        multiple pressures. Support pressure should be negative (because it
        reduces the weight of the beam).
    """

    span: float = 15.0  # (m)
    thickness: float = 1.  # (m)
    stiffness: float = 1.0e10  # (Pa)
    ucs: float = 30_000_000  # (Pa)
    density: float = 3000  # (kg/m**3)
    joint_stiffness: float = 5.0e9  # (Pa/m)
    joint_spacing: float = 1  # (m)
    joint_friction_angle: float = 41  # (degrees)
    rockmass_stiffness: float | None = None
    pressure: float = 0

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

    def solve(self):
        """
        Return a dataframe of:

        [iteration_count, Np, buckling_limit, midspan_displacement, Fm, span]

        After iteratively solving for beam parameters.
        """
        # calc overall stiffness of beam if not provided
        stiffness = self.rockmass_stiffness or self.get_rockmass_stiffness()
        # get effective weight
        weight = self._gravity * self.density + np.sum(self.pressure) / self.thickness
        N_min, N_max = 1, 0
        N = np.arange(0.01, 1.01, 0.01).T
        Fmp = 1e100
        Fm, mm, Np, disp = 0, 0, 0, 0
        for count, i in enumerate(N):
            Zo = self.thickness * (1 - (2 / 3) * i)
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
                    Fm = weight * self.span * self.span / (4 * i * Z)
                    Fav = Fm * ((2 / 3) + i) / 3
                    del_l_prev = del_l
                    del_l = (Fav / stiffness) * L
                    mm += 1
                else:
                    del_l = del_l_prev
                    flag = 1

            if flag == 0 and k > 1:
                N_min = min([i, N_min])
                N_max = max([i, N_max])

            if Fmp > Fm and k > 1 and flag == 0:
                Fmp = Fm
                Np = i
                Zp = Z
                Zop = Zo
                disp = Zop - Zp

        friction_coef = np.tan(np.deg2rad(self.joint_friction_angle))
        buckling_limit = np.min([100 * (1 - (N_max - N_min)), 100])
        if buckling_limit == 200:
            breakpoint()
        results = dict(
            iteration_count=(count + 1) / 100,
            N=Np,
            buckling_limit=buckling_limit,
            buckling_fs=self._maximum_buckling_limit / buckling_limit,
            crushing_fs=self.ucs / Fmp,
            sliding_fs=(Fmp * Np) / (weight * self.span) * friction_coef,
            midspan_displacement=-disp,
            Fm=Fmp,
            span=self.span,
        )
        results.update(self.input_values)
        return pd.Series(results)

    @property
    def input_values(self) -> dict:
        """Return a dict of input values"""
        params = [x for x in list(self.__annotations__) if not x.startswith('_')]
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
                ucs=30_000_000,  # (Pa)
                density=3000,  # (kg/m**3)
                joint_stiffness=5.0e9,  # (Pa)
                joint_spacing=1,  # (m)
            )
            solution = vb.solve()[list(cls._test_columns)]
            out.append(solution)
        return pd.DataFrame(out)


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
        if not (0.6 < self.brittleness_factor < 1.25):
            msg = 'brittleness_factor must be between 0.6 and 1.25'
            raise ValueError(msg)

    @cache
    def get_modified_beam_stiffness(self):
        """Return the modified beam stiffness (eq 7)"""
        num_layers = self.thickness / self.horizontal_joint_spacing
        beam_stiffness = self.get_rockmass_stiffness()
        modified_beam_stiffness = 1.1 / (num_layers ** 2) * beam_stiffness
        return modified_beam_stiffness

    @cache
    def get_modified_thickness(self, mod_stiff=None):
        """Return the modified thickness (EQ 9)"""
        mod_stiff = mod_stiff or self.get_modified_beam_stiffness()
        out = (0.19 * 1 / np.sqrt(mod_stiff) + 1.05) * self.thickness
        return out

    def solve(self):
        """Solve getting new factors of safety."""
        # Solve with modified stiffness to get displacement and buckling.
        vals = dict(self.input_values)
        vals['rockmass_stiffness'] = self.get_modified_beam_stiffness()
        results_1 = DiederichBeam(**vals).solve()
        # now solve with modified thickness to get max stress and crushing limit
        vals = dict(self.input_values)
        vals['thickness'] = self.get_rockmass_stiffness()
        results_2 = DiederichBeam(**vals)
        breakpoint()


def test_with_matlab_results():
    """Ensure matlab results are consistent."""
    results = DiederichBeam()._get_test_output()
    # test against matlab code output
    expected = pd.read_csv("test/matlab_output.csv", names=results.columns)
    assert np.allclose(results.values, expected.values)


def sample_run(cls, run_count=1_000, **kwargs):
    """
    Solve iteratively by randomly sampling values provided in kwargs.
    """
    out = []
    # create an array of possible values for each param.
    variables = {i: np.array(v) for i, v in kwargs.items()}
    for _ in range(run_count):
        # get one value from each parameter, append solution
        inputs = {i: random.sample(v, k=1)[0] for i, v in variables.items()}
        out.append(cls(**inputs).solve())
    return pd.DataFrame(out)


if __name__ == "__main__":
    test_with_matlab_results()
