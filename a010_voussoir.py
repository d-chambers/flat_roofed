"""
Voussoir beam analytical solution.

See Diederichs & Kaiser (1999). 
"""
from dataclasses import dataclass

import pandas as pd
import numpy as np


@dataclass
class VoussoirBeam:
    """A class for analyzing Voussoir beams."""

    span: float = 15.0  # (m)
    thickness: float = 1.  # (m)
    stiffness: float = 1.0e10  # (Pa)
    ucs: float = 30_000_000  # (Pa)
    density: float = 3000  # (kg/m**3)
    joint_stiffness: float = 5.0e9  # (Pa)
    joint_spacing: float = 1  # (m)
    internal_friction = 41  # (degrees)
    overburden = 0

    _gravity: float = 9.81  # (m/s**2)
    _maximum_buckling_limit = 35  # as per detrichs
    _test_columns = (
        "iteration_count",
        "N",
        "buckling_limit",
        "midspan_displacement",
        "Fm",
        "span",
    )

    def solve(self):
        """
        Return a dataframe of:

        [iteration_count, Np, buckling_limit, midspan_displacement, Fm, span]

        After iteratively solving for beam parameters.
        """
        # calc overall stiffness of beam
        stiffness = 1. / (1. / self.stiffness + 1. / (self.joint_stiffness * self.joint_spacing))
        weight = self._gravity * self.density
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

        friction_coef = np.tan(np.deg2rad(self.internal_friction))
        buckling_limit = 100 * (1 - (N_max - N_min))
        results = dict(
            iteration_count=(count + 1) / 100,
            N=Np,
            buckling_limit=buckling_limit,
            buckling_fs=self._maximum_buckling_limit / buckling_limit,
            crushing_fs=self.ucs / Fmp,
            sliding_fs=(Fmp * Np)/(weight * self.span) * friction_coef,
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
            vb = VoussoirBeam(
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




def test_with_matlab_results():
    """Ensure matlab results are consistent."""
    results = VoussoirBeam()._get_test_output()
    # test against matlab code output
    expected = pd.read_csv("test/matlab_output.csv", names=results.columns)
    assert np.allclose(results.values, expected.values)


if __name__ == "__main__":
    test_with_matlab_results()

    vb = VoussoirBeam(
        thickness=1,
        span=15,
        ucs=30e6 / 2,
        stiffness=12e12,
        joint_stiffness=5.0e9,
        joint_spacing=1,
    )
    breakpoint()
    out = vb.solve()
    # out = vb.solve()
    # pass
