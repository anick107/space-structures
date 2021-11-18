"""Stres spectrum generator

It is a simple and user-friendly alternative to the Excel or GUI tools
for the stress spectra generation. Single stress component (S0) spectra
are considered.

Compatible with NASGRO and NASFORM.
"""
from typing import List, Dict
import collections
import math
from scipy import interpolate
import numpy as np

Event = collections.namedtuple("Event", ["neq", "s_min", "s_max", "desc"])


class StressSpectrum:
    """Stress spectrum generator

    Effectively a collector of stress spectrum events with the common
    interface for the safety factors and output.
    """

    def __init__(
        self,
        filename: str,
        prestress: float = 0,
        sf_stress: float = 1.0,
        sf_cycles: float = 4.0,
    ):
        """Inits the stress spectrum

        Keyword arguments:
        filename -- full path to the stress spectrum file (output)
        prestress -- prestress value
        sf_stress -- scaling factor for stress data
        sf_cycles -- scaling factor for number of cycles

        sf_stress is applied only to the alternating stress, and
        is not applied to the prestress.
        """
        self.filename = filename
        self.prestres = prestress
        self.sf_stress = sf_stress
        self.sf_cycles = sf_cycles
        self.events: List[Event] = []

    def append(self, event: Event):
        self.events.append(event)

    def stats(self):
        """Print statistics"""
        s_max = self.sf_stress * max([_.s_max for _ in self.events])
        s_max += self.prestres
        print(f"> Peak max stress: {s_max:4.2f} MPa")

        s_min = self.sf_stress * min([_.s_min for _ in self.events])
        s_min += self.prestres
        print(f"> Peak min stress: {s_min:4.2f} MPa")

        total_n = self.sf_cycles * sum([_.neq for _ in self.events])
        print(f"> Total number of cycles: {math.ceil(total_n)}")

        print("> The following safety factors are included:")
        print(f">> cycles={self.sf_cycles}")
        print(f">> stress={self.sf_stress}")

    def save(self, verbose: bool = True) -> None:
        """Save stress spectrum file in max-min-cycle format

        Keyword arguments:
        verbose -- print the events with the comments

        This format is compatible with NASFORM.
        """
        with open(self.filename, "w") as fh:
            fh.write("max min cycle format\n")
            if verbose:
                print("==============")
            for event in self.events:
                s_min = self.sf_stress * event.s_min + self.prestres
                s_max = self.sf_stress * event.s_max + self.prestres
                neq = math.ceil(self.sf_cycles * event.neq)
                s = f"{s_max:4.2f}  {s_min:4.2f}  {neq}"
                fh.write(s + "\n")

                if verbose:
                    print(s.ljust(26) + " | " + event.desc)


def RandomEvent(
    stress: float,
    fn: float,
    duration: float,
    n: float = 2.0,
    desc: str = "",
) -> Event:
    """Add random vibration event

    Keyword arguments:
    stress -- stress level [MPa]
    fn -- relevant eigen frequency [Hz]
    duration -- duration [seconds]
    n -- Paris' law constrant (n=2 is conservative)
    desc -- short description

    The calculation of the equivalent cycles is performed based on
    NASGRO appendix G.
    """

    table_n = {
        2.0: 0.222,
        3.0: 0.139,
        4.0: 0.099,
        5.0: 0.077,
        6.0: 0.066,
    }
    neq = fn * duration * linear_interpolation(table_n, n)
    neq = math.ceil(neq)
    return Event(neq=neq, s_min=-stress, s_max=stress, desc=desc)


def SineEvent(
    stress_1g: float,
    load: float,
    sweep_rate: float,
    desc: str = "",
) -> Event:
    """Get stress spectrum for sine event

    Keyword arguments:
    stress_1g  -- principle stress for 1g load [MPa]
    load       -- maximum acceleration [g]
    sweep_rate -- sweep rate [oct/min]
    desc       -- description

    Calculation is performed based on NASGO appendix H.

    @todo: one could expand sine sweep as a number of bins
    """
    # typical test is from 20 up to 100 Hz
    df = 100 - 20
    neq = 60 / (sweep_rate * math.log(2)) * df
    neq = math.ceil(neq)
    stress = load * stress_1g
    return Event(neq=neq, s_min=-stress, s_max=stress, desc=desc)


def ThermalEvent(
    stress_1K: float,
    T_max: float,
    T_min: float,
    neq: int = 1,
    T_ref: float = 20,
    desc: str = "",
):
    """Stress spectrum for thermal cycle

    Keyword arguments:
    stress_1K -- principle stress for 1K load [MPa]
    T_max -- max temperature within the cycle
    T_min -- min temperature within the cycle
    neq   -- number of cycles
    T_ref -- room temperature
    """
    s_min = stress_1K * (T_min - T_ref)
    s_max = stress_1K * (T_max - T_ref)
    return Event(neq=neq, s_min=s_min, s_max=s_max, desc=desc)


def linear_interpolation(table: Dict[float, float], argument: float):
    """Linear interpolation based on dict table

    Keyword arguments:
    table -- dictionary[x, y] for the interpolation
    argument -- argument for interpolation (single value)
    """
    args = sorted(table.keys())
    x = np.array(args)
    y = np.array([table[_] for _ in x])
    f = interpolate.interp1d(x, y)
    return f(argument)


if __name__ == "__main__":
    spectrum = StressSpectrum("example.spc", prestress=10.0)
    spectrum.append(
        RandomEvent(
            stress=300 * 1.15,
            fn=200,
            duration=120,
            desc="Random X",
        )
    )
    spectrum.append(
        ThermalEvent(
            stress_1K=10,
            T_max=50,
            T_min=0,
            T_ref=20,
            neq=365,
            desc="Thermoelastic",
        )
    )
    spectrum.stats()
    spectrum.save(verbose=True)
