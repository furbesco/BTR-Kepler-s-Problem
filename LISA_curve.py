"""
LISA sky/polarization-averaged sensitivity curve.

Implements the analytic fit from
    Robson, Cornish & Liu, Class. Quantum Grav. 36, 105011 (2019)
    doi:10.1088/1361-6382/ab1101  (arXiv:1803.01944)
specifically Eq. (13) for the instrument noise combined with Eq. (14)
and Table 1 for the galactic confusion foreground.

Returns the strain power spectral density S_n(f) [1/Hz]; take sqrt for the
amplitude spectral density [1/sqrt(Hz)].
"""

import numpy as np

# 2018 LISA Phase-0 reference design (as adopted in RCL19)
_L = 2.5e9                  # arm length [m]
_C = 299_792_458.0          # speed of light [m/s]
_F_STAR = _C / (2 * np.pi * _L)   # transfer frequency ~ 19.09 mHz

# Galactic confusion fit, Table 1 of RCL19. Amplitude A is fixed at 9e-45.
_SC_A = 9.0e-45
_SC_PARAMS = {
    0.5: dict(alpha=0.133, beta=243.0,  kappa=482.0,  gamma=917.0,  fk=2.58e-3),
    1.0: dict(alpha=0.171, beta=292.0,  kappa=1020.0, gamma=1680.0, fk=2.15e-3),
    2.0: dict(alpha=0.165, beta=299.0,  kappa=611.0,  gamma=1340.0, fk=1.73e-3),
    4.0: dict(alpha=0.138, beta=-221.0, kappa=521.0,  gamma=1680.0, fk=1.13e-3),
}


def P_OMS(f):
    """Single-link optical metrology noise PSD [m^2/Hz], RCL19 Eq. (10)."""
    return (1.5e-11)**2 * (1.0 + (2e-3 / f)**4)


def P_acc(f):
    """Single test-mass acceleration noise PSD [(m/s^2)^2/Hz], RCL19 Eq. (11)."""
    return (3e-15)**2 * (1.0 + (0.4e-3 / f)**2) * (1.0 + (f / 8e-3)**4)


def S_confusion(f, Tobs_yr=4.0):
    """Unresolved galactic-binary foreground, RCL19 Eq. (14) + Table 1."""
    if Tobs_yr not in _SC_PARAMS:
        raise ValueError(f"Tobs_yr must be one of {sorted(_SC_PARAMS)}")
    p = _SC_PARAMS[Tobs_yr]
    return (_SC_A * f**(-7.0/3.0)
            * np.exp(-f**p["alpha"] + p["beta"] * f * np.sin(p["kappa"] * f))
            * (1.0 + np.tanh(p["gamma"] * (p["fk"] - f))))


def S_n(f, Tobs_yr=4.0, include_confusion=True):
    """
    Effective strain PSD S_n(f) [1/Hz]  (RCL19 Eq. 13 + Eq. 14).

    Parameters
    ----------
    f : float or ndarray
        Frequency [Hz]. Valid roughly 1e-5 Hz <= f <= 1 Hz.
    Tobs_yr : {0.5, 1.0, 2.0, 4.0}
        Mission duration used for the confusion-noise fit.
    include_confusion : bool
        If False, returns instrument-only sensitivity (Eq. 13 alone).
    """
    f = np.asarray(f, dtype=float)
    x = f / _F_STAR
    instrument = (10.0 / (3.0 * _L**2)
                  * (P_OMS(f) + 2.0 * (1.0 + np.cos(x)**2) * P_acc(f) / (2*np.pi*f)**4)
                  * (1.0 + 0.6 * x**2))
    if include_confusion:
        return instrument + S_confusion(f, Tobs_yr)
    return instrument


def ASD(f, **kwargs):
    """Amplitude spectral density sqrt(S_n) [1/sqrt(Hz)]."""
    return np.sqrt(S_n(f, **kwargs))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    f = np.geomspace(1e-5, 1.0, 4000)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(f, ASD(f, include_confusion=False), lw=1.2,
              label="Instrument only")
    for T in (0.5, 1.0, 2.0, 4.0):
        ax.loglog(f, ASD(f, Tobs_yr=T), lw=1.2, label=f"+ confusion, {T} yr")

    ax.set_xlabel("Frequency  $f$  [Hz]")
    ax.set_ylabel(r"ASD  $\sqrt{S_n(f)}$  [Hz$^{-1/2}$]")
    ax.set_title("LISA sky-averaged sensitivity (Robson, Cornish & Liu 2019)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_ylim(1e-21, 1e-14)
    fig.tight_layout()
    fig.savefig("lisa_sensitivity.png", dpi=140)
    print("Saved lisa_sensitivity.png")