# water_props.py
"""Single-phase water properties from IAPWS-IF97. T in deg C, P in Pa."""
from iapws import IAPWS97
import config as c


def props(T_C, P_Pa=c.P_SYSTEM):
    """Return dict: rho [kg/m^3], mu [Pa s], k [W/m K], cp [J/kg K], Pr [-].
    Clamps T below saturation (minus 0.5 K) to stay in single-phase liquid region."""
    T_max = T_sat(P_Pa) - 0.5
    if T_C > T_max:
        T_C = T_max
    w = IAPWS97(P=P_Pa / 1e6, T=T_C + 273.15)
    return {"rho": w.rho, "mu": w.mu, "k": w.k,
            "cp": w.cp * 1e3, "pr": w.Prandt}


def T_sat(P_Pa=c.P_SYSTEM):
    """Saturation temperature [deg C] at given pressure."""
    return IAPWS97(P=P_Pa / 1e6, x=0).T - 273.15
