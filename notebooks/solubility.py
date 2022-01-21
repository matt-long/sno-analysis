#! /usr/bin/env python
import numpy as np


T0_Kelvin = 273.15


def _garcia_gordon_polynomial(
    S,
    T,
    A0=0.0,
    A1=0.0,
    A2=0.0,
    A3=0.0,
    A4=0.0,
    A5=0.0,
    B0=0.0,
    B1=0.0,
    B2=0.0,
    B3=0.0,
    C0=0.0,
):

    T_scaled = np.log((298.15 - T) / (273.15 + T))
    return np.exp(
        A0
        + A1 * T_scaled
        + A2 * T_scaled ** 2
        + A3 * T_scaled ** 3
        + A4 * T_scaled ** 4
        + A5 * T_scaled ** 5
        + S * (B0 + B1 * T_scaled + B2 * T_scaled ** 2 + B3 * T_scaled ** 3)
        + C0 * S ** 2
    )


def _umolkg_to_mmolm3(value, rho_ref):
    return value * rho_ref / 1000.0


def Ar(S, T, **kwargs):
    """
    Solubility of Ar in sea water
    INPUT:  (if S and T are not singular they must have same dimensions)
    S = salinity    [PSS]
    T = temperature [degree C]

    conc = solubility of Ar [µmol/kg]

    REFERENCE:
    Roberta Hamme and Steve Emerson, 2004.
    "The solubility of neon, nitrogen and argon in distilled water and seawater."
    Deep-Sea Research I, 51(11), p. 1517-1528.
    """

    # constants from Table 4 of Hamme and Emerson 2004
    return _garcia_gordon_polynomial(
        S,
        T,
        A0=2.79150,
        A1=3.17609,
        A2=4.13116,
        A3=4.90379,
        B0=-6.96233e-3,
        B1=-7.66670e-3,
        B2=-1.16888e-2,
    )


def Ne(S, T, **kwargs):
    """
    Solubility (saturation) of neon (Ne) in sea water
    at 1-atm pressure of air including saturated water vapor

    INPUT:  (if S and T are not singular they must have same dimensions)
    S = salinity    [PSS]
    T = temperature [degree C]

    OUTPUT:
    conc = solubility of Ne  [µmol/kg]

    REFERENCE:
    Roberta Hamme and Steve Emerson, 2004.
    "The solubility of neon, nitrogen and argon in distilled water and seawater."
    Deep-Sea Research I, 51(11), p. 1517-1528.
    """

    # convert from nmol/kg to umol/kg
    return _garcia_gordon_polynomial(
        S, T, A0=2.18156, A1=1.29108, A2=2.12504, B0=-5.94737e-3, B1=-5.13896e-3
    ) / 1000.0


def Xe(S, T, **kwargs):
    """
    Solubility (saturation) of xeon (Xe) in sea water
    at 1-atm pressure of air including saturated water vapor

    INPUT:  (if S and T are not singular they must have same dimensions)
    S = salinity    [PSS]
    T = temperature [degree C]

    OUTPUT:
    conc = solubility of Xe [µmol/kg]

    REFERENCE:
    R. Hamme fit to data of
    D. Wood and R. Caputi (1966) "Solubilities of Kr and Xe in fresh and sea water"
    U.S. Naval Radiological Defense Laboratory, Technical Report USNRDL-TR-988,
    San Francisco, CA, pp. 14.
    """

    return _garcia_gordon_polynomial(
        S, T, A0=-7.48679, A1=5.08637, A2=4.22243, B0=-8.15683e-3, B1=-1.20129e-3
    )



def N2(S, T, **kwargs):
    """
    Solubility (saturation) of nitrogen (N2) in sea water
    at 1-atm pressure of air including saturated water vapor

    INPUT:  (if S and T are not singular they must have same dimensions)
    S = salinity    [PSS]
    T = temperature [degree C]

    OUTPUT:
    conc = solubility of N2  [µmol/kg]

    REFERENCE:
    Roberta Hamme and Steve Emerson, 2004.
    "The solubility of neon, nitrogen and argon in distilled water and seawater."
    Deep-Sea Research I, 51(11), p. 1517-1528.
    """

    return _garcia_gordon_polynomial(
        S,
        T,
        A0=6.42931,
        A1=2.92704,
        A2=4.32531,
        A3=4.69149,
        B0=-7.44129e-3,
        B1=-8.02566e-3,
        B2=-1.46775e-2,
    )


def O2(S, T, **kwargs):
    """
    Solubility of O2 in sea water
    INPUT:
    S = salinity    [PSS]
    T = temperature [degree C]

    conc = solubility of O2 [µmol/kg]

    REFERENCE:
    Hernan E. Garcia and Louis I. Gordon, 1992.
    "Oxygen solubility in seawater: Better fitting equations"
    Limnology and Oceanography, 37, pp. 1307-1312.
    
    Coefficients are in Table 1, using the fit to the Benson & Krause (1984) data.
    """

    # constants from Table 4 of Hamme and Emerson 2004
    return _garcia_gordon_polynomial(
        S,
        T,
        A0=5.80871,
        A1=3.20291,
        A2=4.17887,
        A3=5.10006,
        A4=-9.86643e-2,
        A5=3.80369,
        B0=-7.01577e-3,
        B1=-7.70028e-3,
        B2=-1.13864e-2,
        B3=-9.51519e-3,
        C0=-2.75915e-7,
    )



if __name__ == "__main__":
    """
    Check values:
    - Ne, N2, Ar: Table 4 of Hamme and Emerson (2004)
      https://doi.org/10.1016/j.dsr.2004.06.009
    
    - O2: Table 1 of Garcia and Gordon (1992)
      https://doi.org/10.4319/lo.1992.37.6.1307
    """
    
    S = 35.0
    T = 10.0
    
    check_values = {
        "Ar": {"check": 13.4622, "func": Ar},
        "Ne": {"check": 7.34121 / 1000.0, "func": Ne},
        "N2": {"check": 500.885, "func": N2},
        "O2": {"check": 274.610, "func": O2},
    }

    for gas, d in check_values.items():        
        print(gas)
        check_val = d["check"]
        func_val = d["func"](S, T)
        np.testing.assert_almost_equal(func_val, check_val, )
        diff = check_val - func_val
        print(f"   check: {check_val:0.6f}")
        print(f"    func: {func_val:0.6f}")
        print(f"    diff: {diff:0.6f}")
        print()
    