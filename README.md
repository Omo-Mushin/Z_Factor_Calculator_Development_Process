**Z Factor Calculator**
This software calculates the compressibility factor (Z factor) for gases using various correlations. The Z factor is essential for understanding the behavior of gases under different conditions of temperature and pressure.
**Features**
Dranchuk-Abu-Kassem (DAK) Correlation
Dranchuk-Purvis-Robinson (DPR) Correlation
Hall Yarborough Correlation
Papp Correlation

**Functions**
compute_Z_with_DAK(sp_gravity, temperature, pressure)
Calculates the Z factor using the Dranchuk-Abu-Kassem (DAK) correlation.

**Parameters:**

sp_gravity (float): Specific gravity of the gas.
temperature (float): Temperature in Rankine.
pressure (float): Pressure in psia.
Returns:

Z (float): Compressibility factor.
compute_Z_with_DPR(sp_gravity, temperature, pressure)
Calculates the Z factor using the Dranchuk-Purvis-Robinson (DPR) correlation.

**Parameters:**

sp_gravity (float): Specific gravity of the gas.
temperature (float): Temperature in Rankine.
pressure (float): Pressure in psia.
Returns:

Z (float): Compressibility factor.
**Notes**
Ensure that the specific gravity, temperature, and pressure values are within the valid ranges for the correlations to be accurate.
The valid ranges for the pseudo-reduced pressure (Ppr) and pseudo-reduced temperature (Tpr) are:
Ppr: [0.2, 30)
Tpr: (1, 3.0]
