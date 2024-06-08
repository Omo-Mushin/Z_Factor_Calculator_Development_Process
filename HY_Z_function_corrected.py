import numpy as np


def compute_reduced_rho(Ppr, Tpr):
    if Ppr <= 0 or Tpr <= 0:
        raise ValueError("Ppr and Tpr must be positive.")
    rho_1 = 0.0125 * Ppr * (1 / Tpr) * np.exp(-1.2 * (1 / Tpr - 1) ** 2)
    return rho_1


def calculate_coefficients(t):
    if t <= 0:
        raise ValueError("t must be positive.")
    X1 = 0.06125 * t * np.exp(-1.2 * (1 - t) ** 2)
    X2 = t * (14.76 - 9.76 * t + 4.58 * t ** 2)
    X3 = t * (90.7 - 242.2 * t + 42.4 * t ** 2)
    X4 = 2.18 + 2.82 * t
    return X1, X2, X3, X4


def compute_reduced_rho_function(rho, X1, X2, X3, X4, Ppr):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")
    rho_function = (
                -X1 * Ppr + (rho + rho ** 2 + rho ** 3 - rho ** 4) / (1 - rho) ** 3 - X2 * rho ** 2 + X3 * rho ** X4)
    return rho_function


def compute_derivative_function(rho, X2, X3, X4):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")
    df_rho_function = (1 + 4 * rho + 4 * rho ** 2 - 4 * rho ** 3 + rho ** 4) / (
                1 - rho) ** 4 - 2 * X2 * rho + X3 * X4 * rho ** (X4 - 1)
    return df_rho_function


def compute_effective_reduced_density(Ppr, Tpr, tol=1e-13, verb=False):
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")
    t = 1 / Tpr
    X1, X2, X3, X4 = calculate_coefficients(t)
    rho_1 = compute_reduced_rho(Ppr, Tpr)
    delta = 1
    i = 1  # iterations

    while True:
        abs_rho_func = abs(compute_reduced_rho_function(rho_1, X1, X2, X3, X4, Ppr))
        if np.isnan(abs_rho_func):
            rho_1 += 0.1
        elif abs_rho_func < tol:
            break
        elif np.isnan(rho_1):
            raise ValueError("rho_1 is NA")
        else:
            rho_new = rho_1 - compute_reduced_rho_function(rho_1, X1, X2, X3, X4, Ppr) / compute_derivative_function(
                rho_1, X2, X3, X4)
            delta = abs(rho_1 - rho_new)

            if delta < tol:
                break
            rho_1 = rho_new
            i += 1

    if verb:
        print(f"number of iterations: {i}\n Pseudo-reduced pressure: {Ppr}\n delta: {delta} \
                effective reduced density: {rho_1} \n Pseudo-reduced temperature: {Tpr}")
    return rho_1, X1


def compute_Z_Hall_Yarborough(sp_gravity, temperature, pressure, tol=1e-13, verb=False):
    if sp_gravity <= 0 or temperature <= 0 or pressure <= 0:
        raise ValueError("Specific gravity, temperature, and pressure must be positive.")

    Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
    Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)

    Ppr = round((pressure / Pc), 2)
    Tpr = round((temperature / Tc), 2)

    if Tpr > 1:
        rho_effective, X1 = compute_effective_reduced_density(Ppr, Tpr, tol, verb)
        if rho_effective is None:
            return None
        z = round((X1 * Ppr / rho_effective), 4)
        return z
    else:
        raise ValueError('Pseudo-reduced temperature must be greater than 1')

    if verb:
        print(z)
    return z


# Example
print(compute_Z_Hall_Yarborough(0.854, 610, 2000))
