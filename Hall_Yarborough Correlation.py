# Hall Yarborough Correlation for Calculating Z
import numpy as np


def compute_reduced_rho(Ppr, t):
    rho = 0.0125 * Ppr * t * np.exp(-1.2 * ((1 - t) ** 2))
    return rho


def compute_X_values(Ppr, t):
    X1 = -0.06125 * Ppr * t * np.exp(-1.2 * ((1 - t) ** 2))
    X2 = 14.76 * t - 9.76 * (t ** 2) + 4.58 * (t ** 3)
    X3 = 90.7 * t - 242.2 * (t ** 2) + 42.4 * (t ** 3)
    X4 = 2.18 + 2.82 * t
    return X1, X2, X3, X4


def compute_reduced_rho_function(rho, x1, x2, x3, x4):
    rho_function = x1 + ((rho + rho ** 2 + rho ** 3 + rho ** 4) / ((1 - rho) ** 3)) - (x2 * rho ** 2) + (
                x3 * (rho ** x4))
    return rho_function


def compute_derivative_function(rho, x1, x2, x3, x4):
    df_rho_function = ((1 + 4 * rho + 4 * (rho ** 2) - 4 * (rho ** 3) + rho ** 4) / (1 - rho) ** 4) - (2 * x2 * rho) + (
                x3 * x4 * (rho ** (x4 - 1)))
    return df_rho_function


def compute_effective_reduced_density(Ppr, t, tol=1e-12, runs=200):
    rho = compute_reduced_rho(Ppr, t)
    x1, x2, x3, x4 = compute_X_values(Ppr, t)

    for i in range(runs):
        func = compute_reduced_rho_function(rho, x1, x2, x3, x4)
        df_func = compute_derivative_function(rho, x1, x2, x3, x4)

        if abs(func) <= tol:
            print(f'Found solution after {i} iterations')
            return rho

        if df_func == 0:
            print('Zero derivative. No solution found')
            return None

        try:
            rho = rho - (func / df_func)
        except ZeroDivisionError:
            print('Zero division error encountered')
            return None
        except Exception as e:
            print(f"Error encountered: {e}")
            return None

        if np.isnan(rho) or np.isinf(rho):
            print('Encountered NaN or infinity')
            return None

    print('Exceeded maximum iterations. No solution found')
    print(f'effective reduced density is {rho} after {i} iterations')
    return rho


def compute_Z_with_Hall_Yarborough(sp_gravity, temperature, pressure):
    Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
    Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)

    Ppr = pressure / Pc
    Tpr = temperature / Tc
    t = Tc / temperature

    if Tpr > 1:
        rho_effective = compute_effective_reduced_density(Ppr, t)
        if rho_effective is None:
            return None
        z = round((0.06125 * Ppr * t / rho_effective) * np.exp(1.12 * (1 - t) ** 2), 4)
        return z
    else:
        print('Pseudo-reduced temperature must be greater than 1')
        return None


print(compute_Z_with_Hall_Yarborough(0.854, 610, 2000))

