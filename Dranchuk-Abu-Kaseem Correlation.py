# Computing Z using Dranchuk-Abu-Kassem

import numpy as np

# Coefficients
A1 = 0.3265
A2 = -1.0700
A3 = -0.5339
A4 = 0.01569
A5 = -0.05165
A6 = 0.5475
A7 = -0.7361
A8 = 0.1844
A9 = 0.1056
A10 = 0.6134
A11 = 0.7210

def compute_reduced_density(Tpr, Ppr):
    rho = 0.27 * Ppr / Tpr
    return rho

def compute_rho_function(Tpr, Ppr, rho):
    R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
    R2 = 0.27 * Ppr / Tpr
    R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
    R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
    R5 = A10 / (Tpr ** 3)

    rho_function = R1 * rho - (R2 / rho) + (R3 * (rho ** 2)) - (R4 * (rho ** 5)) + \
                   ((R5 * (1 + (A11 * (rho ** 2)))) * np.exp(-A11 * (rho ** 2))) + 1

    return rho_function

def compute_derivative_function(Tpr, Ppr, rho):
    R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
    R2 = 0.27 * Ppr / Tpr
    R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
    R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
    R5 = A10 / (Tpr ** 3)

    derivative_rho_function = R1 - (R2 / (rho ** 2)) + (2 * R3 * rho) - (5 * R4 * (rho ** 4)) + \
                              (2 * R5 * rho * np.exp(-A11 * (rho ** 2)) * ((1 + 2 * A11 * (rho ** 3)) - ((A11 * rho * 2) * (1 + (A11 * (rho ** 2))))))

    return derivative_rho_function

def compute_effective_reduced_density(Tpr, Ppr, tol=1e-12, runs=200):
    rho = compute_reduced_density(Tpr, Ppr)

    for i in range(runs):
        func = compute_rho_function(Tpr, Ppr, rho)
        df_func = compute_derivative_function(Tpr, Ppr, rho)

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
    return rho

def compute_Z_with_DAK(sp_gravity, temperature, pressure):
    Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
    Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)

    Ppr = pressure / Pc
    Tpr = temperature / Tc

    # print(Ppr)
    # print(Tpr)
    if (Ppr >= 0.2 and Ppr < 30) and (Tpr > 1 and Tpr <= 3.0):
        rho_effective = compute_effective_reduced_density(Tpr, Ppr)
        if rho_effective is None:
            return None
        Z = round((0.27 * Ppr / (rho_effective * Tpr)), 4)
        return Z

print(compute_Z_with_DAK(0.854, 610, 2000))


# import numpy as np
#
# # Coefficients for the DAK correlation
# A1 = 0.3265
# A2 = -1.0700
# A3 = -0.5339
# A4 = 0.01569
# A5 = -0.05165
# A6 = 0.5475
# A7 = -0.7361
# A8 = 0.1844
# A9 = 0.1056
# A10 = 0.6134
# A11 = 0.7210
#
# def compute_rho_function(Tpr, Ppr, rho):
#     R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
#     R2 = 0.27 * Ppr / Tpr
#     R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
#     R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
#     R5 = A10 / (Tpr ** 3)
#
#     rho_function = R1 * rho - (R2 / rho if rho != 0 else 0) + (R3 * (rho ** 2)) - (R4 * (rho ** 5)) + \
#                    ((R5 * (1 + (A11 * (rho ** 2)))) * np.exp(-A11 * (rho ** 2))) + 1
#
#     return rho_function
#
# def compute_derivative_function(Tpr, Ppr, rho):
#     R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
#     R2 = 0.27 * Ppr / Tpr
#     R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
#     R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
#     R5 = A10 / (Tpr ** 3)
#
#     derivative_rho_function = R1 - (R2 / (rho ** 2) if rho != 0 else 0) + (2 * R3 * rho) - (5 * R4 * (rho ** 4)) + \
#                               (2 * R5 * rho * np.exp(-A11 * (rho ** 2)) * ((1 + 2 * A11 * (rho ** 3)) - ((A11 * rho * 2) * (1 + (A11 * (rho ** 2))))))
#
#     return derivative_rho_function
#
# def compute_effective_reduced_density(Tpr, Ppr, tol=0, max_iter=100):
#     rho = 0.27 * Ppr / Tpr  # Initial guess for reduced density
#
#     for i in range(max_iter):
#         func = compute_rho_function(Tpr, Ppr, rho)
#         df_func = compute_derivative_function(Tpr, Ppr, rho)
#
#         if abs(func) <= tol:
#             print(f'Found solution after {i} iterations')
#             return rho
#
#         if df_func == 0:
#             print('Zero derivative. No solution found')
#             return None
#
#         try:
#             rho = rho - (func / df_func)
#         except ZeroDivisionError:
#             print('Zero division error encountered')
#             return None
#         except Exception as e:
#             print(f"Error encountered: {e}")
#             return None
#
#         if np.isnan(rho) or np.isinf(rho):
#             print('Encountered NaN or infinity')
#             return None
#
#     print('Exceeded maximum iterations. No solution found')
#     return rho
#
# def compute_Z_with_DAK(sp_gravity, temperature, pressure):
#     Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
#     Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)
#
#     Ppr = pressure / Pc
#     Tpr = temperature / Tc
#
#     if (Ppr >= 0.2 and Ppr < 30) and (Tpr > 1 and Tpr <= 3.0):
#         rho_effective = compute_effective_reduced_density(Tpr, Ppr)
#         if rho_effective is None:
#             return None
#         Z = round((0.27 * Ppr / (rho_effective * Tpr)), 4)
#         return Z
#
# # Example usage
# sp_gravity = 0.854
# temperature = 610  # in Rankine
# pressure = 2000  # in psia
#
# Z = compute_Z_with_DAK(sp_gravity, temperature, pressure)
# print(f'Compressibility Factor Z: {Z}')
