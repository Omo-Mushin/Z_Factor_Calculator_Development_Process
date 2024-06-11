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

fluid_type = input('Is the fluid a natural gas or condensate: '    )
def compute_reduced_density(Tpr, Ppr):
    if Tpr <= 0 or Ppr <= 0:
        raise ValueError("Tpr and Ppr must be positive.")
    rho = 0.27 * Ppr / Tpr
    return rho


def compute_R_Values(Tpr, Ppr):
    if Tpr <= 0:
        raise ValueError("Tpr must be positive.")

    R1 = A1 + (A2 / Tpr) + (A3 / (Tpr ** 3)) + (A4 / (Tpr ** 4)) + (A5 / (Tpr ** 5))
    R2 = 0.27 * Ppr / Tpr
    R3 = A6 + (A7 / Tpr) + (A8 / (Tpr ** 2))
    R4 = A9 * ((A7 / Tpr) + (A8 / (Tpr ** 2)))
    R5 = A10 / (Tpr ** 3)

    return R1, R2, R3, R4, R5


def compute_rho_function(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    R1, R2, R3, R4, R5 = compute_R_Values(Tpr, Ppr)

    rho_function = R1 * rho - R2 / rho + R3 * rho ** 2 - R4 * rho ** 5 + \
                   R5 * rho ** 2 * (1 + A11 * rho ** 2) * np.exp(-A11 * rho ** 2) + 1

    return rho_function


def compute_derivative_function(Tpr, Ppr, rho):
    if rho < 0 or rho >= 1:
        raise ValueError("rho must be in the range [0, 1).")

    R1, R2, R3, R4, R5 = compute_R_Values(Tpr, Ppr)

    derivative_rho_function = R1 + R2 / rho ** 2 + 2 * R3 * rho - 5 * R4 * rho ** 4 + \
                              2 * R5 * rho * np.exp(-A11 * rho ** 2) * \
                              ((1 + 2 * A11 * rho ** 3) - A11 * rho ** 2 * (1 + A11 * rho ** 2))

    return derivative_rho_function


def compute_effective_reduced_density(Tpr, Ppr, tol=1e-12, runs=200):
    if tol <= 0:
        raise ValueError("Tolerance must be positive.")

    rho = compute_reduced_density(Tpr, Ppr)

    rho_1 = rho
    i = 1

    while i <= runs:
        if abs(compute_rho_function(Tpr, Ppr, rho_1)) < tol:
            return rho_1

        rho_new = rho_1 - compute_rho_function(Tpr, Ppr, rho_1) / compute_derivative_function(Tpr, Ppr, rho_1)
        delta = abs(rho_1 - rho_new)

        if delta < tol:
            return rho_new  # Exit the loop if ideal delta found

        rho_1 = rho_new
        i += 1

        if np.isnan(rho_1) or np.isinf(rho_1):
            raise ValueError('Encountered NaN or infinity in density computation.')

    raise ValueError('Exceeded maximum iterations. No solution found.')
    return rho_1

def compute_Z_with_DAK(sp_gravity, temperature, pressure):
    if sp_gravity <= 0 or temperature <= 0 or pressure <= 0:
        raise ValueError("Specific gravity, temperature, and pressure must be positive.")

    if fluid_type.upper() == 'NATURAL GAS':
        Tc = 168 + 325 * sp_gravity - 12.5 * (sp_gravity ** 2)
        Pc = 677 + 15 * sp_gravity - 37.5 * (sp_gravity ** 2)

    if fluid_type.upper() == 'CONDENSATES':
        Tc = 187 + 330 * sp_gravity - 71.5 * (sp_gravity ** 2)
        Pc = 706 - 51.7 * sp_gravity - 11.1 * (sp_gravity ** 2)

    Ppr = pressure / Pc
    Tpr = temperature / Tc

    if (0.2 <= Ppr < 30) and (1 < Tpr <= 3.0):
        rho_effective = compute_effective_reduced_density(Tpr, Ppr)
        if rho_effective is None:
            return None
        Z = round((0.27 * Ppr / (rho_effective * Tpr)), 4)
        return Z
    else:
        raise ValueError('Ppr must be in the range [0.2, 30) and Tpr in the range (1, 3.0].')

sp_gravity = float(input('Enter specific gravity:   '))
temperature = float(input('Enter Temperature:    '))
pressure = float(input('Enter the pressure:    '))

# Example usage
print(compute_Z_with_DAK(0.854, 610, 2000))
