# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 09:48:45 2024

@author: ieron
"""
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np

# Given parameters
R = 0.01  
LGD = 0.4  
r = 0.03  
T = 1  
N = 4  

# Time points for quarterly premiums
T_i = np.linspace(0, T, N+1)[1:]  
T_i_with_start = np.insert(T_i, 0, 0)  # Helping function to correct for Tmid_i
Tmid_i = (T_i_with_start[:-1] + T_i_with_start[1:]) / 2 

def hazard_rate(u, lambda_1):
    return lambda_1
    
def integral_lambda_u(T_i, lambda_1):
    result, _ = quad(hazard_rate, 0, T_i, args=(lambda_1))
    return result


def CDS_pricing(lambda_1):
    # Calculate T_{i-1} for use in the sums
    T_i_minus_1 = np.append(0, T_i[:-1])  # Prepend 0 for the initial value
    # Calculate Q(tau > T_i) and Q(tau > T_{i-1}) using the survival probability formula
    Q_tau_Ti = np.exp(-np.array([integral_lambda_u(t, lambda_1) for t in T_i]))
    Q_tau_Ti_minus_1 = np.exp(-np.array([integral_lambda_u(t, lambda_1) for t in T_i_minus_1]))
    # Sum_1 Calculation
    Sum_1 = np.sum(np.exp(-r * T_i) * (T_i - T_i_minus_1) * Q_tau_Ti)
    # Sum_2 Calculation
    Sum_2 = np.sum(np.exp(-r * Tmid_i) * (Q_tau_Ti_minus_1 - Q_tau_Ti) * (T_i - T_i_minus_1) / 2)
    # Sum_3 Calculation
    Sum_3 = np.sum(np.exp(-r * Tmid_i) * (Q_tau_Ti_minus_1 - Q_tau_Ti))
    # CDS Pricing Formula Incorporation
    PV_Premium_Leg = R * (Sum_1 + Sum_2)
    PV_Protection_Leg = LGD * Sum_3
    
    return PV_Premium_Leg - PV_Protection_Leg

# Solve for lambda_1
lambda_initial_guess = 0.01  # Initial guess for lambda_1
lambda_1_solution, = fsolve(lambda lambda_1: CDS_pricing(lambda_1) - 0, lambda_initial_guess)

lambda_1_solution

# Since lambda_1 is constant over the period, the integration of a constant over its interval
# results in the constant times the length of the interval. Thus, the average hazard rate
# can be calculated directly in this case.

# Define the hazard rate function, which is constant in this scenario
hazard_rate_function = lambda t: lambda_1_solution

# Integrate the hazard rate function over the interval from 0 to T (1 year)
average_hazard_rate, _ = quad(hazard_rate_function, 0, T)

# Since the hazard rate is constant, the integration over 1 year divided by 1 year
# will just yield the hazard rate itself.
average_hazard_rate


###############################################################################
#%% lambda_2
###############################################################################
# Update parameters
R1 = 0.01  
R3 = 0.011  
LGD = 0.4
r = 0.03
T1 = 1  
T3 = 3  
N = 4 * T3  # Quarterly payments over 3 years

# Time points for quarterly premiums, considering up to 3 years
T_i = np.linspace(0, T3, N+1)[1:]  # Exclude the 0 point to get the end of each quarter
T_i_with_start = np.insert(T_i, 0, 0)  # Helping function to correct for Tmid_i
Tmid_i = (T_i_with_start[:-1] + T_i_with_start[1:]) / 2 


def hazard_rate(u, lambda_1, lambda_2):
    if u <= 1:
        return lambda_1
    else:
        return lambda_2
    
def integral_lambda_u(T_i, lambda_1, lambda_2):
    result, _ = quad(hazard_rate, 0, T_i, args=(lambda_1, lambda_2))
    return result

integral_lambda_u_vectorized = np.vectorize(integral_lambda_u)

def CDS_pricing_lambda_1(lambda_1):
    # Simplified as lambda_2 is not needed for solving lambda_1
    Q_tau_Ti = np.exp(-lambda_1 * T_i[T_i <= 1])
    premium_leg = R1 * np.sum(np.exp(-r * T_i[T_i <= 1]) * (T_i[T_i <= 1] - np.append(0, T_i[T_i <= 1][:-1])) * Q_tau_Ti)
    protection_leg = LGD * np.sum(np.exp(-r * T_i[T_i <= 1]) * (1 - Q_tau_Ti))
    return premium_leg - protection_leg

lambda_1_solution, = fsolve(lambda lambda_1: CDS_pricing_lambda_1(lambda_1) - 0, 0.01)

def CDS_pricing_lambda_2(lambda_2):
    # Ensure lambda_1_solution is defined and accessible in this scope
    global lambda_1_solution
    
    # Apply integral_lambda_u function over T_i array using a list comprehension or similar approach
    Q_tau_Ti = np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2) for t in T_i]))
    
    # The rest of the function remains unchanged
    premium_leg = R3 * np.sum(np.exp(-r * T_i) * (np.diff(np.append(0, T_i)) * Q_tau_Ti))
    protection_leg = LGD * np.sum(np.exp(-r * T_i) * (np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2) for t in np.append(0, T_i[:-1])])) - Q_tau_Ti))
    return premium_leg - protection_leg

lambda_2_solution, = fsolve(lambda lambda_2: CDS_pricing_lambda_2(lambda_2) - 0, 0.01)  # Initial guess

# Average hazard rate for 1Y using lambda_1_solution
average_hazard_rate_1Y = integral_lambda_u(1, lambda_1_solution, 0) / 1

# Average hazard rate for 3Y using lambda_1_solution and lambda_2_solution
average_hazard_rate_3Y = integral_lambda_u(3, lambda_1_solution, lambda_2_solution) / 3
print(lambda_2_solution,average_hazard_rate_3Y)


###############################################################################
#%% lambda_3
###############################################################################
# Given parameters
R1 = 0.01  # 1Y CDS rate in decimal form
R3 = 0.011  # 3Y CDS rate in decimal form
R5 = 0.012  # 5Y CDS rate in decimal form
T1 = 1
T3 = 3
T5 = 5
N = 4 * T5  # Assuming quarterly payments over 5 years

# Time points for quarterly premiums, considering up to 5 years
T_i = np.linspace(0, T5, N+1)[1:]  # Exclude the 0 point to get the end of each quarter
T_i_with_start = np.insert(T_i, 0, 0)  # Helping function to correct for Tmid_i
Tmid_i = (T_i_with_start[:-1] + T_i_with_start[1:]) / 2 


def hazard_rate(u, lambda_1, lambda_2, lambda_3):
    if u <= 1:
        return lambda_1
    elif u <= 3:
        return lambda_2
    else:
        return lambda_3

def integral_lambda_u(T_i, lambda_1, lambda_2, lambda_3):
    result, _ = quad(hazard_rate, 0, T_i, args=(lambda_1, lambda_2, lambda_3))
    return result

def CDS_pricing_lambda_3(lambda_3):
    # Ensure lambda_1_solution and lambda_2_solution are defined and accessible
    Q_tau_Ti = np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3) for t in T_i]))
    
    premium_leg = R5 * np.sum(np.exp(-r * T_i) * (np.diff(np.append(0, T_i)) * Q_tau_Ti))
    protection_leg = LGD * np.sum(np.exp(-r * T_i) * (np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3) for t in np.append(0, T_i[:-1])])) - Q_tau_Ti))
    return premium_leg - protection_leg

lambda_3_solution, = fsolve(lambda lambda_3: CDS_pricing_lambda_3(lambda_3) - 0, 0.01)  # Initial guess


average_hazard_rate_1Y = integral_lambda_u(1, lambda_1_solution, lambda_2_solution, 0) / 1
average_hazard_rate_3Y = integral_lambda_u(3, lambda_1_solution, lambda_2_solution, 0) / 3
average_hazard_rate_5Y = integral_lambda_u(5, lambda_1_solution, lambda_2_solution, lambda_3_solution) / 5
print(lambda_3_solution,average_hazard_rate_5Y)


###############################################################################
#%% lambda_4
###############################################################################
# Given parameters
R1 = 0.01  # CDS rate for 1Y in decimal form (100 bps)
R3 = 0.011  # CDS rate for 3Y in decimal form (110 bps)
R5 = 0.012  # CDS rate for 5Y in decimal form (120 bps)
R7 = 0.012  # CDS rate for 7Y in decimal form (120 bps)
T1, T3, T5, T7 = 1, 3, 5, 7  # Maturities
N = 4 * T7  # Quarterly payments over 7 years

# Time points for quarterly premiums, considering up to 5 years
T_i = np.linspace(0, T7, N+1)[1:]  # Exclude the 0 point to get the end of each quarter
T_i_with_start = np.insert(T_i, 0, 0)  # Helping function to correct for Tmid_i
Tmid_i = (T_i_with_start[:-1] + T_i_with_start[1:]) / 2 

def hazard_rate(u, lambda_1, lambda_2, lambda_3, lambda_4):
    if u <= 1:
        return lambda_1
    elif u <= 3:
        return lambda_2
    elif u <= 5:
        return lambda_3
    else:
        return lambda_4

def integral_lambda_u(T_i, lambda_1, lambda_2, lambda_3, lambda_4):
    result, _ = quad(hazard_rate, 0, T_i, args=(lambda_1, lambda_2, lambda_3, lambda_4))
    return result


def CDS_pricing_lambda_4(lambda_4):
    # Assuming lambda_1_solution, lambda_2_solution, lambda_3_solution are available
    Q_tau_Ti = np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4) for t in T_i]))
    premium_leg = R7 * np.sum(np.exp(-r * T_i) * (np.diff(np.append(0, T_i)) * Q_tau_Ti))
    protection_leg = LGD * np.sum(np.exp(-r * T_i) * (np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4) for t in np.append(0, T_i[:-1])])) - Q_tau_Ti))
    return premium_leg - protection_leg

lambda_4_solution, = fsolve(lambda lambda_4: CDS_pricing_lambda_4(lambda_4) - 0, 0.01)  # Initial guess

average_hazard_rate_1Y = integral_lambda_u(1, lambda_1_solution, 0, 0, 0) / 1
average_hazard_rate_3Y = integral_lambda_u(3, lambda_1_solution, lambda_2_solution, 0, 0) / 3
average_hazard_rate_5Y = integral_lambda_u(5, lambda_1_solution, lambda_2_solution, lambda_3_solution, 0) / 5
average_hazard_rate_7Y = integral_lambda_u(7, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution) / 7
print(lambda_4_solution,average_hazard_rate_7Y)


###############################################################################
#%% lambda_5
###############################################################################
R1 = 0.01   # CDS rate for 1Y in decimal form (100 bps)
R3 = 0.011  # CDS rate for 3Y in decimal form (110 bps)
R5 = 0.012  # CDS rate for 5Y in decimal form (120 bps)
R7 = 0.012  # CDS rate for 7Y in decimal form (120 bps)
R10 = 0.0125  # CDS rate for 10Y in decimal form (125 bps)
T1, T3, T5, T7, T10 = 1, 3, 5, 7, 10  # Maturities
N = 4 * T10  # Quarterly payments over 10 years

# Time points for quarterly premiums, considering up to 5 years
T_i = np.linspace(0, T10, N+1)[1:]  # Exclude the 0 point to get the end of each quarter
T_i_with_start = np.insert(T_i, 0, 0)  # Helping function to correct for Tmid_i
Tmid_i = (T_i_with_start[:-1] + T_i_with_start[1:]) / 2 

def hazard_rate(u, lambda_1, lambda_2, lambda_3, lambda_4, lambda_5):
    if u <= 1:
        return lambda_1
    elif u <= 3:
        return lambda_2
    elif u <= 5:
        return lambda_3
    elif u <= 7:
        return lambda_4
    else:
        return lambda_5

def integral_lambda_u(T_i, lambda_1, lambda_2, lambda_3, lambda_4, lambda_5):
    result, _ = quad(hazard_rate, 0, T_i, args=(lambda_1, lambda_2, lambda_3, lambda_4, lambda_5))
    return result

def CDS_pricing_lambda_5(lambda_5):
    # Assuming lambda_1_solution, ..., lambda_4_solution are available
    Q_tau_Ti = np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5) for t in T_i]))
    premium_leg = R10 * np.sum(np.exp(-r * T_i) * (np.diff(np.append(0, T_i)) * Q_tau_Ti))
    protection_leg = LGD * np.sum(np.exp(-r * T_i) * (np.exp(-np.array([integral_lambda_u(t, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5) for t in np.append(0, T_i[:-1])])) - Q_tau_Ti))
    return premium_leg - protection_leg

lambda_5_solution, = fsolve(lambda lambda_5: CDS_pricing_lambda_5(lambda_5) - 0, 0.01)  # Initial guess

# Average hazard rate calculations for each maturity
average_hazard_rate_1Y = integral_lambda_u(1, lambda_1_solution, 0, 0, 0, 0) / 1
average_hazard_rate_3Y = integral_lambda_u(3, lambda_1_solution, lambda_2_solution, 0, 0, 0) / 3
average_hazard_rate_5Y = integral_lambda_u(5, lambda_1_solution, lambda_2_solution, lambda_3_solution, 0, 0) / 5
average_hazard_rate_7Y = integral_lambda_u(7, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, 0) / 7
average_hazard_rate_10Y = integral_lambda_u(10, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution) / 10
print(lambda_5_solution, average_hazard_rate_10Y)

###############################################################################
###############################################################################
###############################################################################
#%% END

def cumulative_default_probability(T, lambda_1, lambda_2, lambda_3, lambda_4, lambda_5):
    integral, _ = quad(hazard_rate, 0, T, args=(lambda_1, lambda_2, lambda_3, lambda_4, lambda_5))
    Q_tau_T = np.exp(-integral)
    return 1 - Q_tau_T

# Calculate cumulative default probabilities for each maturity
P_1Y = cumulative_default_probability(1, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution)
P_3Y = cumulative_default_probability(3, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution)
P_5Y = cumulative_default_probability(5, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution)
P_7Y = cumulative_default_probability(7, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution)
P_10Y = cumulative_default_probability(10, lambda_1_solution, lambda_2_solution, lambda_3_solution, lambda_4_solution, lambda_5_solution)

print("Cumulative Default Probability (1Y):", P_1Y)
print("Cumulative Default Probability (3Y):", P_3Y)
print("Cumulative Default Probability (5Y):", P_5Y)
print("Cumulative Default Probability (7Y):", P_7Y)
print("Cumulative Default Probability (10Y):", P_10Y)































