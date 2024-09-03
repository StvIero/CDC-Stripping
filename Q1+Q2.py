# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 14:01:58 2024

@author: ieron
"""
from scipy.optimize import fsolve
from scipy.integrate import quad
import numpy as np

# Constants
r = 0.03  # Annual interest rate
LGD = 0.40  # Loss Given Default
R = 0.01  # CDS Rate (1% or 100 bps)
T = 1  # Maturity in years
N = 4  # Quarterly payments
payment_times = np.linspace(0, T, N+1)[1:]  # Excluding 0, payment times at 0.25, 0.5, 0.75, 1

# Function to perform explicit integration of lambda_1 over the period up to time
def integrate_lambda(lambd, time):
    result, _ = quad(lambda u: lambd, 0, time)
    return result

# Function to calculate survival probabilities using explicit integration
def calculate_survival_probabilities(lambd):
    return [np.exp(-integrate_lambda(lambd, t)) for t in payment_times]
'''
For T0.25​: 0.992
For T0.5​: 0.983
For T0.75​: 0.975
For T1​: 0.967
'''
# Function to calculate the value of the CDS contract given a hazard rate lambda
def cds_valuation(lambd):
    survival_probabilities = calculate_survival_probabilities(lambd)
    Q_diff = np.hstack(([1], survival_probabilities[:-1])) - np.array(survival_probabilities)
    mid_times = (payment_times[:-1] + payment_times[1:]) / 2  # Mid points for integral calculation

    # Premium leg calculation
    premium_leg = np.sum(np.exp(-r * payment_times) * (payment_times[1] - payment_times[0]) * np.array(survival_probabilities))
    premium_leg += np.sum(np.exp(-r * mid_times) * Q_diff[1:] * mid_times)
    
    # Protection leg calculation
    protection_leg = np.sum(np.exp(-r * mid_times) * Q_diff[1:])
    
    # CDS value formula
    cds_value = R * premium_leg - LGD * protection_leg
    return cds_value

# Solve for lambda_1
lambda_initial_guess = 0.01  # Initial guess for the hazard rate
lambda_1_solution = fsolve(lambda lambd: cds_valuation(lambd), lambda_initial_guess)[0]

# Display the solution
lambda_1_solution

###############################################################################
# Adjusted parameters for the new setup
R3 = 0.011  # CDS Rate for 3 years (110 bps or 1.1% as a decimal)
T3 = 3  # Maturity in years

# Function to calculate the piece-wise constant hazard rate lambda_u over the period up to time
def integrate_piecewise_lambda(lambda_1, lambda_2, time):
    if time <= 1:
        result, _ = quad(lambda u: lambda_1, 0, time)
    else:
        result_1, _ = quad(lambda u: lambda_1, 0, 1)
        result_2, _ = quad(lambda u: lambda_2, 1, time)
        result = result_1 + result_2
    return result

# Updated function to calculate survival probabilities using explicit integration for the new piece-wise lambda_u
def calculate_survival_probabilities_updated(lambda_1, lambda_2):
    times = np.linspace(0, T3, 3*N+1)[1:]  # Adjust to include all quarterly payment times over 3 years
    return [np.exp(-integrate_piecewise_lambda(lambda_1, lambda_2, t)) for t in times]

# Function to calculate the CDS valuation given lambda_1 and solving for lambda_2
def cds_valuation_updated(lambda_2):
    survival_probabilities = calculate_survival_probabilities_updated(lambda_1_solution, lambda_2)
    Q_diff = np.hstack(([1], survival_probabilities[:-1])) - np.array(survival_probabilities)
    times = np.linspace(0, T3, 3*N+1)[1:]  # Payment times updated for 3 years
    mid_times = (times[:-1] + times[1:]) / 2

    # Premium and protection leg calculations adjusted for 3 years and piece-wise lambda
    premium_leg = np.sum(np.exp(-r * times) * (times[1] - times[0]) * np.array(survival_probabilities))
    premium_leg += np.sum(np.exp(-r * mid_times) * Q_diff[1:] * (times[1] - times[0]) / 2)
    protection_leg = np.sum(np.exp(-r * mid_times) * Q_diff[1:])
    
    # Adjusted CDS value formula for 3 years
    cds_value = R3 * premium_leg - LGD * protection_leg
    return cds_value

# Solve for lambda_2 with lambda_1 known
lambda_2_solution = fsolve(lambda lambda_2: cds_valuation_updated(lambda_2), 0.01)[0]

# Calculate the average hazard rate for the 3-year period
def average_hazard_rate(lambda_1, lambda_2, T3):
    integral, _ = quad(lambda u: lambda_1 if u <= 1 else lambda_2, 0, T3)
    return integral / T3

average_hazard_rate_3Y = average_hazard_rate(lambda_1_solution, lambda_2_solution, T3)

lambda_2_solution, average_hazard_rate_3Y

###############################################################################
# Adjusted parameters for the 5-year setup
R5 = 0.012  # CDS Rate for 5 years (120 bps or 1.2% as a decimal)
T5 = 5  # Maturity in years

# Function to calculate the piece-wise constant hazard rate lambda_u over the period up to time for the 5-year period
def integrate_piecewise_lambda_5Y(lambda_1, lambda_2, lambda_3, time):
    if time <= 1:
        result, _ = quad(lambda u: lambda_1, 0, time)
    elif time <= 3:
        result_1, _ = quad(lambda u: lambda_1, 0, 1)
        result_2, _ = quad(lambda u: lambda_2, 1, time)
        result = result_1 + result_2
    else:
        result_1, _ = quad(lambda u: lambda_1, 0, 1)
        result_2, _ = quad(lambda u: lambda_2, 1, 3)
        result_3, _ = quad(lambda u: lambda_3, 3, time)
        result = result_1 + result_2 + result_3
    return result

# Updated function to calculate survival probabilities for the new piece-wise lambda_u over 5 years
def calculate_survival_probabilities_5Y(lambda_1, lambda_2, lambda_3):
    times = np.linspace(0, T5, 5*N+1)[1:]  # Adjust to include all quarterly payment times over 5 years
    return [np.exp(-integrate_piecewise_lambda_5Y(lambda_1, lambda_2, lambda_3, t)) for t in times]

# Function to calculate the CDS valuation given lambda_1, lambda_2 and solving for lambda_3
def cds_valuation_5Y(lambda_3):
    survival_probabilities = calculate_survival_probabilities_5Y(lambda_1_solution, lambda_2_solution, lambda_3)
    Q_diff = np.hstack(([1], survival_probabilities[:-1])) - np.array(survival_probabilities)
    times = np.linspace(0, T5, 5*N+1)[1:]  # Payment times updated for 5 years
    mid_times = (times[:-1] + times[1:]) / 2

    # Premium and protection leg calculations adjusted for 5 years and piece-wise lambda
    premium_leg = np.sum(np.exp(-r * times) * (times[1] - times[0]) * np.array(survival_probabilities))
    premium_leg += np.sum(np.exp(-r * mid_times) * Q_diff[1:] * (times[1] - times[0]) / 2)
    protection_leg = np.sum(np.exp(-r * mid_times) * Q_diff[1:])
    
    # Adjusted CDS value formula for 5 years
    cds_value = R5 * premium_leg - LGD * protection_leg
    return cds_value

# Solve for lambda_3 with lambda_1 and lambda_2 known
lambda_3_solution = fsolve(lambda lambda_3: cds_valuation_5Y(lambda_3), 0.01)[0]

# Calculate the average hazard rate for the 5-year period
def average_hazard_rate_5Y(lambda_1, lambda_2, lambda_3, T5):
    integral, _ = quad(lambda u: lambda_1 if u <= 1 else lambda_2 if u <= 3 else lambda_3, 0, T5)
    return integral / T5

average_hazard_rate_5Y_solution = average_hazard_rate_5Y(lambda_1_solution, lambda_2_solution, lambda_3_solution, T5)

lambda_3_solution, average_hazard_rate_5Y_solution



















