# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:12:19 2024

@author: ieron
"""

import math

# Given data
maturity_years = [1, 3, 5, 7, 10]
cds_rates = [100, 110, 120, 120, 125]
lgd = 0.4  # Loss Given Default

# Calculate the average hazard rates
average_hazard_rates = [(rate / lgd)/10000 for rate in cds_rates] 


###############################################################################
# Directly use average hazard rate for the first year
forward_hazard_rates = [(0, 1, average_hazard_rates[0])]  

# Calculate forward hazard rates for the remaining intervals
for i in range(1, len(maturity_years)):
    T_i = maturity_years[i]
    T_im1 = maturity_years[i-1]
    lambda_average_Ti = average_hazard_rates[i]
    lambda_average_Tim1 = average_hazard_rates[i-1]
    lambda_forward = (lambda_average_Ti * T_i - lambda_average_Tim1 * T_im1) / (T_i - T_im1)
    forward_hazard_rates.append((T_im1, T_i, lambda_forward))
'''
I use 
forward_lambda(Ti, Ti-1) = (average_lambda(Ti) * Ti - average_lambda(Ti-1) * Ti-1 ) / (Ti-Ti-1)  
to calculate the forward hazard rates
'''


###############################################################################
cumulative_survival_probability = 1
cumulative_default_probabilities = []

for interval in forward_hazard_rates:
    T_im1, T_i, lambda_forward = interval
    interval_length = T_i - T_im1
    # Calculate interval survival probability
    interval_survival_probability = math.exp(-lambda_forward * interval_length)
    # Update cumulative survival probability
    cumulative_survival_probability *= interval_survival_probability
    # Convert to cumulative default probability
    cumulative_default_probability = 1 - cumulative_survival_probability
    cumulative_default_probabilities.append((T_i, cumulative_default_probability))

'''
I use 
P_default(Ti, Ti-1) = 1 - exp[ -forward_lambda * (Ti-Ti-1) ] 
to calculate cumulative probability of default

This gives the cumulative probability of default up to each maturity, 
accounting for the compounding effect of survival (and, by complement, default)
probabilities over time.
'''


###############################################################################
# Print results
print("Average Hazard Rates:")
for maturity, rate in zip(maturity_years, average_hazard_rates):
    print(f"{maturity} years: {rate:.6f}")

print("\nForward Hazard Rates:")
for interval in forward_hazard_rates:
    print(f"Between {interval[0]} and {interval[1]} years: {interval[2]:.6f}")

print("\nCumulative Default Probabilities:")
for maturity, cum_prob in cumulative_default_probabilities:
    print(f"Up to {maturity} years: {cum_prob:.6f}")