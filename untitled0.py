# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:33:33 2023

@author: ISA
"""

import matplotlib.pyplot as plt

# Data
years = [2017, 2018, 2019, 2020, 2021]
values = [8.4, 7.6, 7.7, 7.3, 7.8]

# Additional data
additional_years = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
additional_values = [13.4, 14.1, 13.7, 12.3, 10.9, 9.6, 9.3]

# Combine old and new data
all_years = additional_years + years
all_values = additional_values + values

# High-resolution figure
plt.figure(figsize=(12,7), dpi=150)
plt.plot(all_years, all_values, marker='o', linestyle='-', color='blue')

plt.title('Incidence (cases per 100,000)')
plt.xlabel('Year')
plt.ylabel('Value')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xticks(all_years)  # Display all years on x-axis
plt.ylim(0, max(all_values) + 1)
plt.tight_layout()

# Show plot
plt.show()



#############

import pandas as pd

# Load the uploaded CSV file
data = pd.read_csv('Book1.csv')
# Convert the 'UK Born number of notifications' and 'Non-UK Born number of notifications' columns to integers
data['UK Born number of notifications'] = data['UK Born number of notifications'].str.replace(',', '').astype(int)
data['Non-UK Born number of notifications'] = data['Non-UK Born number of notifications'].str.replace(',', '').astype(int)

# Plotting the data
plt.figure(figsize=(12,7), dpi=150)
plt.plot(data['Year'], data['UK Born number of notifications'], marker='o', linestyle='-', color='blue', label='UK Born')
plt.plot(data['Year'], data['Non-UK Born number of notifications'], marker='o', linestyle='-', color='red', label='Non-UK Born')

plt.title('TB Notifications by Place of Birth')
plt.xlabel('Year')
plt.ylabel('Number of Notifications')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xticks(data['Year'])  # Display all years on x-axis
plt.ylim(0, max(data['Non-UK Born number of notifications'].max(), data['UK Born number of notifications'].max()) + 1000)
plt.legend()
plt.tight_layout()

# Show plot
plt.show()