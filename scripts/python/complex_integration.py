import numpy as np
import sympy as sp
import scipy.integrate as spi

# Define the symbolic variables for phi and theta
phi, theta = sp.symbols('phi theta')

# Define the associated Legendre polynomial P(1, 0, cos(theta))
legendre_polynomial = sp.assoc_legendre(1, 0, sp.cos(theta))

# Define your complex-valued function of phi and theta (modify this according to your specific function)
# For example, let's use a complex function f(phi, theta) = exp(i*phi) * sin(theta)
function_of_phi_theta = sp.exp(sp.I * phi) * sp.sin(theta)

# Define the product of the Legendre polynomial and the function
product_function = legendre_polynomial * function_of_phi_theta

# Convert the SymPy expression to a callable Python function
function_lambda = sp.lambdify((phi, theta), product_function, 'numpy')

print(product_function)
print(function_lambda(np.pi/2,np.pi/2))

# Define the range of phi and theta values
phi_values = np.linspace(0, 2 * np.pi, 100)  # Adjust the number of points as needed
theta_values = np.linspace(0, np.pi, 50)  # Adjust the number of points as needed

# Create a grid of phi and theta values
phi_mesh, theta_mesh = np.meshgrid(phi_values, theta_values)

# Calculate the Jacobian factor (r^2 * sin(theta))
jacobian = np.sin(theta_mesh)

# Calculate the product function values for each phi and theta value in the grid
product_values = function_lambda(phi_mesh, theta_mesh)

# Perform the double integral over the sphere
result = spi.simps(spi.simps(product_values * jacobian, dx=phi_values[1] - phi_values[0], axis=1), dx=theta_values[1] - theta_values[0])

print("Integral result:", result)