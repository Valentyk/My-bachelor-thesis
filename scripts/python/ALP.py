import numpy as np
from sympy import assoc_legendre as alp
from sympy.abc import x, theta
import sympy as sp
import scipy.integrate as spi
import math

from sphere_plot import plot

#######################################################################################

def normalize(l):
    Norm = []
    for l in range(l+1):
        new_line_of_norm = []
        for m in range(-l,l+1):
            function_expr = alp(l,m,sp.cos(phi_sym))*sp.exp(sp.I*m*theta_sym)*alp(l,m,sp.cos(phi_sym))*sp.exp(-sp.I*m*theta_sym)   # eigen function*complex conjugate eigen function

            function_lambda = sp.lambdify((phi_sym, theta_sym), function_expr, 'numpy')  # Convert the symbolic function to a callable Python function
            function_values = function_lambda(phi_mesh, theta_mesh) # Calculate the function values for each phi and theta value in the grid
            
            jacobian = np.sin(phi_mesh)

            result = spi.simps(spi.simps(function_values * jacobian, dx=phi_values[1] - phi_values[0], axis=1), dx=theta_values[1] - theta_values[0]) # Integrate over the sphere
            new_line_of_norm.append(result**(1/2))

            #print(f"Integral result ‖P({l},{m})‖**2:", result)
            #print(f"Normalization of P({l}, {m}) polynomial ‖P({l},{m})‖:", result**(1/2))
        Norm.append(new_line_of_norm)
    
    #print(Norm)
    
    return(Norm)

def eigen_components(L, g):
    components = []
    for l in range(L+1):
        new_line_components = []
        for m in range(-l,l+1):
            eigen_func = alp(l,m,sp.cos(phi_sym))*sp.exp(sp.I*m*theta_sym)/norm[l][m+l]
            product_function = g*eigen_func
            function_lambda = sp.lambdify((phi_sym, theta_sym), product_function, "numpy")

            #print(product_function)

            function_values = function_lambda(phi_mesh, theta_mesh)
        
            if (l == 0) and (m == 0):
                print(product_function)
                #result = 0
                result = (spi.simps(spi.simps(function_values * jacobian, dx=phi_values[1] - phi_values[0], axis=1), dx=theta_values[1] - theta_values[0]))
                print(result)
            else:
                result = (spi.simps(spi.simps(function_values * jacobian, dx=phi_values[1] - phi_values[0], axis=1), dx=theta_values[1] - theta_values[0]))/(l*(l+1))

            new_line_components.append(result)
            print(f" \"Eigen\" components are (l = {l}, m = {m}) =", result)
        components.append(new_line_components)
    return(components)

def approx_func(L, g):
    components = eigen_components(L, g)
    v = 0
    for l in range(L+1):
        for m in range(-l,l+1):
            v += components[l][m+l]*alp(l,m,sp.cos(phi_sym))*sp.exp(sp.I*m*theta_sym)/(norm[l][m+l])
    return(v)

def norm_control(L):
    control = []
    for l in range(L+1):
        control_new_line = []
        for m in range(-l,l+1):
            function_expr = alp(l,m,sp.cos(phi_sym))*sp.exp(sp.I*m*theta_sym)*alp(l,m,sp.cos(phi_sym))*sp.exp(-sp.I*m*theta_sym)/(norm[l][m+l]**2)
            function_lambda = sp.lambdify((phi_sym, theta_sym), function_expr, 'numpy')  
            function_values = function_lambda(phi_mesh, theta_mesh) 
            jacobian = np.sin(phi_mesh)

            result = spi.simps(spi.simps(function_values * jacobian, dx=phi_values[1] - phi_values[0], axis=1), dx=theta_values[1] - theta_values[0])   
            print(result, norm[l][m+l])
            control_new_line.append(result)
        control.append(control_new_line)
    return(control)

#######################################################################################

N_integration = 100
N_vykres = 100
l = 5

phi_values = np.linspace(0, np.pi, int(N_integration/2))  # Adjust the number of points as needed
theta_values = np.linspace(0, 2*np.pi, int(N_integration))  # Adjust the number of points as needed
phi_mesh, theta_mesh = np.meshgrid(phi_values, theta_values)

jacobian = np.sin(phi_mesh)

phi_sym, theta_sym = sp.symbols('phi theta')    # Define the symbolic variables for phi and theta

norm = normalize(l)
#norm_control(l)

#g = sp.cos(phi_sym)
g = theta_sym-sp.pi
#g = 1

v = approx_func(l, g)
#print(v)

phi_values = np.linspace(0, np.pi, int(N_vykres/2))  # Adjust the number of points as needed
theta_values = np.linspace(0, 2*np.pi, int(N_vykres))  # Adjust the number of points as needed
phi_mesh, theta_mesh = np.meshgrid(phi_values, theta_values)

x = np.cos(phi_mesh)
y = np.sin(phi_mesh)*np.cos(theta_mesh)
z = np.sin(phi_mesh)*np.sin(theta_mesh)

v_lambda = sp.lambdify((phi_sym, theta_sym), v, "numpy")
g_lambda = sp.lambdify((phi_sym, theta_sym), g, "numpy")

v_values = v_lambda(phi_mesh, theta_mesh)
g_values = g_lambda(phi_mesh, theta_mesh)

v_real = np.real(v_values)
g_real = np.real(g_values)

plot(v_real, phi_sym, theta_sym, x, y, z)
plot(g_real, phi_sym, theta_sym, x, y, z)

alp(0,0,sp.cos(phi_sym))