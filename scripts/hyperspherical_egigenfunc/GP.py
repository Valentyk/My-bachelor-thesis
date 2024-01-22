import sympy as sp
import numpy as np
from  sympy.functions.special.polynomials import gegenbauer as gp
from sympy import assoc_legendre as alp
from sympy.polys.polytools import pdiv
from sympy.abc import x, nu, l, m
from sympy import assoc_legendre as alp
from sympy.functions.special.polynomials import legendre as lp
from  scipy.integrate import tplquad

#######################################################################################

def agp(n,l,nu):
    poly = gp(n,nu,x)
    result = (-1)**l*sp.sqrt((1-x**2)**l)*sp.diff(poly,x,l)
    return(result.subs(x,sp.cos(phi_sym)))

def eigen_func(n,l,m):
    eigen_func_unnorm = agp(n,l,1)*alp(l,m,sp.cos(theta_sym))*sp.exp(sp.I*m*zeta_sym)
    return(eigen_func_unnorm/norm(eigen_func_unnorm, phi, theta, zeta))

def norm(func, phi, theta, zeta):

    func_sqrt = ((sp.conjugate(func))*func)
    jacobian = (sp.sin(phi_sym))**2*sp.sin(theta_sym)

    function_real = sp.lambdify((phi_sym, theta_sym, zeta_sym), (func_sqrt*jacobian).as_real_imag()[0], "numpy")
    function_imag = sp.lambdify((phi_sym, theta_sym, zeta_sym), (func_sqrt*jacobian).as_real_imag()[1], "numpy")

    result_real, error_real = tplquad(function_real, np.min(zeta), np.max(zeta), np.min(theta), np.max(theta), np.min(phi), np.max(phi))
    result_complex, error_complex = tplquad(function_imag, np.min(zeta), np.max(zeta), np.min(theta), np.max(theta), np.min(phi), np.max(phi))
    result = result_real+result_complex*1j

    return(result**(1/2))

def projection_comp(eigen_func, projected_function):
    projection_func = eigen_func*projected_function

    function_real = sp.lambdify((phi_sym, theta_sym, zeta_sym), projection_func.as_real_imag()[0], "numpy")
    function_imag = sp.lambdify((phi_sym, theta_sym, zeta_sym), projection_func.as_real_imag()[1], "numpy")

    result_real, error_real = tplquad(function_real, np.min(zeta), np.max(zeta), np.min(theta), np.max(theta), np.min(phi), np.max(phi))
    result_complex, error_complex = tplquad(function_imag, np.min(zeta), np.max(zeta), np.min(theta), np.max(theta), np.min(phi), np.max(phi))
    result = result_real+result_complex*1j

    return(result)

#######################################################################################

phi_sym, theta_sym, zeta_sym = sp.symbols('phi theta zeta')    # Define the symbolic variables for phi and theta
N = 60
phi = np.linspace(0,np.pi,N)
theta = np.linspace(0,np.pi,N)
zeta = np.linspace(0,2*np.pi,2*N)

g_func = phi_sym-sp.pi/2
v_func = 0

max_n = 3

for n in range(max_n+1):
    for l in range(n+1):
        for m in range(-l,l+1):
            v_func += projection_comp(eigen_func(n,l,m), g_func)
        print("l = ", l, "done")

g_lambda = sp.lambdify((phi_sym, theta_sym, zeta_sym), g_func, "numpy")
v_lambda = sp.lambdify((phi_sym, theta_sym, zeta_sym), v_func, "numpy")

print(g_lambda(np.pi/2,np.pi/2,np.pi/2)-v_lambda(np.pi/2,np.pi/2,np.pi/2))






















# phi_mesh, theta_mesh, zeta_mesh = np.meshgrid(phi, theta, zeta)

# print(norm(eigen_func(5,1,1), phi, theta, zeta))

# print(eigen_func(2,1,1)(0,0,0))

#print(norm(phi_sym+1j*theta_sym, phi, theta, zeta))

# test = sp.lambdify(x,agp(3,3,1/2,x))

# for i in range(5):
#     for j in range(0,i+1):
#         test = sp.lambdify(x,agp(i,j,1/2,x))
#         if test(0.2) - alp(i,j,0.2) == 0:
#             a = 2
#             #print(agp(i,j,0.5,x), "kontrola: ", alp(i,j,x), "l =", i,"m =", j, "\n")
#         else:
#             print(agp(i,j,0.5,x), "kontrola: ", alp(i,j,x), "l =", i,"m =", j, "\n")
#             print(test(0.2), alp(i,j,0.2))

# print(agp(3,-2,1/2,x))
# print(alp(3,-2,x))
#print(lp(2,x))

