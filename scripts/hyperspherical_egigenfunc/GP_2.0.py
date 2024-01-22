import numpy as np
from scipy.integrate import tplquad
from scipy.special import gegenbauer
from scipy.special import lpmv
from scipy.integrate import simpson

#######################################################################################

def agp(n,l,alpha,x):
    polynomial = gegenbauer(n,alpha)
    return((np.poly1d([-1,0,1])(x))**(l/2)*np.polyder(polynomial, l)(x))

def eigen_func(zeta,theta,phi,n,l,m):
    result = agp(n,l,1,np.cos(phi))*lpmv(m,l,np.cos(theta))*np.exp(-1j*m*zeta)
    return(result)

def square_eigen_jacobian(zeta,theta,phi,n,l,m):
    return(eigen_func(phi,theta,zeta,n,l,m)*np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*np.sin(phi)**2*np.sin(theta))

def norm(n,l,m):
    result, error = tplquad(square_eigen_jacobian, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
    return(np.sqrt(result))

def integrand_components_real(zeta,theta,phi,n,l,m):
    result = np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*func_g(zeta,theta,phi)*np.sin(phi)**2*np.sin(theta)
    return(np.real(result))

def integrand_components_imag(zeta,theta,phi,n,l,m):
    result = np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*func_g(zeta,theta,phi)*np.sin(phi)**2*np.sin(theta)
    return(np.imag(result))

def components(n,l,m):
    a = norm(n,l,m)
    result_real, error_real = tplquad(integrand_components_real, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
    result_imag, error_imag = tplquad(integrand_components_imag, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
    return((result_real+result_imag)/a)

def norm_simps(n,l,m):
    integrand = eigen_func(zeta,theta,phi,n,l,m)*np.conjugate(eigen_func(zeta,theta,phi,n,l,m))*np.sin(phi)**2*np.sin(theta)
    result = simpson(simpson(simpson(integrand,zeta_values, axis = 2), theta_values, axis = 1),phi_values, axis = 0)
    return(np.sqrt(result))

def components_simps(n,l,m):
    function_values = np.conjugate(eigen_func(zeta,theta,phi,n,l,m))*func_g(zeta,theta,phi)
    jacobian = (np.sin(phi))**2*np.sin(theta)
    result = simpson(simpson(simpson(function_values*jacobian,zeta_values, axis = 2), theta_values, axis = 1),phi_values, axis = 0)
    return(result)

#######################################################################################

def func_g(zeta,theta,phi):
    return(np.cos(phi-np.pi/2)*np.cos(theta-np.pi/2))

N_integration = 100

phi_values = np.linspace(0,np.pi,int(N_integration/2))
theta_values = np.linspace(0,np.pi,int(N_integration/2))
zeta_values = np.linspace(-np.pi,np.pi,int(N_integration))
phi,theta,zeta = np.meshgrid(phi_values, theta_values, zeta_values)

N = 6
g_val = 0
i = 0

# print(components_simps(0,0,0))

for n in range(N+1):
    for l in range(n+1):
        for m in range(-l,l+1):
            comp = components_simps(n,l,m)
            g_val += comp*eigen_func(0,np.pi/2,np.pi/2,n,l,m)/norm_simps(n,l,m)**2
            i += 1
            print("done", i)
            print("Comp",n,l,m,"=",comp)
            # print(gegenbauer(n,1))
            # print(np.polyder(gegenbauer(n,1), l))
            # print(n,l,m,"-",eigen_func(0,0,0,n,l,m))

print(func_g(0,np.pi/2,np.pi/2), g_val)