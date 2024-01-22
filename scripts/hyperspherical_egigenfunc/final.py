import plotly.graph_objects as go
import numpy as np
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

def R3_R4(x_values, y_values, z_values):
    X = (2*x_values)/(x_values**2+y_values**2+z_values**2+1) 
    Y = (2*y_values)/(x_values**2+y_values**2+z_values**2+1) 
    Z = (2*z_values)/(x_values**2+y_values**2+z_values**2+1) 
    W = (x_values**2+y_values**2+z_values**2-1)/(x_values**2+y_values**2+z_values**2+1) 
    return(X,Y,Z,W)

def R4_S3(X,Y,Z,W):
    zeta = np.arctan2(W,Z)
    theta = np.arctan2(np.sqrt(Z**2+W**2),Y)
    phi = np.arctan2(np.sqrt(W**2+Z**2+Y**2),X)
    return(phi, theta, zeta)

#######################################################################################

def func_g(zeta,theta,phi):
    return(np.sin(phi)*np.exp(theta))

N_integration = 100

phi_values = np.linspace(0,np.pi,int(N_integration/2))
theta_values = np.linspace(0,np.pi,int(N_integration/2))
zeta_values = np.linspace(-np.pi,np.pi,int(N_integration))
phi,theta,zeta = np.meshgrid(phi_values, theta_values, zeta_values)

N = 7
comp = []
norma = []
i = 0 

for n in range(N+1):
    for l in range(n+1):
        for m in range(-l,l+1):
            comp.append(components_simps(n,l,m))
            norma.append(norm_simps(n,l,m))

g_val = 0
i = 0

N_fig = 100

x_values = np.linspace(-2,2,N_fig)
y_values = np.linspace(-2,2,N_fig)
z_values = np.linspace(-2,2,N_fig)

x,y,z = np.meshgrid(x_values, y_values, z_values)

X,Y,Z,W = R3_R4(x,y,z)
phi, theta, zeta = R4_S3(X,Y,Z,W)


for n in range(N+1):
    for l in range(n+1):
        for m in range(-l,l+1):
            g_val += comp[i]*eigen_func(zeta,theta,phi,n,l,m)/norma[i]**2
            i += 1

func_values = func_g(zeta,theta,phi)

fig = go.Figure(data = go.Volume(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value = func_values.flatten(),

    isomin = np.min(func_values),
    isomax = np.max(func_values),
    surface_count = 7,
    opacity = 0.5,
    #surface_fill = 0.9,
    colorscale = "plasma",
    caps = dict(x_show = False, y_show = False, z_show = False)
))

fig.show()

func_values = np.real(g_val)

fig = go.Figure(data = go.Volume(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value = func_values.flatten(),

    isomin = np.min(func_values),
    isomax = np.max(func_values),
    surface_count = 7,
    opacity = 0.5,
    #surface_fill = 0.9,
    colorscale = "plasma",
    caps = dict(x_show = False, y_show = False, z_show = False)
))

fig.show()