import plotly.graph_objects as go
import numpy as np
import numpy as np
from scipy.special import gegenbauer
from scipy.special import lpmv
from scipy.integrate import simpson
# from scipy.constants import epsilon_0
# from scipy.integrate import tplquad
# from plotly.subplots import make_subplots

#######################################################################################
# FUNCTIONS USED IN THE CODE:

def agp(n,l,alpha,x):                                                           # Defining the asssociated Gegenbauer polynomials   
    polynomial = gegenbauer(n,alpha)
    return((np.poly1d([-1,0,1])(x))**(l/2)*np.polyder(polynomial, l)(x))

def eigen_func(zeta,theta,phi,n,l,m):                                           # Defining the hyperspherical harmonics (eigenfunciton of the Laplace operator on S^3)
    result = agp(n,l,1,np.cos(phi))*lpmv(m,l,np.cos(theta))*np.exp(-1j*m*zeta)
    return(result)

def components_simps(n,l,m):                                                    # Calculation of the Fourier coefficient for the input function
    function_values = np.conjugate(eigen_func(zeta,theta,phi,n,l,m))*func_g(zeta,theta,phi)
    jacobian = (np.sin(phi))**2*np.sin(theta)
    result = simpson(simpson(simpson(function_values*jacobian,zeta_values, axis = 2), theta_values, axis = 1),phi_values, axis = 0)
    return(result)

def R3_R4(x_values, y_values, z_values):                                        # Inverse stereographic projection
    X = (2*x_values)/(x_values**2+y_values**2+z_values**2+1) 
    Y = (2*y_values)/(x_values**2+y_values**2+z_values**2+1) 
    Z = (2*z_values)/(x_values**2+y_values**2+z_values**2+1) 
    W = (x_values**2+y_values**2+z_values**2-1)/(x_values**2+y_values**2+z_values**2+1) 
    return(X,Y,Z,W)

def R4_S3(X,Y,Z,W):                                                             # Inverse hyperspherical map
    zeta = np.arctan2(W,Z)
    theta = np.arctan2(np.sqrt(Z**2+W**2),Y)
    phi = np.arctan2(np.sqrt(W**2+Z**2+Y**2),X)
    return(phi, theta, zeta)

def norm_analytical(n,l,m):                                                      # Nomralization constant of the hyperspherical harmonic
    return(np.sqrt(((2*l+1)*(n+1)*np.math.factorial(l-m)*np.math.factorial(n-l))/(2*np.pi**2*np.math.factorial(l+m)*np.math.factorial(n+l+1))))

#######################################################################################
#OLD FUNCTIONS GRAVEYARD:

# def square_eigen_jacobian(zeta,theta,phi,n,l,m):                                
    # return(eigen_func(phi,theta,zeta,n,l,m)*np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*np.sin(phi)**2*np.sin(theta))

# def norm(n,l,m):
#     result, error = tplquad(square_eigen_jacobian, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
#     return(np.sqrt(result))

# def integrand_components_real(zeta,theta,phi,n,l,m):
#     result = np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*func_g(zeta,theta,phi)*np.sin(phi)**2*np.sin(theta)
#     return(np.real(result))

# def integrand_components_imag(zeta,theta,phi,n,l,m):
#     result = np.conjugate(eigen_func(phi,theta,zeta,n,l,m))*func_g(zeta,theta,phi)*np.sin(phi)**2*np.sin(theta)
#     return(np.imag(result))

# def components(n,l,m):
#     a = norm(n,l,m)
#     result_real, error_real = tplquad(integrand_components_real, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
#     result_imag, error_imag = tplquad(integrand_components_imag, 0,np.pi, 0,np.pi, 0,2*np.pi, args=(n,l,m))
#     return((result_real+result_imag)/a)

# def norm_simps(n,l,m):
#     integrand = eigen_func(zeta,theta,phi,n,l,m)*np.conjugate(eigen_func(zeta,theta,phi,n,l,m))*np.sin(phi)**2*np.sin(theta)
#     result = simpson(simpson(simpson(integrand,zeta_values, axis = 2), theta_values, axis = 1),phi_values, axis = 0)
#     return(np.sqrt(result))

#######################################################################################
# START OF THE CODE

stereographic_projection = True                                                 # Turn on/off the stereographic projection of the func_g (see below!)
approx_fourier_series = False                                                    # Turn on/off the approximation of the function func_g by the Fourier series and stereographic projection of the approximatiion
solve_poisson = False                                                            # Turn on/off the solving of the Poisson equation (func_g serves as the charge density function) by the Fourier series and eigenvalues (more on this in the Thesis) 

n,l,m = 2,2,-2                                                                # Defiing coefficients of hyperspherical harmonics 

def func_g(zeta,theta,phi):                                                     # Function for stereographic projection or solving the poisson equation  
    return(norm_analytical(n,l,m)*eigen_func(zeta,theta,phi,n,l,m))           # Option for the hyperspherical harmonic
    # return(np.cos(phi)*np.cos(theta)*np.cos(zeta))
    # return((np.sin(phi))*np.sin(theta)*(np.cos(zeta))**2)
    # return((np.cos(phi)+2*np.cos(theta))*np.sin(2*phi)*np.sin(zeta))

N_integration = 100                                                             # Number of steps for integration

phi_values = np.linspace(0,np.pi,int(N_integration/2))                          # Definition of the coordinate meshgrid for integration and plotting
theta_values = np.linspace(0,np.pi,int(N_integration/2))
zeta_values = np.linspace(-np.pi,np.pi,int(N_integration))
phi,theta,zeta = np.meshgrid(phi_values, theta_values, zeta_values)

N = 7                                                                           # Limit of the number of terms (N - cap on the coefficient n). We used N = 7 in the thesis
i = 0                                                                           # Terms counter variable

g_val = 0                                                                       # Defiing the aproximation function - Fourier series  
l_val = 0                                                                       # Defiing the solution of Poisson equaiton                                                                         
N_fig = 120                                                                     # Grid number definig the "resolution" of the Figures

x_values = np.linspace(-2,2,N_fig)                                              # Defiing the grid for the plotting
y_values = np.linspace(-2,2,N_fig)
z_values = np.linspace(-2,2,N_fig)
x,y,z = np.meshgrid(x_values, y_values, z_values)                               

X,Y,Z,W = R3_R4(x,y,z)                                                          # Inverse stereographic projection
phi_g, theta_g, zeta_g = R4_S3(X,Y,Z,W)                                         # Inverse hyperspherical map - gives each point in the plot its hyperpherical coordinate phi, theta and zeta


if approx_fourier_series == True or solve_poisson == True:                      # Calculating the Fourier coefficients if it is needed by our setup of the script

    for n in range(N+1):                                                        # Here we used the variable N for limiting the number of calculated Fourier coefficients
        for l in range(n+1):
            for m in range(-l,l+1):
                comp = components_simps(n,l,m)

                if n == 0:                                                      # Information about the mean value of the function on the hypersphere - we need it to be 0 if we want to solve the Poisson equation
                    print("Mean value of the function times the volume of the hypersphere:", comp)                                          
                    if comp > 1e-4:
                        print("The mean value of the function is not zero")
                        # quit()                                                # We can end the sript here

                norma = norm_analytical(n,l,m)                                  # Not exactly norm but normalization constant :)
                a = comp*eigen_func(zeta_g,theta_g,phi_g,n,l,m)*norma**2        # Fourier coefficient time the corresponding hyperspherical harmonics
                
                if approx_fourier_series == True:
                    g_val += a                                                  # Adding the term to the approximation

                if solve_poisson == True:    
                    if n != 0:                                                  # Adding the term to the solution of the Poisson equation
                        l_val += a/(n*(n+2))                                       

                i += 1
                print("n:",n,"l:",l,"m:",m,"number of terms:",i,"   ",end="\r") # Information about the progress of calculating the Fourier coefficients


name = ("")                                                                     # If we do not use the hyperspherical harmonics the naming of the output files is general
# name = (str(n)+","+str(l)+","+str(m))                                         # If the hyperspherical harmonic is used as an function func_g we can name the output files by its coefficients

camera = dict(                                                                  # Position of the "camera" for the final plots
    eye=dict(x=0., y=0.0, z=2.0)
    # eye=dict(x=1.25, y=1.25, z=1.25)
)

layout = go.Layout(                                                             # Modification of the image margins
    margin=go.layout.Margin(
        l=0,                                                                    # Left margin
        r=0,                                                                    # Right margin
        b=5,                                                                    # Bottom margin
        t=0,                                                                    # Top margin
        )
    )


if stereographic_projection == True:                                            # Block of the stereographic projection of the func_g

    func_values_g = np.real(func_g(zeta_g,theta_g,phi_g))                       # Evaluation of the func_g for each point of the plot grid 

    fig = go.Figure(layout=layout, data = go.Volume(                            # Making the plot
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value = func_values_g.flatten(),
        isomin = np.min(func_values_g),
        isomax = np.max(func_values_g),
        # isomin = -1,
        # isomax = 1,
        surface_count = 10,
        opacity = 0.5,
        #surface_fill = 0.9,
        colorscale = "plasma",
        caps = dict(x_show = False, y_show = False, z_show = False)
    ))

    fig.update_layout(scene_camera=camera)

    fig.show()                                                                                                                                # Open the interactive plot
    # fig.write_html(f"/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/inter_orig_{name}.html")                   # Save the interactive plot (THE SIZE OF THE PLOT FILE MAY EXCEED 200 MB!) 
    # fig.write_image(f"/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/fig_orig_{name}.png", scale = 5)          # Save an image of the plot


if approx_fourier_series == True:                                               # Block of the approximation of the func_g and its stereographic projection 
    
    func_values = np.real(g_val)                                                # Giving each point of the plot grid its approximation value

    fig = go.Figure(layout=layout, data = go.Volume(                            # Making the plot
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value = func_values.flatten(),
        isomin = np.min(func_values),                                        
        isomax = np.max(func_values),                                         
        # isomin = -1,
        # isomax = 1, 
        surface_count = 10,
        opacity = 0.5,
        #surface_fill = 0.9,
        colorscale = "plasma",
        caps = dict(x_show = False, y_show = False, z_show = False)
    ))

    fig.update_layout(scene_camera=camera)

    fig.show()                                                                                                                              # Open the interactive plot
    # fig.write_html("/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/inter_aprox.html")                        # Save the interactive plot (THE SIZE OF THE PLOT FILE MAY EXCEED 200 MB!)
    fig.write_image("/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/fig_aprox.png", scale = 5)               # Save an image of the plot

if solve_poisson == True:                                                      # Block for solving the Poisson equations (func_g is the charge density function) and stereographic projection of the solution 

    func_values_l = np.real(l_val)                                             # Giving each point of the plot grid its solution value

    fig = go.Figure(layout=layout, data = go.Volume(                           # Making the plot
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value = func_values_l.flatten(),
        isomin = np.min(func_values_l),
        isomax = np.max(func_values_l),
        # isomin = -1,
        # isomax = 1, 
        surface_count = 10,
        opacity = 0.5,
        #surface_fill = 0.9,
        colorscale = "plasma",
        caps = dict(x_show = False, y_show = False, z_show = False)
    ))

    fig.update_layout(scene_camera=camera)

    fig.show()                                                                                                                              # Open the interactive plot
    # fig.write_html("/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/inter_sol.html")                          # Save the interactive plot (THE SIZE OF THE PLOT FILE MAY EXCEED 200 MB!)
    fig.write_image("/home/michalvalentik/OneDrive/MUNI/Bakalářka/vizualizace_v_praci/zatím_mimo/fig_sol.png", scale = 5)                 # Save an image of the plot   
