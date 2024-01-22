import plotly.graph_objects as go
import numpy as np
import plotly.express as px

###########################################################################

N = 100

x_values = np.linspace(-2,2,N)
y_values = np.linspace(-2,2,N)
z_values = np.linspace(-2,2,N)

x,y,z = np.meshgrid(x_values, y_values, z_values)

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

X,Y,Z,W = R3_R4(x,y,z)
phi, theta, zeta = R4_S3(X,Y,Z,W)

func_values = zeta #np.cos(zeta)

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


print("W:", np.min(W), np.max(W), "Z:", np.min(Z), np.max(Z))