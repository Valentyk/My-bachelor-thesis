import plotly.graph_objects as go
import numpy as np
import plotly.express as px

#plot the 3D function with isosurfaces

# N_plot = 50

# x = np.linspace(0,1,N_plot)
# y = np.linspace(0,1,N_plot)
# z = np.linspace(0,1,N_plot)

# xv, yv, zv, = np.meshgrid(x,y,z)
# value = xv**2+yv**2+zv**2

# fig = go.Figure(data = go.Isosurface(
#     x=xv.flatten(),
#     y=yv.flatten(),
#     z=zv.flatten(),
#     value = value.flatten(),

#     isomin = 0,
#     isomax = 1,
#     surface_count = 5,
#     opacity = 0.5,
#     surface_fill = 0.5,
#     caps = dict(x_show = False, y_show = False, z_show = False)
# ))

# fig.show()

#Define function of phi, theta and zeta

N = 100

###########################################################################

x_values = np.linspace(-2,2,N)
y_values = np.linspace(-2,2,N)
z_values = np.linspace(-2,2,N)

x, y, z = np.meshgrid(x_values, y_values, z_values)

phi = np.arccos((2*x)/(x**2+y**2+z**2+1))
theta = np.arccos((2*y)/((x**2+y**2+z**2+1)*np.sin(phi)))
zeta = np.arccos((2*z)/((x**2+y**2+z**2+1)*np.sin(phi)*np.sin(theta)))

###########################################################################

# phi_values = np.linspace(0,np.pi,int(N/2))
# theta_values = np.linspace(0,np.pi,int(N/2))
# zeta_values = np.linspace(0,2*np.pi,int(N))

# phi, theta, zeta = np.meshgrid(phi_values, theta_values, zeta_values)

###########################################################################

# STEREOGRAPHICAL PROJECTION

# xv = ((np.cos(phi))/(1-np.sin(phi)*np.sin(theta)*np.sin(zeta)))
# yv = ((np.sin(phi)*np.cos(theta))/(1-np.sin(phi)*np.sin(theta)*np.sin(zeta)))
# zv = ((np.sin(phi)*np.sin(theta)*np.cos(zeta))/(1-np.sin(phi)*np.sin(theta)*np.sin(zeta)))

func_values = zeta

fig = go.Figure(data = go.Volume(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value = func_values.flatten(),

    isomin = np.min(func_values),
    isomax = np.max(func_values),
    surface_count = 7,
    opacity = 0.5,
   #  surface_fill = 0.5,
    caps = dict(x_show = False, y_show = False, z_show = False)
))

# trace = go.Scatter3d(
#    x = xv.flatten(), y = yv.flatten(), z = zv.flatten(), mode = 'markers', marker = dict(
#       size = 1,
#       color = zeta.flatten(), # set color to an array/list of desired values
#       )
#    )
# layout = go.Layout(title = '3D Scatter plot')
# fig = go.Figure(data = [trace], layout = layout)

fig.show()