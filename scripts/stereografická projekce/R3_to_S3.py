import plotly.graph_objects as go
import numpy as np

N = 60

###########################################################################

x_values = np.linspace(-5,5,N)
y_values = np.linspace(-5,5,N)
z_values = np.linspace(-5,5,N)

x, y, z = np.meshgrid(x_values, y_values, z_values)

phi = np.arccos((2*x)/(x**2+y**2+z**2+1))
theta = np.arccos((2*y)/((x**2+y**2+z**2+1)*np.sin(phi)))
zeta = np.arccos((2*z)/((x**2+y**2+z**2+1)*np.sin(phi)*np.sin(theta)))
helper = (x**2+y**2+z**2-1)/((x**2+y**2+z**2+1)*np.sin(phi)*np.sin(theta))

grid_shape = phi.shape

for i in range(grid_shape[0]):
    for j in range(grid_shape[1]):
        for k in range(grid_shape[2]):
            if helper[i,j,k] < 0:
                zeta[i,j,k] = 2*np.pi - zeta[i,j,k]


# func_values = phi
func_values = theta**2-np.exp(zeta)
# func_values = np.sin(zeta)
# func_values = zeta


fig = go.Figure(data = go.Volume(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value = func_values.flatten(),

    isomin = np.min(func_values),
    isomax = np.max(func_values),
    surface_count = 7,
    opacity = 0.5,
    # surface_fill = 0.5,
    caps = dict(x_show = False, y_show = False, z_show = False)
))

fig.show()