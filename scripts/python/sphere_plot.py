import plotly.graph_objects as go
import numpy as np
from sympy import assoc_legendre as alp
import sympy as sp

def plot(function_values, phi_sym, theta_sym, x, y, z):
    # Create a 3D surface plot with surfacecolor
    fig = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=function_values, colorscale='Viridis')])

    # Set axis labels
    fig.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'))

    # Set the title of the plot
    fig.update_layout(title=r'Function $f = f(\varphi, \vartheta)$y on a Sphere')

    # Show the plot
    fig.show()

if __name__ == "__main__":

    phi_sym, theta_sym = sp.symbols('phi theta')
    function_expr = theta_sym
    N = 50

    print(function_expr)   

    plot(function_expr, phi_sym, theta_sym, N)