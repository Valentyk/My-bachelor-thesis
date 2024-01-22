import sympy as sp
from sympy import solve
from sympy.abc import x,y,z,phi,theta,zeta

#result = solve([(sp.cos(phi))/(1-sp.sin(phi)*sp.sin(theta)*sp.sin(zeta)), ((sp.sin(phi)*sp.cos(theta)))/(1-sp.sin(phi)*sp.sin(theta)*sp.sin(zeta)),(sp.sin(phi)*sp.sin(theta)*sp.cos(zeta))/(1-sp.sin(phi)*sp.sin(theta)*sp.sin(zeta))],[x,y,z], set = True)
result = solve([(2*x**2)/(x**2+y**2+z**2+1), (2*y**2)/(x**2+y**2+z**2+1), (2*z**2)/(x**2+y**2+z**2+1)], [sp.cos(phi), sp.sin(phi)*sp.cos(theta), sp.sin(phi)*sp.sin(theta)*sp.cos(zeta)], set = True)
print(result)