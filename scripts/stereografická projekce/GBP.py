import scipy as sp
import numpy as np
from sympy.functions.special.polynomials import gegenbauer as gbp
from sympy import assoc_legendre as alp
from sympy.abc import x, n, l

print(gbp(1,1,x))