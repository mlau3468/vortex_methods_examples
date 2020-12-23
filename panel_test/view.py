import pymesh
from dust.tools.DUST_utils import *

from numpy import genfromtxt
ee = genfromtxt('ee.csv', delimiter=',')
print(ee)
rr = genfromtxt('rr.csv', delimiter=',')
print(rr)
view_ee_rr(ee,rr)