# This file is part of gILC.
#
# gILC - Generic Iterative Learning Control for Nonlinear Systems
# Copyright (C) 2012 Marnix Volckaert, KU Leuven
#
# gILC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gILC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gILC. If not, see <http://www.gnu.org/licenses/>.
#

from plant import *
import random
from casadi import *

class invpen(plant):

    def __init__(self):
        self.nx = 2
        self.ny = 1
        self.noisestd = 0.0

    # Define System Dynamics
    #=======================
    def f(self,xk,uk):
    #---------------
        I = 0.006
        m = 0.5
        L = 0.3
        c = 0.1
        Ts = 1/500.0
        #------
        xkp1 = vertcat(xk[0] + Ts*xk[1],\
                xk[1] + (Ts/(I+m*L**2))*(-c*xk[1] - 9.81*m*L*sin(xk[0]) + uk[0]))
        return xkp1
    #---------------
    def h(self,xk,uk):
        #---------------
        y = xk[0]
        y += random.gauss(0.0,self.noisestd)
        return y
    #=======================
