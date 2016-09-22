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

from numpy import *
from casadi import *

class plant:

        def simout(self,xinit,u):
            N = u.shape[1]
            x = DM.zeros(self.nx, N)
            y = DM.zeros(self.ny, N)
            x[:,0] = xinit

            for k in range(N-1):
                x[:,k+1] = self.f(x[:,k],u[:,k])
                y[:,k] = self.h(x[:,k],u[:,k])
            y[:,N-1] = self.h(x[:,N-1],u[:,k])

            return y.full()

        def simstates(self,xinit,u):
            N = u.shape[1]
            x = DM.zeros(self.nx, N)
            x[:,0] = xinit

            for k in range(N-1):
                x[:,k+1] = self.f(x[:,k],u[:,k])

            return x.full()

        def sim(self,xinit,u):
            y = self.simout(xinit,u)
            states = self.simstates(xinit,u)

            return (states,y)

        def simssout(self,u,periods):
            ut = tile(u,(1,periods))    # Repeat the input for 'periods' times

            N = u.shape[1]
            Nt = periods*N
            xt = zeros((len(self.f(zeros((self.nx)),u[:,0])),Nt))
            yt = zeros((len(self.h(zeros((self.nx)),u[:,0])),Nt))
            y = zeros((len(self.h(zeros((self.nx)),u[:,0])),N))

            for k in range(Nt-1):
                xt[:,k+1] = self.f(xt[:,k],ut[:,k])
                yt[:,k] = self.h(xt[:,k],ut[:,k])
            yt[:,Nt-1] = self.h(xt[:,Nt-1],ut[:,k])

            for k in range(N):
                y[:,k] = yt[:,Nt-N+k]

            return y
