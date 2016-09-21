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

"""gILC - an open source tool for model based iterative learning control

gilcClass comprises wrappers for the estimation and control step. Nlpproblem
defines a generic frame for a nonlinear problem to be solved by utilizing
optimization.
"""

from numpy import *
from casadi import *

class steps:
    pass

class gilcClass:

    def __init__(self):
        self.est = steps()
        self.con = steps()

        # Default: no minimized and constrained outputs
        #----------------------------------------------
        self.est.hzc = None
        self.est.hzm = None
        self.est.Qy = None
        self.est.xguess = None
        self.con.hzm = None
        self.con.hzc = None
        self.con.Qy = None
        self.con.xguess = None
        self.end = 'dummy'
        #----------------------------------------------

    def estimation(self,u_i,y_i):

        eststep = nlpproblem(1,self.nu,self.na,self.nx,self.nxa,self.ny)
        eststep.setSystemdynamics(self.f,self.h)
        if self.est.hzm != None:
            eststep.setMinimizedoutputs(self.est.hzm,self.est.rz,self.est.Qz)
        if self.est.hzc != None:
            eststep.setConstrainedoutputs(self.est.hzc,self.est.zmin,self.est.zmax)
        if self.est.Qy != None:
            eststep.setWeights(self.est.Qy)
        if self.est.xguess != None:
            eststep.setInitialguess(self.est.xguess,self.est.alphaguess)
        if self.end == "periodic":
            eststep.setReference(y_i,self.end)
        else:
            eststep.setInitstate(self.xinit)
            eststep.setReference(y_i,'dummy')
        eststep.setExogin(u_i)
        self.alpha = eststep.solve()
        return (self.alpha,eststep.xopt)

    def control(self,alpha):

        constep = nlpproblem(2,self.nu,self.na,self.nx,self.nxa,self.ny)
        constep.setSystemdynamics(self.f,self.h)
        if self.con.hzm != None:
            constep.setMinimizedoutputs(self.con.hzm,self.con.rz,self.con.Qz)
        if self.con.hzc != None:
            constep.setConstrainedoutputs(self.con.hzc,self.con.zmin,self.con.zmax)
        if self.con.Qy != None:
            constep.setWeights(self.con.Qy)
        if self.con.xguess != None:
            constep.setInitialguess(self.con.xguess,self.con.uguess)
        if self.end == 'periodic':
            constep.setReference(self.r,self.end)
        else:
            constep.setInitstate(self.xinit)
            constep.setReference(self.r,'dummy')
        constep.setExogin(alpha)
        self.u_ip1 = constep.solve()

        return (self.u_ip1,constep.xopt)

    def update(self,u_i,y_i):

        alpha,xalpha = self.estimation(u_i,y_i)
        u,xu = self.control(alpha)
        return (alpha,u)


class nlpproblem:

    def __init__(self,step,nu,na,nx,nxa,ny):

        self.vinitProvided = False
        self.QyProvided = False
        self.Constrainedoutputs = False
        self.Minimizedoutputs = False

        self.step = step
        self.nx = nx
        self.nxa = nxa
        self.ny = ny

        if step == 1:
            self.nuc = na    # number of control inputs
            self.nw = nu     # number of exogeneous inputs
        elif step == 2:
            self.nuc = nu
            self.nw = na

        # Create the symbolic variables to be used in f,h,hy,hz
        #------------------------------------------------------
        self.xk = SX.sym("xk",self.nx)
        self.uck = SX.sym("uck",self.nuc)
        self.wk = SX.sym("wk",self.nw)
        #------------------------------------------------------

    def setSystemdynamics(self,f,h):

        # Initialise dynamics functions
        #------------------------------
        if self.step == 1:
            self.ffcn = Function('ffcn', [self.xk,self.uck,self.wk], [f(self.xk,self.wk,self.uck)])
            self.hfcn = Function('hfcn', [self.xk,self.uck,self.wk], [h(self.xk,self.wk,self.uck)])
        elif self.step == 2:
            self.ffcn = Function('ffcn', [self.xk,self.uck,self.wk], [f(self.xk,self.uck,self.wk)])
            self.hfcn = Function('hfcn', [self.xk,self.uck,self.wk], [h(self.xk,self.uck,self.wk)])
        #------------------------------

    def setReference(self,r,end):
        self.N = r.shape[1]
        self.r = concatenate(numpy.transpose(r))
        self.end = end

    def setInitstate(self,xinit):
        self.xinit = xinit

    def setExogin(self,w):
        self.w = w

    def setConstrainedoutputs(self,hzc,zmin,zmax):
        self.Constrainedoutputs = True
        self.nz = zmin.shape[0]
        self.zmin = zmin
        self.zmax = zmax

        # Initialise Constrained output function
        #---------------------------------------
        if self.step == 1:
            self.hzcfcn = Function('hzcfcn', [self.xk,self.uck,self.wk], hzc(self.xk,self.wk,self.uck))
        if self.step == 2:
            self.hzcfcn = Function('hzcfcn', [self.xk,self.uck,self.wk], hzc(self.xk,self.uck,self.wk))
        #---------------------------------------

    def setMinimizedoutputs(self,hzm,rz,Qz):
        self.Minimizedoutputs = True

        # Initialise Constrained output function
        #---------------------------------------
        if self.step == 1:
            self.hzmfcn = Function('hzmfcn', [self.xk,self.uck,self.wk], hzm(self.xk,self.wk,self.uck))
        if self.step == 2:
            self.hzmfcn = Function('hzmfcn', [self.xk,self.uck,self.wk], hzm(self.xk,self.uck,self.wk))
        #---------------------------------------
        # Set the reference
        #------------------
        # self.rz = concatenate(rz)
        self.rz = concatenate(rz.transpose())
        #------------------
        # Set weight
        #-----------
        self.Qz = diag(concatenate(Qz.transpose()))
        #-----------

    def setInitialguess(self,xguess,ucguess):
        self.vinit = concatenate(transpose(vstack((xguess,ucguess))))
        self.vinitProvided = True

    def setWeights(self,Qy):
        self.QyProvided = True
        self.Qy = diag(concatenate(transpose(Qy)))

    def solve(self):

        # Check if an initial guess provided
        #-----------------------------------
        if self.vinitProvided == False:
            self.vinit = zeros((self.nx+self.nuc)*self.N)
        #-----------------------------------

        # Check if residual weights provided
        #-----------------------------------
        if self.QyProvided == False:
            self.Qy = diag(numpy.concatenate(numpy.transpose(1e8*ones((self.ny,self.N)))))
        #-----------------------------------

        # Define optimization variables
        #------------------------------
        x = SX.sym("x",self.nx,self.N)
        uc = SX.sym("u",self.nuc,self.N)
        #------
        v = []
        for k in range(self.N):
            v = vertcat(v, x[:,k])
            v = vertcat(v, uc[:,k])

        #------------------------------

        # Define residual vector
        #-----------------------
        ybar = []
        for k in range(0,self.N):
            ybar = vertcat(ybar, self.hfcn(x[:,k],uc[:,k],self.w[:,k]))
        ybar = self.r - ybar

        #-----------------------
        # Objective function
        #-------------------
        J = mtimes(transpose(ybar),mtimes(self.Qy,ybar))    # Residual of main tracking error

        #------
        if self.Minimizedoutputs:
            zmbar = []
            for k in range(0,self.N):
                zmbar = vertcat(zmbar, self.hzmfcn(x[:,k],uc[:,k],self.w[:,k]))
            zmbar = self.rz - zmbar
            J += mtimes(zmbar.T,mtimes(self.Qz, zmbar))

        # Constraint function
        #--------------------
        g = []
        gmin = []
        gmax = []
        for k in range(0,self.N-1):
            xkp1 = self.ffcn(x[:,k],uc[:,k],self.w[:,k])    # evaluate the state k+1
            g = vertcat(g, x[:,k+1]-xkp1)                   # add the constraint
            gmin = vertcat(gmin, zeros(self.nx))
            gmax = vertcat(gmax, zeros(self.nx))
        xkp1 = self.ffcn(x[:,self.N-1],uc[:,self.N-1],self.w[:,self.N-1])    # State at N
        #--- Periodic signal or end in rest ---
        if self.end == 'periodic':
            g = vertcat(g, x[:,0] - xkp1)
            gmin = vertcat(gmin, zeros(self.nx))
            gmax = vertcat(gmax, zeros(self.nx))
        else:
            g = vertcat(g, x[0:self.nx-self.nxa,0]-self.xinit[0:self.nx-self.nxa])       # Initial state
            gmin = vertcat(gmin, zeros(self.nx-self.nxa))
            gmax = vertcat(gmax, zeros(self.nx-self.nxa))
        #--- Constrained outputs ---
        if self.Constrainedoutputs:
            for k in range(0,self.N):
                g = vertcat(g, self.hzcfcn(x[:,k],uc[:,k],self.w[:,k]))
                gmin = vertcat(gmin, self.zmin[:,k])
                gmax = vertcat(gmax, self.zmax[:,k])
        #--------------------

        # Solve problem using IPOPT
        #--------------------------
        prob = {'f': J, 'x': v, 'g': g}
        nlpopt = {'verbose': False, 'ipopt': {'tol': 1e-10}}
        solver = nlpsol('solver', 'ipopt', prob, nlpopt)
        arg = {'x0': self.vinit, 'lbg': gmin, 'ubg': gmax, 'lbx': -inf*ones(v.size()), 'ubx': inf*ones(v.size())}
        sol = solver.call(arg)

        #--------------------------
        # Get the solution
        #-----------------
        vopt = sol['x']
        self.xopt = zeros((self.nx,self.N))
        ucopt = zeros((self.nuc,self.N))
        vk = 0    # Counter for which time point
        for k in range(self.N):
            for i in range(self.nx):
                self.xopt[i,k] = vopt[vk+i]
            vk += self.nx
            for i in range(self.nuc):
                ucopt[i,k] = vopt[vk+i]
            vk += self.nuc
        #-----------------
        return ucopt

    def returnStates(self):
        return self.xopt
