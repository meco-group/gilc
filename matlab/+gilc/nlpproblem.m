% This file is part of gILC.
%
% gILC - Generic Iterative Learning Control for Nonlinear Systems
% Copyright (C) 2012 Marnix Volckaert, KU Leuven
%               2016 Armin Steinhauser, KU Leuven
%
% gILC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% gILC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with gILC. If not, see <http://www.gnu.org/licenses/>.
%

classdef nlpproblem
  %NLPPROBLEM Definition of a generic NLP problem.

  properties
    vinitProvided = false;
    QyProvided = false;
    Constrainedoutputs = false;
    Minimizedoutputs = false;
    step, nx, nxa, ny, nuc, nw, xk, uck, wk, ffcn, hfcn, xinit, N, r;
    mode_end, w, vinit, Qy, xopt, nz, zmin, zmax, hzcfcn, hzmfcn, rz, Qz;
  end

  methods
    function obj = nlpproblem(step, nu, na, nx, nxa, ny)
      import casadi.*

      obj.step = step;
      obj.nx = nx;
      obj.nxa = nxa;
      obj.ny = ny;
      if step == 1
        obj.nuc = na;
        obj.nw = nu;
      elseif step == 2
        obj.nuc = nu;
        obj.nw = na;
      end
      obj.xk = SX.sym('xk',obj.nx);
      obj.uck = SX.sym('uck',obj.nuc);
      obj.wk = SX.sym('wk',obj.nw);
    end

    function this = setSystemdynamics(this, f, h)
      import casadi.*

      if this.step == 1
        this.ffcn = Function('ffcn',{this.xk, this.uck, this.wk},{vertcat(f(this.xk, this.wk, this.uck))});
        this.hfcn = Function('hfcn',{this.xk, this.uck, this.wk},{vertcat(h(this.xk, this.wk, this.uck))});
      elseif this.step == 2
        this.ffcn = Function('ffcn',{this.xk, this.uck, this.wk},{vertcat(f(this.xk, this.uck, this.wk))});
        this.hfcn = Function('hfcn',{this.xk, this.uck, this.wk},{vertcat(h(this.xk, this.uck, this.wk))});
      end
    end

    function this = setReference(this, r, mode_end)
      this.N = numel(r);
      this.r = r;
      this.mode_end = mode_end;
    end

    function this = setInitstate(this, xinit)
      this.xinit = xinit;
    end

    function this = setExogin(this, w)
      this.w = w;
    end

    function this = setConstrainedoutputs(this, hzc, zmin, zmax)
      import casadi.*

      this.Constrainedoutputs = true;
      this.nz = size(zmin, 1);
      this.zmin = zmin;
      this.zmax = zmax;
      if this.step == 1
        this.hzcfcn = Function('hzcfcn',{this.xk,this.uck,this.wk},{vertcat(hzc(this.xk,this.wk,this.uck))});
      elseif this.step == 2
        this.hzcfcn = Function('hzcfcn',{this.xk,this.uck,this.wk},{vertcat(hzc(this.xk,this.uck,this.wk))});
      end
    end

    function this = setMinimizedoutputs(this, hzm, rz, Qz)
      import casadi.*

      this.Minimizedoutputs = true;
      % Initialise Constrained output function
      %---------------------------------------
      if this.step == 1
        this.hzmfcn = Function('hzmfcn',{this.xk,this.uck,this.wk},{vertcat(hzm(this.xk,this.wk,this.uck))});
      elseif this.step == 2
        this.hzmfcn = Function('hzmfcn',{this.xk,this.uck,this.wk},{vertcat(hzm(this.xk,this.uck,this.wk))});
      end
      %---------------------------------------
      % Set the reference
      %------------------
      this.rz = rz;
      %------------------
      % Set weight
      %-----------
      this.Qz = diag(Qz);
    end

    function this = setInitialguess(this, xguess, ucguess)
      import casadi.*

      this.vinit = vertcat(xguess, ucguess);
      this.vinitProvided = true;
    end

    function this = setWeights(this, Qy)
      this.QyProvided = true;
      this.Qy = diag(Qy);
    end

    function [this, ucopt] = solve(this)
      import casadi.*

      % Check if an initial guess provided
      %-----------------------------------
      if this.vinitProvided == false
        this.vinit = zeros((this.nx+this.nuc)*this.N, 1);
      end

      %-----------------------------------
      % Check if residual weights provided
      %-----------------------------------
      if this.QyProvided == false
        this.Qy = 1e8*ones(this.ny*this.N, 1);
      end

      %-----------------------------------
      % Define optimization variables
      %------------------------------
      x = SX.sym('x', this.nx, this.N);
      uc = SX.sym('u', this.nuc, this.N);

      % optimization variables
      v = vertcat(x, uc);
      v = v(:);

      %------------------------------
      % Define residual vector
      %-----------------------
      hfcn_map = this.hfcn.map('hfcn_map','serial',this.N);
      ybar = (this.r - hfcn_map(x,uc,this.w)).';

      %-----------------------
      % Objective function
      %-------------------
      J = (ybar.') * (this.Qy.*ybar);    % Residual of main tracking error
      %------
      if this.Minimizedoutputs
        hzmfcn_map = this.hzmfcn.map('hzmfcn_map','serial',this.N);
        zmbar = (this.rz - hfcn_map(x,uc,this.w)).';
        J = J + (zmbar.') * (this.Qz * zmbar);
      end

      % Constraint function
      %--------------------
      xkp1_map = this.ffcn.map('xkp1_map','serial',this.N);
      xkp1_all = xkp1_map(x,uc,this.w);
      xk_gaps = (x(:,2:this.N)-xkp1_all(:,1:this.N-1));
      g = xk_gaps(:);                          % add the constraint
      gmin = zeros((this.N-1)*this.nx, 1);
      gmax = zeros((this.N-1)*this.nx, 1);

      %--- Periodic signal or end in rest ---
      if strcmp(this.mode_end, 'periodic')
        g = vertcat(g, x(:,1) - xkp1_all(:,this.N));
        gmin = vertcat(gmin, zeros(this.nx, 1));
        gmax = vertcat(gmax, zeros(this.nx, 1));
      else
        g = vertcat(g, x(1:this.nx-this.nxa, 1)-this.xinit(1:this.nx-this.nxa));       % Initial state
        gmin = vertcat(gmin, zeros(this.nx-this.nxa, 1));
        gmax = vertcat(gmax, zeros(this.nx-this.nxa, 1));
      end

      %--- Constrained outputs ---
      if this.Constrainedoutputs
        hzcfcn_map = this.hzcfcn.map('hzcfcn_map','serial',this.N);
        g = vertcat(g, hzcfcn_map(x,uc,this.w).');
        gmin = vertcat(gmin, this.zmin.');
        gmax = vertcat(gmax, this.zmax.');
      end

      % Solve problem using IPOPT
      %--------------------------
      prob = struct('f', J, 'x', v, 'g', g);
      solver = nlpsol('solver', 'ipopt', prob);
      arg = struct('x0', this.vinit, 'lbg', gmin, 'ubg', gmax);
      sol = solver.call(arg);
      x_sol = full(sol.x);
      %--------------------------

      % Get the solution
      xopt = x_sol(logical(kron(ones(this.N,1),[ones(this.nx,1); zeros(this.nuc)])));
      xopt = reshape(xopt, [this.nx, this.N]);
      ucopt = x_sol(logical(kron(ones(this.N,1),[zeros(this.nx,1); ones(this.nuc)])));
      ucopt = reshape(ucopt, [this.nuc, this.N]);
      this.xopt = xopt;
    end

    function [xopt] = returnStates(this)
      xopt = this.xopt;
    end

  end

end
