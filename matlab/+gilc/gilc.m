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

classdef gilc
  %GILC An (almost literal) translation of gILC's Python implementation.

  properties
    % Default: no minimized and constrained outputs
    est = struct('hzc', NaN, 'hzm', NaN, 'Qy', NaN, 'xguess', NaN);
    con = struct('hzm', NaN, 'hzc', NaN, 'Qy', NaN, 'xguess', NaN);
    mode_end, r, nu, na, nx, nxa, ny, xinit, f, h, u_ip1, alpha;
  end

  methods
    function [this, alpha, xopt] = estimation(this, u_i, y_i)
      eststep = gilc.nlpproblem(1, this.nu, this.na, this.nx, this.nxa, this.ny);
      eststep = eststep.setSystemdynamics(this.f, this.h);
      if isa(this.est.hzm, 'function_handle')
        eststep = eststep.setMinimizedoutputs(this.est.hzm, this.est.rz, this.est.Qz);
      end
      if isa(this.est.hzc, 'function_handle')
        eststep = eststep.setConstrainedoutputs(this.est.hzc, this.est.zmin, this.est.zmax);
      end
      if ~isnan(this.est.Qy)
        eststep = eststep.setWeights(this.est.Qy);
      end
      if ~isnan(this.est.xguess)
        eststep = eststep.setInitialguess(this.est.xguess, this.est.alphaguess);
      end
      if strcmp(this.mode_end, 'periodic')
        eststep = eststep.setReference(y_i, this.mode_end);
      else
        eststep = eststep.setInitstate(this.xinit);
        eststep = eststep.setReference(y_i, 'dummy');
      end
      eststep = eststep.setExogin(u_i);
      [eststep, alpha] = eststep.solve();
      xopt = eststep.xopt;
    end

    function [this, u_ip1, xopt] = control(this, alpha)
      constep = gilc.nlpproblem(2, this.nu, this.na, this.nx, this.nxa, this.ny);
      constep = constep.setSystemdynamics(this.f, this.h);
      if isa(this.con.hzm, 'function_handle')
        constep = constep.setMinimizedoutputs(this.con.hzm, this.con.rz, this.con.Qz);
      end
      if isa(this.con.hzc, 'function_handle')
        constep = constep.setConstrainedoutputs(this.con.hzc, this.con.zmin, this.con.zmax);
      end
      if ~isnan(this.con.Qy)
        constep = constep.setWeights(this.con.Qy);
      end
      if ~isnan(this.con.xguess)
        constep = constep.setInitialguess(this.con.xguess, this.con.uguess);
      end
      if strcmp(this.mode_end, 'periodic')
        constep = constep.setReference(this.r, this.mode_end);
      else
        constep = constep.setInitstate(this.xinit);
        constep = constep.setReference(this.r, 'dummy');
      end
      constep = constep.setExogin(alpha);
      [constep, u_ip1] = constep.solve();
      xopt = constep.xopt;
    end

  end

end
