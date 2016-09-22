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

function [fdyn, hdyn, simstep] = pendulum_dynamics(params)
  import casadi.*

  xk = SX.sym('xk',2);
  uk = SX.sym('uk',1);
  ak = SX.sym('ak',1);

  % the discrete state equations; integrator and dynamics
  xkp1 = vertcat( xk(1) + params.Ts*xk(2), ...
                  xk(2) + (params.Ts/(params.I+params.m*params.L^2))*(-params.c*xk(2) - 9.81*params.m*params.L*sin(xk(1)) + (uk(1)+ak(1))) ...
                );
  fdyn = Function('fdyn',{xk, uk, ak}, {xkp1});

  % the output; simply the first state
  y = xk(1);
  hdyn = Function('fdyn',{xk, uk, ak}, {y});

  % a compact form for simulating
  simstep = Function('simstep',{xk, uk, ak}, {xkp1, y});

end
