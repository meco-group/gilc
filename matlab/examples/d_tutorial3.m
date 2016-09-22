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

clear all
clear global
clc

%% Imports
import casadi.*
addpath('../')
import gilc.*

%% Setup
% plant; "the real plant"
plant = struct();
% ... parameters
plant.params.I = 0.006;
plant.params.m = 0.5;
plant.params.L = 0.3;
plant.params.c = 0.1;
plant.params.Ts = 1/500.0;
% ... dynamics
[plant.fdyn, plant.hdyn, plant.simstep] = pendulum_dynamics(plant.params);

% model; "the theoretical model"
model = struct();
% ... parameters incl. deviations from plant
model.params = plant.params;
model.params.m = 0.495;
model.params.c = 0.105;
% ... dynamics
[model.fdyn, model.hdyn, model.simstep] = pendulum_dynamics(model.params);

% construct the gILC object
ilc = gilc();

% read the reference
ref = csvread('reference.txt').';

%% ILC settings
trials = 5;                            % Number of trials in the simulation
N = 1000;                              % Number of evaluation points
fs = 500;                              % Sampling frequency
time = linspace(0,(N-1)/fs,N);         % Resulting time vector

ilc.r = interp1(linspace(0,1,numel(ref)),ref,linspace(0,1,N));
ilc.nu = 1;
ilc.na = 1;
ilc.nx = 2;
ilc.nxa = 0;
ilc.ny = 1;
ilc.xinit = zeros(2,1);
ilc.f = model.fdyn;
ilc.h = model.hdyn;

% add constraints on input signal
ilc.con.hzc  = @(x,u,a) (u);           % define a constrained output function
ilc.con.zmin = -2*ones(1,N);
ilc.con.zmax = +2*ones(1,N);

% add penalty to control signal
ilc.con.hzm = @(x,u,a) (u);            % define a minimized output function
ilc.con.rz  = zeros(1,N);
ilc.con.Qz  = 1e2*ones(1,N);

% Storage allocation
uSave = zeros(trials,N);               % Initialize saved inputs
ySave = zeros(trials,N);               % Initialize saved outputs
eSave = zeros(trials,N);               % Initialize saved tracking errors
normeSave= zeros(trials,1);            % Initialize saved 2-norm of tracking error
aSave = zeros(trials,N);               % Initialize saved correction signals

%% Running the ILC
a_i = zeros(ilc.na, N);                % Initialize the correction signal for the first trial
[ilc, u_i, ~] = ilc.control(a_i);      % Solve the control step for the first trial
for i = 1:1:trials
  x_i = NaN(ilc.nx, N);
  y_i = NaN(ilc.ny, N);
  x_in = ilc.xinit;
  for n = 1:1:N;
    [x_in, y_in] = plant.simstep(x_in, u_i(:,n), 0);     % Simulate the system
    x_i(:,n) = full(x_in);
    y_i(:,n) = full(y_in);
  end
  %------
  uSave(i,:) = u_i;                    % Save the applied input signal
  ySave(i,:) = y_i;                    % Save the measured output signal
  eSave(i,:) = ilc.r - y_i;            % Save the tracking error
  normeSave(i) = norm(ilc.r - y_i);    % Save the 2-norm of the tracking error
  aSave(i,:) = a_i;                    % Save the correction signal
  %------
  [ilc, a_i, ~] = ilc.estimation(u_i,y_i);  % Solve the estimation problem
  [ilc, u_i, ~] = ilc.control(a_i);         % Solve the control problem
end

%% Plotting
figure(1)
subplot(3,2,1)
  plot(time,uSave(1,:),'k--')
  hold on
  plot(time,uSave(trials,:),'k')
  grid on
  ylabel('Input [Nm]')
  xlabel('Time [s]')
  legend('Trial 1','Trial 5')
%------
subplot(3,2,3)
  plot(time,aSave(2,:),'k--')
  hold on
  plot(time,aSave(trials,:),'k')
  grid on
  ylabel('Correction signal [Nm]')
  xlabel('Time [s]')
  legend('Trial 2','Trial 5')
%------
subplot(3,2,5)
  plot((1:1:trials), normeSave, 'ko-')
  grid on
  xlim([1, trials])
  xlabel('Iteration number [-]')
  ylabel('2-norm of tracking error [rad]')
%------
subplot(5,2,2)
  plot(time, (180/pi)*ilc.r, 'k--')
  hold on
  plot(time, (180/pi)*ySave(1,:),'k')
  grid on
  ylabel('Output [deg]')
  legend('Reference','Trial 1')
%------
subplot(5,2,4)
  plot(time, (180/pi)*ilc.r,'k--')
  hold on
  plot(time, (180/pi)*ySave(2,:),'k')
  grid on
  ylabel('Output [deg]')
  legend('Reference','Trial 2')
%------
subplot(5,2,6)
  plot(time, (180/pi)*ilc.r,'k--')
  hold on
  plot(time, (180/pi)*ySave(3,:),'k')
  grid on
  ylabel('Output [deg]')
  legend('Reference','Trial 3')
%------
subplot(5,2,8)
  plot(time, (180/pi)*ilc.r,'k--')
  hold on
  plot(time, (180/pi)*ySave(4,:),'k')
  grid on
  ylabel('Output [deg]')
  legend('Reference','Trial 4')
%------
subplot(5,2,10)
  plot(time, (180/pi)*ilc.r,'k--')
  hold on
  plot(time, (180/pi)*ySave(5,:),'k')
  grid on
  xlabel('Time [s]')
  ylabel('Output [deg]')
  legend('Reference','Trial 5')
