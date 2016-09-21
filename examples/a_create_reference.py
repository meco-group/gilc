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
import matplotlib.pyplot as plt

# Create the reference signal
#----------------------------
newpos = 180                                        # New position [deg]
twait = 190                                         # Wait time before and after the motion [number of samples]
tac = 30                                            # Time to accelerate [number of samples]
tvel = 500                                          # Time for constant angular velocity [number of samples]
thold = 500                                         # Time to hold the new position [number of samples]
N = 2*twait + 4*tac+2*tvel+thold                    # Total number of samples in the signal
print("N = ")
print(N)
fs = 500.0                                          # Sample time [Hz]
time = linspace(0,(N-1)/fs,N)                       # Time vector
#------
v = zeros(N)                                        # Velocity profile
v[twait:twait+tac] = linspace(0,tac,tac)
v[twait+tac:twait+tac+tvel] = tac
v[twait+tac+tvel:twait + 2.0*tac+tvel] = linspace(tac,0,tac)
v[twait+ 2*tac+tvel+thold:twait + 3.0*tac+tvel+thold] = linspace(0,-tac,tac)
v[twait+ 3*tac+tvel+thold:twait+3.0*tac+2.0*tvel+thold] = -tac
v[twait+ 3*tac+2.0*tvel+thold:twait+4.0*tac+2.0*tvel+thold] = linspace(-tac,0,tac)
ref = cumsum(v)                                     # Cumulative sum of velocity profile
ref = newpos*(pi/180)*ref/amax(ref)                 # Scale to reach newpos
#----------------------------

# Write the reference signal to file
#-----------------------------------
with open('reference.txt','w') as f:
    for k in range(N):
        f.write(str(ref[k]))
        f.write('\n')
#-----------------------------------

# Plotting
#---------
plt.figure(1)
plt.clf()
plt.plot(time,(180/pi)*ref.T,'k')
plt.grid()
plt.ylim([(180/pi)*min(ref)-20,(180/pi)*max(ref)+20])
plt.ylabel('Reference output [deg]')
plt.xlabel('Time [s]')
plt.show()
#---------
