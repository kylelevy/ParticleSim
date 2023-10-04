# MATLAB Particle Simulator
This is a simple particle simulator in the MATLAB environment. It was created as an exercise for my Kinematics and Dynamics course and uses the Euler's Method (Momentum Update) to simulate the trajectory of particles under the influence of gravity.

$$F_{net} = \frac{dp}{dt}$$

$$dp = F_{net} * dt$$

Using this property, we can extract the following relations:

$$P_0 = P_1 + F_{net} * dt$$

$$X_0 = X_1 + \frac{P_0}{m} * dt$$

Using these equations, we can simulate the trajectory of these particles by calculating their updated position per incriment time ($dt$) in a loop. 