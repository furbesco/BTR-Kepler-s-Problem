# Analytical 1PN framework

In this repository, I have simulated a simple Kepler's problem using two bodies orbiting each other and treating it as a one body system for ease of calculation with 2 degrees of freedom, taking one of the bodies as fixed and the second one as the rotating one. 

In this code, the semi-major axis and eccentricity are constants. Therefore, the calculation yields the following varying parameters:

1. Time
2. Position X
3. Position Y
4. Mean anomaly
5. Kepler's equation
6. Radius
7. True anomaly 

As a next step, the 1PN approximation was implemented analytically in code and used to simulate the orbital evolution with precession to it. There, the precession comes from angle $\Phi$ being shifted by $2\pi + k$ each orbital revolution. 
