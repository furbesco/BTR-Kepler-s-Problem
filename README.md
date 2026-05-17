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

## How to use the parameters!

The user has the choice of whether they want to run a specific number of orbits vs specific time in months vs time in years. To do that, go to the configuration file. There, plug in the value you want to run, and the other two parameters set to -1 so that they don't work.

## Running the files

To run the files in the correct order, it is recommended to first run the Keplerian file (which simulates the simple Newtonian problem) or the 1PN file (which incorporates the first-order relativistic approximation). Next, the post-processing file is run.

### 2PN Approximation

The 2PN approximation is currently under making. As of now, it still does not fucntion fully properly, therefore the file is still undergoing changes to it. 
