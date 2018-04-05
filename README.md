# Flight Dynamics - Autopilot
Program for calculation and graphical presentation of Lateral and Longitudinal Dynamics  - Autopilot Characteristics

## Introduction
Given the physical characteristics of the aircraft (mass, moments of inertia), the flight parameters (speed and altitide)
and the partial derivatives that describe the behaviour of the aircraft, the program can calculate and draw the response curves
of the aircraft depending of the external factors acting on the aircraft during flight.     
These external factor are represented as force and momentum acting on the aircraft and changing it's velocities and position in
the aircraft's relative Cartesian coordinate system. The center of relative coordinate system coincides with the aircraft's center of mass.     

The shape of the curves can give the description of aircrafts behaviour in flight, whether it will stabilize itself or will it 
deviate from it's trajectory and change altitude or attitude, and the time frame in which these changes will occur.

Based on these curves the autopilot can move the control suffaces of the aircraft and keep it flying inside the flight envelope,
or by the predefined flight parameters. e.g. maintaining constant altitude or speed of descent.

## Installation
Dependencies for the application are python libraries : 
  - scipy
  - numpy
  - matplotlib

On Linux (in terminal):   
    `sudo apt-get install python-scipy python-numpy python-matplotlib`
    
## Run
To run the application execute following command in the terminal:       
    `python Autopilot_Dynamics.py`
