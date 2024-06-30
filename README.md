# NISAR ADCS Model

This is a project that models the orbit and ADCS system of NASA and ISRO's joint NISAR mission

# Objective
Develop an ADCS for NISAR to achieve 0.075° pointing accuracy for synthetic aperture radar (SAR) operations in LEO.

# Tools
MATLAB, Simulink

# Key Tasks
- Identified satellite specifications and physical characteristics.
- Modeled dynamics and kinematics, analyzed stability, and introduced environmental perturbations.
- Implemented sensor models and MEKF for attitude determination and simulated controller-in-the-loop operations.
- Designed and validated reaction wheel desaturation using magnetorquers.

# Achievements:
- Stabilized the satellite within pointing accuracy constraints despite disturbance torques.
- Demonstrated effective reaction wheel desaturation for extended operations.

# Project layout
- ```doc``` is a folder that contains different reports developed throughout the project development process.
- ```res``` is a folder that contains key properties of the satellite model.
- ```https://tinyurl.com/nisaradcs``` the link where the CAD model of the satellite is stored
- ```src``` is a folder where all relevant source code is located
- ```src/model.slx``` is the simulink file that the ADCS model is stored
- ```src/modelVars.m``` is where key workspace variables needed to run model.slx are stored
- ```src/psX.m``` are the files that generate the plots seen in the report (where X ranges from 1 to 10)

