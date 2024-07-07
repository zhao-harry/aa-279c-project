# Satellite ADCS Project: NISAR

This project follows the modeling of the orbital and rotational mechanics of a spacecraft and introduces the design of a conventional attitude determination and control system. We choose the NASA/ISRO NISAR mission as reference.

# Objective
Develop an ADCS for NISAR to achieve 0.075° pointing accuracy for Earth-pointing synthetic aperture radar (SAR) in LEO.

# Tools
MATLAB, Simulink

# Key Tasks
- Identified satellite specifications and physical characteristics
- Modeled dynamics and kinematics, analyzed stability, and introduced environmental perturbations
- Implemented sensor models and MEKF for attitude determination and simulated controller-in-the-loop operations
- Designed and validated reaction wheel desaturation using magnetorquers

# Achievements:
- Stabilized the satellite within pointing accuracy constraints despite disturbance torques
- Demonstrated effective reaction wheel desaturation for extended operations

# Project layout
- ```doc```: problem set reports developed over the course of the project
- ```res```: figures generated from experiments and simulations
- ```https://tinyurl.com/nisaradcs```: CAD model of the satellite
- ```src```: source code
- ```src/model.slx```: Simulink model of the satellite
- ```src/modelVars.m``` important workspace variables, such as those representing physical characteristics of the satellite and ADCS parameters
- ```src/psX.m``` problem set files which run experiments and generate figures
