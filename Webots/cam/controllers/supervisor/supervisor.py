"""supervisor controller."""

from controller import Robot, Camera, RangeFinder
import numpy as np
import math
supervisor = Supervisor()
# create the Robot instance.
robot_node = supervisor.getFromDef("framework1")

# get the time step of the current world.
timestep = int(robot.getBasicTimeStep())

camera = robot_node.getDevice('camera')
camera.enable(4*timestep)
#setter opp rangefinderen
rangeFinder = robot_node.getDevice('3d')
rangeFinder.enable(4*timestep)
#setter opp keyboard detection
keyboard.enable(timestep)
keyboard = getKeyboard()

# Main loop:
# - perform simulation steps until Webots is stopping the controller
while robot.step(timestep) != -1:
    # Read the sensors:
    # Enter here functions to read sensor data, like:
    #  val = ds.getValue()

    # Process sensor data here.

    # Enter here functions to send actuator commands, like:
    #  motor.setPosition(10.0)
    pass

# Enter here exit cleanup code.
