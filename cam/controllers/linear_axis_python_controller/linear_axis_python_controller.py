
from controller import Robot, Motor
from math import pi, sin

TIME_STEP = 32

robot = Robot()
motor = robot.getDevice("linear")

F = 2.0   # frequency 2 Hz
t = 0.0   # elapsed simulation time

while robot.step(TIME_STEP) != -1:
    position = sin(t * 2.0 * pi * F)
    motor.setPosition(position)
    t += TIME_STEP / 1000.0