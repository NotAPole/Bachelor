
"""lidar_test_controller controller."""

# You may need to import some classes of the controller module. Ex:
#  from controller import Robot, Motor, DistanceSensor
from controller import Robot
    
def run_robot(robot): 
    timestep = 32 
    max_speed = 6.28 
    
    left_motor = robot.getDevice('motor_1')
    right_motor = robot.getDevice('motor_2')

    left_motor.setPosition(float('inf'))
    right_motor.setPosition(float('inf'))
    
    left_motor.setVelocity(0.0)
    right_motor.setVelocity(0.0)
    
    
    
    lidar = robot.getLidar('lidar')
    lidar.enable(timestep)
    lidar.enablePointCloud()
    # You should insert a getDevice-like function in order to get the
    # instance of a device of the robot. Something like:
    #  motor = robot.getDevice('motorname')
    #  ds = robot.getDevice('dsname')
    #  ds.enable(timestep)
    
    # Main loop:
    # - perform simulation steps until Webots is stopping the controller
    while robot.step(timestep) != -1:
        # Read the sensors:
        # Enter here functions to read sensor data, like:
        #  val = ds.getValue()
        range_image = lidar.getRangeImage()
        print('{}'.format(range_image))
        # Process sensor data here.
    
        # Enter here functions to send actuator commands, like:
        #  motor.setPosition(10.0)
            
        left_motor.setVelocity(max_speed*0.25)
        right_motor.setVelocity(max_speed*0.25)

# Enter here exit cleanup code.

if __name__== "__main__":

    my_robot = Robot()
    run_robot(my_robot)