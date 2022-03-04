"""lidar_test_controller controller."""

# You may need to import some classes of the controller module. Ex:
#  from controller import Robot, Motor, DistanceSensor
from controller import Robot



    
def run_robot(robot): 
    timestep = 32 
    max_speed = 6.28 
    
    left_motor1 = robot.getDevice('front_left_wheel_joint')
    left_motor2 = robot.getDevice('back_left_wheel_joint')

    right_motor1 = robot.getDevice('front_right_wheel_joint')
    right_motor2 = robot.getDevice('back_right_wheel_joint')

    left_motor1.setPosition(float('inf'))
    left_motor2.setPosition(float('inf'))
    right_motor1.setPosition(float('inf'))
    right_motor2.setPosition(float('inf'))
    
    left_motor1.setVelocity(0.0)
    left_motor2.setVelocity(0.0)
    right_motor1.setVelocity(0.0)
    right_motor2.setVelocity(0.0)
    
    lidar = robot.getDevice('lidar(1)')
    lidar.enable(timestep)
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
        rangeimage = lidar.getRangeImage()
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
