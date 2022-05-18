
"""lidar_test_controller controller."""

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

    
    # Main loop:
    while robot.step(timestep) != -1:

        range_image = lidar.getRangeImage()
        print('{}'.format(range_image))

            
        left_motor.setVelocity(max_speed*0.25)
        right_motor.setVelocity(max_speed*0.25)


if __name__== "__main__":

    my_robot = Robot()
    run_robot(my_robot)