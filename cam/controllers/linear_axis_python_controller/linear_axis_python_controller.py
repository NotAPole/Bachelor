
from controller import Robot, Motor, Camera, RangeFinder
from math import pi, sin

class platform (Robot):
    timestep = 32

    def __init__(self):
        super(platform, self).__init__()
        
        self.motorAngle = self.getDevice("motor")
        self.motorAngle.setPosition(float('inf'))
        self.motorAngle.setVelocity(0.0)

        self.pSensor = self.getDevice("ps")
        self.pSensor.enable(self.timestep)
        
        self.vSensor = self.getDevice("vs")
        self.vSensor.enable(self.timestep)
        
        self.motor = self.getDevice("linear")
        #setter opp kamera
        self.camera = self.getDevice('camera1')
        self.camera.enable(4*self.timestep)
        #setter opp rangefinderen
        self.rangeFinder = self.getDevice('3d1')
        self.rangeFinder.enable(4*self.timestep)
        #setter opp keyboard detection
        self.keyboard.enable(self.timestep)
        self.keyboard = self.getKeyboard()
        #lar deg gjøre forskjellige ting med hver robot
        self.ID = self.getName()
        if self.ID == 'metal_pole':
            print('The linear axis identified')
            self.filename = ['cam1_dist.txt','cam1_rgb.txt']
        
    def run(self):
        while True:
            #2d array med alle dybder
            self.depth = self.rangeFinder.getRangeImageArray()
            #array med RGB info for siste bilde
            #formen er array[x][y][n] hvor n = 0 er rød, 1 er grøn og 2 er blå
            self.image = self.camera.getImageArray()
            
            k = self.keyboard.getKey()

            
            if k == ord('Q'):
                file = open(self.filename[0], "w+")
                content = str(self.depth)
                file.write(content)
                file.close()
                file = open(self.filename[1], "w+")
                content = str(self.image)
                file.write(content)
                file.close()
            if self.step(self.timestep) == -1:
                break
    def camera_movement(self):

  
        speed=1
        vertical_value = 0 
        horizontal_value = 0
        position = 0.0
        constant = 1
        while self.step(self.timestep) != -1:

            n = self.keyboard.getKey()  
            
            # Verical controle
            self.motor.setPosition(position)
            vertical_value = self.vSensor.getValue()
            print("Vertical position", vertical_value)
            
            # Angle controle 
            self.motorAngle.setVelocity(speed)
            horizontal_value = self.pSensor.getValue()
            print("Horizontal angle", horizontal_value)   
            
            if n == ord('Q'):
                print("Going down")
                position -= constant
                speed -= 1

         

            if n == ord('O'):
                print("Going up")
                position += constant
                speed += 1
   
            
            pass   
 
            
controller = platform()
controller.camera_movement()
#controller.run()