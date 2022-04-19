
from controller import Robot, Motor, Camera, RangeFinder
from math import pi, sin

class platform (Robot):
    timestep = 32

    def __init__(self):
        super(platform, self).__init__()
        
        #m=self.getDevice("motor")
        self.motorAngle = self.getDevice("motor")
        self.motorAngle.setPosition(float('inf'))
        self.motorAngle.setVelocity(0.0)

        #pSensor=self.getDevice("ps")
        #pSensor.enable(timestep)
        self.pSensor = self.getDevice("ps")
        self.pSensor.enable(self.timestep)
        
        self.motor = self.getDevice("linear")
        #setter opp kamera
        self.camera = self.getDevice('camera')
        self.camera.enable(4*self.timestep)
        #setter opp rangefinderen
        self.rangeFinder = self.getDevice('3d')
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
        F = 2.0   # frequency 2 Hz
        t = 0.0   # elapsed simulation time
        speed=1
        k=0 
        while self.step(self.timestep) != -1:
             # Verical controle
             position = sin(t * 2.0 * pi * F)
             self.motor.setPosition(position)
             t += self.timestep / 3000.0
             # Angle controle 
             self.motorAngle.setVelocity(speed)
             k= self.pSensor.getValue()
             print(k)    
            
 
             
    def hingeJoint(self):  
        speed=1
        k=0 
        while self.step(self.timestep) != -1:
            self.motorAngle.setVelocity(speed)
            k= self.pSensor.getValue()
            print(k)    
            
controller = platform()
controller.camera_movement()
#controller.hingeJoint()
controller.run()
