
from controller import Robot, Motor, Camera, RangeFinder
from math import pi, sin

#from controller import Robot, Camera, RangeFinder

#TIME_STEP = 32

#robot = Robot()
#motor = robot.getDevice("linear")



#while robot.step(TIME_STEP) != -1:
#    position = sin(t * 2.0 * pi * F)
#    motor.setPosition(position)
#    t += TIME_STEP / 1000.0
    
   

#lagger en klasse som brukes for alle robotene, navnet på child komponentene må vere likt    
class platform (Robot):
    timestep = 32
    #F = 2.0   # frequency 2 Hz
    #t = 0.0   # elapsed simulation time
    #position = sin(t * 2.0 * pi * F)
    def __init__(self):
        super(platform, self).__init__()
        
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
            
            #while self.step(self.timestep) != -1:
            #    position = sin(t * 2.0 * pi * F)  
            #    self.motor.setPosition(position)
            #    t += self.timestep / 1000.0
            
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
    
    def test(self):
        F = 2.0   # frequency 2 Hz
        t = 0.0   # elapsed simulation time
        while self.step(self.timestep) != -1:
             position = sin(t * 2.0 * pi * F)
             self.motor.setPosition(position)
             t += self.timestep / 3000.0
            
controller = platform()
controller.test()
controller.run()
