"""camera controller."""

from controller import Robot, Camera, RangeFinder
import numpy as np
import math

#lagger en klasse som brukes for alle robotene, navnet på child komponentene må vere likt    
class platform (Robot):
    timestep = 32
    def __init__(self):
        super(platform, self).__init__()
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
        if self.ID == 'framework2':
            print('robot2 identified')
            self.filename = ['cam2_dist.txt','cam2_rgb.txt']
            image_width = self.rangeFinder.getWidth()
            image_height = self.rangeFinder.getHeight()
            focal_length = 0.5 * image_width * (1 / math.tan(0.5 * self.rangeFinder.getFov()))
            k_matrix = np.array([
                [focal_length, 0, image_width / 2],
                [0, focal_length, image_height / 2],
                [0, 0, 0]
            ])
            file = open('k_matrix.txt', "w+")
            for row in k_matrix:
                np.savetxt(file, row)
            file.close()
            
        elif self.ID == 'framework1':
            print('robot1 identified')
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
                #content = str(self.depth)
                #file.write(content)
                for row in self.depth:
                    np.savetxt(file, row)
                file.close()
                file = open(self.filename[1], "w+")
                #content = str(self.image)
                #file.write(content)
                for row in self.image:
                    np.savetxt(file, row)
                file.close()
                

            if self.step(self.timestep) == -1:
                break

controller = platform()
controller.run()
            