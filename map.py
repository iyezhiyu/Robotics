#!/usr/bin/env python 

# Some suitable functions and data structures for drawing a map and particles

import time
#import random
import math
import numpy as np
from numpy import *
# Functions to generate some dummy particles data:
def calcX():
    return random.gauss(80,3) + 70*(math.sin(t)); # in cm

def calcY():
    return random.gauss(70,3) + 70*(math.sin(t)); # in cm

def calcW():
    return random.random();

def calcTheta():
    return random.randint(0,360);

# A Canvas class for drawing a map and particles:
#     - it takes care of a proper scaling and coordinate transformation between
#      the map frame of reference (in cm) and the display (in pixels)
class Canvas:
    def __init__(self,map_size=210):
        self.map_size    = map_size;    # in cm;
        self.canvas_size = 768;         # in pixels;
        self.margin      = 0.05*map_size;
        self.scale       = self.canvas_size/(map_size+2*self.margin);

    def drawLine(self,line):
        x1 = self.__screenX(line[0]);
        y1 = self.__screenY(line[1]);
        x2 = self.__screenX(line[2]);
        y2 = self.__screenY(line[3]);
        print "drawLine:" + str((x1,y1,x2,y2))

    def drawParticles(self,data):
        display = [(self.__screenX(d[0]),self.__screenY(d[1])) + tuple(d[2:]) for d in data];
        print "drawParticles:" + str(display);

    def __screenX(self,x):
        return (x + self.margin)*self.scale

    def __screenY(self,y):
        return (self.map_size + self.margin - y)*self.scale

# A Map class containing walls
class Map:
    def __init__(self):
        self.walls = [];

    def add_wall(self,wall):
        self.walls.append(wall);

    def clear(self):
        self.walls = [];

    def draw(self):
        for wall in self.walls:
            canvas.drawLine(wall);

# Simple Particles set
class Particles:
    def __init__(self):
        self.n = 100;    
        self.data = [];

    def update(self):
        self.data = [[calcX(), calcY(), calcTheta(), calcW()] for i in range(self.n)];
        
    def start_point(self,x,y,theta):
        self.data = [[x,y,theta,0.01] for i in range(self.n)]
    
    def draw(self):
        canvas.drawParticles(self.data);
        
    def set_w(self,i,w):
        self.data[i][3] = w
    
    def mean_x(self):
        sum_x = 0.
        #print(len(self.data))
        #print(len(self.data[0]))
        #print(self.data)
        for i in range(self.n):
            sum_x = sum_x + self.data[i][0]*self.data[i][3]
            
        return sum_x
    
    def mean_y(self):
        sum_y = 0.
        for i in range(self.n):
            sum_y = sum_y + self.data[i][1]*self.data[i][3]
            
        return sum_y

    def mean_theta(self):
        sum_theta = 0.
        for i in range(self.n):
            '''
            while self.data[i][2]>np.pi:
                self.data[i][2] = self.data[i][2] - 2*np.pi
            
            while self.data[i][2]<= -np.pi:
                self.data[i][2] = self.data[i][2] + 2*np.pi
            '''
            sum_theta = sum_theta + self.data[i][2]*self.data[i][3]
            
            while sum_theta >np.pi:
                sum_theta = sum_theta -2*np.pi
            while sum_theta<=-np.pi:
                sum_theta = sum_theta + 2*np.pi
                
        print sum_theta   
        return sum_theta
    
#   def __getitem__(self, item):
#       return self.Particles[item]
    
    def __getitem__(self, item):
        return self.item

    
    def set_one_data(self,x,y,theta,w):
        
        a = [x,y,theta,w]
        self.data.append(a)
        
canvas = Canvas();
        
        
'''
canvas = Canvas();    # global canvas we are going to draw on



mymap = Map();
# Definitions of walls
# a: O to A
# b: A to B
# c: C to D
# d: D to E
# e: E to F
# f: F to G
# g: G to H
# h: H to O
mymap.add_wall((0,0,0,168));        # a
mymap.add_wall((0,168,84,168));     # b
mymap.add_wall((84,126,84,210));    # c
mymap.add_wall((84,210,168,210));   # d
mymap.add_wall((168,210,168,84));   # e
mymap.add_wall((168,84,210,84));    # f
mymap.add_wall((210,84,210,0));     # g
mymap.add_wall((210,0,0,0));        # h
mymap.draw();

particles = Particles();

t = 0;
while True:
    particles.update();
    particles.draw();
    t += 0.05;
    time.sleep(0.05);
'''
