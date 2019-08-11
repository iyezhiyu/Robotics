import brickpi
#from brickpi import * 
import time
#import numpy as np
#from numpy import *
import sys
#import random
import math
from map import *


interface=brickpi.Interface()
interface.initialize()
motors = [3,0]

port = 0

interface.motorEnable(motors[0])
interface.motorEnable(motors[1])

k_left_i = 500
k_right_i = 500


motorParams0 = interface.MotorAngleControllerParameters()
motorParams0.maxRotationAcceleration = 6.0
motorParams0.maxRotationSpeed = 12.0
motorParams0.feedForwardGain = 255/20.0
motorParams0.minPWM = 17.882
motorParams0.pidParameters.minOutput = -255
motorParams0.pidParameters.maxOutput = 255
motorParams0.pidParameters.k_p = 500
motorParams0.pidParameters.k_i = k_right_i#400# k_right_i 350 10
motorParams0.pidParameters.k_d = 1



motorParams1 = interface.MotorAngleControllerParameters()

motorParams1.maxRotationAcceleration = 6.0
motorParams1.maxRotationSpeed = 12.0
motorParams1.feedForwardGain = 255/20.0
motorParams1.minPWM = 18.0
motorParams1.pidParameters.minOutput = -255
motorParams1.pidParameters.maxOutput = 255
motorParams1.pidParameters.k_p = 500
motorParams1.pidParameters.k_i = k_left_i#30#10 #k_left_i
motorParams1.pidParameters.k_d = 1


interface.setMotorAngleControllerParameters(motors[0],motorParams0)
interface.setMotorAngleControllerParameters(motors[1],motorParams1)
#interface.setMotorRotationSpeedReferences(motors,[6,6])
#angle = 14.5 #float(input("Enter a angle to rotate (in radians): "))

#interface.increaseMotorAngleReferences(motors,[angle,angle])
#interface.startLogging("kil_"+str(int(motorParams0.pidParameters.k_i))+"_File.log")
#while True:
#while not interface.motorAngleReferencesReached(motors) :

#    motorAngles = interface.getMotorAngles(motors)
#    if motorAngles :
#        print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]    
#    time.sleep(0.01)




#left90deg()

#time.sleep(1)

#right90deg()

#interface.stopLogging()

#interface.stopLogging()
#interface.terminate()


def forward40():
    count = 0
    angle = -11.5
    interface.increaseMotorAngleReferences(motors,[angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1
    
    while not interface.motorAngleReferencesReached(motors) :
        motorAngles = interface.getMotorAngles(motors)
        
        if motorAngles :
            pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
#        count = count + 1
        if count % 20 == 0:
            if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                break
            old_angle1 = motorAngles[0][0]
            old_angle2 = motorAngles[1][0]

        count = count + 1

        time.sleep(0.05)
                #if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                #        break
                #old_angle1 = motorAngles[0][0]
                #old_angle2 = motorAngles[1][0]

    print "Destination reached!"

def backward40():
    count = 0
    angle = 11.5
    interface.increaseMotorAngleReferences(motors,[angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1

    while not interface.motorAngleReferencesReached(motors):
            motorAngles = interface.getMotorAngles(motors)

            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
           # count = count + 1
            if count % 20 == 0:
                    if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                            break
                    old_angle1 = motorAngles[0][0]
                    old_angle2 = motorAngles[1][0]

            count = count + 1

            time.sleep(0.05)

#                time.sleep(0.2)
#                if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
#                        break
#                old_angle1 = motorAngles[0][0]
#                old_angle2 = motorAngles[1][0]

    print "Destination reached!"

def backwardcm(a):
    count = 0
    angle = a*11.5/40
    interface.increaseMotorAngleReferences(motors,[angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1

    while not interface.motorAngleReferencesReached(motors) :
            motorAngles = interface.getMotorAngles(motors)

            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
           # count = count + 1
            if count % 20 == 0:
                    if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                            break
                    old_angle1 = motorAngles[0][0]
                    old_angle2 = motorAngles[1][0]

            count = count + 1

            time.sleep(0.05)

#                time.sleep(0.2)
#                if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
#                        break
#                old_angle1 = motorAngles[0][0]
#                old_angle2 = motorAngles[1][0]

    print "Destination reached!"

def forwardcm(a):
    count = 0
    angle = -a*11.5/40
    interface.increaseMotorAngleReferences(motors,[angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1

    while not interface.motorAngleReferencesReached(motors) :
            motorAngles = interface.getMotorAngles(motors)

            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
           # count = count + 1
            if count % 20 == 0:
                    if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                            break
                    old_angle1 = motorAngles[0][0]
                    old_angle2 = motorAngles[1][0]

            count = count + 1

            time.sleep(0.05)

#                time.sleep(0.2)
#                if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
#                        break
#                old_angle1 = motorAngles[0][0]
#                old_angle2 = motorAngles[1][0]

    print "Destination reached!"


def right90deg():
    angle = 2.88    
    interface.increaseMotorAngleReferences(motors,[-angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1
    while not interface.motorAngleReferencesReached(motors) :

        motorAngles = interface.getMotorAngles(motors)
        if motorAngles :
            pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]
        
        if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
            break        
        
        old_angle1 = motorAngles[0][0]
        old_angle2 = motorAngles[1][0]
        time.sleep(0.2)
    
    print "Rotation reached!"

    



def left90deg():
    angle = 2.88
    interface.increaseMotorAngleReferences(motors,[angle,-angle])
    old_angle1 = -0.1
    old_angle2 = 0.1
    while not interface.motorAngleReferencesReached(motors) :

            motorAngles = interface.getMotorAngles(motors)
            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]

            if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                    break

            old_angle1 = motorAngles[0][0]
            old_angle2 = motorAngles[1][0]
            time.sleep(0.2)

    print "Rotation reached!"


def leftTurndeg(deg=0.):
    angle = deg*2.88/90
    interface.increaseMotorAngleReferences(motors,[angle,-angle])
    old_angle1 = -0.1
    old_angle2 = 0.1
    while not interface.motorAngleReferencesReached(motors) :

            motorAngles = interface.getMotorAngles(motors)
            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]

            if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                    break

            old_angle1 = motorAngles[0][0]
            old_angle2 = motorAngles[1][0]
            time.sleep(0.2)

    print "Rotation finished!"


def rightTurndeg(deg=0.):
    angle = deg*2.88/90
    interface.increaseMotorAngleReferences(motors,[-angle,angle])
    old_angle1 = -0.1
    old_angle2 = 0.1
    while not interface.motorAngleReferencesReached(motors) :

            motorAngles = interface.getMotorAngles(motors)
            if motorAngles :
                    pass#print "Motor angles: ", motorAngles[0][0], ", ", motorAngles[1][0]

            if motorAngles[0][0]==old_angle1 and motorAngles[1][0]==old_angle2:
                    break

            old_angle1 = motorAngles[0][0]
            old_angle2 = motorAngles[1][0]
            time.sleep(0.2)

    print "Rotation finished!"


def touch():
    touch_port1 = 1
    touch_port2 = 0
    speed = 6
    interface.sensorEnable(touch_port1,brickpi.SensorType.SENSOR_TOUCH)
    interface.sensorEnable(touch_port2,brickpi.SensorType.SENSOR_TOUCH)

    #interface.setMotorRotationSpeedReferences(motors,[speed,speed])
    
    while True:
        
            #interface.setMotorRotationSpeedReferences(motors,[speed,speed])

        result1 = interface.getSensorValue(touch_port1)[0]
        result2 = interface.getSensorValue(touch_port2)[0]
        print result1
        print result2
                
        if not result1 and result2:
            print "touched 1"
            time.sleep(1)
            #backwardcm(10)
            #left90deg()
            continue
        if not result2 and result1:
            print "touched 2"
            time.sleep(1)
            #backwardcm(10)
            #right90deg()
            continue
        if  result1 and  result2:
            print "touched 1 and 2"
            time.sleep(1)
            #backwardcm(10)
            continue
                    
        print "run!"
        time.sleep(1)
    
    
def us():
    ref = 38
    port = 0
    interface.sensorEnable(port,brickpi.SensorType.SENSOR_ULTRASONIC)
    vc = 6.0    
    Kp = 0.1
    lastReading = []
    while True:
        
        usReading = float(interface.getSensorValue(port)[0])
        '''
        interface.setMotorRotationSpeedReferences(motors,[vl,vr])
        lastReading.append(usReading)
        if(len(lastReading)>11):
            del lastReading[0]
        xx = median(lastReading)

        if(xx>100):
            continue
        e = ref - xx
         
        vl = vc - 1.0/2 * Kp * e
        vr = vc + 1.0/2 * Kp * e
        #interface.setMotorRotationSpeedReferences(motors,[vr,vl])
        '''
        if usReading:
            print usReading #,'      ',xx, '      ',e,'            ',vl,'          ',vr
        else:
            print 'Failed us reading'

        time.sleep(0.02)

def q5_3():
    forward40()
    time.sleep(1)
    backward40()
    time.sleep(1)
    left90deg()
    time.sleep(1)
    right90deg()
    interface.terminate()

def q5_4():

    interface.startLogging("kil_"+str(int(motorParams0.pidParameters.k_i))+"_File.log")

    forward40()
    time.sleep(0.5)
    left90deg()
    time.sleep(0.5)
    forward40()
    time.sleep(0.5)
    left90deg()
    time.sleep(0.5)
    forward40()
    time.sleep(0.5)
    left90deg()
    time.sleep(0.5)
    forward40()
    time.sleep(0.5)
    left90deg()
    time.sleep(0.5)
    interface.stopLogging()

def getRandomE():
    return np.random.normal(0, 1)

def getRandomF():
    return np.random.normal(0, 0.05)

def getRandomTheta(): 
    return np.random.normal(0,0.05)
    
    
def dist(last_angle1,last_angle2):
    
    motorAngles = interface.getMotorAngles(motors)
    angle1 = motorAngles[0][0]
    angle2 = motorAngles[1][0]
    angle = 0.5*(angle1+angle2)
    
    return (angle - 0.5*(last_angle1+last_angle2))
    
    
def straightLine(data,dis,w):
    
    ttheta = data[2]+getRandomF()
    x = data[0] + (dis+getRandomE())*cos(data[2])
    y = data[1] + (dis+getRandomE())*sin(data[2])
    
    return [x,y,ttheta,w]

def rotation(point,alpha,beta,w):
    #if alpha>=0:
    return [point[0],point[1],point[2]+beta+getRandomTheta(),w]
   # else:
     #   return [point[0],point[1],point[2]-beta+getRandomTheta(),w]
    
    
    
def plot_forward(_position,distance):
    motorAngles = interface.getMotorAngles(motors)
    
    last_angle1 = motorAngles[0][0]
    last_angle2 = motorAngles[1][0]   

    forwardcm(20)
    
    for i in range(100):
        _position[i] = straightLine(_position[i][0],_position[i][1],_position[i][2],distance)
    


    
    return _position


'''
def plot_rotation(_position,alpha):
    
    
    motorAngles = interface.getMotorAngles(motors)
    last_angle1 = motorAngles[0][0]
    last_angle2 = motorAngles[1][0]   

    #left90deg()
    if alpha > 0:
        
    
    leftTurndeg(100)
    
    distance = dist(last_angle1,last_angle2)
    
    for i in range(100):
        _position[i] = rotation(_position[i][0],_position[i][1],_position[i][2],1.57)


    print "drawParticles:" + str(_position)
    
    
    return _position
'''

'''
def cw2():

    points = []
    for i in range(100):
        points.append((300,700,0))
    
        
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_rotation(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_rotation(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_rotation(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    points = plot_forward(points)
    
'''    
'''    
def navigateToWaypoint(X,Y,lp):
    
    dx = X - lp[0]
    dy = Y - lp[1]        
    alpha = math.atan2(dy,dx)
    beta = alpha - lp[2]
    
    if beta > np.pi:
        beta = beta - np.pi * 2.
    if beta < -np.pi:
        beta = beta + np.pi * 2.

    dist = np.sqrt(dx*dx+dy*dy)

    if alpha >= 0:
        leftTurndeg(beta*180/(np.pi))
    else:
        rightTurndeg(beta*180/(np.pi))
    
    while dist > 20:
        forwardcm(20)
        dist = dist - 20
    forwardcm(dist)
    
    return X,Y,alpha

def navi():
    
    lp = (0,0,0)
    
    while True:
        X = input()
        Y = input()
        
        print(lp)

'''
        
def init_map():
    mymap = Map()
    mymap.add_wall((0,0,0,168));        # a
    mymap.add_wall((0,168,84,168));     # b
    mymap.add_wall((84,126,84,210));    # c
    mymap.add_wall((84,210,168,210));   # d
    mymap.add_wall((168,210,168,84));   # e
    mymap.add_wall((168,84,210,84));    # f
    mymap.add_wall((210,84,210,0));     # g
    mymap.add_wall((210,0,0,0));        # h
    mymap.draw();
    
    return mymap


def cal_m(x,y,theta,B):
    
    mm = np.inf
    num = -1
    flag = 0
    for i in range(8):
        #print('******************')
        #print(B)
        
        ax = B[i][0]
        ay = B[i][1]
        bx = B[i][2]
        by = B[i][3]

        m = ((by-ay)*(ax-x)-(bx-ax)*(ay-y))/((by-ay)*cos(theta)-(bx-ax)*sin(theta))
        
        
        v1x = ax - x
        v1y = ay - y
        
        v2x = bx - x
        v2y = by - y
        
        if(v1x*v2x+v1y*v2y <= 0) and m<mm:
            
            mm = m
            num = i
    
    
    ax = B[num][0]
    ay = B[num][1]
    bx = B[num][2]
    by = B[num][3]   
        
    beta = np.arccos((cos(theta)*(ay-by)+sin(theta)*(bx-ax))/(np.sqrt((ay-by)**2+(bx-ax)**2)))
    
    if beta > 0.5 or beta < -0.5:
        #print('invalid measurement')
        flag = 0
        # if beta>30 ?
    
    return mm,num,flag#############################################################
    
    
        
def calculate_likelihood(x,y,theta,z,mymap):
    
        
    sigma = 1
    K = 0.01
    
    m = cal_m(x,y,theta,mymap.walls)
        
        
    l = exp(-(z-m[0])**2/(2*sigma**2)) + K        
        
    return l,m[1]

def resample(particles):
    
    mx = particles.mean_x()
    my = particles.mean_y()
    mtheta = particles.mean_theta()
    z = float(interface.getSensorValue(port)[0])
    print("zzzzzzzzzzzzzzzzz")
    print(z)
    #z = 255
    flag = 0#cal_m(mx,my,mtheta,mymap.walls)[2]
    
    lp = [mx,my,mtheta]
    
    if(flag):
        return particles.data,lp
        
    if z > 150 and z < 10:
        return particles.data,lp
        
    z = z-6
    
    #*********normalize**************
    print('start normalize')
    print(particles.data)
    sum_w = 0.
    for i in range(100):
        x = particles.data[i][0]
        y = particles.data[i][1]
        theta = particles.data[i][2]
        w = particles.data[i][3]


        ll = calculate_likelihood(x,y,theta,z,mymap)
        particles.set_w(i,ll[0]*w)
        sum_w = sum_w + ll[0]*w
    
    
    cdf = [0]
    sam = [[0]*100]
    sam = sam[0]
    pb = []


    for i in range(100):

        w = particles.data[i][3] / sum_w
        particles.set_w(i,w)
        cdf.append(cdf[i]+w)
        pb.append(w)
    
    
    print('normalize finished')
    print(particles.data)
    #************resampling*************
    print('start reampling')
    
    for i in range(100):
        t = random.uniform(0,1)
        l = 0
        r = 100
        
        while True:
            m = int((l+r)/2)
            if t>cdf[m] and r-l > 1:
                l = m
                continue
            if t<cdf[m] and r-l > 1:
                r = m
                continue
            if r-l == 1:
                sam[l] = sam[l] + 1
                break
                
            if t==cdf[m]:
                sam[m] = sam[m] + 1
                break
        
    new_particles = Particles();
    index = 0
    
    print('samsamsamsamsamsam')
    print(sam)
    
    for i in range(100):
        
        while sam[i]>0:
            sam[i] = sam[i]-1
            x = particles.data[i][0]
            y = particles.data[i][1]
            theta = particles.data[i][2]
            w = 0.01
            
            new_particles.set_one_data(x,y,theta,w)
            index = index+1
  
    
    mx = new_particles.mean_x()
    my = new_particles.mean_y()
    mtheta = new_particles.mean_theta()
        
    lp = [mx,my,mtheta]
    
    #print('half way')
    #**********************************************
    
    print('resampling finished')
    print(new_particles.data)
    #print(lp)
    #print('------------------------')
            
    return new_particles.data,lp
    
def rotae_particles(particles,alpha,beta):
    print(alpha,beta)
    if alpha>= 0:
        leftTurndeg(beta*180/(np.pi))
    else:
        rightTurndeg(-beta*180/(np.pi))
        
    for i in range(100):
        w = particles[i][3]
        particles[i] = rotation(particles[i],alpha,beta,w)
    
    
    return particles


def straight_particles(particles,dist):
    
    forwardcm(dist)
    for i in range(100):
        w = particles[i][3]
        particles[i] = straightLine(particles[i],dist,w)
        
    return particles
        
        
   

def cw3(_id,_theta):
    
    
    interface.sensorEnable(port,brickpi.SensorType.SENSOR_ULTRASONIC)
    canvas = Canvas()
    mymap = init_map()
    
    particles = Particles(); 
    
    
    #particles.update();
    
    particles.draw()
    
    
    XX = [84,180,180,138,138,138]
    YY = [30,30,54,54,168,30]
    
    particles.start_point(XX[_id],YY[_id],_theta)
    lp = [XX[_id],YY[_id],_theta]
    for _ in range(0,6):
        
        
        #*******************************************************
        #navigateToWaypoint(XX[N],YY[N],lp)
        
        
        if _id ==5:
            _id=-1;
            
        X = XX[_id+1]
        Y = YY[_id+1]
        dx = X - lp[0]
        dy = Y - lp[1]
        print('@@@@@',_id+1,X,Y,dx,dy)
        alpha = math.atan2(dy,dx)
        print(alpha)
        print(lp[2])
        beta = alpha - lp[2]
        _id = _id + 1
        
        #beta = beta/1.05
        
        print('rotation angle should be:   ', beta)
        
        
        #print('alpha')
        #print(alpha)
        #print('beta')
        #print(beta)

        if beta > np.pi:
            beta = beta - np.pi * 2.
        if beta <= -np.pi:
            beta = beta + np.pi * 2.

        dist = np.sqrt(dx*dx+dy*dy)
        
        particles.data = rotae_particles(particles.data,alpha,beta)
        
        re = resample(particles)
        particles.data = re[0]
        lp = re[1]
        
        
        print('lp is')
        print(lp)
        particles.draw()

        print('go!__________________________________________________')

        while dist > 50:
            straight_particles(particles.data,50)
            particles.draw()#!!!
            
            re = resample(particles)
            particles.data = re[0]
            lp = re[1]
            print('lp is')
            print(lp)
            particles.draw()
            #*********************************************************
            dx = X - lp[0]
            dy = Y - lp[1]        
            alpha = math.atan2(dy,dx)
            beta = alpha - lp[2]
            
            #beta = beta/1.05

            if beta > np.pi:
                beta = beta - np.pi * 2.
            if beta <= -np.pi:
                beta = beta + np.pi * 2.
            
            print('alpha is')
            print(alpha, beta)
            dist = np.sqrt(dx*dx+dy*dy)
            #print('dist is')
            #print(dist)

            particles.data = rotae_particles(particles.data,alpha,beta)
            particles.draw()#!!!
            re = resample(particles)
            particles.data = re[0]
            lp = re[1]
            print('lp is')
            print(lp)
            particles.draw()
            #*********************************************************        
        
        #print('******************************************************')
        #print(dist)
        straight_particles(particles.data,dist)#!!!!
        particles.draw()#!!!
        re = resample(particles)
        particles.data = re[0]
        lp = re[1]
        print('lp is')
        print(lp)
        particles.draw()
        
        print('finish!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #***************************************************
       
        
        
        
#canvas = Canvas();    # global canvas we are going to draw on

mymap = init_map();
    
#print(mymap.walls)
    
#mymap = init_map()

#particles = Particles();

#cw3()

#navi()


from map import *
from numpy import genfromtxt



motor = 2

interface.motorEnable(motor)

k_left_i = 500
k_right_i = 500

motorParams2 = interface.MotorAngleControllerParameters()
motorParams2.maxRotationAcceleration = 2.0
motorParams2.maxRotationSpeed = 8.0
motorParams2.feedForwardGain = 255/20.0
motorParams2.minPWM = 15#17.882
motorParams2.pidParameters.minOutput = -255
motorParams2.pidParameters.maxOutput = 255
motorParams2.pidParameters.k_p = 100
motorParams2.pidParameters.k_i = k_right_i#400# k_right_i 350 10
motorParams2.pidParameters.k_d = 1

interface.setMotorAngleControllerParameters(motor,motorParams2)


def rotate(degree):
    angle = degree * (np.pi / 180)
    oldAngle = interface.getMotorAngles([motor])[0][0]
    print(oldAngle)
    
    interface.increaseMotorAngleReferences([motor],[-angle])
    old_angle = 5
    
    
    port = 0
    interface.sensorEnable(port, brickpi.SensorType.SENSOR_ULTRASONIC)
    
    reads = []
    a=0
    while not interface.motorAngleReferencesReached([motor]) :
        
        
        
        #reading = float(interface.getSensorValue(port)[0])
        #print(len(reads))

        motorAngle = interface.getMotorAngles([motor])[0][0]
        t = abs(motorAngle - oldAngle)*180./np.pi
        
        
        if t-5*a>=5.:
            a=a+1
            reading = float(interface.getSensorValue(port)[0])
            if reading >=160.:
                reading = 255.
            reads.append(reading)
        

            
        #if motorAngle[0][0]==old_angle:
         #   break        
        
        #old_angle = motorAngle[0][0]
        
        
    #for i in range(5):
    #    reading = float(interface.getSensorValue(port)[0])
    #    reads.append(reading)
    
    
    
    print "Destination reached!"
    return reads[:72]


def get_72reads():
    reads = []
    for i in range(72):
        read = rotate(5)
        reads.append(read)
    return reads
    

def saveSig(reads,i,name='sig'):
    with open(name+str(i)+".csv", "w") as outfile:
        for entries in reads:
            outfile.write(str(entries))
            outfile.write("\n")
            
def readSig():
    for i in range(1,6):   
        hist = []
        my_data = np.array(genfromtxt('sig'+str(i)+'.csv', delimiter=','))
        for j in range(85):
            t = np.where((my_data > 3*j) & (my_data <= 3*(j+1)))
            #print('t:',t)
            hist.append(len(t[0]))
        saveSig(hist,i,'hist')
    
    
def judgeSig():
    my_data = []
    for i in range(1,6):
        d = np.array(genfromtxt('hist'+str(i)+'.csv',delimiter=','))
        my_data.append(d)
    #freq
    cur_hist = []
    us_cur = np.array(rotate(370))
    #print('us:',us_cur)
    for j in range(85):
        t = np.where((us_cur > 3*j) & (us_cur <= 3*(j+1)))
        cur_hist.append(len(t[0]))
    saveSig(cur_hist,'9','test')
    cur_hist = np.array(cur_hist)
    print('cur_hist:',cur_hist)
    #find sig
    cos = []
    for d in my_data:
        t1 = np.sqrt(np.sum(d**2))
        t2 = np.sqrt(np.sum(cur_hist**2))
        cos.append(d.dot(cur_hist)/(t1*t2))
    cur_pos = np.argmax(cos)
    print('cos:',cos)
    print('sig position:',cur_pos+1)
    return cur_pos, us_cur
  

def judgeOrien():   
    cur_pos, us_cur = judgeSig()
    # rotate(-370)
    saveSig(us_cur,'10','test')
    true_sig = np.array(genfromtxt('sig'+str(cur_pos+1)+'.csv',delimiter=','))
    var = []
    a = list(us_cur)
    for i in range(len(us_cur)-1):
        b = a[1:]
        b.append(a[0])
        a = b
        #print(b)
        b = np.array(b)
        t1 = np.sqrt(np.sum(b**2,axis=0))
        t2 = np.sqrt(np.sum(true_sig**2,axis=0))
        var.append(b.dot(true_sig)/(t1*t2))
        #t = np.sum((b - true_sig)**2,axis=0)
        #var.append(t)
    ind = np.argmax(var)
    print('var:',var)
    print('orientation:',ind*5)
    
    ang = ind*5
    if ang>180:
        theta = (ang-360)*np.pi/180.
    else:
        theta = (ang)*np.pi/180.
        
    
    cw3(cur_pos,theta)
    
    return ind*5
    
        
        
judgeOrien()  

#judgeSig()
#readSig()    
#reads = rotate(370)
#reads = get_72reads()
#saveSig(reads,5)
#time.sleep(0.2)
#rotate(-90)
#time.sleep(0.2)
#rotate(45)
#cw2()
    
#us()
#q5_4()

#forwardcm(200)
#left90deg()
#right90deg()

#backwardcm(100)
#time.sleep(1)
#forward40()
#while 1:
        
#left90deg()
#time.sleep(1)
#left90deg()
#time.sleep(1)
#forwardcm(40)
interface.terminate()
