# Matthew Lisondra
# ME8135 Assignment 3 Q2b

# Import the necessary modules
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.stats import multivariate_normal
import math
import pygame
from pygame.locals import *
import random

# -----------------------------------------------------------------------------------------------------

def convertToPolar(xytheta):
    currentXYTHETA = np.copy([xytheta[0]*1000,xytheta[1]*1000,xytheta[2]])
    
    rho = math.dist(currentXYTHETA[:2], CENTER)
    theta = np.arctan2((currentXYTHETA[1]-CENTER[1]),(currentXYTHETA[0]-CENTER[0])) # should have removed from beginning
    
    polarx = np.array([rho, theta])
    
    return polarx

N = 100
particles = []
particle_weights = np.ones(N)/N # should add to 1
particle_distances = np.ones(N)/N # this is set to not allow 1/0. also, weights will be 1/particle_distances and then normalized using some python function
Zx = np.array([0.,0.,0.])
d=10. # OFFSET FROM CENTER STARTING 270 DEG COUNTER CLOCKWISE
CENTER = np.array([640.,360.])
PFx = np.array([CENTER[0]/1000., CENTER[1]/1000. + (d*10/1000. + d*5/1000.), 0.])

T = 1/8 # 1/8 seconds
r = 0.1 # in m
L = 0.3 #  L is the distance between the wheel, known as wheelbase, and is 0.3m.
Gt = np.identity(3)
Ht = np.array([[1.,0.,0.],[0.,2.,0.],[0.,0.,1.]]) # Has to be made float
rsigmar2 = 0.1
rsigmab2 = 0.01
# transpose_Gt = np.transpose(Gt)
# transpose_Ht = np.transpose(Ht)

# wsigmax2 = 0.01
# wsigmay2 = 0.1
# Rt = np.array([[wsigmax2,0.,0.],[0.,wsigmax2,0.],[0.,0.,wsigmay2]])*T # no need to square again

# Qt = np.array([[rsigmax2,0,0],[0,rsigmay2,0],[0,0,0]]) # no need to square again

def GT(x):
    omegaw = np.random.normal(0,0.1) / 10. # scale
    omegaphi = np.random.normal(0,0.01) / 10. # scale
    epsilonvect = np.array([T*omegaw,T*omegaw,T*omegaphi])
    
    ur = 0.5
    ul = 0.5
    
    currentXY = np.copy([x[0]*1000,x[1]*1000])
    if math.dist(currentXY, CENTER) < d*10:
        ur = 1.
        ul = 0.1
    if math.dist(currentXY, CENTER) > (d+1)*10:
        ur = 0.1
        ul = 1.
        
    # print(math.dist(currentXY, CENTER))
    ut = np.array([(ur+ul)/2,(ur-ul)])
    
    THETAM = np.array([[T*r*np.cos(x[2]),0.],[T*r*np.sin(x[2]),0.],[0.,T*r/L]])
    
    return np.matmul(Gt, x) + np.matmul(THETAM, ut) + epsilonvect # noise is the epsilonvect

# def KF_prediction():
#     global Px
#     firstP = np.matmul(Gt, Px)
#     P_bar = np.matmul(firstP, transpose_Gt) + Rt
#     Px = P_bar
#     # print(Px)

# def zt():
#     global GTx
#     global Zx
#     rx = np.random.normal(0, rsigmax2)
#     ry = np.random.normal(0, rsigmay2)
#     delt = np.array([rx,ry,0.])
#     Zxnew = Ht.dot(GTx) + delt
#     Zx = Zxnew
#     # print(Zxnew)

def zt():
    global Zx
    global PFx
    global CENTER
    
    rr = np.random.normal(0, rsigmar2) # range measurement noise
    rb = np.random.normal(0, rsigmab2) # bearing measurement noise
    delt = np.array([rr,rb])
    Zxnew = convertToPolar(PFx) + delt
    Zx = Zxnew
    
    # print(Zx)
    # print([CENTER[0] + Zx[0]*np.cos(Zx[1]), CENTER[1] + Zx[0]*np.sin(Zx[1])])
    
# def zt():
#     global Zx
#     global PFx
    
#     rx = np.random.normal(0, rsigmax2)
#     ry = np.random.normal(0, rsigmax2)
#     delt = np.array([rx,ry,0.])
#     Zxnew = Ht.dot(PFx) + delt
#     Zx = Zxnew
    
#     print(Zx)
    
def update_weights():
    global particles
    global particle_weights
    global particle_distances
    global Zx
    global CENTER
    
    ZxINTOXY = [CENTER[0] + Zx[0]*np.cos(Zx[1]), CENTER[1] + Zx[0]*np.sin(Zx[1])]
    
    for i in range(N):
        particle_distances[i] = math.dist([particles[i][0],particles[i][1]], ZxINTOXY)
        particle_weights[i] = 1/particle_distances[i]
        
    vector_norm = np.linalg.norm(particle_weights)
    particle_weightsnew = np.copy(particle_weights)/vector_norm
    correction = np.sum(particle_weightsnew) # this is to make the particle_weights a probability that sums to 1
    particle_weights = particle_weightsnew / correction

    print(particle_weights)
    print(sum(particle_weights))

def pickPF():
    global particles
    global particle_weights
    global PFx
    
    # print(PFx)
    
    randomchoice = random.choices(particles, weights=particle_weights, k=1)
    PFxnew = [randomchoice[0][0],randomchoice[0][1], randomchoice[0][2]]
    PFx = PFxnew
    
    # print(PFx)
    
    # print(PFxnew)

# def KF_measurement():
#     global GTx
#     global Px
#     global Zx
#     global KFx
    
#     part1INV = np.matmul(Ht, Px)
#     part2INV = np.matmul(part1INV, transpose_Ht) + Qt
#     invPart = np.linalg.inv(part2INV)
    
#     part1Kt = np.matmul(Px, transpose_Ht)
#     Kt = np.matmul(part1Kt, invPart)
    
#     muPart = Zx - Ht.dot(GTx)
#     munew = GTx + np.matmul(Kt, muPart)
#     KFx = munew
#     # print(KFx)
    
#     PPart = np.identity(3) - np.matmul(Kt, Ht)
#     Pnew = np.matmul(PPart, Px)
#     Px = Pnew
#     # print(Px)

# -----------------------------------------------------------------------------------------------------

for i in range(N):
    particles.append(GT(PFx))

pygame.init()
display_width = 1280
display_height = 720
gameDisplay = pygame.display.set_mode((display_width,display_height))
pygame.display.set_caption('Matthew Lisondra ME8135 Assignment 3 Q2b')

black = (0,0,0)
white = (255,255,255)

clock = pygame.time.Clock()
crashed = False

PFpoints=[]

imp = pygame.image.load("INFO-TEXT.png").convert()

t = 0
while not crashed:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            crashed = True

    gameDisplay.fill(white)
    WHITE=(255,255,255)
    BLUE=(0,0,255)
    RED=(255,0,0)
    GREEN=(0,255,0)
    PINK = (255,192,203)
    
    gameDisplay.blit(imp, (600, 0))
    pygame.display.flip()

    for i in range(N):
        particles[i] = GT(PFx)
        pygame.draw.circle(gameDisplay, RED, ([particles[i][0]*1000,particles[i][1]*1000]),3)
    
    if (t%8)==0:
        zt()
        update_weights()
        pickPF()
        
        # print([Zx[0]*1000,Zx[1]*1000])
        ZxINTOXY2 = [CENTER[0] + Zx[0]*np.cos(Zx[1]), CENTER[1] + Zx[0]*np.sin(Zx[1])]
        print(Zx)
        print(ZxINTOXY2)
        pygame.draw.circle(gameDisplay, PINK, ([ZxINTOXY2[0],ZxINTOXY2[1]]),10) # scaled, Ht[2,2] row column is element = 2 which scales too large
        
        PFpoints.append([PFx[0]*1000,PFx[1]*1000])

    # PF DRAWING
    PFpoints.append([PFx[0]*1000,PFx[1]*1000])
    pygame.draw.lines(gameDisplay,GREEN,False,PFpoints,5)

    # for i in range(N):
    #     temp_particle = np.copy([particles[i][0][0],particles[i][1][0]])
    #     particles[i] = GT(temp_particle)
    
    pygame.draw.circle(gameDisplay, RED, (CENTER[0],CENTER[1]),5)
    
    pygame.display.update()
    clock.tick(8) 
    t += 1

