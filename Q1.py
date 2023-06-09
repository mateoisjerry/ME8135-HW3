# Matthew Lisondra
# ME8135 Assignment 3 Q1

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

N = 100
particles = []
particle_weights = np.ones(N)/N # should add to 1
particle_distances = np.ones(N)/N # this is set to not allow 1/0. also, weights will be 1/particle_distances and then normalized using some python function
Zx = np.array([0,0])
PFx = np.array([0,0])

T = 1/8 # 1/8 seconds
r = 0.1 # in m
At = np.identity(2)
Bt = np.zeros((2,2)) + (T*r/2)
ut = np.array([1.,0.1])
Ct = np.array([[1.,0.],[0.,2.]]) # Has to be made float
rsigmax2 = 0.05
rsigmay2 = 0.075

def GT(x):
    omegax = np.random.normal(0,0.1)
    omegay = np.random.normal(0,0.15)
    epsilonvect = np.array([[omegax*T],[omegay*T]])
    
    return np.matmul(At, x) + np.matmul(Bt, ut) + epsilonvect

def zt():
    global Zx
    global PFx
    
    rx = np.random.normal(0, rsigmax2)
    ry = np.random.normal(0, rsigmay2)
    delt = np.array([rx,ry])
    Zxnew = Ct.dot(PFx) + delt
    Zx = Zxnew
    
    # print(Zxnew)

def update_weights():
    global particles
    global particle_weights
    global particle_distances
    global Zx
    
    if np.size(Zx) != 2:
        Zx = Zx[0]
    
    for i in range(N):
        particle_distances[i] = math.dist([particles[i][0][0],particles[i][1][0]], Zx)
        particle_weights[i] = 1/particle_distances[i]
        
    vector_norm = np.linalg.norm(particle_weights)
    particle_weightsnew = np.copy(particle_weights)/vector_norm
    correction = np.sum(particle_weightsnew) # this is to make the particle_weights a probability that sums to 1
    particle_weights = particle_weightsnew / correction

    # print(particle_weights)
    # print(sum(particle_weights))

def pickPF():
    global particles
    global particle_weights
    global PFx
    
    # print(PFx)
    
    randomchoice = random.choices(particles, weights=particle_weights, k=1)
    PFxnew = [randomchoice[0][0][0],randomchoice[0][1][0]]
    PFx = PFxnew
    
    # print(PFx)
    # print(PFxnew)
    
# -----------------------------------------------------------------------------------------------------

for i in range(N):
    particles.append(GT(PFx))

pygame.init()
display_width = 1280
display_height = 720
gameDisplay = pygame.display.set_mode((display_width,display_height))
pygame.display.set_caption('Matthew Lisondra ME8135 Assignment 3 Q1')

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
        pygame.draw.circle(gameDisplay, RED, ([particles[i][0][0]*1000,particles[i][1][0]*1000]),3)
        # try:
        #     pygame.draw.circle(gameDisplay, RED, ([particles[i][0][0][0]*1000,particles[i][0][1][0]*1000]),3)
        # except:
        #     pygame.draw.circle(gameDisplay, RED, ([particles[i][0][0]*1000,particles[i][1][0]*1000]),3)
        
    if (t%8)==0:
        zt()
        update_weights()
        pickPF()
        
        pygame.draw.circle(gameDisplay, PINK, ([Zx[0]*1000,Zx[1]*1000]),10)
        
        PFpoints.append([PFx[0]*1000+50,PFx[1]*1000+50])

    # PF DRAWING
    pygame.draw.polygon(gameDisplay, GREEN,
                        [[PFx[0]*1000+50,PFx[1]*1000+50],[PFx[0]*1000+40,PFx[1]*1000+35] ,
                        [PFx[0]*1000+40,PFx[1]*1000+65]])
    PFpoints.append([PFx[0]*1000+50,PFx[1]*1000+50])
    pygame.draw.lines(gameDisplay,GREEN,False,PFpoints,5)

    # for i in range(N):
    #     temp_particle = np.copy([particles[i][0][0],particles[i][1][0]])
    #     particles[i] = GT(temp_particle)
    
    pygame.display.update()
    clock.tick(8) 
    t += 1

