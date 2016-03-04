# Read the instructions file before implementing

import sys, pygame
from pygame.locals import *
from pygame.constants import *
from OpenGL.GL import *
from OpenGL.GLU import *
from numpy import *
# IMPORT OBJECT LOADER
from load_room_obj import *

pygame.init()
viewport = (800,600)
hx = viewport[0]/2
hy = viewport[1]/2
srf = pygame.display.set_mode(viewport, OPENGL | DOUBLEBUF)
glEnable(GL_DEPTH_TEST)
obj = OBJ(sys.argv[1], swapyz=True)
clock = pygame.time.Clock()
glMatrixMode(GL_PROJECTION)
glLoadIdentity()
width, height = viewport
gluPerspective(60.0, width/float(height), 1, 100.0)
glEnable(GL_DEPTH_TEST)
glMatrixMode(GL_MODELVIEW)
rx, ry = (0,0)
tx, ty = (0,0)
zpos = 0
rotate = move = False
entry = False
while 1:
    clock.tick(30)
    for e in pygame.event.get():
        if e.type == QUIT or e.type == K_ESCAPE:
            sys.exit()
        elif e.type == KEYDOWN:
            rx, ry = (0,0)
            mat = glGetDoublev(GL_MODELVIEW_MATRIX)
            print mat
            x = mat[3][0]
            z = mat[3][2]    
            if x <= 7 and x >= -5 and z <= 10 and z >= -2:
                #Code for the entry gate begins
                if z == -2:
                    if x >= -1 and x <= 2:
                        # print "."
                        if e.key == K_UP:
                            ty += 1
                            entry = True
                        if e.key == K_DOWN:
                            ty -= 1
                            entry = False
                    else:
                        # print ".."
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                #Code for the exit from gate begins
                elif z == -1:
                    if x >= -1 and x <= 2:
                        # print "..."
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                    else:
                        # print "...."
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1            
                elif z >= 0 and z <= 1:
                    if x == 1 or x == 5 or x == -5:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                    else:
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1                    
                elif z == 2:
                    if x == 1 or x == 5 or x == -5:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1        
                    else:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                elif z == 4:
                    if x >= -3 and x <= 3:
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                    else:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                elif z == 5 or z == 6:
                    if x == -4:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_RIGHT:
                            tx += 1
                    elif x == 4:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                    else:
                        if e.key == K_UP:
                            ty += 1
                        if e.key == K_DOWN:
                            ty -= 1
                        if e.key == K_LEFT:
                            tx -= 1
                        if e.key == K_RIGHT:
                            tx += 1                                
                else:
                    # print "....."
                    if e.key == K_UP:
                        ty += 1
                    if e.key == K_DOWN:
                        ty -= 1
                    if e.key == K_LEFT:
                        tx -= 1
                    if e.key == K_RIGHT:
                        tx += 1                                                                             
            else:
                if e.key == K_UP:
                    ty += 1
                if e.key == K_DOWN:
                    ty -= 1
                if e.key == K_LEFT:
                    tx -= 1
                if e.key == K_RIGHT:
                    tx += 1        
        elif e.type == MOUSEBUTTONDOWN:
            if e.button == 1: rotate = True
            # elif e.button == 3: move = True
        elif e.type == MOUSEBUTTONUP:
            rx, ry = (0,0)
            if e.button == 1: rotate = False
        elif e.type == MOUSEMOTION:
            i, j = e.rel
            if rotate:
                rx += i
                ry += j
                print rx
                print ry
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()
        # RENDER OBJECT at -3, 3, -9        
    if entry == True:
        if tx < -5:
            tx = -5
        elif tx > 5:
            tx = 5
        elif ty > 19:
            ty = 19
    elif entry == False:
        if ty >= 8 and ty <= 21:
            if tx == 7:
                # print ".."
                tx = 8
            elif tx == -8:
                tx = -9

                # print ""
        elif ty == 22 and tx <= 7 and tx >= -9:
            ty = 23
                    

    print tx
    print ty
    if rotate:
        print rx/10.0
        print ry/10.0

        gluLookAt(-tx, 10-ty, 0.0, -tx+rx/10,-100, -ry, 0.0, 0.0, 1.0)
    else:
        gluLookAt(0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)            
        glTranslate(tx, ty, - zpos)        
    # glRotate(ry, 1, 0, 0)
    # glRotate(rx, 0, 1, 0)
    glCallList(obj.gl_list)
    pygame.display.flip()