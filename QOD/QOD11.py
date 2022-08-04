import cv2
from tkinter import *
import turtle
import numpy as np

def bad():

    img=cv2.imread("turtle.jpg", 2)
    ret,nimg=cv2.threshold(img,127,160,cv2.THRESH_BINARY)
    h,w=img.shape[0],img.shape[1]

    scn=turtle.Screen()
    scn.screensize(w,h)
    t=turtle.Turtle()
    scn.tracer(0)

    for i in range(int(h/2),-int(h/2),-1):
        t.penup()
        t.goto(-int(w/2),i)
        
        for j in range(-int(w/2),int(w/2),1):
            ii=int(h/2-i-1)
            jj=int(w/2-j-1)
            if nimg[ii][jj]==0:
                t.pendown()
                t.forward(1)
            else:
                t.penup()
                t.forward(1)

        t.hideturtle()
        scn.update()
            

    res=turtle.getscreen()
    res.getcanvas().postscript(file="QOD11Result.png")
    turtle.done()

def good():
##    img=cv2.imread("turtle3.jpg", 2)
##    ret,nimg=cv2.threshold(img,127,160,cv2.THRESH_BINARY)
##    h,w=img.shape[0],img.shape[1]
##
    scn=turtle.Screen()
    scn.screensize(w,h)
    t=turtle.Turtle()
    scn.tracer(0)
##
##    for i in range(int(h/2),-int(h/2),-1):
##        t.penup()
##        t.goto(-int(w/2),i)
##        
##        for j in range(-int(w/2),int(w/2),1):
##            ii=int(h/2-i-1)
##            jj=int(w/2-j-1)
##            if nimg[ii][jj]==0:
##                t.pendown()
##                t.forward(1)
##            else:
##                t.penup()
##                t.forward(1)
##
##        t.hideturtle()
##        scn.update()
    
    counter = 10
    while(counter>0):
        turtle.listen()
        turtle.onclick(saveCoord)
        counter-=1

    f=open("qodcoords.txt",'a')
    for i in range(len(info)):
        f.write(info[0],info[1])
    f.close()

    res=turtle.getscreen()


def generate():
    scn=turtle.Screen()
    scn.screensize(100,100)
    scn.bgpic('turtle3.gif')
    t=turtle.Turtle()
    t.penup()
    t.goto(200,0)
    t.pendown()
    

    info=[]
    f=open("qodcoords.txt",'a')

    def saveCoord(x, y):
        t.pendown()
        t.goto(x,y)
        f.write(str(x)+" "+str(y)+"\n")

    scn.onclick(saveCoord)
    scn.mainloop()
    f.close()

#generate()

def draw():
    scn=turtle.Screen()
    scn.screensize(100,100)
    t=turtle.Turtle()
    t.penup()
    
    f=np.loadtxt("qodcoords.txt",delimiter=" ")
    t.goto(f[0][0],f[0][1])
    t.pendown()
    for pos in f:
        t.goto(pos[0],pos[1])

    t.hideturtle()
        

draw()

