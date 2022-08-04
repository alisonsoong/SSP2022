import math
import matplotlib.pyplot as plt
import numpy as np


def noAir():
    delT=0.02
    times=np.arange(0,2,delT)
    v=8.0
    theta=50.0*math.pi/180
    g=9.8

    info="dt="+str(delT)+"\nv="+str(v)+"\ntheta="+str(theta)

    v0x=v*math.cos(theta)
    v0y=v*math.sin(theta)

    plt.figure(figsize=(6,6))
    cntr=0

    
    
    while(True):
        t=times[cntr]
        cntr+=1
        x=v0x*t
        y=v0y*t-0.5*g*t*t
        if y<0: break
        plt.plot(x,y,'bo')

        plt.ylim(0,3)
        plt.xlim(0,10)


        plt.xlabel('x')
        plt.ylabel('y')
        plt.text(0.5,2.5,info)
        
        plt.pause(delT)

    results="Results:\nTime:"+str(cntr*delT)+"\nRange:"+str(x)
    print(results)
    plt.text(0.5,2.05,results)
    plt.show()
        
noAir()

def airRes(v):
    delT=0.04
    dt=delT
    theta=35.0*math.pi/180#45.0*math.pi/180

    a=0.023#0.032
    c=0.2#0.04
    p=1.29
    b=c*p*a
    g=9.8
    m=0.42

    info="dt="+str(delT)+"\nv="+str(v)+"\ntheta="+str(theta)

    vx2=v*math.cos(theta)
    vy2=v*math.sin(theta)
    vx,vy=vx2,vy2
    v=math.sqrt(vx2*vx2+vy2*vy2)
    
    ax,ay=(-b*vx2*v/m),(-b*vy2*v/m-g)

    plt.figure(figsize=((6,6)))
    cnt=0
    res=0
    x2,y2=0,0
    t=0
    x,y=0,0
    flag=False
    
    while(True):
        t+=dt

        cnt+=1

        x=vx*t
        y=vy*t-0.5*g*t*t
        
        v=math.sqrt(vx2*vx2+vy2*vy2)
        ax,ay=(-b*vx2*v/m),(-b*vy2*v/m-g)
        vx2+=ax*dt
        vy2+=ay*dt
        x2+=vx2*dt
        y2+=vy2*dt
        
        if y2>=0: plt.plot(x2,y2,'ro')
        elif not flag:
            res=(t,x2)
            flag=True
        if y>=0: plt.plot(x,y,'bo')
        if y<0 and y2<0: break
        
        plt.ylim(0,60)
        plt.xlim(0,300)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.text(10,40,info)
        
        plt.pause(0.000001)

    results="Results:\nTime:"+str(res[0])+"\nRange:"+str(res[1])
    print(results)
    plt.text(10,35,results)
    plt.show()

airRes(55)#37.4)
