import math
import matplotlib.pyplot as plt
import numpy as np

def keplers(M,v1,v2,r1,r2):
    dt=0.05
    plt.figure(figsize=(7,7))
    G=1

    t=0
    #print(v1,v2,r1,r2)
    r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
    r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])
        
    a1=np.array([-G*M/(r1mag**3)*r1[i] for i in range(2)])
    a2=np.array([-G*M/(r2mag**3)*r2[i] for i in range(2)])

    r1all=np.array([])
    r2all=np.array([])

    r1start=np.copy(r1)
    r2start=np.copy(r2)

    maxMag1=0
    maxMag2=0
    minMag1=1E10
    minMag2=1E10

    t1=0
    t2=0
    
    while(True):
        t+=dt
        
        r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
        r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])

        maxMag1=max(maxMag1,r1mag)
        maxMag2=max(maxMag2,r2mag)
        minMag1=min(minMag1,r1mag)
        minMag2=min(minMag2,r2mag)

        a1=np.array([-G*M/(r1mag**3)*r1[i] for i in range(2)])
        a2=np.array([-G*M/(r2mag**3)*r2[i] for i in range(2)])
        
        prevA1=np.copy(a1)
        prevA2=np.copy(a2)

        r1=np.array([r1[i]+v1[i]*dt+1/2*a1[i]*dt*dt for i in range(2)])
        r2=np.array([r2[i]+v2[i]*dt+1/2*a2[i]*dt*dt for i in range(2)])
        
        r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
        r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])

        a1=np.array([-G*M/(r1mag**3)*r1[i] for i in range(2)])
        a2=np.array([-G*M/(r2mag**3)*r2[i] for i in range(2)])
        
        v1=np.array([v1[i]+1/2*(a1[i]+prevA1[i])*dt for i in range(2)])
        v2=np.array([v2[i]+1/2*(a2[i]+prevA2[i])*dt for i in range(2)])
        #print(r1,r2)
        

        if abs(r1[0]-r1start[0])<0.5 and abs(r1[1]-r1start[1])<0.5 and t1==0 and t>0.3: t1=t
        if abs(r2[0]-r2start[0])<0.5 and abs(r2[1]-r2start[1])<0.5 and t2==0 and t>0.3: t2=t

        if t1!=0 and t2!=0: break

        #print(r1)
        #np.append(r1all,np.copy(r1))
        #np.append(r2all,np.copy(r2))

        #plt.cla()
        plt.plot(r1[0],r1[1],'bo')
        plt.plot(r2[0],r2[1],'ro')
        plt.plot(0,0,'yo')
        #print(r1all,r2all)
        #plt.plot(r1all[0,:],r1all[1,:])
        #plt.plot(r2all[0,:],r2all[1,:])
        plt.ylim(-30,30)
        plt.xlim(-30,30)

        #if t>20: break

        plt.xlabel('x')
        plt.ylabel('y')
        #plt.text(0.5,2.5,info)
        
        plt.pause(0.0001)

    print("C1:",(t1**2)/(1/2*(maxMag1+minMag1))**3)
    print("C2:",(t2**2)/(1/2*(maxMag2+minMag2))**3)

    plt.show()

M=1000
v1=np.array([10,0])
v2=np.array([0,10])
r1=np.array([8,7])
r2=np.array([-10,9])
keplers(M,v1,v2,r1,r2)
