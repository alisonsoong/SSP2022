import math
import matplotlib.pyplot as plt
import numpy as np

def threeBodySim(m,v1,v2,v3,r1,r2,r3):
    dt=0.05
    plt.figure(figsize=(7,7))
    G=1.6

    m1,m2,m3=m[0],m[1],m[2]

    t=0

    # magnitudes/distances
    r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
    r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])
    r3mag=math.sqrt(r3[0]*r3[0]+r3[1]*r3[1])

    dist12=[r2[0]-r1[0],r2[1]-r1[1]]
    dist21=[r1[0]-r2[0],r1[1]-r2[1]]
    magDist12=math.sqrt(dist12[0]**2+dist12[1]**2)

    dist13=[r3[0]-r1[0],r3[1]-r1[1]]
    dist31=[r1[0]-r3[0],r1[1]-r3[1]]
    magDist13=math.sqrt(dist13[0]**2+dist13[1]**2)

    dist23=[r3[0]-r2[0],r3[1]-r2[1]]
    dist32=[r2[0]-r3[0],r2[1]-r3[1]]
    magDist23=math.sqrt(dist23[0]**2+dist23[1]**2)

    # update accelerations
    a1=np.array([G*m2/(magDist12**3)*dist12[i] + G*m3/(magDist13**3)*dist13[i] for i in range(2)])
    a2=np.array([G*m1/(magDist12**3)*dist21[i] + G*m3/(magDist23**3)*dist23[i] for i in range(2)])
    a3=np.array([G*m1/(magDist13**3)*dist31[i] + G*m2/(magDist23**3)*dist32[i] for i in range(2)])
    
    r1start=np.copy(r1)
    r2start=np.copy(r2)
    r3start=np.copy(r3)
    
    while(True):
        t+=dt
        
        # magnitudes/distances
        r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
        r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])
        r3mag=math.sqrt(r3[0]*r3[0]+r3[1]*r3[1])

        dist12=[r2[0]-r1[0],r2[1]-r1[1]]
        dist21=[r1[0]-r2[0],r1[1]-r2[1]]
        magDist12=math.sqrt(dist12[0]**2+dist12[1]**2)

        dist13=[r3[0]-r1[0],r3[1]-r1[1]]
        dist31=[r1[0]-r3[0],r1[1]-r3[1]]
        magDist13=math.sqrt(dist13[0]**2+dist13[1]**2)

        dist23=[r3[0]-r2[0],r3[1]-r2[1]]
        dist32=[r2[0]-r3[0],r2[1]-r3[1]]
        magDist23=math.sqrt(dist23[0]**2+dist23[1]**2)

        # update accelerations
        a1=np.array([G*m2/(magDist12**3)*dist12[i] + G*m3/(magDist13**3)*dist13[i] for i in range(2)])
        a2=np.array([G*m1/(magDist12**3)*dist21[i] + G*m3/(magDist23**3)*dist23[i] for i in range(2)])
        a3=np.array([G*m1/(magDist13**3)*dist31[i] + G*m2/(magDist23**3)*dist32[i] for i in range(2)])
    
        prevA1=np.copy(a1)
        prevA2=np.copy(a2)
        prevA3=np.copy(a3)

        r1=np.array([r1[i]+v1[i]*dt+1/2*a1[i]*dt*dt for i in range(2)])
        r2=np.array([r2[i]+v2[i]*dt+1/2*a2[i]*dt*dt for i in range(2)])
        r3=np.array([r3[i]+v3[i]*dt+1/2*a3[i]*dt*dt for i in range(2)])
        
        # magnitudes/distances
        r1mag=math.sqrt(r1[0]*r1[0]+r1[1]*r1[1])
        r2mag=math.sqrt(r2[0]*r2[0]+r2[1]*r2[1])
        r3mag=math.sqrt(r3[0]*r3[0]+r3[1]*r3[1])

        dist12=[r2[0]-r1[0],r2[1]-r1[1]]
        dist21=[r1[0]-r2[0],r1[1]-r2[1]]
        magDist12=math.sqrt(dist12[0]**2+dist12[1]**2)

        dist13=[r3[0]-r1[0],r3[1]-r1[1]]
        dist31=[r1[0]-r3[0],r1[1]-r3[1]]
        magDist13=math.sqrt(dist13[0]**2+dist13[1]**2)

        dist23=[r3[0]-r2[0],r3[1]-r2[1]]
        dist32=[r2[0]-r3[0],r2[1]-r3[1]]
        magDist23=math.sqrt(dist23[0]**2+dist23[1]**2)

        # update accelerations
        a1=np.array([G*m2/(magDist12**3)*dist12[i] + G*m3/(magDist13**3)*dist13[i] for i in range(2)])
        a2=np.array([G*m1/(magDist12**3)*dist21[i] + G*m3/(magDist23**3)*dist23[i] for i in range(2)])
        a3=np.array([G*m1/(magDist13**3)*dist31[i] + G*m2/(magDist23**3)*dist32[i] for i in range(2)])
    
        v1=np.array([v1[i]+1/2*(a1[i]+prevA1[i])*dt for i in range(2)])
        v2=np.array([v2[i]+1/2*(a2[i]+prevA2[i])*dt for i in range(2)])
        v3=np.array([v3[i]+1/2*(a3[i]+prevA3[i])*dt for i in range(2)])


        plt.plot(r1[0],r1[1],'bo',markersize=2)
        plt.plot(r2[0],r2[1],'ro',markersize=2)
        plt.plot(r3[0],r3[1],'yo',markersize=2)
        plt.ylim(-50,50)
        plt.xlim(-50,50)

        #print ke + pe
        ke1=1/2*m1*(v1[0]**2+v1[1]**2)
        ke2=1/2*m2*(v2[0]**2+v2[1]**2)
        ke3=1/2*m3*(v3[0]**2+v3[1]**2)
        pe1=(G*m1*m2)/(magDist12)+(G*m1*m3)/(magDist13)
        pe2=(G*m2*m1)/(magDist12)+(G*m2*m3)/(magDist23)
        pe3=(G*m1*m3)/(magDist13)+(G*m2*m3)/(magDist23)

        print("total energy:", ke1+ke2+ke3+pe1+pe2+pe3)
    
        


        plt.xlabel('x')
        plt.ylabel('y')

        
        plt.pause(0.0001)

    plt.show()

M=[1000,10,0.1]
v1=np.array([0,0])
v2=np.array([0,8.94])
v3=np.array([0.1,12.94])
r1=np.array([0,0])
r2=np.array([20,0])
r3=np.array([19,0])
threeBodySim(M,v1,v2,v3,r1,r2,r3)
