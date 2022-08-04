import math
def formatInfo():
    f=open("SoongInput.txt",'r')
    info=f.readlines()
    res=[]
    print(len(info))
    for i in range(len(info)//4):
        rline=i*4+1
        vline=i*4+2
        l=info[rline].split("=")
        pos=[float(l[1][:-2]),float(l[2][:-2]),float(l[3])]
        l=info[vline].split("=")
        vel=[float(l[1][:-3]),float(l[2][:-3]),float(l[3])]
        res.append([pos,vel])  
    
    f.close()
    # return [[[x,y,z],[vx,vy,vz]]...]
    return res

def cross(v1:list, v2:list)->list:
    '''Returns the cross product of two 3D vectors'''
    return [(v1[1]*v2[2] - v1[2]*v2[1]),-(v1[0]*v2[2] - v1[2]*v2[0]),(v1[0]*v2[1] - v1[1]*v2[0])]

def angMoment(val):
    info=formatInfo()
    try:
        pos=info[val][0]
        vel=info[val][1]
        print(pos,vel)
        res=cross(pos,vel)
        return [round(res[0]/(2*math.pi)*365.2568983,6), round(res[1]/(2*math.pi)*365.2568983,6), round(res[2]/(2*math.pi)*365.2568983,6)]
    except: raise Exception ("Out of range")
    
    
print(angMoment(24))
