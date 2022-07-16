import numpy as np
import math
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

# angles/radians (declination and right ascension)

def decToRad(dec:float)->float:
    """ Given decimal value return in radians
        Args:
            dec (float): decimal value
        Returns:
            float: value in radians
    """
    return float(dec)/180*math.pi

def HMStoDeg(h:float,m:float,s:float,convert=False)->float:
    """ Given HMS, return in degrees or radians
        Args:
            h (float): hours
            m (float): minutes
            s (float): seconds
            convert (bool): True to convert to radians
        Returns:
            float: value in degrees or radians
    """
    try: h,m,s=float(h),float(m),float(s)
    except: raise Exception("Values must be floats")
    res = h*15 + m/60*15 + s/3600*15
    res%=360
    return np.deg2rad(res) if convert else res

def DMStoDeg(d:float,arcm:float,arcs:float,convert=False)->float:
    """ Given DMS, return in degrees or radians
        Args:
            d (float): degrees
            m (float): minutes
            s (float): seconds
            convert (bool): True to convert to radians
        Returns:
            float: value in degrees or radians
    """
    return np.deg2rad(math.copysign(abs(d)+arcm/60+arcs/3600,d)) if convert else math.copysign(abs(d)+arcm/60+arcs/3600,d)

def RAdecimalToHMS(dec:float)->tuple:
    """ Converts from decimal to HMS (right ascension)
        Args:
            dec (float): decimal value
        Returns:
            tuple: (hours, minutes, seconds)
    """
    dec%=360
    h=dec/15
    m=(dec%15)/15*60
    s=m%1*60
    h//=1
    m//=1
    return (h,m,s)

def DECdecimalToDMS(dec:float)->tuple:
    """ Converts from decimal to DMS (declination)
        Args:
            dec (float): decimal value
        Returns:
            tuple: (degrees, minutes, seconds)
    """
    d=math.trunc(dec)
    dec=abs(dec)-abs(d)
    m=math.trunc(dec*60)
    dec-=m/60
    s=dec*3600
    return (d,m,s)

def getAngle(sin:float,cos:float)->float:
    """ Returns the angle (in radians) in the correct quadrant given sin and cos
        Args:
            sin (float): sin value
            cos (float): cos value
        Returns:
            float: resulting angle in radians
    """
    return math.atan2(sin,cos) % (math.pi*2)

# vector stuff

def cross(v1:list, v2:list)->list:
    """ Returns the cross product of two 3D vectors
        Args:
            v1 (list): vector1
            v2 (list): vector2
        Returns:
            list: resulting cross product
    """
    return [(v1[1]*v2[2] - v1[2]*v2[1]),-(v1[0]*v2[2] - v1[2]*v2[0]),(v1[0]*v2[1] - v1[1]*v2[0])]

def dot(v1:list,v2:list)->float:
    """ Returns the dot product of two vectors
        Args:
            v1 (list): vector1
            v2 (list): vector2
        Returns:
            float: resulting dot product
    """
    return sum([v1[i]*v2[i] for i in range(len(v1))])

def det2x2(m)->float:
    """ Returns the determinant of a 2x2 matrix
        Args:
            m (list): matrix
        Returns:
            float: resulting determinant
    """
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]

def getMag(vec:list)->float:
    """ Returns the magnitude given a vector
        Args:
            vec (list): vector
        Returns:
            float: the magnitude of the vector
    """
    return math.sqrt(sum([vec[i]**2 for i in range(len(vec))]))
    
def rotZX(v:list,alpha:float,beta:float)->list:
    """ Rotates a vector around z axis with alpha, then rotates around x axis with beta'''
        Args:
            v (list): original vector
            alpha (float): angle to rotate around z axis
            beta (float): angle to rotate around x axis
        Returns:
            list: rotated vector
    """
    z=np.array([[np.cos(alpha),-np.sin(alpha),0],
                [np.sin(alpha),np.cos(alpha),0],
                [0,0,1]])
    x=np.array([[1,0,0],
                [0,np.cos(beta),-np.sin(beta)],
                [0,np.sin(beta),np.cos(beta)]])
    return np.dot(np.matmul(x,z),v)

# assorted other functions

def newton(func,der,init:float,err:float)->float: #function, derivative, inital_guess, error_tolerance
    """ Performs the Newton-Rhapson method
        Args:
            func (function): function
            der (function): derivative of function
            init (float): initial guess
            err (float): tolerance
        Returns:
            float: the answer
    """
    prev=-1e10
    while abs(init-prev)>=err:
        prev,init=init,init-func(init)/der(init)
    return init

def error(acc:float,exp:float)->float:
    """ Returns the error given accepted and experimental values
        Args:
            acc (float): accepted value
            exp (float): experimental value
        Returns:
            float: the error in percent
    """
    return math.copysign(abs(acc-exp)/acc*100,1)


# OD code
class ODElements:
    '''Class that represents and calculates orbital elements'''
    
    def __init__(self, pos:list, vel:list, time:float):
        """ Initializes ODElements class
            Args:
                pos (list): position of asteroid
                vel (list): velocity of asteroid (non Gaussian)
                time (float): given time in julian days
            Returns:
                None
        """
        # constants
        self.k = 0.0172020989484 # au^(3/2)/day
        
        self.time=time
        self.vel=vel
        self.vel=np.array([self.vel[i]/(2*math.pi)*365.2568983 for i in range(3)]) # convert to Gaussian
        self.pos=pos
        self.angMoment=self.getAngMoment()
        
        self.mu=1
        self.a=self.getSemiMajor() # semi-major axis
        self.e= self.getEcc() # eccentricity
        self.i=self.getInc() # inclination
        self.o=self.getLongAsc() # longitude of ascending node
        self.v=self.getTrueAnom() # true anomaly
        self.w=self.getArgPer() # argument of perihelion
        self.M=self.getMeanAnomaly() # mean anomaly
        self.T=self.getPeriT() # time of perihelion passage T
      
    def getInfo(self)->list:
        """ Returns info from given day from initialization
            Args:
                None
            Returns:
                list: all info
        """
        return self.info
    
    def getPos(self)->list:
        """ Returns the position of the asteroid
            Args:
                None
            Returns:
                list: the position
        """
        return self.pos
    
    def getVel(self)->list:
        """ Returns the velocity of the asteroid
            Args:
                None
            Returns:
                list: the velocity
        """
        return self.vel
    
    def getPosMag(self)->float:
        """ Returns the magnitude of the position vector
            Args:
                None
            Returns:
                float: magnitude of position
        """
        return math.sqrt(sum([self.pos[i]**2 for i in range(3)]))
    
    def getVelMag(self)->float:
        """ Returns the magnitude of the velocity vector
            Args:
                None
            Returns:
                float: magnitude of velocity
        """
        return math.sqrt(sum([self.vel[i]**2 for i in range(3)]))
    
    def getAngMoment(self)->list:
        """ Calculates and returns the specific angular momentum
            Args:
                None
            Returns:
                list: specific angular momentum components
        """
        pos,vel=self.getPos(),self.getVel()
        res=cross(pos,vel)
        return np.array([res[0], res[1], res[2]])
    
    def getAngMomentMag(self)->float:
        """ Returns the magnitude of the specific angular momentum
            Args:
                None
            Returns:
                float: the magnitude of specific angular momentum
        """
        return math.sqrt(sum([self.angMoment[i]**2 for i in range(3)]))
    
    def getSemiMajor(self)->float:
        """ Calculates and returns the semi major axis using vis-viva
            Args:
                None
            Returns:
                float: the semi major axis
        """
        return 1/(2/self.getPosMag() - dot(self.vel,self.vel)/self.mu)
    
    def getEcc(self)->float:
        """ Calculates and returns the eccentricity
            Args:
                None
            Returns:
                float: eccentricity
        """
        return math.sqrt(1-dot(self.angMoment, self.angMoment)/(self.mu*self.a))
    
    def getInc(self, rad:bool=False)->float:
        """ Calculates and returns the inclination in degrees
            Args:
                rad (bool): True if return in radians
            Returns:
                float: the inclination in degrees or radians
        """
        return math.acos(self.angMoment[2]/self.getAngMomentMag()) if rad else np.rad2deg(math.acos(self.angMoment[2]/self.getAngMomentMag()))
    
    def getLongAsc(self, rad:bool=False):
        """ Calculates and returns the longitude of ascending node in degrees
            Args:
                rad (bool): True if return in radians
            Returns:
                float: the longitude of ascending node in degrees or radians
        """
        s=self.angMoment[0]/(self.getAngMomentMag()*math.sin(np.deg2rad(self.i)))
        c=-self.angMoment[1]/(self.getAngMomentMag()*math.sin(np.deg2rad(self.i)))
        return getAngle(s,c) if rad else np.rad2deg(getAngle(s,c))
    
    def getArgPer(self, rad:bool=False)->float:
        """ Calculates and returns the argument of perihelion in degrees
            Args:
                rad (bool): True if return in radians
            Returns:
                float: the longitude of ascending node in degrees or radians
        """
        x,y,z=self.pos[0],self.pos[1],self.pos[2]
        s=z/(self.getPosMag()*math.sin(np.deg2rad(self.i)))
        c=(x*math.cos(np.deg2rad(self.o)) + y*math.sin(np.deg2rad(self.o)))/self.getPosMag()
        U=np.rad2deg(getAngle(s,c))
        return np.deg2rad((U-self.v)%360) if rad else (U-self.v)%360
    
    def getTrueAnom(self, rad:bool=False)->float:
        """ Calculates and returns the true anomaly
            Args:
                rad (bool): True if return in radians
            Returns:
                float: the true anomaly in degrees or radians
        """
        c=1/self.e*(self.a*(1-self.e**2)/self.getPosMag() - 1)
        s=self.a*(1-self.e**2)/(self.getAngMomentMag() * self.e)*dot(self.pos,self.vel)/self.getPosMag()
        return getAngle(s,c) if rad else np.rad2deg(getAngle(s,c))
    
    def getPeriT(self)->float:
        """ Calculates and returns the time of perihelion
            Args:
                None
            Returns:
                float: the time of perihelion
        """
        n=self.k/(self.a**(3/2))
        return self.time-np.deg2rad(self.M)/n
    
    def getMeanAnomaly(self,rad:bool=False)->float:
        """ Calculates and returns the mean anomaly
            Args:
                rad (bool): True if return in radians
            Returns:
                float: the mean anomaly in degrees or radians
        """
        s=(self.getPosMag()*math.sin(np.deg2rad(self.v)))/(self.a*math.sqrt(1-self.e**2))
        c=(self.e+math.cos(np.deg2rad(self.v)))/(1+self.e*math.cos(np.deg2rad(self.v)))
        E=getAngle(s,c)
        return E-self.e*math.sin(E) if rad else np.rad2deg(E-self.e*math.sin(E))
    
    def printError(self, results:list):
        """ Prints everything
            Args:
                None
            Returns:
                None
        """
        # EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR,
        print("Semi-major axis:", results[9], self.a, error(results[9],self.a))
        print("Eccentricity:", results[0], self.e, error(results[0],self.e))
        print("Inclination:",results[2],self.i, error(results[2],self.i))
        print("Longitude of Ascending Node:",results[3],self.o, error(results[3],self.o))
        print("True anomaly:",results[8],self.v,error(results[8],self.v))
        print("Argument of perihelion:",results[4],self.w,error(results[4],self.w))
        print("Time of Perihelion Passage T:",results[5],self.T,error(results[5],self.T))
        print("Mean Anomaly:",results[7],self.M,error(results[7],self.M))
        #print(od.a, od.e, od.i, od.o, od.v, od.w)
        
    def getElements(self):
        """ Returns all orbital elements
            Args:
                rad (bool): True if return in radians
            Returns:
                floats: a,e,i,o,v,w,T,M
        """
        return self.a, self.e, self.i, self.o, self.v, self.w, self.T, self.M
    
    

