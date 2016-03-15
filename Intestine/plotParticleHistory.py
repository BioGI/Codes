import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc_file
from pylab import loadtxt
rc_file('rcFile')

vcf = 6e-2
tcf = 2e-3
xcf = 1.2e-4
zcf = 1.2e-4
L = 0.012

def plotParticleHistory():
    par1 = loadtxt('particle1-history.dat')

    xA,yA,zA = computeAnalyticalParticleHistory(par1[0,2],par1[0,3],par1[0,4],par1[0,1],par1[:,1])
    
    fig = plt.figure()
    ax = plt.axes([0.19,0.18,0.96-0.19,0.92-0.18])
    plt.plot(par1[:,1],par1[:,4],'r+-',label='Dual Lattice')
    plt.plot(par1[:,1],zA,label='Analyt. Soln.')
    plt.legend(loc=0)    
    plt.xlabel('Time (s)')
    plt.ylabel('z location (m)')
    plt.savefig('particle1_zLocation.png')

    plt.xlim(0,0.31)
    plt.savefig('particle1_zLocationZoom.png')
    plt.show()

    
    fig = plt.figure()
    ax = plt.axes([0.19,0.18,0.96-0.19,0.92-0.18])
    plt.plot(par1[:,1],par1[:,3],'r+-',label='Dual Lattice')
    plt.plot(par1[:,1],yA,label='Analyt. Soln.')
    plt.legend(loc=0)        
    plt.axhline(0.00072)
    plt.axhline(0.00069)
    plt.axhline(-0.00069)
    plt.axhline(-0.00072)
    plt.xlabel('Time (s)')
    plt.ylabel('y location (m)')
    plt.savefig('particle1_yLocation.png')

    plt.xlim(0,0.31)    
    plt.ylim(-0.0008,0.0008)
    plt.savefig('particle1_yLocationZoom.png')    
    plt.show()

def computeAnalyticalParticleHistory(x0,y0,z0,t0,t):
    """Computes the analytical solution to the particle trajectory at a discrete set of time points 't', given the initial location of the particle (x0,y0,z0) at time t0"""
    nt = np.size(t)
    x = np.empty(nt)
    y = np.empty(nt)
    z = np.empty(nt)

    x[0] = x0
    y[0] = y0
    z[0] = z0
    for i in range(1,nt):
        x[i],y[i],z[i] = computeOneTimeStep(x[i-1],y[i-1],z[i-1],t[i]-t[i-1])

    return np.array([x,y,z])

        
        
def computeOneTimeStep(x0,y0,z0,dt):
        """Computes one time step of the particle initially at location (x0,y0,z0) using second order Runge-Kutta time stepping over a time-step dt."""

        #Calculate intermediate time step
        r = np.sqrt(x0**2 + y0**2)
        v0 = np.array([0,v(r),w(r)])
        xi = x0 + v0[0]*dt
        yi = y0 + v0[1]*dt
        zi = z0 + v0[2]*dt
        r = np.sqrt(xi**2 + yi**2)
        vi = np.array([0,v(r),w(r)])
        
        #Calculate final loc of particle
        xn = x0 + 0.5*(v0[0]+vi[0])*dt
        yn = y0 + 0.5*(v0[1]+vi[1])*dt
        zn = z0 + 0.5*(v0[2]+vi[2])*dt

        if(zn > L):
            zn = zn - L

        return(np.array([xn,yn,zn]))
    

        
def w(r):
    """Returns the w velocity as a function of radial location"""
    return (1.0 - r/0.005)*vcf

def v(r):
    """Returns the v velocity as a function of radial location"""
    return -0.1*vcf

if __name__=="__main__":
    plotParticleHistory()
    

