import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_file
rc_file('rcFile')
from glob import glob

nProcs = len(glob('scalar*dat'))
gridRatio = 4
tcf = 2.0000000000000000E-003
Tmix = 2.5

scalarData = []
for i in sorted(glob('scalar*dat')):
    scalarData.append(np.loadtxt(i,skiprows=2))
time = scalarData[0][:589,0]

scalarAbsorbed = np.zeros(np.size(time))
scalarTotal = np.zeros(np.size(time))
for i in range(nProcs):
    scalarAbsorbed = scalarAbsorbed + scalarData[i][:589,1]
    scalarTotal = scalarTotal + scalarData[i][:589,5]    

time = time*tcf


def analyzeScalarAbsorbedTime():
    global time

    fig = plt.figure()
    ax = plt.axes([0.13,0.18,0.97-0.13,0.92-0.18])
    plt.plot(time[10:], scalarAbsorbed[10:],label='Absorbed')
#    plt.plot(time, scalarTotal,label='Total')
#    plt.plot(time, scalarAbsorbed/scalarTotal,label='Absorbed')
#    plt.xlim(0,20)
    plt.legend(loc=0)
    plt.xlabel('Iterations')
    plt.ylabel('Scalar ($\\times 10^{-5}$ units)')
    plt.savefig('scalarAbsorbed.png')
    plt.show()
    plt.close(fig)

#    ratio = scalarAbsorbed[2000:589]/scalarTotal[2000:589]

#    print "Absorption over 90% after ", (np.where(ratio > 0.9)[0][0] + 2000)*tcf, " s"
    
    # fig = plt.figure()
    # ax = plt.axes([0.13,0.18,0.97-0.13,0.95-0.18])
    # plt.plot(time[2000:589], scalarAbsorbed[2000:589]/scalarTotal[2000:589],label='Absorbed')
    # plt.xlim(0,20)
    # plt.xlabel('Iterations')
    # plt.ylabel('Scalar Absorbed (units)')
    # plt.savefig('scalarAbsorbedRatio.png')

def compareScalarDomainVsReleased():
    global time

    sr = np.loadtxt('drugReleaseScalar.dat')  # Scalar Released
    for i in range(1,np.size(sr)):
        sr[i] = sr[i]+sr[i-1]
    
    fig = plt.figure()
    ax = plt.axes([0.13,0.18,0.97-0.13,0.92-0.18])
    plt.plot(time, scalarTotal,'+-', label='Domain')
    print np.shape(time)
#    print np.shape(sr[(gridRatio-1)::gridRatio])
#    plt.plot(time, sr[(gridRatio-1)::gridRatio], '+-', label='Released')
    print np.shape(sr)
    plt.plot(time, sr, '+-', label='Released')        
#    plt.xlim(0,20)
    plt.legend(loc=0)
    plt.xlabel('Time (s)')
    plt.ylabel('Scalar (mol)')
    plt.savefig('scalarDomainVsReleased.png')
    plt.show()
    plt.close(fig)


    
if __name__=="__main__":
#    analyzeScalarAbsorbedTime()
    compareScalarDomainVsReleased()
                                                    

    

    
