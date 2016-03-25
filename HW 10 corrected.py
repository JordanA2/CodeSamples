import numpy as N
from matplotlib import pylab as P
#mags is an array that contains the array of coordinates for each of the magnets
mags = [[-.250,.433],[-.250,-.433],[.500,000]]
R = .2
C = .5
d = .25
dt = .1
tmax = 100

def dxdt(v):
    """diff eq for position, takes 1x2 array for velocity as input"""
    return N.array([v[0],v[1]]) 
    
def dvdt(x,v):
    """diff eq given for velocity, takes 1x2 array for velocity and position
    as input"""
    xv = -R*v[0]-C*x[0] + sum([((mags[s][0]-x[0])/(((mags[s][0]-x[0])**2
    +(mags[s][1]-x[1])**2 + d**2)**(3./2.))) for s in range(len(mags))])
    yv = -R*v[1]-C*x[1] + sum([((mags[s][1]-x[1])/(((mags[s][0]-x[0])**2
    +(mags[s][1]-x[1])**2 + d**2)**(3./2.))) for s in range(len(mags))])
    return N.array([xv,yv])

def RungKuta(dxdt,dvdt,x0,v0):
    """Prints the implicit graph of x(t) and y(t) described by the function
given in the homework. Takes as inputs dxdt: a diff eq for x and y; dvdt: a diff
eq for Vx and Vy; x0: an array of initial positions in a 1x2 array; v0: an array
of initial velocities in a 1x2 array"""
    t = N.arange(0,tmax,dt)
    x = list(N.arange(0,tmax,dt))
    v = list(N.arange(0,tmax,dt))
    plotx = list(N.arange(0,tmax,dt))
    ploty = list(N.arange(0,tmax,dt))
    idx = 0
    x[idx] = x0
    v[idx] = v0
    for s in t[0:(len(t)-1)]:
        plotx[idx] = x[idx][0]
        ploty[idx] = x[idx][1]
        idx = idx + 1
        K1 = dt*dvdt(x[idx-1],v[idx-1])
        Ka1 = dt*dxdt(v[idx-1])
        K2 = dt*dvdt(x[idx-1]+0.5*Ka1,v[idx-1]+0.5*K1)
        Ka2 = dt*dxdt(v[idx-1]+0.5*K1)
        K3 = dt*dvdt(x[idx-1]+0.5*Ka2,v[idx-1]+0.5*K2)
        Ka3 = dt*dxdt(v[idx-1]+0.5*K2)
        K4 = dt*dvdt(x[idx-1]+Ka3,v[idx-1]+K3)
        Ka4 = dt*dxdt(v[idx-1]+K3)
        v[idx] = v[idx-1] + (1./6.)*(K1+2*K2+2*K3+K4)
        x[idx] = x[idx-1] + (1./6.)*(Ka1+2*Ka2+2*Ka3+Ka4)
    plotx[-1] = plotx[-2]
    ploty[-1] = ploty[-2]
    return plotx, ploty

        
def endmag(x0,v0):
    """takes 1x2 arrays of initial position and velocity, returns as an array
    the location of the magnet the pendulum approaches"""
    a,b = RungKuta(dxdt,dvdt,x0,v0)
    c = [((a[-1]-mags[s][0])**2+(b[-1]-mags[s][1])**2)**(0.5) for s in \
    range(len(mags))]
    for s in range(len(c)):
        if c[s] == min(c):
            return ([mags[s][0],mags[s][1]])
            
def maggraph(x0,v0):
    """test function, please ignore"""
    a,b = RungKuta(dxdt,dvdt,x0,v0)
    P.plot(a,b,'c')
    P.show()
    
def maggraph3():
    """produces the graph for question 3"""
    a,b = RungKuta(dxdt,dvdt,[-1,1],[0,0])
    c,d = RungKuta(dxdt,dvdt,[-.75,-.75],[0,0])
    e,f = RungKuta(dxdt,dvdt,[.7,.1],[0,0])
    P.plot(a,b,'c',c,d,'r',e,f,'g',-.25,-.433,'+r',.500,0,'+g')
    P.show()
    
def maggraph4():
    """produces the graph for question 4"""
    a,b = RungKuta(dxdt,dvdt,[-1,.3],[0,0])
    c,d = RungKuta(dxdt,dvdt,[.5,0],[0,0])
    e,f = RungKuta(dxdt,dvdt,[-1,-.3],[0,0])
    g,h = RungKuta(dxdt,dvdt,[-.3,0],[0,0])
    i,j = RungKuta(dxdt,dvdt,[-.3,-.001],[0,0])
    k,l = RungKuta(dxdt,dvdt,[-.3,.001],[0,0]) 
    P.subplot(2,1,1)
    P.plot(a,b,'c',c,d,'r',e,f,'g')
    P.subplot(2,1,2)
    P.plot(g,h,'m',i,j,'b',k,l,'y')
    P.show() 

if __name__ == '__main__':
    a, b = RungKuta(dxdt, dvdt, [-.24,.4], [3,3])
    maggraph3()
    maggraph4()
    P.plot(a, b, 'c', -.25, -.433, '+b', -.25, .433, '+b', .5, 0, '+b')
    P.show()

