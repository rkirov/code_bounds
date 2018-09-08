''' Functions from Bras-Amoros O'Sullivan's paper.

'''
from math import *

#q0=a=4
#q=q0**2
#g=(a-1)*a/2
#m=q0+1
#K=2*g-2

def lamb(i,a):
    if(i<=(a-1)*a/2):
        col=0
        count=0
        while(count<i-col):
            count+=(col+1)
            col+=1
        return col*a+i-count
    else:
        return i+(a-1)*a/2
    
def r(t,a):
    '''Reduncancy 
    t - errors corrected 
    '''
    TAFLOOR=(2*t)//(a+1)
    if t<=a/2:
        return t*(2*t+1)
    if t<a*(TAFLOOR+1)/2:
        return (a**2-a)/2+(a+1)*TAFLOOR
    else:
        return (a**2-a)/2+2*t
    
def deltaxt(x,t):
    return (x+int((x**2+4*x-8*t)**0.5)+1)%2

def rtilda(t,a):
    TAFLOOR=(2*t)//(a+1)
    ROOTFLOOR=lambda x:int((x**2+4*x-8*t)**0.5)
    CEILING=int(ceil(2*(2*t+1)**0.5-2))
    if t<=a/2:
        return t*(2*t+1)-sum([ROOTFLOOR(x)+deltaxt(x,t) for x in range(CEILING,2*t)])
    if t<a*(TAFLOOR+1)/2:
        return (a**2-a)/2+(a+1)*TAFLOOR-sum([ROOTFLOOR(x)+deltaxt(x,t) for x in range(CEILING,a-1+TAFLOOR)])
    if t<=a*(a+1)/2:
        return (a**2-a)/2+2*t-sum([ROOTFLOOR(x)+deltaxt(x,t) for x in range(CEILING,a+TAFLOOR)])
    return (a**2-a)/2+2*t
    
def v(i,a):
    '''that is bound for L((k+1)P)/L(kP) for the i-th non-empty coset.'''
    l=lamb(i,a)
    x=l/a
    y=l-a*x
    if(-a+x<=y<=a-1):
        return (x-y+1)*(y+1)
    else:
        return l-a*(a-1)+1

def d(i,a):
    l=lamb(i,a)
    x=l/a
    y=l-a*x
    if(x<a):
        if(y==x):
            return x+2
        else:
            return x+1
    if(-a+x<=y<a-1):
        return a*(x-a+2)
    else:
        return l-a*(a-1)+2
