#Hodg and Hux

#A. L. Hodgkin and A. F. Huxley, "A Quantitative Description
#of Membrane Current and Its Application to Conduction and
#Excitation in Nerve," Journal of Physiology, 117, 1952 pp. 500–544.

#Adapted from Anderson, Britt

#C(dV(t)/dt) = I(t) − [(gNa*(m^3) * h* (V(t) − ENa) + gKn^4(V(t)−EK)+gL(V(t)−EL)]


import matplotlib.pyplot as plt
import math as m

vInit = 0.0
dt = 0.01

#Constraints from Gerstner and Kistler (2002)
eNa = 115       #Sodium reversal
gNa = 120       #Sodium Conductance
eK = -12          #Potassium reversal
gK = 36           #Potassium Conductance
eL = 10.6        #Leak voltage
gL = 0.3          #Leak Voltage

#Updating Rule
def upd(x, dltaX):
    return(x + dltaX * dt)

def mnh0(a,b):
    return(a/(a+b))

#AB
def alphaM(v):
    return ((2.5 - 0.1*v) / (m.exp (2.5 - 0.1*v) - 1))

def betaM(v):
    return (4 * m.exp((-1) * v / 18))

def alphaH(v):
    return ( 0.07 * m.exp(( -1) * v / 20))

def betaH(v):
    return ( 1.0 / (m.exp(3 - (0.1) * v) +1))

def alphaN(v):   
    return ((0.1 - 0.01 * v) / (m.exp(1 - (0.1 * v)) -1))

def betaN(v):
    return 0.125 / m.exp((-1) *  v / 80.0)

am0 = alphaM(0)
bm0 = betaM(0)
an0 = alphaN(0)
bn0 = betaN(0)
ah0 = alphaH(0)
bh0 = betaH(0)

m0 = mnh0(am0,bm0)
n0 = mnh0(an0,bn0)
h0 = mnh0(ah0,bh0)

#Membrane Activity
#  Sodium 
def INa(m, h, V):
    return gNa * (m**3) * h * (V - eNa)
#  Potassium
def IK(n, V):
    return gK  * (n**4) * (V - eK)
#  Leak
def IL(V):
    return gL  * (V - eL)


#def iStim(t): 
#    return 10*(t>100) - 10*(t>200) + 35*(t>300)

#Updates values
def newS(V,m,n,h,t):
    if (t < 6.0) or (t > 7.0):
        iStim = 0.0
    else:
        iStim = 20.0

    dv = iStim - (INa(m, h, V) + IK(n, V)  + IL(V))
    dm = alphaM(V) * (1.0 - m) - betaM(V) * m 
    dn = alphaN(V) * (1.0 - n) - betaN(V) * n
    dh = alphaH(V) * (1.0 - h) - betaH(V) * h

    vp = upd(V, dv)
    tp = t + dt
    mp = upd(m, dm)
    np = upd(n, dn)
    hp = upd(h, dh)

    return vp, mp, np, hp, tp


vs = []
ms = []
ns = []
hs = []
ts = []

a,b,c,d,e = newS(vInit ,m0,n0,h0,0.0)
vs.append(a)
ms.append(b)
ns.append(c)
hs.append(d)
ts.append(e)

for i in range(2, 3000):
    a,b,c,d,e = newS(vs[-1], ms[-1], ns[-1], hs[-1], ts[-1])
    vs.append(a)
    ms.append(b)
    ns.append(c)
    hs.append(d)
    ts.append(e)

plt.plot(ts, vs)
plt.show()









