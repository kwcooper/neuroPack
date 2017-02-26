#Hodg and Hux
#Adapted from Anderson, Britt

import matplotlib.pyplot as plt
import math as m

vInit = 0.0
dt = 0.01

#Constraints from Gerstner and Kistler (2002)
eNa = 115 #Sodium reversal
gNa = 120 #Sodium Conductance
eK = -12 #Potassium reversal
gK = 36 #Potassium Conductance
eL = 10.6 #Leak voltage
gL = 0.3 #Leak Voltage


def upd(x, dltaX):
    return(x + dltaX * dt)

def mnh0(a,b):
    return(a/(a+b))

#AB
def alphaM(v):
    return 0.1*(v+40.0)/(1.0 - m.exp(-(v+40.0) / 10.0))

def betaM(v):
    return 4.0*m.exp(-(v+65.0) / 18.0)

def alphaH(v):
    return 0.07*m.exp(-(v+65.0) / 20.0)

def betaH(v):
    return 1.0/(1.0 + m.exp(-(v+35.0) / 10.0))

def alphaN(v):
    return 0.01*(v+55.0)/(1.0 - m.exp(-(v+55.0) / 10.0))

def betaN(v):
    return 0.125*m.exp(-(v+65) / 80.0)

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
    return gNa * m**3 * h * (V - eNa)
#  Potassium
def IK(n, V):
    return gK  * n**4 * (V - eK)
#  Leak
def IL(V):
    return gL  * (V - el)


#def iStim(t): 
#    return 10*(t>100) - 10*(t>200) + 35*(t>300)

def newS(V,m,n,h,t):
    if (t < 5.0) or (t > 6.0):
        iStim = 0.0
    else:
        iStim = 20.0

    dv = iStim - (INa(m, h, V) + IK(n, V)  + IL(V))
    dm = alphaM(V) * (1.0-m) - betaM(V) * m 
    dn = alphaN(V) * (1.0-n) - betaN(V) * n
    dh = alphaH(V) * (1.0-h) - betaH(V) * h
    










