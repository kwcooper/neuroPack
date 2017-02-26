#compNeuro
#Anderson, Britt


#Newton's method

#dy/dx = C

'''To “solve” a differential equation requires us to find
another equation that will make our DE true.
Restated, if our DE involves dy/dx we have to find
some function y = . . . that when we take its derivative will
give us the right-hand side of the dy/dx = . . . equation we started with. '''

#y = dy/dx = C


'''Is our model too simple? We can only deter-
mine thatafter implementing it. If we fail to get behavior
that matches empirical observations we will have to increase the
model’s complexity, but until then we will simplify without mercy. This eases our work,
and it makes our model understandable. This penchant for the whole- sale adoption of
simplifying assumptions is the basis for many jokes at the expense of the stereotypical
theoretical physicist. A poultry farmer approaches his neighbor, the vacationing physicist,
to help the farmer figure out the most efficient way to house his hens:
“That’s trivial,” replies the physicist. “First, we assume a spherical chicken. . .”'''

'''
import pylab as pyl

dt = 0.05
p = -5.0
sp = 5.0

acc = [p*sp]
vel = [0.0]
s = [sp]
t = [0.0]
for i in range (1 ,100):
    acc . append ( s[−1]*p)
    vel.append(vel[-1] + acc[-1]*dt)
    s.append(s[−1] + vel[−1]*dt)
    t . append ( dt*i )

dp = pyl.plot(t,s)
pyl . show ( )
'''

'''Gerstner and Kistler (2002): Spiking Neuron Models.
This is an excellent textbook that has an online
version http://lcn.epfl.ch/ gerstner/BUCH.html.'''





