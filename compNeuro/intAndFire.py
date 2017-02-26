#integrate and fire model, a simplification of the HH model

#  j(dV(t)/dt) = RI(t) - V(t)

#where,
#dV/dt = rate of change of voltage over time
#j (traditional tau) is membrane constant. resistance and capasistance
# -(V(t)) allows for neuronal self correction
# I refers to current, or intensity

#  the voltage in the future will be a sum of whatever current is being
# injected into our cell minus some portion of the current voltage.

#additionally, 
#Charges, (Q)
#I = dQ/dt
#Capacatiance, (C) = Q/V


import matplotlib.pyplot as plt

r = 1
c = 1
j = r * c
dt = 0.05
t = 0
v = 0
thresh = 5
i = []
tData = []
vData = []

#current pulse
for z in range(0,40):
    num = 10
    i.append(num)

#return to 0
for z in range(40,75):
    num = 0
    i.append(num)

#voltage Calc

for ii in range(0,75):
    dvdt = (1/j) * (r*i[ii] - v)
    v = v + dvdt * dt
    if v > thresh:
        v = 0
    t = t + dt
    tData.append(t)
    vData.append(v)

plt.plot(tData, vData)
plt.axis([0,t,-1,7])
plt.xlabel('time')
plt.ylabel('V')
plt.show()
        




















