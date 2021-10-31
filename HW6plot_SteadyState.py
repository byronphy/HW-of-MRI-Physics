import numpy as np
import matplotlib.pyplot as plt

TR = 52 #ms
TE = 14 #ms
T1 = 500 #ms
T2 = 300 #ms
T22 = 200 #ms

# rGE

def R(theta):
    matrix_R = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    return matrix_R

def b(t):
    vector_b = np.array([[0],[1-np.exp(-t/T1)]])
    return vector_b

def Gamma(t):
    star_Gamma = np.array([[np.exp(-t/T22),0],[0,np.exp(-t/T1)]])
    return star_Gamma

#rGE
def signal(theta):
    M1 = np.array([[np.exp(-TE / T22) * np.sin(theta)], [1 + (np.cos(theta) - 1) * np.exp(-TE / T1)]])
    signallist = [np.exp(-TE / T22) * np.sin(theta)]
    A = np.dot(np.dot(Gamma(TE), R(theta)), Gamma(TR - TE))
    B = np.dot(np.dot(Gamma(TE), R(theta)), b(TR - TE)) + b(TE)
    for i in range(n - 1):
        M1 = np.dot(A, M1) + B
        signallist.append(M1[0][0])
    return signallist

#aGE
def age_signal(theta):
    M1 = np.array([[np.exp(-TE / T22) * np.sin(theta)], [1 + (np.cos(theta) - 1) * np.exp(-TE / T1)]])
    age_signallist = [np.exp(-TE / T22) * np.sin(theta)]
    angle = theta
    for i in range(n - 1):
        angle = -angle
        A = np.dot(np.dot(Gamma(TE), R(angle)), Gamma(TR - TE))
        B = np.dot(np.dot(Gamma(TE), R(angle)), b(TR - TE)) + b(TE)
        M1 = np.dot(A, M1) + B
        age_signallist.append(abs(M1[0][0]))
    return age_signallist

#rSE
def rse_signal():
    M1 = np.array([[np.exp(-TE / T22)], [1 - np.exp(-TE / T1)]])
    rse_signallist = [np.exp(-TE / T2)]
    nostar_Gamma = np.array([[np.exp(-TE / T2), 0], [0, np.exp(-TE / T1)]])
    maR = np.array([[0, 1], [-1, 0]])
    A = np.dot(np.dot(nostar_Gamma, maR), Gamma(TR - TE))
    B = np.dot(np.dot(nostar_Gamma, maR), b(TR - TE)) + b(TE)
    for i in range(n - 1):
        M1 = np.dot(A, M1) + B
        rse_signallist.append(M1[0][0])
    return rse_signallist

#steady state
def steady(theta):
    Ms0 = np.sin(theta) * (1 - np.exp(-TR / T1)) * np.exp(-TE / T22) / (
                1 - (np.exp(-TR / T1) + np.exp(-TR / T22)) * np.cos(theta) + np.exp(-TR / T22 - TR / T1))
    Ms = np.sin(theta) * (1 - np.exp(-TR / T1)) * np.exp(-TE / T22) / (
                1 - (np.exp(-TR / T1) - np.exp(-TR / T22)) * np.cos(theta) - np.exp(-TR / T22 - TR / T1))
    Mss = (1 - np.exp(-TR/T1))* np.exp(-TE / T2) / (1 + np.exp(-(-TE+TR)/T22-TE/T2-TR/T1))
    return [Ms0,Ms,Mss]


n = 60
theta1 = 22.5 # angle, not rad
theta2 = 45
theta3 = 67.5
theta4 = 90

#check
callist = [signal(theta1/180*np.pi)[n-1],age_signal(theta1/180*np.pi)[n-1],rse_signal()[n-11]]
sslist = steady(theta1/180*np.pi)
print(callist)
print(sslist)
print((sslist[0]-callist[0])/callist[0])
print((sslist[1]-callist[1])/callist[1])
print((sslist[2]-callist[2])/callist[2])

x = list(range(1,n+1))

#plt.plot(x,age_signal(theta1/180*np.pi),'o-',color = 'silver',label = '%.1f$^\circ$'%theta1)
#plt.plot(x,age_signal(theta2/180*np.pi),'^-',color = 'darkgrey',label = '%.1f$^\circ$'%theta2)
#plt.plot(x,age_signal(theta3/180*np.pi),'s-',color = 'dimgrey',label = '%.1f$^\circ$'%theta3)
#plt.plot(x,age_signal(theta4/180*np.pi),'*-',color = 'black',label = '%.1f$^\circ$'%theta4)

plt.plot(x,rse_signal(theta4/180*np.pi),'*-',color = 'black',label = 'rSE')

plt.xlabel('n')
plt.xticks(np.arange(0,n,5))
plt.ylabel('M_y(TE)/M_0')
plt.legend()
plt.show()
