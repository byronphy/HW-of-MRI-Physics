import numpy as np
import matplotlib.pyplot as plt

G_z = 0.08 #T/m
gamma = 2*np.pi*42.58e6 #rad/(sT)
z_0 = 0.001 #m
tau_m = 0.005 #s
omega_0 = gamma*G_z*z_0

t_min = 0
t_max = tau_m
N = 1000
t = np.linspace(t_min, t_max, N) #s
#y = np.sin(omega_0*(t))/(np.pi*gamma*(t)) #center
#y = np.sin(omega_0*(t-tau_m/2))/(np.pi*gamma*(t-tau_m/2)) #case1
#y = np.sin(omega_0/2*(t-tau_m/2))/(np.pi*gamma*(t-tau_m/2)) #case2
#y = 2*np.sin(omega_0/2*(t-tau_m/2))/(np.pi*gamma*(t-tau_m/2))*np.cos(3*omega_0/2*(t-tau_m/2)) #case3
#y = -2*np.sin(omega_0/2*(t-tau_m/2))/(np.pi*gamma*(t-tau_m/2))*np.sin(3*omega_0/2*(t-tau_m/2)) #case4

"""
#signal of rf pulse

plt.plot(t, y, label="$H(t)$", color="blue", linewidth=2)

plt.xlabel("t/s")
plt.ylabel("H(t)/T")

plt.title("rf pulse")

plt.xlim(t_min, t_max)
#plt.ylim(-0.0001, 0.0001)

plt.legend()

plt.show()
"""


#ifft of rf pulse

tt = np.linspace(-tau_m/2, tau_m/2, N) #s
#yy = np.fft.fft(np.sin(omega_0*(tt))/(np.pi*gamma*(tt))) #case1
yy = np.fft.fft(np.sin(omega_0/2*(tt))/(np.pi*gamma*(tt))) #case2, and xx should plus 0.0015
#yy = np.fft.fft(2*np.sin(omega_0/2*(tt))/(np.pi*gamma*(tt))*np.cos(3*omega_0/2*(tt))) #case3

yyf = abs(yy) #case1,2,3

"""
#case4
yy = np.fft.fft(-2*np.sin(omega_0/2*(tt))/(np.pi*gamma*(tt))*np.sin(3*omega_0/2*(tt))) # case4
c = 1
yyf = []
for yyy in yy:
    yyyf = abs(yyy)
    yyyphase = np.arctan(np.imag(yyy)/np.real(yyy))
    if yyyphase>0:
        c = -1
        yyyf = c*yyyf
    yyf.append(yyyf)
"""

xx = np.fft.fftfreq(N,tau_m/N)/gamma/G_z*2*np.pi #z-axis
#print(yyf)

plt.plot(np.fft.fftshift(xx),np.fft.fftshift(yyf), label="$iFFT$", color="blue", linewidth=2)
plt.title("iFFT of rf pulse")
plt.xlim(-0.003, 0.003)
plt.xlabel("z_0/m")
plt.ylabel("Relative Intensity")
plt.show()



