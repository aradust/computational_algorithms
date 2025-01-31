import numpy as np
import matplotlib.pyplot as plt
from scipy. special import logsumexp
def TR(x,t):
    if x<=t:
        return 1
    if x>-2*t or x<t:
        return -3*x+1
    if x>=-2*t:
        return -2
def R(N):
 N_t=N*2
 h=1/N
 tau=1/N_t
 X=np.arange(0,1+h,h)
 T=np.arange(0,1+tau,tau)
 mu=1
 co_1=logsumexp(tau/(4*h))
 co_2=logsumexp((tau/(2*h)))
 co_3=logsumexp(mu*tau/2)
 co_4=logsumexp((1-mu*tau))
 V=np.zeros((N_t+1,N+1))
 L=np.zeros((N_t+1,N+1))
 C_1=np.zeros(N+1)
 C_2=np.zeros(N+1)
 C_3=np.zeros(N+1)
 L[0][0]=1
 for i in range(1,N):
     L[0][i]=-3*X[i]+1
 L[0][N]=-2
 n=0
 while n<N_t:
  for i in range(0,N):
      V[n][i]=logsumexp(co_1*((L[n][i])**2)-co_1*((L[n][i+1])**2)+0.5*L[n][i+1]+0.5*L[n][i])
  L[n+1][0]=1
  for i in range(0,N-1):
      L[n+1][i+1]=logsumexp(co_2*((V[n][i])**2)-co_2*((V[n][i+1])**2)+co_3*L[n][i]+co_4*L[n][i+1]+co_3*L[n][i+2])
  L[n+1][N]=-2
  for i in range(0,N+1):
      C_1[i]=L[n][i]
      C_2[i]=TR(X[i],T[n])
  plt.plot(X,C_1 )
  plt.show()
  n=n+1
 #K=np.zeros(N+1)
 #for i in range(0,N+1):
   #  K[i]=TR(X[i],T[2])
 #plt.plot(X,K,'r',X,C_3,'b')
R(10)
plt.show()