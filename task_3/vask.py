import numpy as np
import matplotlib.pyplot as plt
def toch(t, x):
    if (t < x):
        return 0
    if (t >= x):
        return 1
def R(N):
 N_t=N*2
 mu=0
 tau=1/N_t
 h=1/N
 X=np.arange(0,1+h,h)
 T=np.arange(0,1+tau,tau)
 A=-1/(2*h)
 C=3/(2*tau)
 B=-A
 co_3=2/tau
 co_4=-1/(2*tau)
 co_5=((2/tau)-mu)
 co_6=mu/2
 u=np.zeros(N+1)
 co_1=(0.5-(tau/(2*h)))
 co_2=(0.5+(tau/(2*h)))
 u_1=co_2
 u_2=co_1
 L=np.zeros((N_t+1,N+1))
 k=np.zeros(N+1)
 v=np.zeros(N+1)
 res=np.zeros(N+1)
 u[0]=1
 u[1]=u_1*tau/h
 L[0][0]=1
 for i in range(1,N+1):
     L[0][i]=0
 for i in range(2,N-1):
      u[i]=0
 u[N-1]=-(tau/h)*u_2
 u[N]=0
 for i in range(0,N+1):
     L[1][i]=u[i]
 n=0
 while n<N_t-1:
  k[0]=0
  v[0]=1
  for i in range(1, N):
      F = co_5 * L[n+1][i]+ co_4*L[n][i]+co_6*L[n+1][i-1] +co_6*L[n+1][i+1]
      k[i] = -B / (C + A * k[i - 1])
      v[i] = -(A * v[i - 1] - F) / (C + A * k[i - 1])
  res[N]=0
  for i in range(N-1,0,-1):
      res[i] = k[i] * res[i + 1] + v[i]
  res[0]=1
  n=n+1
  for i in range(0,N+1):
      L[n+1][i]=res[i]
  for i in range(0, N + 1):
     res[i] = L[n][i]
  plt.plot(X, res)
  plt.show()
R(50)

plt.show()