import matplotlib.pyplot as plt
import numpy as np
def toch(x,t): # точное решение
    return np.cos(x + t)
def f(x,t): # правая часть
    return -np.sin(x + t) + np.cos(x + t)
N=500
T=500
a=0
b=np.pi
h=(b-a)/N
k=np.zeros(N+1)
v=np.zeros(N+1)
Setka=np.arange(a,b+h,h)
tau=1/(4*N)
sgm=0.5
ko_1=((1-sgm)/h**2)
ko_2=((1/tau)-(2*(1-sgm)/(h**2)))
A=-sgm/(h**2)
B=A
C=(1/(tau))+(2*sgm/(h**2))
res=np.zeros(N+1)
TR=np.zeros(N+1)
n=0
n_1=1/(h**2)
n_2=-2/(h**2)
L=np.zeros(N+1)
for i in range(0,N+1):
    L[i]=toch(Setka[i],n*tau)
while n<=T:
   k[0] = 0
   v[0] = toch(Setka[0],n*tau)
   for i in range(1,N):
       #F=ko_1*toch(Setka[i-1],n*tau)+ko_2*toch(Setka[i],n*tau)+ko_1*toch(Setka[i+1],n*tau)+f(Setka[i],n*tau)
       F=ko_1*L[i-1]+ko_2*L[i]+ko_1*L[i+1]+f(Setka[i],n*tau)
       k[i]=-B/(C+A*k[i-1])
       v[i]=-(A*v[i-1]-F)/(C+A*k[i-1])
   k[N]=0
   v[N]=toch(Setka[N],n*tau)
   res[N]=(v[N]+k[N]*v[N-1])/(1-k[N]*k[N-1])
   for i in range(N-1,-1,-1):
       res[i]=k[i]*res[i+1]+v[i]
   L=res
   for i in range(0,N+1):
       TR[i]=toch(Setka[i],n*tau)
   n = n + 1
plt.plot(Setka,res,'r--',Setka,TR,'b*')
plt.show()
# h + прогонка
