import matplotlib.pyplot as plt
import numpy as np
def toch(x,t): # точное решение
    return np.cos(x + t)
def f(x,t): # правая часть
    return -np.sin(x + t) + np.cos(x + t)
def w(a,t,alpha):
    return np.cos(a + t)+np.sin(a+t)/alpha
def V(b,t,betta):
    return np.cos(b + t)-np.sin(b+t)/betta
def Gauss(N,h,alpha,betta,n,tau,sgm,ko_1,ko_2,a,b):
    Setka = np.arange(a, b + h, h)
    A_1=((-1.5)/h)-alpha
    C_1=2/h
    B_1=(-0.5)/h
    F_1=-alpha*w(a,n,alpha)
    A=-sgm/(h**2)
    C=((1/tau)+2*sgm/(h**2))
    B=A
    A_2=((1.5)/h)+betta
    C_2=-2/h
    B_2=(0.5)/h
    F_2=betta*V(b,n,betta)
    G_1=ko_1 * toch(Setka[0], n) + ko_2 * toch(Setka[1], n) + ko_1 * toch(Setka[2], n) + f(Setka[1], n)
    G_2 = ko_1 * toch(Setka[N-2], n ) + ko_2 * toch(Setka[N-1], n) + ko_1 * toch(Setka[N], n) + f(Setka[N-1],  n)
    co_1=(-B / B_1)
    co_2=(-A / B_2)
    C_1=C_1*co_1+C
    A_1=A_1*co_1+A
    F_1=F_1*co_1+G_1
    C_2=C_2*co_2+C
    A_2=A_2*co_2+B
    F_2=F_2*co_2+G_2
    a=[F_1/A_1,-C_1/A_1,A_2,F_2,C_2]
    return a
def RES(N):
 t=10
 a=0
 b=4
 h=(b-a)/N
 k=np.zeros(N+1)
 v=np.zeros(N+1)
 Setka=np.arange(a,b+h,h)
 alpha=1
 betta=1
 tau=0.4*h**2
 T=np.arange(0,t+tau,tau)
 sgm=0.5
 ko_1=((1-sgm)/h**2)
 ko_2=((1/tau)-(2*(1-sgm)/(h**2)))
 A=-sgm/(h**2)
 B=A
 C=(1/(tau))+(2*sgm/(h**2))
 res=np.zeros(N+1)
 TR=np.zeros(N+1)
 n=T[1]
 L=np.zeros(N+1)
 GAU=Gauss(N, h, alpha, betta, n, tau, sgm, ko_1, ko_2, a, b)
 list_2=[]
 Nev=0
 n=0
 for i in range(0,N+1):
     L[i]=toch(Setka[i],T[0])
 while n!=len(T)-1:
    k[0] = GAU[1]
    v[0] = GAU[0]
    for i in range(1,N):
        F=ko_1*L[i-1]+ko_2*L[i]+ko_1*L[i+1]+f(Setka[i],T[n])
        k[i]=-B/(C+A*k[i-1])
        v[i]=-(A*v[i-1]-F)/(C+A*k[i-1])
    res[N]=(GAU[3]-GAU[4]*v[N-1])/(GAU[2]+GAU[4]*k[N-1])
    for i in range(N-1,-1,-1):
        res[i]=k[i]*res[i+1]+v[i]
    L=res
    list_2.append(n)
    n = n + 1
    if n!=len(T)-1:
     GAU = Gauss(N, h, alpha, betta, T[n+1], tau, sgm, ko_1, ko_2, a, b)
 for i in range(0, N + 1):
     TR[i] = toch(Setka[i], T[n])
 Nev = np.max(np.abs(L - TR))
 plt.plot(Setka,res,'r*')
 return Nev

N=200
K=np.zeros(N+1)
TR=np.zeros(N+1)
a = 0
b = 4
h = (b - a) / N
Setka=np.arange(a,b+h,h)
for i in range(0,N+1):
    TR[i]=toch(Setka[i],10 )
plt.plot(Setka,TR,'b--')

RES(20)
RES(40)
a=RES(80)
RES(100)
RES(120)
b=RES(160)
#RES(200)
plt.show()
print("Порядок: ")
print(np.log2(a / b))