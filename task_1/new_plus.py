import matplotlib.pyplot as plt
import numpy as np
def toch(x,t): # точное решение
    return np.cos(x + t)
def f(x, t, tau, h):  # правая часть
    return (-np.sin(x + t) + np.cos(x + t)) + (tau / 2) * (-np.cos(x + t) - np.sin(x + t)) + ((h ** 2) / 12) * (np.sin(x + t) - np.cos(x + t))
# правая часть
def w(a,t,alpha):
    return np.cos(a + t)+np.sin(a+t)/alpha
def V(b,t,betta):
    return np.cos(b + t)-np.sin(b+t)/betta
def Gauss(N,h,alpha,betta,n,tau,sgm,ko_1,ko_2,a,b):
    Setka = np.arange(a, b + h, h)
    A_1=(-25/(12*h))-alpha
    C_1=4/h
    B_1=-3/h
    D_1=4/(3*h)
    E_1=-1/(4*h)
    A_2=(25/(12*h))+betta
    C_2=-C_1
    B_2=-B_1
    D_2=-D_1
    E_2=-E_1
    A = -sgm / (h ** 2)
    C = ((1 / tau) + 2 * sgm / (h ** 2))
    B = A
    G_1 = ko_1 * toch(Setka[0],n) + ko_2 * toch(Setka[1],n) + ko_1 * toch(Setka[2],n) + f(Setka[1],n,tau,h)
    G_2 = ko_1 * toch(Setka[1],n) + ko_2 * toch(Setka[2],n) + ko_1 * toch(Setka[3],n) + f(Setka[2],n,tau,h)
    G_3 = ko_1 * toch(Setka[2],n) + ko_2 * toch(Setka[3],n) + ko_1 * toch(Setka[4],n) + f(Setka[3],n, tau, h)
    G_4 = ko_1 * toch(Setka[N-4],n) + ko_2 * toch(Setka[N-3],n) + ko_1 * toch(Setka[N-2],n) + f(Setka[N-3], n, tau, h)
    G_5 = ko_1 * toch(Setka[N-3],n) + ko_2 * toch(Setka[N-2],n) + ko_1 * toch(Setka[N-1],n) + f(Setka[N-2],n, tau, h)
    G_6 = ko_1 * toch(Setka[N-2],n) + ko_2 * toch(Setka[N-1],n) + ko_1 * toch(Setka[N],n) + f(Setka[N-1],n, tau, h)
    F_1 = -alpha * w(a, n, alpha)
    F_2 = betta * V(b, n, betta)
    co_1=(-B/E_1)
    D_1=D_1*co_1+C
    co_2 = (-B / D_1)
    B_1=B_1*co_1+A
    C_1=C_1*co_1
    A_1=A_1*co_1
    F_1=F_1*co_1+G_3
    B_1 = B_1 * co_2 + C
    co_3 = (-B / B_1)
    C_1 = C_1 * co_2+A
    A_1 = A_1 * co_2
    F_1 = F_1 * co_2 + G_2
    C_1 = C_1 * co_3 + C
    A_1 = A_1 * co_3 + A
    F_1 = F_1 * co_3 + G_1

    co_1 = (-A / E_2)
    D_2 = D_2 * co_1 + C
    co_2 = (-A / D_2)
    B_2 = B_2 * co_1 + B
    C_2 = C_2 * co_1
    A_2 = A_2 * co_1
    F_2 = F_2 * co_1 + G_4
    B_2 = B_2 * co_2 + C
    co_3 = (-A / B_2)
    C_2 = C_2 * co_2 + B
    A_2 = A_2 * co_2
    F_2 = F_2 * co_2 + G_5
    C_2 = C_2 * co_3 + C
    A_2 = A_2 * co_3 + B
    F_2 = F_2 * co_3 + G_6
    a = [F_1 / A_1, -C_1 / A_1, A_2, F_2, C_2]
    return a
def RES(N):
 t=1
 a=0
 b=10
 h=(b-a)/N
 k=np.zeros(N+1)
 v=np.zeros(N+1)
 Setka=np.arange(a,b+h,h)
 alpha=1
 betta=1
 tau=h**2
 T=np.arange(0,t+tau,tau)
 sgm=0.5-((h**2)/(12*tau))
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
 n=0
 for i in range(0,N+1):
     L[i]=toch(Setka[i],T[0])
 while n!=len(T)-1:
    k[0] = GAU[1]
    v[0] = GAU[0]
    for i in range(1,N):
        F=ko_1*L[i-1]+ko_2*L[i]+ko_1*L[i+1]+f(Setka[i],T[n],tau,h)
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
N=800
TR=np.zeros(N+1)
a = 0
b = 10
h = (b - a) / N
Setka=np.arange(a,b+h,h)
for i in range(0,N+1):
    TR[i]=toch(Setka[i],1   )
plt.plot(Setka,TR,'b--')

RES(20)
RES(40)
a=RES(80)
RES(100)
RES(140)
RES(180)
RES(200)
#b=RES(320)
#RES(800)
plt.show()
print("Порядок: ")
print(np.log2(a / b))