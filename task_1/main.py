import matplotlib.pyplot as plt
import numpy as np
def toch(x,y): # точное решение
    z = np.cos(x + y)
    return z
def f(x,y): # правая часть
    z = -np.sin(x + y) + np.cos(x + y)
    return z
def A(a,m): # создание матрицы
    B = np.zeros((a,a))
    for i in range(1,a-1):
        for j in range(0,a):
            if i==j:
                B[i,j-1] = m[0]
                B[i,j]=m[1]
                B[i,j+1]=m[2]
    return B
def F(a,b,h,N,sgm,n,tau):
    Z=np.arange(a,b+h,h)
    y=np.zeros(N+1)
    y[0]=0
    for i in range(1,N):
        y[i]=-(1-sgm)*toch(Z[i+1],n*tau)-((h**2/tau)-2*(1-sgm))*toch(Z[i],n*tau)-(1-sgm)*toch(Z[i-1],n*tau)-h**2*f(Z[i],n*tau)
    y[N]=0
    return y

def Metod_Progonki(A,B):
    #известные константы
    k1 = -A[0,1]
    m1 = B[0]
    k2 = -A[A.shape[0] - 1, A.shape[1] - 2]
    m2 = B[B.shape[0] - 1]
    alfa = k1
    beta = m1
    #поиск альф и бет
    c = 2
    a = 0
    b = 1
    alf = [alfa]
    bet = [beta]
    for i in range(1, A.shape[0] - 1):
        beta = (B[i] - A[i,a] * beta) / (A[i,a] * alfa + A[i,b])
        alfa = -A[i,c] / (A[i,a] * alfa + A[i,b])
        a += 1
        b += 1
        c += 1
        alf.append(alfa)
        bet.append(beta)
    #расчет игриков
    y = (k2 * beta + m2) / (1 - k2 * alfa)
    otv = [y]
    for i in range(len(alf) - 1, -1, -1):
        y = alf[i] * y + bet[i]
        otv.append(y)
    #переворачиваем значения в списке
    otvet = []
    for i in reversed(otv):
        otvet.append(i)
    return otvet
# параметры
alpha=1
betta=1
sgm=0.5
a=0
b=2*np.pi
N=29
h=(b-a)/N
n=0
tau=1/(4*N)
t=n*tau
A_1=sgm
B_1=A_1
C_1=-((h**2/tau)+2*sgm)
r=np.zeros(N+1)
r=F(a,b,h,N,sgm,n,tau)
r[0]=1
r[N]=1
m=[0,0,0]
m[0]=A_1
m[1]=C_1
m[2]=B_1
Matrix=A(N+1,m)
Matrix[0,0]=1
Matrix[N,N]=1
x = Metod_Progonki(Matrix, r)
print(x)
t_1=np.arange(a,b+h,h)
u=np.zeros(N+1)
for i_1 in range(0,N+1):
    u[i_1]=toch(t_1[i_1],0)
plt.plot(t_1,x,'b*',t_1,u,'r--')
plt.show()