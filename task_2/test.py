import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
def toch(x,y):
    return np.cos(x+y)*np.sin(x*y)
def f(x,y):
    return -np.cos(x+y)*np.sin(x*y)*(2+x**2+y**2)-2*np.sin(x+y)*np.cos(x*y)*(x+y)
def fp(a,b,c,d,N_x,N_y,h_1,h_2):
    A=np.zeros((N_y+1,N_x+1))
    B=np.arange(a,b+h_1,h_1)
    C=np.arange(c,d+h_2,h_2)
    for i in range(0,N_y+1):
        for j in range(0,N_x+1):
            A[i][j]=toch(B[j],C[i])
    return A


# параметры
def RES(N_x,N_y,par):
 Z_1=np.zeros(N_x+1)
 a=(np.pi)/2
 b=3*(np.pi)/2
 c=(np.pi)/2
 d=3*(np.pi)/2
 h_1=(b-a)/N_x
 h_2=(d-c)/N_y
 tau=10*h_1
 A_1=-1/(h_1**2)
 C_1=(1/(tau))+(2/(h_1**2))
 B_1=A_1
 A_2=-1/(h_2**2)
 C_2=(1/(tau))+(2/(h_2**2))
 B_2=A_2
 TR=fp(a, b, c, d, N_x, N_y, h_1, h_2)
 Setka_j=np.arange(a,b+h_1,h_1)
 Setka_m=np.arange(c,d+h_2,h_2)
 ko_1=(1/(h_2**2))
 ko_2=((1/tau)-(2/(h_2**2)))
 ko_3=(2/(h_2**2))
 ko_4=(1/tau)
 n_1=-A_1
 n_2=(-2/(h_1**2))
 n_3=n_1
 n_4=-A_2
 n_5=(-2/(h_2**2))
 n_6=n_4
 eps =2*10**(-3)
 list_1=[0]
 list_2=[0]
#прогонка для n-1/2

 Matrix_1=np.zeros((N_x+1,N_y+1))
 Matrix_2=np.zeros((N_x+1,N_y+1))
 k=np.zeros(N_y+1)
 v=np.zeros(N_y+1)
 L=np.zeros((N_x+1,N_y+1))
 y=0
 Nev_max=0
 Nev=0
 while  y<=100:
     if(y==0):
         for i in range(0,N_y+1):
             for j in range(0,N_x+1):
                 L[i][j]=toch(Setka_j[j],Setka_m[i])
     for j in range(1, N_x):
         k[0] = 0
         v[0] = L[0][j]
         for m in range(1, N_y):
             G = ko_1 * L[m][j-1] + ko_2 * L[m][j] + ko_1 * L[m][j+1] - f(Setka_j[j], Setka_m[m])
             k[m] = -B_1 / (C_1 + A_1 * k[m - 1])
             v[m] = -(A_1 * v[m - 1] - G) / (C_1 + A_1 * k[m - 1])
         k[N_y] = 0
         v[N_y] = L[N_y][j]
         Matrix_1[N_y][j] = (v[N_y] + k[N_y] * v[N_y - 1]) / (1 - k[N_y] * k[N_y - 1])
         for m in range(N_y - 1, -1, -1):
             Matrix_1[m][j] = k[m] * Matrix_1[m + 1][j] + v[m]

     for j in range(0, N_x + 1):
         Matrix_2[0][j] = L[0][j]
     for m in range(1, N_y):
         k[0] = 0
         v[0] = L[m][0]
         for j in range(1, N_x):
             G = ko_4 * Matrix_1[m][j] - ko_1 * L[m][j-1] + ko_3 * L[m][j] - ko_1 * L[m][j+1]
             k[j] = -B_2 / (C_2 + A_2 * k[j - 1])
             v[j] = -(A_2 * v[j - 1] - G) / (C_2 + A_2 * k[j - 1])
         k[N_x] = 0
         v[N_x] = L[m,N_x]
         Matrix_2[m][N_x] = (v[N_x] + k[N_x] * v[N_x - 1]) / (1 - k[N_x] * k[N_x - 1])
         for j in range(N_x - 1, -1, -1):
             Matrix_2[m][j] = k[j] * Matrix_2[m][j + 1] + v[j]
     for j in range(0, N_x + 1):
         Matrix_2[N_y][j] = L[N_y][j]
     for m in range(1,N_y):
         for j in range(1,N_x):
             Nev=np.abs(n_1*Matrix_2[m-1][j]+n_2*Matrix_2[m][j]+n_3*Matrix_2[m+1][j]+n_4*Matrix_2[m][j-1]+n_5*Matrix_2[m][j]+n_6*Matrix_2[m][j+1]-f(Setka_j[j],Setka_m[m]))
             if(Nev> Nev_max):
                 Nev_max=Nev
     Nev_max=Nev
     y = y + 1
     L=Matrix_2
 for i in range(0,N_x+1):
     Z_1[i]=Matrix_2[par][i]
 plt.plot(Setka_m,Z_1,'r*')
Z_2=np.zeros(21)
par=5
a=(np.pi)/2
b=3*(np.pi)/2
c=(np.pi)/2
d=3*(np.pi)/2
h_1=(b-a)/20
h_2=(d-c)/20
h=(b-a)/10
Setka=np.arange(a,b+h_1,h_1)
TR=fp(a,b,c,d,20,20,h_1,h_2)
for i in range(0,21):
    Z_2[i]=TR[21][i]
RES(20,20,par)
a=RES(10,10,par)
for i in range(0,11):
    a[i]=TR[par][i]
Setka_1=np.arange(a,b+h,h)
plt.plot(Setka,Z_2,'b--',Setka_1,)
plt.show()
#относительная ошибка
#порядок
#разные N_x и N_y