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
a=(np.pi)/2
b=3*(np.pi)/2
c=(np.pi)/2
d=3*(np.pi)/2
N_x=50
N_y=N_x
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
#прогонка для n-1/2

Matrix_1=np.zeros((N_x+1,N_y+1))
Matrix_2=np.zeros((N_x+1,N_y+1))
k=np.zeros(N_y+1)
v=np.zeros(N_y+1)
L=np.zeros((N_x,N_y))
for j in range(1,N_x):
    k[0] = 0
    v[0] = toch(Setka_j[j], Setka_m[0])
    for m in range(1,N_y):
        G = ko_1 * toch(Setka_j[j-1], Setka_m[m]) + ko_2 * toch(Setka_j[j], Setka_m[m]) + ko_1 * toch(Setka_j[j+1],Setka_m[m]) - f(Setka_j[j], Setka_m[m])
        k[m]=-B_1/(C_1+A_1*k[m-1])
        v[m]=-(A_1*v[m-1]-G)/(C_1+A_1*k[m-1])
    k[N_y] = 0
    v[N_y] = toch(Setka_j[j], Setka_m[N_y])
    Matrix_1[N_y][j]=(v[N_y]+k[N_y]*v[N_y-1])/(1-k[N_y]*k[N_y-1])
    for m in range(N_y-1,-1,-1):
        Matrix_1[m][j]=k[m]*Matrix_1[m+1][j]+v[m]
for j in range(0, N_x + 1):
    Matrix_2[0][j] = toch(Setka_j[j], Setka_m[0])
for m in range(1, N_y):
    k[0] = 0
    v[0] = toch(Setka_j[0], Setka_m[m])
    for j in range(1, N_x):
        G =ko_4*Matrix_1[m][j]-ko_1*toch(Setka_j[j-1],Setka_m[m])+ko_3*toch(Setka_j[j],Setka_m[m])-ko_1*toch(Setka_j[j+1],Setka_m[m])
        k[j] = -B_2 / (C_2 + A_2 * k[j - 1])
        v[j] = -(A_2 * v[j - 1] - G) / (C_2 + A_2 * k[j - 1])
    k[N_x] = 0
    v[N_x] = toch(Setka_j[N_x], Setka_m[m])
    Matrix_2[m][N_x] = (v[N_x] + k[N_x] * v[N_x - 1]) / (1 - k[N_x] * k[N_x - 1])
    for j in range(N_x - 1, -1, -1):
        Matrix_2[m][j] = k[j] * Matrix_2[m][j+1] + v[j]
for j in range(0, N_x + 1):
    Matrix_2[N_y][j] = toch(Setka_j[j], Setka_m[N_y])
H =np.zeros(N_x+1)
H_1=np.zeros(N_x+1)
for i in range(0,N_x+1):
    H_1[i]=Matrix_2[5][i]
for i in range(0,N_x+1):
    H[i]=TR[5][i]
plt.plot(Setka_j, H, 'r*',Setka_j,H_1,'b--')
plt.show()



