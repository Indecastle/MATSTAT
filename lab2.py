import math
import numpy as np
import sympy as sm
from itertools import groupby
import matplotlib.pyplot as plt
from tabulate import tabulate

x, y = sm.symbols('x, y', real = True)


def get_aba(m,d, yx):
    b = m - 3**0.5 * d
    a = 2*m - b
    #a, b = -1, 1
    if a > b:
        a, b = b, a
    print("a, b: ", a, b)
    print('f(a), f(b): ', yx(a), yx(b))
    return a, b


def get_XY(a, b, n, check=False):
    if check:
        Ei = lambda i: i/(n-1)
    else:
        Ei = lambda _: np.random.uniform()
    Xi = lambda i: a + Ei(i)*(b-a)
    X = [Xi(i) for i in range(n)]
    #print('X:', X)
    #Y = [0.04, 0.06, 0.06, 0.08, 0.08, 0.09, 0.09, 0.12, 0.12, 0.13, 0.13, 0.14, 0.17, 0.19, 0.2, 0.2, 0.22, 0.2, 0.23, 0.24, 0.26, 0.27, 0.29, 0.29, 0.32, 0.32, 0.32, 0.33, 0.38, 0.38, 0.38, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.47, 0.51, 0.52, 0.53, 0.53, 0.53, 0.54, 0.54, 0.55, 0.55, 0.59, 0.6, 0.6, 0.6, 0.61, 0.65, 0.65, 0.7, 0.73, 0.74, 0.75, 0.75, 0.75, 0.76, 0.8, 0.81, 0.81, 0.82, 0.82, 0.86, 0.86, 0.86, 0.86, 0.86, 0.87, 0.88, 0.88, 0.88, 0.89, 0.91, 0.92, 0.92, 0.93, 0.94, 0.94, 0.94, 0.94, 0.96, 0.96, 0.97, 0.98, 0.98, 0.99, 0.99, 0.99, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    #Y = [-6.237, -6.229, -5.779, -5.139, -4.950, -4.919, -4.636, -4.560, -4.530, -4.526, -4.523, -4.511, -4.409, -4.336, -4.259, -4.055, -4.044, -4.006, -3.972, -3.944, -3.829, -3.794, -3.716, -3.542, -3.541, -3.431, -3.406, -3.384, -3.307, -3.181, -3.148, -3.124, -3.116, -2.892, -2.785, -2.734, -2.711, -2.637, -2.633, -2.428, -2.381, -2.339, -2.276, -2.222, -2.167, -2.111, -2.034, -1.958, -1.854, -1.803, -1.774, -1.755, -1.745, -1.713, -1.709, -1.566, -1.548, -1.480 -1.448, -1.353, -1.266, -1.229, -1.179, -1.130, -1.102, -1.060, -1.046, -1.035, -0.969, -0.960, -0.903, -0.885, -0.866, -0.865, -0.774, -0.721, -0.688, -0.673, -0.662, -0.626, -0.543 ,-0.445, -0.241, -0.174, -0.131, 0.115, 0.205, 0.355, 0.577, 0.591, 0.795, 0.986, 1.068, 1.099, 1.195, 1.540, 2.008, 2.160, 2.534, 2.848]
    #Y = [1,1,0,0,5,0,1,2,2,1,0,0,1,1,3,4,4,3,2,1,1,0,0,0,2,1,3,3,4,3,7,7,6,7,5,4,1,0,0,0]
    Y = [yx(x) for x in X]
    Y.sort()
    #print('Y:', Y)
    return X, Y

def func1(Y, n):
    #M = int(math.sqrt(n))
    h = (Y[-1]-Y[0]) / M
    A = [Y[0] + (i-1)*h for i in range(1, M+1)]
    B = A[1:]+[Y[-1]]
    def vi(A, B):
        counter = 0
        for y in Y:
            if A==y or B==y:
                counter += 0.5
            if B <= y:
                break
            if A < y:
                counter += 1
        return counter
    v = [vi(A[i], B[i]) for i in range(M)]
    f = [v[i]/(n*h) for i in range(M)]
    
    XX = [A[0]]
    YY = [0]
    for i in range(M):
        XX.append(A[i])
        XX.append(B[i])
        YY.append(f[i])
        YY.append(f[i])
    XX.append(B[-1])
    YY.append(0)
    #print(M, h)
    #print('Y:', list(map(lambda x: round(x, 3), Y)))
    #print('A:', list(map(lambda x: round(x, 3), A)))
    #print('B:', list(map(lambda x: round(x, 3), B)))  
    #print('v:', v)
    #print('f:', list(map(lambda x: round(x, 3), f)))
    #print('XX:', list(map(lambda x: round(x, 3), XX)))
    #print('YY:', list(map(lambda x: round(x, 3), YY)))
    table = {'A':A, 'B':B, 'h':[h]*M, 'n':v, 'f':f}
    return XX, YY, table

def func2(Y, n):
    #M = int(math.sqrt(n))
    v = int(n/M)
    A = [Y[0]] + [(Y[v*i-1]+Y[v*i])/2 for i in range(1, M)]
    B = A[1:]+[Y[-1]]
    h = [B[i] - A[i] for i in range(M)]
    f = [v/(n*h[i]) for i in range(M)]
    XX = [A[0]]
    YY = [0]
    for i in range(M):
        XX.append(A[i])
        XX.append(B[i])
        YY.append(f[i])
        YY.append(f[i])
    XX.append(B[-1])
    YY.append(0)
    #print(M, v)
    #print(Y[39], Y[40], (Y[41]+Y[40])/2)
    #print('A:', list(map(lambda x: round(x, 3), A)))
    #print('B:', list(map(lambda x: round(x, 3), B)))  
    #print('h:', list(map(lambda x: round(x, 3), h)))
    #print('f:', list(map(lambda x: round(x, 3), f)))
    #print('XX:', list(map(lambda x: round(x, 3), XX)))
    #print('YY:', list(map(lambda x: round(x, 3), YY)))
    table = {'A':A, 'B':B, 'h':h, 'n':[v]*M, 'f':f}
    return XX, YY, table

def func3(Y, N, n):
    #group =  [[key, len(list(group))] for key, group in groupby(Y)]
    N = list(map(int, N))
    group = [[yi, ni] for yi,ni in zip(Y, N)]
    group.sort(key = lambda x: x[0])
    group = [[group[0][0] - 0.5, 0]] + group
    #print("groups:", group)
    XX = [group[0][0]]
    YY = [0]
    for i in range(1, len(group)):
        XX.append(group[i][0])
        XX.append(group[i][0])
        YY.append(YY[2*(i-1)])
        YY.append(YY[2*(i-1)] + group[i][1]/n)
    XX.append(group[-1][0] + 0.5)
    YY.append(YY[-1])
    return XX, YY


def get_XY_teor(fy, fx, fa, func=False):
    expr = sm.Eq(y, fy)
    xx = sm.solve(expr,x)
    xx0 = xx[0]
    #print('xx:', xx)
    if func:
        Gy = fx * (xx0 - xx0.subs({y:fa}))
    else:
        diff_xx = sm.diff(xx[0], y, 1)
        Gy = fx * diff_xx
    #=====
    #Gy = 0.5*sm.real_root(y,3) - 0.5*sm.real_root(-1,3)
    #Gy = 0.5*sm.sqrt(y)
    #=====
    #print('Gy:', Gy)
    Gy = sm.lambdify(y, Gy)
    return Gy

def plots(X1,Y1, X2,Y2, X3,Y3,   Y, table1, table2):
    #M = int(math.sqrt(n))
    
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        nrows=3, ncols=2,
        figsize=(8, 8)
    )
    AX = [ax1, ax2, ax3, ax4, ax5, ax6]
    ax1.plot(X1,Y1, label='эмпирическая плотность')
    ax1.plot(X3,Y3, label='теоретическая плотность')

    ax2.plot(X2,Y2, label='эмпирическая плотность')
    ax2.plot(X3,Y3, label='теоретическая плотность')

    pol_x1 = [(a+b)/2 for a, b in zip(table1['A'], table1['B'])]
    ax3.plot(pol_x1, table1['f'])
    ax3.plot(X3,Y3, label='теоретическая плотность')
    pol_x2 = [(a+b)/2 for a, b in zip(table2['A'], table2['B'])]
    ax4.plot(pol_x2, table2['f'])
    ax4.plot(X3,Y3, label='теоретическая плотность')

    Gy = get_XY_teor(fi2, fx, a2, func=True)
    vals = np.linspace(a2, yx(b), 10000)
    yyy = [Gy(x) for x in vals]
    
    XX, YY = func3(pol_x1, table1['n'], n)
    ax5.plot(XX, YY, label='эмпирическая')
    ax5.plot(vals, yyy, label='теоретическая')
    XX, YY = func3(pol_x2, table2['n'], n)
    ax6.plot(XX, YY, label='эмпирическая')
    ax6.plot(vals, yyy, label='теоретическая')

    for ax in AX:
        ax.grid()
    ax1.set_title('равноинтервальный метод')
    ax2.set_title('равновероятный метод')
    ax3.set_title('полигон, равноинтервальный метод')
    ax4.set_title('полигон, равновероятный метод')
    ax5.set_title('функция, равноинтервальный метод')
    ax6.set_title('функция, равновероятный метод')
    for ax in AX:
        ax.set_xlabel('Y')
        ax.set_ylabel('f*y')
    ax5.set_ylabel('Fy')
    ax6.set_ylabel('Fy')
    ax1.legend(loc='upper center', shadow=True, fontsize='small', frameon=False)
    ax2.legend(loc='upper center', shadow=True, fontsize='small', frameon=True)
    ax5.legend(loc='upper center', shadow=True, fontsize='small', frameon=False)
    ax6.legend(loc='upper center', shadow=True, fontsize='medium', frameon=False)
    fig.tight_layout()
    #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()

fi = sm.real_root(x,3)
fi2 = x**(1/3)
yx = sm.lambdify(x, fi)
yx2 = sm.lambdify(x, fi2)

n = 200
if n <= 100:
    M = int(math.sqrt(n))
else:
    M = int(2*math.log10(n))
m = -1
d = 3**0.5
a, b = get_aba(m,d, yx)
a2 = yx(a)
fx = 1/(b-a)

X, Y = get_XY(a, b, n, check=False)
XX, YY, table1 = func1(Y, n)
XX2, YY2, table2 = func2(Y, n)
print('равноинтервальный метод:')
print(tabulate(table1, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))
print('равновероятный метод:')
print(tabulate(table2, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))


Gy = get_XY_teor(fi2, fx, a2)
vals = np.linspace(a2, yx(b), 100000)
yyy = [Gy(x) for x in vals]

plots(XX, YY, XX2, YY2, vals, yyy,   Y, table1, table2)



#group2 = group[1:]
#table = [['X']+[t[0] for t in group2], ['n']+[t[1] for t in group2]]
#print(tabulate(table, tablefmt='fancy_grid'))

#rnd = lambda x: round(x,2)
#yyy = [Gy(x) for x in XX]
#print(tabulate([['Xi']+list(map(rnd,XX)),['Y1i']+list(map(rnd,YY)), ['Y2i']+list(map(rnd,yyy))]))
#print(tabulate([list(map(rnd,XX)),list(map(rnd,yyy)) ], tablefmt='fancy_grid'))

















