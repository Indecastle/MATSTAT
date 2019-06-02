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


def get_XY(a, b, n):
    #Ei = lambda i: i/(n-1)
    Ei = lambda: np.random.uniform()
    Xi = lambda i: a + Ei()*(b-a)
    X = [Xi(i) for i in range(n)]
    #print('X:', X)
    Y = [yx(x) for x in X]
    Y.sort()
    #print('Y:', Y)
    return X, Y


def func1(X,Y):
    group =  [[key, len(list(group))] for key, group in groupby(Y)]
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
    return XX, YY, group


def get_XY_teor(fy, fx, fa):
    expr = sm.Eq(y, fy)
    xx = sm.solve(expr,x)
    xx0 = xx[0]
    print('xx:', xx)
    #diff_xx = sm.diff(xx[0], y, 1)
    #print('diff_xx:', diff_xx)
    #Gy = fx * sm.integrate(diff_xx, (y, fa, y))
    Gy = fx * (xx0 - xx0.subs({y:fa}))
    #=====
    #Gy = 0.5*sm.real_root(y,3) - 0.5*sm.real_root(-1,3)
    #Gy = 0.5*sm.sqrt(y)
    #=====
    print('Gy:', Gy)
    Gy = sm.lambdify(y, Gy)
    return Gy


fi = sm.real_root(x,3)
fi2 = x**(1/3)
yx = sm.lambdify(x, fi)
yx2 = sm.lambdify(x, fi2)

n = 500
m = -1
d = 3**0.5
a, b = get_aba(m,d, yx)
a2 = yx(a)
fx = 1/(b-a)

X, Y = get_XY(a, b, n)
XX, YY, group = func1(X,Y)
plt.plot(XX, YY, label='Адамс')

Gy = get_XY_teor(fi2, fx, a2)

vals = np.linspace(a2, yx(b), 100000)
yyy = [Gy(x) for x in vals]

plt.plot(vals, yyy)
plt.grid()
plt.show()



group2 = group[1:]
table = [['X']+[t[0] for t in group2], ['n']+[t[1] for t in group2]]
print(tabulate(table, tablefmt='fancy_grid'))

#rnd = lambda x: round(x,2)
#yyy = [Gy(x) for x in XX]
#print(tabulate([['Xi']+list(map(rnd,XX)),['Y1i']+list(map(rnd,YY)), ['Y2i']+list(map(rnd,yyy))]))
#print(tabulate([list(map(rnd,XX)),list(map(rnd,yyy)) ], tablefmt='fancy_grid'))

















