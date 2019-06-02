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


def get_XY(a, b, yx, n, check=False):
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


def func2(Y, M, n):
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


def func1(Y, n):
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


def func2_1(table, Gy, n):
    Xak = 11.34  #  a = 0.01, k = 3
    FA = [Gy(A) for A in table['A']]
    FB = [Gy(B) for B in table['B']]
    pi = [fb-fa for fa,fb in zip(FA, FB)]
    pi2 = [v/n for v in table['n']]
    Xi = [n*(p - p2)**2/p for p,p2 in zip(pi, pi2)]
    Xi2 = sum(Xi)
    #test = sum([(v - n*p)**2/(n*p) for p,v in zip(pi, table['n'])])
    #print("X^2:", Xi2)
    if Xi2 < Xak:
        print(f'X^2({round(Xi2, 3)}) < Xak({Xak}) - гипотеза не протеворечит имеющимся данным')
    else:
        print(f'X^2({round(Xi2, 3)}) > Xak({Xak}) - гипотеза протеворечит')
    table3 = {'FA':FA, 'FB':FB, 'pi':pi, 'pi*':pi2, 'Xu':Xi}
    return table3


def func2_2(XX, YY, Gy, n):
    diff = []
    maxdiff = [0, 0, 0]
    for i in range(1, len(XX)-1):
        df = abs(Gy(XX[i])-YY[i])
        diff.append(df)
        if df > maxdiff[0]:
            maxdiff = [df, XX[i], YY[i]]
    fact = math.sqrt(n) * maxdiff[0]
    crit = 1.64  #  a = 0.01
    print("max_diff:", maxdiff[0])
    if fact < crit:
        print(f'fact({round(fact, 3)}) < crit({crit}) - гипотеза не протеворечит имеющимся данным')
    else:
        print(f'fact({round(fact, 3)}) > crit({crit}) - гипотеза протеворечит')
    return maxdiff


def func2_3(Y, Gy, n):
    I = [i for i in range(1,n+1)]
    Fx = [Gy(y) for y in Y]
    Fnx = [(i-0.5)/n for i in range(n)]
    D = [(f1-f2)**2 for f1,f2 in zip(Fnx, Fx)]
    factM = 1/(12*n) + sum(D)
    kritM = 0.744
    if factM < kritM:
        print(f'factM({round(factM, 3)}) < kritM({kritM}) - гипотеза не протеворечит имеющимся данным')
    else:
        print(f'factM({round(factM, 3)}) > kritM({kritM}) - гипотеза протеворечит')
    table = {'i':I, 'x':Y, 'Fx':Fx, 'Fnx':Fnx, 'D':D}
    return table


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


def plots(X1, Y1, gy,   X2, Y2, Gy, maxdiff,   yx, a2, b):
    #M = int(math.sqrt(n))
    vals = np.linspace(a2, yx(b), 100000)
    
    fig, (ax1, ax2) = plt.subplots(
        nrows=1, ncols=2,
        figsize=(8, 4)
    )
    AX = [ax1, ax2]
    ax1.plot(X1,Y1, label='эмпирическая плотность')
    yyy = [gy(x) for x in vals]
    ax1.plot(vals, yyy, label='теоретическая плотность')

    ax2.plot(X2,Y2, label='эмпирическая функция')
    yyy = [Gy(x) for x in vals]
    ax2.plot(vals, yyy, label='теоретическая функция')
    ax2.plot([0], [0], 'r', label=f"max_diff: {round(maxdiff[0], 3)}")  #  only for legend
    arrowprops = { 'arrowstyle': '<->', 'color':'red'} #'connectionstyle':'arc3,rad=-0.5'
    ax2.annotate('', # f'diff: {round(maxdiff[0],3)}'
                 xy=(maxdiff[1], maxdiff[2]),
                 xytext = (maxdiff[1], Gy(maxdiff[1])),
                 arrowprops = arrowprops,
                 ha='center', va="center", size=15)

    ax1.set_title('равновероятный метод')
    ax2.set_title('Эмпирическая функция')
    for ax in AX:
        ax.grid()
        ax.set_xlabel('Y')
        ax.set_ylabel('f*y')
    ax1.legend(loc='upper center', shadow=True, fontsize='medium', frameon=False)
    ax2.legend(loc='upper center', shadow=True, fontsize='small', frameon=True)
    fig.tight_layout()
    #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()


def main():
    fi = sm.real_root(x,3)
    fi2 = x**(1/3)
    yx = sm.lambdify(x, fi)
    yx2 = sm.lambdify(x, fi2)

    n = 200
    if n <= 100:
        M = int(math.sqrt(n))
    else:
        M = int(4*math.log10(n))
    m = -1
    d = 3**0.5
    a, b = get_aba(m,d, yx)
    a2 = yx(a)
    fx = 1/(b-a)

    X, Y = get_XY(a, b, yx, n, check=False)

    XX, YY, table = func2(Y, M, n)
    print('равновероятный метод:')
    print(tabulate(table, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))

    gy = get_XY_teor(fi2, fx, a2)
    Gy = get_XY_teor(fi2, fx, a2, func=True)

    print('Хи-квадрат Пирсона:')
    ##  Хи-квадрат Пирсона
    table22 = func2_1(table, Gy, n)
    print(tabulate(table22, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))

    ##  Колмогорова
    print("Колмогорова:")
    n = 30
    X, Y = get_XY(a, b, yx, n, check=False)
    XX2, YY2, _ = func1(Y, n)
    maxdiff = func2_2(XX2, YY2, Gy, n)

    ##  Мизеса
    print("Мизеса:")
    n = 50
    X, Y = get_XY(a, b, yx, n, check=False)
    table23 = func2_3(Y, Gy, n)
    print(tabulate(table23, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))

    plots(XX, YY, gy,   XX2, YY2, Gy, maxdiff, yx, a2, b)






if __name__ == "__main__":
    main()












