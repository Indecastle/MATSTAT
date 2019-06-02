import math
import numpy as np
import sympy as sm
from itertools import groupby
import matplotlib.pyplot as plt
from tabulate import tabulate

#from table_t import *

from scipy.stats import t
from scipy.stats import chi2
from scipy.special import erf
from scipy.stats import norm


x, y = sm.symbols('x, y', real = True)

n_list = np.arange(10, 150, 1)
n_list2 = [20, 30, 50, 70, 100, 150]
a_t = np.arange(0.001, 1, 0.001)
a_t2 = [0.80, 0.90, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999]


def get_aba(m,d, yx):
    b = m - 3**0.5 * d
    a = 2*m - b
    #a, b = -2, 2
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
    Y = [yx(x) for x in X]
    Y.sort()
    #print('Y:', Y)
    return X, Y


def tochn_m(X):
  return sum(X) / len(X)


def tochn_d(X, xs):
  n = len(X)
  return 1/(n-1) * sum([(x-xs)**2 for x in X])


"""
def get_inter_m(xs, ds, tyk, n, inter=False):
    if inter:
        return (xs - ds*tyk/(n-1), xs + ds*tyk/(n-1))
    else:
        return (xs + ds*tyk/(n-1)) - (xs - ds*tyk/(n-1))

def get_inter_d(ds2, tyk, n, inter=False):
    if inter:
        return (n*ds2/tyk[0], n*ds2/tyk[1])
    else:
        return n*ds2/tyk[0] - n*ds2/tyk[1]
"""


def print_inter_m():
    pass


def get_inter_m(md, md_teor, a_list, n_list, inter=False, teor=False):
    inter_m = []
    inter_m_teor = []
    bol = callable(md)
    if not bol:
        m = md[0]
        ds = md[1]**0.5
    for a,n in zip(a_list, n_list):
        if bol:
            _,Y = md(n)
            m = tochn_m(Y)
            ds2 = tochn_d(Y, m)
            ds = math.sqrt(ds2)
            
        u = ds * norm.ppf(1-a/2.0) / n**0.5
        u2 = md_teor[1]**0.5 * t.ppf(1 - a / 2, n - 1) / n**0.5
            
        if inter:
            inter_m.append( (m - u, m + u) )
            inter_m_teor.append( (m - u2, m + u2) )
        else:
            inter_m.append( 2 * u )
            inter_m_teor.append( 2 * u2 )
    return inter_m, inter_m_teor


def get_inter_d(md, md_teor, a_list, n_list, inter=False):
    inter_d = []
    inter_d_teor = []
    bol = callable(md)
    if not bol:
        m = md[0]
        ds2 = md[1]
        ds2_teor = md_teor[1]
    for a,n in zip(a_list, n_list):
        if bol:
            _,Y = md(n)
            m = tochn_m(Y)
            ds2 = tochn_d(Y, m)
            ds2_teor = tochn_d(Y, md_teor[0])
        left  = chi2.isf((1 + (1 - a)) / 2, n - 1)
        right = chi2.isf((1 - (1 - a)) / 2, n - 1)
        if inter:
            inter_d.append( ((n-1)*ds2/left, (n-1)*ds2/right) )
            inter_d_teor.append( ((n-1)*ds2_teor/left, (n-1)*ds2_teor/right) )
        else:
            inter_d.append( (n-1)*ds2/left - (n-1)*ds2/right )
            inter_d_teor.append( (n-1)*ds2_teor/left - (n-1)*ds2_teor/right )
    return inter_d, inter_d_teor


def get_XY_teor(fy, fx, fa, func=False):
    expr = sm.Eq(y, fy)
    xx = sm.solve(expr,x)
    xx0 = xx[0]
    #print('xx:', xx)
    if func:
        Gy = fx * (xx0 - xx0.subs({y:fa}))
    else:
        diff_xx = sm.diff(xx0, y, 1)
        Gy = fx * diff_xx
    #=====
    #Gy = 0.5*sm.real_root(y,3) - 0.5*sm.real_root(-1,3)
    #Gy = 0.5*sm.sqrt(y)
    #Gy = 0.125/sm.sqrt(y) * 2
    #=====
    #print('Gy:', Gy)
    #Gy = sm.lambdify(y, Gy)
    return Gy


def get_md_teor(gy, fa, fb, D=False):
    if D:
        m_teor = sm.integrate(y*gy, (y, fa, fb))
        result = sm.integrate(y*y*gy, (y, fa, fb)) - m_teor**2
    else:
        result = sm.integrate(y*gy, (y, fa, fb))
    #print("privet",sm.integrate(y*y*gy, y).subs({y:36}))
    return result


def plots(inter_m, inter_m_teor, a_t, inter2_m, inter2_m_teor, inter_d, inter_d_teor, inter2_d, inter2_d_teor):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows=2, ncols=2,
        figsize=(8, 8)
    )
    AX = [ax1, ax2, ax3, ax4]

    ax1.plot(1-a_t, inter_m, 'b')
    ax1.plot(1-a_t, inter_m_teor, 'r')
    ax1.plot([], [], 'b', label='факт.')
    ax1.plot([], [], 'r', label='теор.')

    ax3.plot(n_list, inter2_m, 'b')
    ax3.plot(n_list, inter2_m_teor, 'r')
    ax3.plot([], [], 'b', label='факт.')
    ax3.plot([], [], 'r', label='теор.')

    ax2.plot(1-a_t, inter_d, 'b')
    ax2.plot(1-a_t, inter_d_teor, 'r')
    ax2.plot([], [], 'b', label='факт.')
    ax2.plot([], [], 'r', label='теор.')

    ax4.plot(n_list, inter2_d, 'b')
    ax4.plot(n_list, inter2_d_teor, 'r')
    ax4.plot([], [], 'b', label='факт.')
    ax4.plot([], [], 'r', label='теор.')
    
    ax1.set_title('для мат. ожидания от уровня. знач.')
    ax3.set_title('для мат. ожидания от объема с дов. знач. %.2f' % (1-0.01))
    ax2.set_title('для дисперсии от уровня. знач.')
    ax4.set_title('для дисперсии от объема с дов. знач %.2f' % (1-0.01))
    
    ax1.set_xlabel('a'); ax1.set_ylabel('величины дов. интервалов');
    ax3.set_xlabel('n'); ax3.set_ylabel('величины дов. интервалов');
    ax2.set_xlabel('a'); ax3.set_ylabel('величины дов. интервалов');
    ax4.set_xlabel('n'); ax3.set_ylabel('величины дов. интервалов');
    
    ax1.legend(loc='upper center', shadow=True, fontsize='medium', frameon=False)
    ax3.legend(loc='upper center', shadow=True, fontsize='small', frameon=True)
    ax2.legend(loc='upper center', shadow=True, fontsize='medium', frameon=False)
    ax4.legend(loc='upper center', shadow=True, fontsize='small', frameon=True)
    #fig.tight_layout()
    #fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()



def main():
    fi = sm.real_root(x,3)
    fi2 = x**(1/3)
    yx = sm.lambdify(x, fi)
    yx2 = sm.lambdify(x, fi2)

    n = 20
    m = -1
    d = 3**0.5
    a, b = get_aba(m,d, yx)
    a2, b2 = yx(a), yx(b)
    fx = 1/(b-a)

    X, Y = get_XY(a, b, yx, n, check=False)
    func = lambda n: get_XY(a, b, yx, n, check=False)

    xs = tochn_m(Y)
    ds2 = tochn_d(Y, xs)
    #print(xs, ds2, ds)
    gy = get_XY_teor(fi2, fx, a2)
    xs_teor = get_md_teor(gy, a2, b2, D=False)
    ds2_teor = get_md_teor(gy, a2, b2, D=True)
    print("ds и ds_teor:", ds2, ds2_teor)

    print_inter_m()
    inter_m, inter_m_teor = get_inter_m((xs, ds2), (xs_teor, ds2_teor), a_t, [n]*len(a_t))
    #print(list(map(lambda x: round(x,3), inter_m)))
    #inter_m_teor = get_inter_m(xs, ds_teor, a_t, [n]*len(a_t), teor=True)
    #print(list(map(lambda x: [round(x[0], 3), round(x[1], 3)], inter_m_teor)))

    inter2_m, inter2_m_teor = get_inter_m(func, (xs_teor, ds2_teor), [0.01]*len(n_list), n_list)
    #inter2_m_teor = get_inter_m(func, ds_teor, [0.01]*len(n_list), n_list)
    

    #inter_d = [get_inter_d(ds2, a, n) for a in xi_t]
    inter_d, inter_d_teor = get_inter_d((xs, ds2), (xs_teor, ds2_teor), a_t, [n]*len(a_t))
    #inter_d_teor = get_inter_d(ds2_teor, a_t, [n]*len(a_t))
    #print(list(map(lambda x: round(x,3), inter_d)))

    inter2_d, inter2_d_teor = get_inter_d(func, (xs_teor, ds2_teor), [0.01]*len(n_list), n_list)
    #inter2_d_teor = get_inter_d(ds2_teor, [0.01]*len(n_list), n_list)


    #inter_m_2 = [get_inter_m(xs, ds, t, n, inter=True) for t in table_t[n]]
    #table1 = {"":table_t[20]}
    #print(tabulate(table1, headers='keys', floatfmt=".2f", tablefmt='fancy_grid'))





    plots(inter_m, inter_m_teor, a_t, inter2_m, inter2_m_teor, inter_d, inter_d_teor, inter2_d, inter2_d_teor)



if __name__ == "__main__":
    main()
