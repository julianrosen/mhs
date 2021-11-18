from .sage_compatibility import Int, Rat, binomial as bi, Number
from mhs import MhsAlg, aperybp
from .lifts import P, Hp
from .global_vars import max_weight
from sympy import var, Symbol, poly, S as SS
from .qsym import qsym
from .rational import MHS, MHS_approx, red, v
from .misc import psc, psc_

import sympy.polys.partfrac as pf
import sys

Inf = None

unit = (0,Inf,0,())
# (factors of p, location of pole, order (order is non-neg, and if 0, location must be Inf)), H_n(_)

def weight_s(L,small=False):
    # small is True if we can only plug in values smaller than p
    if small:
        if L[1] is Inf or L[1] >= 0:
            return L[0]
        else:
            return L[0] - L[2]
    else:
        if L[1] is Inf:
            return L[0] - maxz(L[3])
        elif L[1] > 0:
            return L[0] - L[2] - maxz(L[3])
        else:
            return L[0] - maxz((L[2],)+L[3])

def mult(L1,L2):
    var('x')
    q = (qsym(L1[3])*qsym(L2[3])).data
    #D = my_part(un_part(L1[1:3]) * un_part(L2[1:3]))
    D = mult_part(L1[1:3],L2[1:3])
    T = {}
    for a in q:
        for d in D:
            T[(L1[0]+L2[0],d[0],d[1],a)] = q[a] * D[d]
    return T

class SumAlg():
    def __init__(self, s=0, err=99999):
        if isinstance(s,Number):
            self.data = {}
            if s != 0:
                self.data[unit] = s
            self.err = err
            self.small = False
            self.desc = str(s)
        elif isinstance(s,SumAlg):
            self.data = dict(s.data)
            self.err = min(s.err, err)
            self.small = s.small
            self.desc = s.desc
        elif isinstance(s,dict):
            self.data = dict(s)
            self.err = err
            self.small = False
            self.desc = r'\text{A sequence depending on $p$}'
        elif isinstance(s,type(unit)):
            if s[2] == 0:
                self.data = {(s[0],Inf,0,s[3]):Rat(1)}
            else:
                self.data = {s: Rat(1)}
            self.err = err
            self.small = False
            self.desc = r'\text{A sequence depending on $p$}'
        else:
            print((type(s)))
            raise ValueError('Cannot construct instance from this data type')

            
    def __repr__(self):
        if self.desc is not None:
            self.disp()
        return r'\text{A sequence depending on $p$}'
    
    def disp(self):
        """Displays a human-readable description, if one is available"""
        if self.desc is None:
            print(r'\text{A sequence depending on $p$}')
        else:
            display(Math(self.desc))
        return None
    
    def clean(self):
        D = dict(self.data)
        for x in D:
            if D[x] == 0 or weight_s(x,self.small) >= self.err:
                self.data.pop(x)
        return self

    def __add__(self, s):
        T = SumAlg(s)
        T.small = self.small or T.small
        T.err = min(T.err, self.err)
        for x in self.data:
            if self.data[x] != 0 and weight_s(x,T.small) < T.err:
                if x in T.data:
                    T.data[x] += self.data[x]
                else:
                    T.data[x] = self.data[x]
        return T

    def __mul__(self, s):
        if isinstance(s,SumAlg):
            T = SumAlg()
            T.small = self.small or s.small
            T.err = min(self.err+s.v(),s.err+self.v())
            for x in self.data:
                for y in s.data:
                    c = self.data[x]*s.data[y]
                    if c != 0 and weight_s(x,T.small)+weight_s(y,T.small)<T.err:
                        T += c*SumAlg(mult(x,y))
            return T
        else:
            try:
                s = Rat(s)
            except:
                try:
                    s = Rat(int(s))
                except:
                    raise ValueError(s)
            if s == 0:
                return SumAlg()
            T = SumAlg(self)
            T.clean()
            for x in T.data:
                T.data[x] *= s
            return T

    def __sub__(self,s):
        return self + Rat(-1)*s

    def __radd__(self,s):
        return self + s

    def __rmul__(self,s):
        return self*s

    def __rsub__(self,s):
        return Rat(-1)*self + s

    def __pow__(self,n):
        if n == 0:
            return SumAlg(1)
        return self*self**(n-1)

    def v(self):
        err = self.err
        self.clean()
        for x in self.data:
            err = min(err,weight_s(x,self.small))
        return err
    
    def disp(self):
        print((self.data))
        return None
    
    def old__repr__(self):
        self.clean()
        return str(self.data)
        T = ""
        p = Symbol('p', commutative=False)
        n = Symbol('n', commutative=False)
        
        for x in self.data:
            T += str(self.data[x]) + "\t:\t"
            s = p**x[0]
            if x[1] is Inf:
                s *= n**x[2]
            else:
                s *= (n-x[1])**(-x[2])
            T += latex(s)
            if x[3] != ():
                T += "H_n" + tstr(x[3])
        return T
    
    def inv(self):
        self.clean()
        vv = self.v()
        q = 0
        c = 1
        for x in self.data:
            if weight_s(x,self.small) == vv:
                q += 1
                if q >= 2 or x[3] != ():
                    raise ValueError('Element not invertible')
                if x[1] is Inf:
                    c = Rat(1)/self.data[x] * PN(-x[0],-x[2])
                else:
                    c = Rat(1)/self.data[x] * PN(-x[0],0)
                    c *= (PN(0,1) - x[1]) ** x[2]
        y = SumAlg(self) * SumAlg(c)
        T = SumAlg(1) - SumAlg(y)
        vv = T.v()
        if vv < 1:
            print(("c: ", c.data))
            print(("self: ", self.data))
            print(("y: ", y.data))
            print(("T: ", T.data))
            print("Element is not invertible")
            return None
        A = SumAlg(1)
        B = SumAlg(1)
        A.err = T.err
        B.err = T.err
        for n in range(min(A.err // vv + 1, max_weight()+1)):
            B *= T
            A += B
        A.err = min(A.err, max_weight())
        return A * c
    
    def is_inv(self):
        self.clean()
        vv = self.v()
        q = 0
        for x in self.data:
            if weight_s(x,self.small) == vv:
                q += 1
                if q >= 2 or x[2] != ():
                    return False
        return True
    
    def __div__(self, other):
        if isinstance(other,Number):
            return self*(Rat(1)/other)
        else:
            return self * other.inv()
    
    def __rdiv__(self, other):
        return self.inv() * other
    
    def ev(self,n,p):
        return sum(self.data[s]*Rat(n)**(s[1])*\
                   Rat(p)**(s[0])*\
                   MHS_approx(p, 14, n-1,*(s[2])) for s in self.data)
    
    def aerr(self,err):
        self.err = min(self.err,err)
        return None
    
    def e_p(self):
        self.clean()
        T = MhsAlg()
        for t in self.data:
            c = self.data[t]
            if t[1] is Inf:
                T += c * P(t[0] + t[2]) * Hp(*t[3])
            else:
                T += c * P(t[0]) * (P(1) - t[1])**(-t[2]) * Hp(*t[3])
        T.aerr(self.err)
        return T
    
    def e(self,n,p=None,err=None):
        self.clean()
        if err is None:
            E = self.err
        else:
            E = min(self.err,err)
        T = 0
        if p is None:
            for x in self.data:
                if x[0] != 0:
                    raise ValueError("This quantity depends on p")
            for x in self.data:
                c = self.data[x]
                if x[1] is None:
                    T += c * Rat(n)**x[2] * MHS(n-1,*x[3])
                else:
                    T += c * Rat(n-x[1])**(-x[2]) * MHS(n-1,*x[3])
        elif p is not None:
            if E > 9000:
                for x in self.data:
                    c = self.data[x]
                    P = p**x[0]
                    if x[1] is None:
                        T += c * P * Rat(n)**x[2] * MHS(n-1,*x[3])
                    else:
                        T += c * P * Rat(n-x[1])**(-x[2]) * MHS(n-1,*x[3])
            else: # E <= 9000
                for x in self.data:
                    c = self.data[x]
                    P = p**x[0]
                    if x[1] is None:
                        T += c * P * Rat(n)**x[2] * MHS_approx(p,E-x[0],n-1,*x[3])
                    else:
                        T += c * P * Rat(n-x[1])**(-x[2]) * MHS_approx(p,E-x[0],n-1,*x[3])
        return T
    
    def add_one(self):
        T = SumAlg()
        for y in self.data:
            D = add(y[1:3])
            for d in D:
                c = D[d] * self.data[y]
                T += SumAlg({(y[0],d[0],d[1],y[3]):c})
                if len(y[3]) >= 1:
                    T += SumAlg({(y[0],d[0],d[1],y[3][1:]):c}) * SumAlg((0,0,y[3][0],()))
        T.err = self.err
        return T
    
    def add_k(self,k):
        T = SumAlg(self)
        if k > 0:
            for _ in range(k):
                T = T.add_one()
            return T
        else:
            return self.sub_k(-k)
    
    def sub_one_old(self):
        T = SumAlg()
        for y in self.data:
            D = add(y[1:3],-1)
            for d in D:
                c = D[d] * self.data[y]
                T += SumAlg({(y[0],d[0],d[1],y[3]):c})
                if len(y[3]) >= 1:
                    T -= SumAlg((0,Inf,0,y[3][1:])).sub_one() * SumAlg({(y[0],d[0],d[1],()):c})\
                        * SumAlg((0,1,y[3][0],()))
        return T
    
    def sub_one(self):
        T = SumAlg(0)
        for y in self.data:
            D = add(y[1:3],-1)
            Q = SumAlg({(y[0],w[0],w[1],()):self.data[y] * D[w] for w in D})
            s = y[3]
            QQ = SumAlg({(0,1,sum(s[:i]),s[i:]):(-1)**i for i in range(len(s)+1)})
            T += Q * QQ
        T.err = self.err
        return T
    
    def sub_k(self,k):
        T = SumAlg(self)
        if k >= 0:
            for _ in range(k):
                T = T.sub_one()
            return T
        else:
            return self.add_k(-k)
    
    def sum(self, i=1, w=0, alert=False):
        T = SumAlg(0)
        assert w <= 1 or not self.small
        if alert:
            print(("There are %i entries" % len(self.data)))
            A = 0
        for x in self.data:
            if alert:
                print(("Starting entry %i" %A))
                print(("x: ", x))
                print("")
                sys.stdout.flush()
                A += 1
            if weight_s(x,self.small) >= T.err:
                continue
            c = self.data[x]
            if x[1] is Inf or x[1] <= 0:
                if x[1] is Inf:
                    Q = H(-x[2],*x[3])
                else:
                    Q = H_geq(1-x[1],x[2],*x[3]).add_k(-x[1])
                if i >= 2:
                    Q -= Q.e(i)        
                if i <= 0:
                    Q -= Q.e(i)
                if w != 0:
                    Q = Q.add_k(w)
                T += c * PN(x[0],0) * Q
            else: # x[1] > 0
                Q = SumAlg(0)
                Q.data[x] = c
                Q = Q.add_k(x[1])
                T += Q.sum(i-x[1], w-x[1])
        T.aerr(self.err)
        T.small = self.small
        return T


def Sm(*s):
    """Sum_{n-1 >= n_1 > ... > n_k >= 1} n_1^{-s_1} ... n_k ^ {-s_k}"""
    if len(s) == 1 and isinstance(s[0], tuple):
        ss = tuple(s[0])
    else:
        ss = tuple(s)
    if len(ss) == 0:
        return SumAlg(1)
    elif len(ss) == 1:
        if ss[0] <= 0:
            T = sum(SumAlg({(0,Inf,j,()):psc_(-ss[0], j)}) for j in range(2 - ss[0]))
            return T
        else:
            return SumAlg({(0,Inf,0,ss):1})
    if ss[0] <= 0:
        return sum((SumAlg({(0,Inf,j,()):psc_(-ss[0], j)}) * Sm(ss[1:])\
                    if psc_(-ss[0], j) != 0 else 0) - (psc(-ss[0], j) *\
                    Sm((ss[1] - j,) + ss[2:]) if psc(-ss[0], j) != 0 else 0)\
                   for j in range(2 - ss[0]))
    if ss[-1] <= 0:
        return sum(psc_(-ss[-1], j) * Sm(ss[:-2] + (ss[-2] - j,))\
                   if psc_(-ss[-1], j) != 0 else 0\
                   for j in range(2 - ss[-1]))
    for i in range(1, len(ss) - 1):
        if ss[i] <= 0:
            return sum((psc_(-ss[i], j) * Sm(ss[:i - 1] + (ss[i - 1] - j,)\
                        + ss[i + 1:]) if psc_(-ss[i], j) != 0 else 0)\
                       - (psc(-ss[i], j) * Sm(ss[:i] + (ss[i + 1] - j,) + ss[i+2:]\
                        ) if psc(-ss[i], j) != 0 else 0) for j in range(2 - ss[i]))
    for q in ss:
        assert q>0
    return SumAlg({(0,Inf,0,ss):1})

def S(*s):
    """Sum_{n >= n_1 > ... > n_k >= 1} n_1^{-s_1} ... n_k ^ {-s_k}"""
    if len(s) == 0:
        return SumAlg(1)
    return Sm(*s) + PN(0,-s[0])*Sm(*s[1:])

def PN(a=0, b=0):
    if b >= 0:
        return SumAlg({(a,Inf,b,()):1})
    else:
        return SumAlg({(a,0,-b,()):1})

def ff(x, n):
    if n == 0:
        return 1
    return x * ff(x-1,n-1)

def roots_among_integers(f):
    var('x')
    Q = poly(f,x).all_roots()
    for r in Q:
        if not r.is_integer:
            return False
    return True

Inf = None
def my_part(F):
    x = var('x')
    f = SS(F)
    if f.is_Atom:
        if f is x:
            return {(Inf,1):1}
        else:
            return {(Inf,0):f}
    D = {}
    Q = pf.apart_list(f)
    for i,x in enumerate(Q[1].all_coeffs()[::-1]):
        if x != 0:
            D[(Inf,i)] = x/Q[0]
    for T in Q[2]:
        for n in T[0].all_roots():
            D[(n,T[3])] = T[1](n)/Q[0]
    return D

def un_part(D):
    var('x')
    T = 0
    if isinstance(D,dict):
        for a,b in D:
            if a is Inf:
                T += D[(a,b)] * x**b
            else:
                T += D[(a,b)] / (x-a)**b
        return T
    elif isinstance(D,(list,tuple)):
        a,b = D[0],D[1]
        if a is Inf:
            T += x**b
        else:
            T += 1 / (x-a)**b
        return T

def add(D,k=1):
    if isinstance(D,tuple):
        return add({D:1},k)
    elif isinstance(D,dict):
        T = {}
        for x in D:
            if x[0] is Inf:
                for j in range(x[1]+1):
                    add_entry(T,(Inf,j),bi(x[1],j)*k**(x[1]-j)*D[x])
            else:
                add_entry(T,(x[0]-k,x[1]),D[x])
        return T

def from_rat(f):
    n,d = f.as_numer_denom()
    if roots_among_integers(d):
        D = my_part(f)
        T = SumAlg()
        for a,b in D:
            T.data[(0,a,b,())] = D[(a,b)]
        return T
    else:
        print(("n,d: ",n,d))
        raise ValueError('Pole at non-integral value')
        
def maxz(L):
    return max(L) if len(L) > 0 else 0



def H(*s):
    return Sm(*s)

def H_eq(k,*s):
    assert len(s)>=2
    var('x')
    T = SumAlg()
    D = my_part(x**(-s[0]) * (x-k)**(-s[1]))
    for i in D:
        c = D[i]
        assert i[0] in (None, 0, k)
        if i[0] is None:
            T += c * H_geq(k+1,-i[1],*s[2:])
        elif i[0] == 0:
            T += c * H_geq(k+1,i[1],*s[2:])
        else:
            T += c * H(i[1],*s[2:]).sub_k(k)
    return T

def H_geq(k,*s):
    assert len(s) >= 1
    T = H(*s)
    if len(s) == 1:
        return T - MHS(k-1,*s)
    else:
        for i in range(1,k):
            T -= H_eq(i,*s)
    return T

def BIN(a,b,k=1):
    """ The binomial coefficient {kp+n+a \choose n + b}"""
    if a == b:
        T = sum(k**i*PN(i,0)*H(*i*(1,)).add_k(a+1) for i in range(max_weight() + 2))
        T.err = max_weight()
        return T
    elif a > b:
        T = BIN(a,a,k)
        for i in range(b+1,a+1):
            T *= PN(0,1) + i
        for i in range(1,a-b+1):
            T /= k*PN(1,0) + i
        T.err = max_weight()
        return T
    elif a < b:
        T = BIN(a,a,k)
        T.err += b - a
        for i in range(a+1,b+1):
            T *= SumAlg((0,-i,1,()))
        for i in range(b-a):
            T *= k*PN(1,0) - i
        T.err = max_weight()
        return T
    
def nn(i,p=1):
    if p > 0:
        return SumAlg((0,i,p,()))
    elif p == 0:
        return SumAlg(1)
    elif p < 0:
        return (PN(0,1) - i) ** (-p)

def add_entry(D,x,c=1):
    if x in D:
        D[x] += c
    else:
        D[x] = c
    return None

def mult_part(L1,L2):
    a1,b1 = L1
    a2,b2 = L2
    if a1 is Inf and a2 is Inf:
        return {(Inf,b1+b2):1}
    elif a1 is Inf:
        D = {(a2,b2-j):bi(b1,j)*a2**(b1-j) for j in range(0,min(b1+1,b2))}
        for j in range(b2,b1+1):
            for i in range(j-b2+1):
                add_entry(D,(Inf,i),bi(b1,j)*a2**(b1-j)*bi(j-b2,i)*(-a2)**(j-b2-i))
        return D
    elif a2 is Inf:
        return mult_part(L2,L1)
    elif a1 == a2:
        return {(a1,b1+b2):1}
    else:
        D = {}
        for j in range(b2):
            D[(a2,b2-j)] = bi(-b1,j) * Rat(a2-a1)**(-b1-j)
        for j in range(b1):
            D[(a1,b1-j)] = bi(-b2,j) * Rat(a1-a2)**(-b2-j)
        return D
            
        
def test(a,b):
    diff = 99999
    fail = False
    for p in [17,19,23]:
        diff = min(diff,v(p,a.e(p)-b(p)) - a.err)
        if diff < 0:
            fail = True
            break
    ds = str(diff)
    if diff >= 0:
        ds = '+' + ds
    print(("Supposed to be:\t %i"%a.err))
    print(("Difference:\t %s"%ds))
    if fail:
        print("Failed!!!")
    else:
        print("Passed")

def BINN(a,b,k=1):
    """The binomial coefficient (-1)^n {kp+a\choose n+b}"""
    return Rat(-1)**(b%2) * BIN(b-a-1,b,-k)

def z3_altp():
    a = sum(PN(t-1,-2-t) for t in range(max_weight()))
    a.small = True
    b = sum((-1)**(t+1)*SumAlg((2*t,Inf,0,t*(2,))) for t in range(1,max_weight()//2))
    b.small = True
    c = sum(b**t for t in range(max_weight()//2))
    T = Hp(3) - (a*c).sum().e_p()/2
    T.desc = r"\frac{5}{2}\sum_{n=1}^{p-1}\frac{(-1)^{n-1}}{n^3 {2n\choose n}}"
    return T

def aperyap():
    a = sum(PN(t-1,-2-t) for t in range(max_weight()))
    a.small = True
    b = sum((-1)**(t+1)*SumAlg((2*t,Inf,0,t*(2,))) for t in range(1,max_weight()//2))
    b.small = True
    c = sum(b**t for t in range(max_weight()//2))
    d = (a*c).sum(1,1)*(Rat(1)/2)
    x = BIN(-1,0)
    y = BINN(-1,0)
    x.small = True
    y.small = True
    T = Hp(3)*aperybp() - (d*x**2*y**2).sum().e_p()
    T.desc = r"Apery number a_{p-1}"
    return T

def aperyp():
    return aperyap() / aperybp()
