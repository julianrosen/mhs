from sympy import Sum, var, poly
from .sage_compatibility import Int, Rat, Number, prod, latex, factorial, binomial
from .qsym import qsym
from .global_vars import max_weight, num, mzv_alg, mono_to_basis, mhs_to_mono, mzv_mono, mzv_to_mono, weights, mhs_basis
from .rational import MHS_approx, red
from .misc import psc, psc_, D_add, comps, comps_i, lists, wt
from .tex import error_tex, tex_write, tstr, monomial, un_texify, tex_write_single, fix_pow, fix_prod
from .misc import D_add, wt, bn, decomps, S_min_v
from IPython.display import display, Math


def _unit():
    return (0,tuple(0 for _ in range(len(mzv_alg()))))

def _weight(x):
    return x[0]

def _mul(x,y):
    return {(x[0]+y[0],tuple(x[1][i]+y[1][i] for i in range(len(mzv_alg())))):1}

class MhsAlg():
    def __init__(self, s=0, err=99999):
        if isinstance(s,Number):
            self.data = {}
            if s != 0:
                self.data[_unit()] = s
            self.err = err
            self.desc = latex(s)
        elif isinstance(s,MhsAlg):
            self.data = dict(s.data)
            self.err = min(s.err, err)
            self.desc = s.desc
        elif isinstance(s,dict):
            self.data = dict(s)
            self.err = err
            self.desc = None
        elif isinstance(s,type(_unit())):
            self.data = {s: 1}
            self.err = err
            self.desc = None
        elif isinstance(s,qsym):
            T = sum(s.data[q] * Hp(*q) for q in s.data)
            self.data = T.data
            self.err = T.err
            self.desc = None
        else:
            raise ValueError("Cannot create element of MHS algebra from data type %s"%str(type(s)))
    
    def __repr__(self):
        if self.desc is not None:
            self.disp()
        return '<An element of the MHS algebra>'
    
    def disp(self):
        """Displays a human-readable description, if one is available"""
        if self.desc is None:
            print("<An element of the MHS algebra>")
        else:
            display(Math(self.desc))
        return None
               
    
    def mzv(self, err=None, reverse=False,p='p',mot=False,tex=True):
        """Displays an element in the MZV basis"""
        if err is not None:
            T = MhsAlg(self)
            T.aerr(err)
            return T.mzv(reverse=reverse,p=p,mot=mot)
        self.clean()
        A = tex_write_single(self.data,
                         lambda x:fix_prod(fix_pow(p,x[0]-wt(x[1])),monomial([r"\zeta_p"+tstr(q) for q in mzv_alg()],x[1])),
                             lambda x:(x[0],wt(x[1])),
                         err=self.err,
                         var='p')
        if tex:
            return Math(A)
        else:
            return A
        
    def mzv_old(self, err=None, reverse=False,p='p',mot=False,tex=True):
        """Displays an element in the MZV basis"""
        if err is not None:
            T = MhsAlg(self)
            T.aerr(err)
            return T.mzv(reverse=reverse,p=p,mot=mot)
        self.clean()
        if reverse:
            A = tex_write({(x[0]-wt(x[1]),x[1]):self.data[x] for x in self.data},
                          lambda x:monomial([r"\zeta"+(r"^{\mathfrak{dr}}" if mot else "_p")+tstr(q) for q in mzv_alg()],x[1]),
                         lambda x:(wt(x[1]),) + x[1],
                         lambda x:monomial([p],[_weight(x)]),
                         lambda x:_weight(x),
                         err=self.err,
                         first=not mot,
                         var=r'\mathbb{L}' if mot else 'p')
        else:
            A = tex_write({(x[0]-wt(x[1]),x[1]):self.data[x] for x in self.data},
                         lambda x:monomial([p],[_weight(x)]),
                         lambda x:_weight(x),
                         lambda x:monomial([r"\zeta"+(r"^{\mathfrak{dr}}" if mot else "_p")+tstr(q) for q in mzv_alg()],x[1]),
                         lambda x:(wt(x[1]),) + x[1],
                         err=self.err,
                         first=not mot,
                         var=r'\mathbb{L}' if mot else 'p')
        if tex:
            return Math(A)
        else:
            return un_texify(A)
    
    def mhs(self,err=None):
        """Displays an element in the MHS basis"""
        if err is not None:
            T = MhsAlg(self)
            T.aerr(err)
            return T.mhs()
        self.clean()
        D = {}
        e = self.err
        for x in self.data:
            if wt(x[1]) > max_weight():
                e = min(e,_weight(x))
            else:
                i = mzv_mono().index(x[1])
                for j in range(num()):
                    D_add(D,(x[0]-wt(x[1]),j),self.data[x]*mono_to_basis()[i][j])
                    if i != 0:
                        e = min(e,x[0] + max_weight() + 1)
        K = list(D.keys())
        for x in K:
            if x[0]+weights()[x[1]] >= e:
                D.pop(x)
        T = tex_write(D,
             lambda x:monomial(['p'],[x[0]+weights()[x[1]]]),
             lambda x:x[0]+weights()[x[1]],
             lambda x:"H_{p-1}"+tstr(mhs_basis()[x[1]]) if x[1] != 0 else "1",
             lambda x:(weights()[x[1]],) + (x[0],),
             e,
             True)
        return Math(T)
    
                    
    
    def clean(self):
        K = list(self.data.keys())
        for x in K:
            if self.data[x] == 0 or x[0] >= self.err:
                self.data.pop(x)

    def __add__(self, s):
        if isinstance(s,Number):
            return self + MhsAlg(s)
        elif isinstance(s,MhsAlg):
            T = MhsAlg(s)
            T.err = min(T.err, self.err)
            for x in self.data:
                if self.data[x] != 0 and _weight(x) < T.err:
                    D_add(T.data,x,self.data[x])
            T.desc = None
            return T

    def __mul__(self, s):
        if isinstance(s,Number):
            if s == 0:
                return MhsAlg()
            T = MhsAlg(self)
            T.clean()
            for x in T.data:
                T.data[x] *= s
            T.desc = None
            return T
        elif isinstance(s,MhsAlg):
            T = MhsAlg()
            T.err = min(self.err+s.v(),s.err+self.v())
            for x in self.data:
                for y in s.data:
                    if self.data[x] != 0 and s.data[y] != 0 and _weight(x)+_weight(y)<T.err:
                        T += self.data[x]*s.data[y]*MhsAlg(_mul(x,y))
            T.desc = None
            return T
        print(("Mul just ran off the bottom, ", type(s)))
        raise ValueError("Oops")

    def __sub__(self,s):
        return self + (-1)*s

    def __radd__(self,s):
        return self + s

    def __rmul__(self,s):
        return self*s

    def __rsub__(self,s):
        return (-1)*self + s

    def __pow__(self,n):
        if n == 0:
            return MhsAlg(1)
        elif n > 0:
            return self*self**(n-1)
        elif n < 0:
            return self.inv() ** (-n)

    def v(self):
        """Computes a lower bound for valuation (conjecturally the liminf)"""
        err = self.err
        self.clean()
        for x in self.data:
            err = min(err,x[0])
        return err

    def inv(self):
        """Find the multiplicative inverse, if it exists"""
        self.clean()
        VV = self.v()
        u = len(mzv_alg()) * (0,)
        for x in self.data:
            if x[0] == VV and x[1] != mzv_mono()[0]:
                print("Element is not invertible")
                return None
        c = self.data[(VV,u)]
        T = MhsAlg(1) - self*P(-VV)/c
        A = MhsAlg(1)
        B = MhsAlg(1)
        A.err = self.err
        B.err = self.err
        for n in range(min(A.err + 1,15)):
            B *= T
            A += B
            if B.data == {}:
                break
        if B.data != {}:
            A.aerr(16)
        return A/c*P(-VV)
    
    def __div__(self,s):
        if isinstance(s,Number):
            return self * (Rat(1)/s)
        if isinstance(s,MhsAlg):
            return self*s.inv()
    
    def __rdiv__(self,s):
        return s * self.inv()
    
    def __neg__(self):
        return (-1)*self
    
    def aerr(self,err):
        self.err = min(self.err, err)
        return self
    
    def add(self,s):
        """Adds s to self, leaving self in place"""
        if isinstance(s,Number):
            D_add(self.data,_unit(),s)
        elif isisntance(s,MhsAlg):
            for x in s.data:
                D_add(self.data,x,s.data[x])
        return None
    
    def Gm(self, r):
        """Acts on an element by an element r in G_m(Q)"""
        T = MhsAlg(self)
        for x in T.data:
            T.data[x] *= r ** (weight(x) - x[0])
        T.desc = None
        return T
    
    def Gmp(self):
        """Acts on an element by (p) in G_m(Q_{p\to\infty})"""
        T = MhsAlg()
        for x in self.data:
            T.data[_weight(x), x[1]] = self.data[x]
        T.desc = None
        return T
    
    def new_e(self, p):
        mzv = [0 for i in range(len(mzv_alg()))]
    
    def e(self, p):
        """Evaluate at a prime p"""
        S = 0
        self.clean()
        k = self.err
        for i in range(num()):
            L = mzv_mono()[i]
            T = 0
            v = 99999
            for x in self.data:
                if x[1] == L:
                    T = red(T + self.data[x] * Rat(p)**x[0], p, k)
                    v = min(v, x[0])
            if T != 0:
                Q = mono_to_basis()[i]
                S += T*sum(Q[j]*Rat(p)**(weights()[j]-wt(L))*MHS_approx(p,k - v - weights()[j] + wt(L),p-1,*mhs_basis()[j]) for j in range(num()) if Q[j]!=0)
        return red(S, p, k)

def P(r):
    """p^r"""
    T = MhsAlg((r,_unit()[1]))
    T.desc = "p^r"
    return T

def Hp(*s):
    """Multiple harmonic sum H_{p-1}(s)"""
    T = P(-sum(s))*hp(*s)
    T.desc = "H_{p-1}" + tstr(s)
    return T
