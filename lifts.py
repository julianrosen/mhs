from mhs_alg import MhsAlg, P
from sympy import Sum, var, poly
from sage_compatibility import Int, Rat, Number, prod, latex, factorial, binomial, bernoulli
from qsym import qsym
from global_vars import max_weight, num, mzv_alg, mono_to_basis, mhs_to_mono, mzv_mono, mzv_to_mono, weights, mhs_basis
from rational import MHS_approx, red
from misc import psc, psc_, D_add, comps, comps_i, lists
from tex import error_tex, tex_write, tstr, monomial
from misc import D_add, wt, bn, decomps, S_min_v

def hp(*s):
    """Weighted multiple harmonic sum p^{|s|}H_{p-1}(s)"""
    T = MhsAlg()
    if len(s) == 1 and isinstance(s[0],(list,tuple)):
        s = s[0]
    T.desc = ("p^{%i} H_{p-1}"%sum(s)) + tstr(s)
    if sum(s) > max_weight():
        T.err = sum(s)
        return T
    M = mhs_to_mono()[s]
    for i in range(num()):
        if M[i] != 0:
            T.data[(wt(mzv_mono()[i]),mzv_mono()[i])] = M[i]
    T.err = max_weight() + 1
    return T
    
def Hp(*s):
    """Multiple harmonic sum H_{p-1}(s)"""
    T = P(-sum(s))*hp(*s)
    T.desc = "H_{p-1}" + tstr(s)
    return T

def binp(k, r):
    """Binomial coefficient {kp \\choose rp}, as element of the MhsAlg algebra"""
    T = binomial(k, r) * prod(sum(i ** j * hp(* j * (1, )) for j in range(num())) for i in range(k-r,k)) / prod(sum(i ** j * hp(* j * (1, )) for j in range(num())) for i in range(r))
    T.desc = "{%ip \choose %ip}" % (k, r)
    return T

def bin_linp(a, b, c, d):
    """Binomial coefficient {ap + b \\choose cp + d}, as element of the MhsAlg algebra"""
    if a < c or (a == c and b < d):
        return MhsAlg(0)
    T = binp(a, c)
    if b >= 0:
        for i in range(1, b + 1):
            T *= (a*P(1) + i) / ((a-c)*P(1) + i)
    else:
        for i in range(-b):
            T *= ((a-c)*P(1) - i) / (a*P(1) - i)
    if d >= 0:
        for i in range(d):
            T *= ((a-c)*P(1) + b - i) / (c*P(1) + i + 1)
    else:
        for i in range(-d):
            T *= (c*P(1) - i) / ((a-c)*P(1) + b + 1 + i)
    T.desc = "{%ip " % a + ("+ %i" % b if b > 0 else "- %i" % -b if b < 0 else "") + " \choose %ip " % c + ("+ %i" % d if d > 0 else "- %i" % -d if d < 0 else "") + "}"
    return T


def aperybp():
    """Apery number b_{p-1}, as element of the MhsAlg algebra"""
    a = sum((-1) ** i * qsym(i * (2,)) for i in range(max_weight() / 2))
    a = a * a
    a = 1 + a.lcat(4) - 2*a.lcat(3) + a.lcat(2)
    T = sum(a.data[s] * hp(*s) for s in a.data)
    T.desc = r"\sum_{n=0}^{p-1}{p-1\choose n}^2{p+n-1\choose n}^2"
    T.aerr(max_weight()-(max_weight()%2)+2)
    return T

def zetap(n):
    """p-adic zeta value zeta(n), as element of the MhsAlg algebra"""
    """Wrong normalization, don't use"""
    if n%2 == 0:
        return MhsAlg(0)
    if n > max_weight():
        return MhsAlg(err=n, desc="p-adic zeta values z(%i)" % n)
    T = sum(Rat((-1) ** (s + n + n)) / (n-1) * binomial(s - 1, n - 2) * bernoulli(s + 1 - n) * hp(s) for s in range(n-1, max_weight()+1))
    T.desc = "\zeta_p(%i)" % n
    return T

def mzvp(*s):
    """p-adic multizeta value zeta(s), as element of the MhsAlg algebra"""
    if len(s) == 1 and isinstance(s[0], (tuple, list)):
        s = s[0]
    if sum(s) > max_weight():
        return MhsAlg(err=0)
    T = MhsAlg()
    x = mzv_to_mono()[s]
    for i in range(num()):
        if x[i] != 0:
            T.data[(0,mzv_mono()[i])] = x[i]
    T.desc = "\zeta_p" + tstr(s)
    return T

def Spr(r, s, err=max_weight() + 1):
    """Sum_{ap^r >= n_1 > ... > n_k >= 1m p \\nmid n_i} 1 / (n_1^{s_1} ... n_k^{s_k}),
    as element of the MhsAlg algebra"""
    if r <= 0:
        raise ValueError('First argument must be positive')
        return None
    elif len(s) == 0:
        T = MhsAlg(1)
    elif r == 1:
        T = Spr1(s, err)
    else: #r >= 2
        T = sumSr(r, res_pr(len(s)), s, err)
    T.desc = (r"S_{p^%i}^{(p)}"%r) + tstr(s)
    return T

def Spp(r, s, err=max_weight() + 1):
    """# Sum_{p^r >= n_1 > ... > n_k >= 1} 1 / (n_1^{s_1} ... n_k^{s_k}), 
    as element of the MhsAlg algebra"""
    s = tuple(s)
    T = None
    if err <= S_min_v(r, s):
        T = MhsAlg()
        T.err = err
    elif r <= 0:
        raise ValueError('First argument must be positive')
        return None
    elif len(s) == 0:
        T = MhsAlg(1)
    elif r == 1:
        T = Spr1(s, err) + P(-s[0]) * Spr1(s[1:], err + s[0])
    else:
        assert r >= 2
        if len(s) == 0:
            T = MhsAlg(1)
        elif len(s) == 1:
            if s[0] <= 0:
                T = sum(psc(-s[0], j) * P(r * j) for j in range(min(2 - s[0], err)))
                if not isinstance(T, MhsAlg):
                    T = MhsAlg(T)
                T.aerr(err)
            else:
                T = sumSr(r, res_all(len(s)), s, err)
        else:
            assert len(s) >= 2
            if s[0] <= 0:
                T = MhsAlg(0)
                for j in range(2 - s[0]):
                    if psc(-s[0], j) != 0:
                        T += psc(-s[0], j) * (P(r*j) * Spp(r, s[1:], err - r * j) - Spp(r, (s[1] - j,) + s[2:], err))
                    T.aerr(err)
            elif s[-1] <= 0:
                return sum(psc_(-s[-1], j) * Spp(r, s[:-2] + (s[-2] - j,), err) if psc_(-s[-1], j) != 0 else 0 for j in range(2 - s[-1])).aerr(err)
            else:
                for i in range(1, len(s) - 1):
                    if s[i] <= 0:
                        T = sum((psc_(-s[i], j) * Spp(r, s[:i - 1] + (s[i - 1] - j,) + s[i + 1:], err) if psc_(-s[i], j) != 0 else 0) - (psc(-s[i], j) * Spp(r, s[:i] + (s[i + 1] - j,) + s[i+2:], err) if psc(-s[i], j) != 0 else 0) for j in range(2 - s[i])).aerr(err)
                        break
                if T is None:
                    T = sumSr(r, res_all(len(s)), s, err)
    T.desc = (r"H_{p^%i}"%r) + tstr(s)
    return T


def Spr1(s, err=99999):
    """Sum_{p-1 >= n_1 > ... > n_k >= 1} n_1^{-s_1} ... n_k ^ {-s_k}.
    as element of the MhsAlg algebra"""
    ss = tuple(s)
    if len(ss) == 0:
        return MhsAlg(1)
    elif len(ss) == 1:
        if s[0] <= 0:
            T = sum(psc_(-ss[0], j) * P(j) for j in range(min(2 - ss[0], err)))
            if not isinstance(T, MhsAlg):
                T = MhsAlg(T)
            T.aerr(err)
            return T
        else:
            return (P(-ss[0]) * hp(*ss)).aerr(err)
    if ss[0] <= 0:
        return sum((psc_(-ss[0], j) * P(j) * Spr1(ss[1:], err - j) if psc_(-ss[0], j) != 0 else 0) - (psc(-ss[0], j) * Spr1((ss[1] - j,) + ss[2:], err) if psc(-ss[0], j) != 0 else 0) for j in range(2 - ss[0])).aerr(err)
    if ss[-1] <= 0:
        return sum(psc_(-ss[-1], j) * Spr1(ss[:-2] + (ss[-2] - j,), err) if psc_(-ss[-1], j) != 0 else 0 for j in range(2 - ss[-1])).aerr(err)
    for i in range(1, len(ss) - 1):
        if ss[i] <= 0:
            return sum((psc_(-ss[i], j) * Spr1(ss[:i - 1] + (ss[i - 1] - j,) + ss[i + 1:], err) if psc_(-ss[i], j) != 0 else 0) - (psc(-ss[i], j) * Spr1(ss[:i] + (ss[i + 1] - j,) + ss[i+2:], err) if psc(-ss[i], j) != 0 else 0) for j in range(2 - ss[i])).aerr(err)
    return (P(-sum(ss)) * hp(*ss)).aerr(err)



def Srst(r, A, J, s, err=max_weight() + 2):# Should check default err
    """Arbitrary restricted power sum, as element of the MhsAlg algebra.

    This is a sum of (a_1*p - j_1)^{-s_1} * ... * (a_k*p - j_k)^{-s_k},
    for p^{r-1} >= a_1 >= ... >= a_k >= 1 and
    0 <= j_1, ..., j_k <= p-1 satisfying constraints indicated by arguments.
    Arguments:
        r is the power of p
        A is a binary list of length k-1 saying which a's are equal
            A 0 means two a's are equal, 1 means inequality is strict
        J is a list of lists saying how the j's are ordered,
            The first element is the list of j indices where j is 0
            The next element is the list of j indices taking the next smaller value
            And so on
        s is a composition
    """
    if len(s) == 0:
        return 1
    j = 0
    Li = []
    for i in range(len(s)):
        if i in J[0]:
            Li.append(0)
        else:
            Li.append(j)
            j += 1
    no_p = j
    Sv = S_min_v(r - 1, s)
    po = sum(s[j] for j in J[0])
    TT = MhsAlg()
    for L in lists(no_p, err - Sv):     # Might not need this many L
        c = 1
        for i, a in enumerate(s):
            if i not in J[0]:
                c *= binomial(a + L[Li[i]] - 1, L[Li[i]]) * (-1) ** (a % 2)
        T = [a if i in J[0] else -L[Li[i]] for i, a in enumerate(s)]
        AN = [T[0]]
        for i, a in enumerate(A):
            if a == 0:
                AN[-1] += T[i + 1]
            else:
                AN.append(T[i + 1])
        JN = []
        for x in J[1:]:
            JN = [sum(s[i] + L[Li[i]] for i in x)] + JN
        e = err + po - sum(L)
        Sv = S_min_v(r - 1, AN)
        Q = c * Spr1(JN, max(e - Sv, 0))
        WW = Spp(r - 1, AN, Q.err + Sv + 1)
        Q *= WW
        Q *= P(-po + sum(L))
        TT += Q
        err = min(err, TT.err)
    TT.aerr(err)
    return TT

def curp(r, k, err=max_weight() + 1):
    """Compute Sum__{n_1 + ... + n_k = p^r, p \\nmid n_i} (n_1...n_k)^{-1},
    as an element of the MhsAlg algebra"""
    if r <= 0:
        raise ValueError("First argument must be positive")
    if r == 1:
        return factorial(k) * P(-1) * hp(*(1 for i in range(k - 1)))
    T = factorial(k) * P(-r)* sumSr(r, res_cur(k), [1] * (k - 1), err + r)
    T.desc = r"\sum_{\substack{n_1"
    for i in range(2, k + 1):
        T.desc += " + n_%i" % i
    T.desc += " = p^%i\\\\p\\nmid "%r
    for i in range(1,k+1):
        T.desc += "n_{%i}"%i
    T.desc += r"}}\frac{1}{"
    for i in range(1, k + 1):
        T.desc += "n_%i" % i
    T.desc += "}"
    return T

# Generator for restrictions to give unrestricted harmonic sum
def res_all(k):
    for s in comps(k):
        for y in bn(k - 1):
            for J in decomps(s, range(k)):
                B = True
                for i, A in enumerate(y):
                    if A == 0 and find_index(J, i) >= find_index(J, i + 1):
                        B = False
                        break
                if B:
                    yield y, J
                    yield y, [[]] + J

def find_index(L, i):
    """Find the smallest j such that i is in L[j]"""
    for a, q in enumerate(L):
        if i in q:
            return a

# Generator for restrictions to give curious harmonic sum
def res_cur(k):
    for A, J in res_all(k - 1):
        if 0 not in J[0] and k - 2 not in J[0]:
            B = True
            for i in range(k - 2):
                if find_index(J, i) == find_index(J, i + 1):
                    B = False
                    break
            if B:
                yield A, J

# Generator for restrictions to give p-restricted harmonic sum
def res_pr(k):
    for s in comps(k):
        for y in bn(k - 1):
            for J in decomps(s, range(k)):
                B = True
                for i, A in enumerate(y):
                    if A == 0 and find_index(J, i) >= find_index(J, i + 1):
                        B = False
                        break
                if B:
                    yield y, [[]] + J

def sumSr(r, it, s, err=99999):
    T = MhsAlg(0, err=err)
    for A, J in it:
        T += Srst(r, A, J, s, T.err)
    return T

def hpm(r, *s):
    if len(s) == 0:
        T = MhsAlg(1)
    else:
        T = P(r*sum(s))*Spp(r,s) - hpm(r,*s[1:])
    #T = MhsAlg(sum((-1)**i * P(r*sum(s[i:])) * Spp(r,s[i:]) for i in range(len(s) + 1)))
    #T = MhsAlg(sum((-1)**i * P(r*sum(s[i:])) * Spp(r,s[i:]) for i in range(2)))
    T.desc = "p^blah * H_{p^%i-1}" % r + tstr(s)
    return T

def sum_equal(L, q):
    return sum(sum(a * q.coef[s] * (Spr1((-i,) + s) + (Spr1((s[0] - i,) + s[1:]) if s != () else 0)) for i, a in enumerate(L)) for s in q.coef)

def sum_not(L, q):
    return sum(sum(a * q.coef[s] * Spr1((-i,) + s) for i, a in enumerate(L)) for s in q.coef)

def sun_sump(a, b):
    return sum_equal([1 if i == a else 0 for i in range(a+1)], qsym((1,)) ** b)

def harm_halfp(s):
    if s == 1:
        raise ValueError("Argument must be at least 2")
    T = zetap(s) + sum(binomial(-s, j) * Rat(1-2**(s+j)) / 2**j * zetap(s+j) for j in range(max_weight()-s+1))
    T.desc = "p^{%i} \sum_{n=1}^{\\frac{p-1}{2}} \\frac{1}{n^{%i}}"%(s,s)
    return T

def altp(s):
    if s == 1:
        raise ValueError("Argument must be at least 2")
    T = Rat(1) / 2**(s-1) * harm_halfp(s) - hp(s)
    T.desc = "p^{%i} \\sum_{n=1}^{p-1} \\frac{(-1)^n}{n^{%i}}"%(s,s)
    return T

def hrp(r, *s):
    """ p-restricted multiple harmonic sum with upper bound rp"""
    if r == 1:
        return hp(*s)
    k = len(s)
    T = MhsAlg()
    for i in range(k+1):
        for L in lists(i,max_weight() - sum(s)):
            T += hrp(r-1,*s[i:]) * hp(*[s[j] + L[j] for j in range(len(L))]) *\
                    prod(binomial(-s[j], L[j]) for j in range(i))
    T.desc = "p^{%i}H_{%ip}^{(p)}"%(sum(s),r) + tstr(s)
    return T

def S_poly_pr(L,*s):
    T = MhsAlg()
    p = var('p')
    T.desc = "H^{(p)}_{" + latex(sum(L[i]*p**i for i in range(len(L)))) + "}" + tstr(s)
    return T

def power_sum_poly(L, s, res=None):
    """Power sum with non-negative powers, lower limit is 1, upper limit is a polynomial in p, res gives inequalities (1 means inequality, 0 means equality)"""
    if len(s) == 0:
        return MhsAlg(1)
    if res is not None:
        S = [s[0]]
        for i in range(len(res)):
            if res[i] == 1:
                S.append(s[i+1])
            else:
                S[-1] += s[i+1]
        return power_sum_poly(L,S)
    N = len(s)
    VARS = []
    var('p')
    for i in range(N):
        VARS.append(var('n%i'%i))
    B = 1
    for i in range(N-1,0,-1):
        B *= VARS[i] ** s[i]
        B = Sum(B,(VARS[i],1,VARS[i-1]-1)).doit()
    B *= VARS[0] ** s[0]
    B = Sum(B,(VARS[0],1,sum(L[i]*p**i for i in range(len(L))))).doit()
    b = poly(B,p).all_coeffs()
    b.reverse()
    T = sum(b[i] * P(i) for i in range(len(b)))
    #T.desc = r"\sum_{" + latex(sum(L[i]*p**i for i in range(len(L)))) + "\geq "
    #for i in range(1,N):
    #    T.desc += "n_{%i}"%i
    #    T.desc += ">" if res[i-1] == 1 else r"\geq "
    #T.desc += "n_{%i}\geq 1}"%N
    #for i in range(1,N+1):
    #    if s[i-1] > 0:
    #        T.desc += "n_{%i}"%i
    #        if s[i-1] > 1:
    #            T.desc += "^{%i}"%s[i-1]
    return T

def ordered_tuples(a,b,k):
    if k == 0:
        yield ()
    else:
        for x in range(a,b+1):
            for t in ordered_tuples(a,x-1,k-1):
                yield (x,) + t

def H_poly_pr(L, s, err=None):
    """p-restricted MHS, upper limit is a polynomial in p"""
    K = len(s)
    if K == 0:
        return MhsAlg(1)
    T = MhsAlg()
    if err is None:
        err = max_weight() + 1 - sum(s)
    if L[0] > 0:
        Q = sum(L[i] * P(i) for i in range(1,len(L)))
        for i in range(K+1):
            for t in ordered_tuples(1,L[0],i):
                T += prod((Q+Rat(t[j]))**(-s[j]) for j in range(i))\
                    * H_poly_pr([0]+L[1:],s[i:],err)
        T.err = err
    elif L[0] < 0:
        Q = sum(L[i] * P(i) for i in range(1,len(L)))
        T = H_poly_pr([0]+L[1:],s,err)
        for i in range(1,K+1):
            for t in ordered_tuples(L[0]+1, -1, i):
                T -= prod((Q+Rat(t[j]))**(-s[j]) for j in range(i))\
                    * H_poly_pr(L, s[i:], err)
        T.err = err
    else:
        T = MhsAlg()
        for Q in lists(K, err - 1):
            c = prod((-1)**((Q[i]+s[i])%2) * binomial(-s[i],Q[i]) for i in range(K))
            for res in res_pr(K):
                if c != 0:
                    T += c * P(sum(Q))\
                        * power_sum_poly(L[1:],Q,res=res[0])\
                        * Spr1([sum(s[i]+Q[i] for i in q) for q in res[1][:0:-1]])
        T.err = err
    var('p')
    f = sum(L[i]*p**i for i in range(len(L)))
    T.desc = "H_{" + latex(f) + "}^{(p)}" + tstr(s)
    return T

def schurp(*t):
    """Schur function indexed by t, evaluated at 1/1,...,1/(p-1)"""
    tt = list(t)
    tt.sort(key=lambda x:-x)
    s = SymmetricFunctions(QQ).schur()
    e = SymmetricFunctions(QQ).elementary()
    L = list(e(s(tt)))
    T = MhsAlg()
    for q in L:
        a,b = q
        T += Hp(*a)*b
    return T