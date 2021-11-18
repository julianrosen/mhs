from .sage_compatibility import binomial, bernoulli
from .global_vars import mzv_alg, mzv_mono

# Coefficient of x^d in Sum_{n=1}^x n^k
def psc(k, d):
    if d <= 0 or d >= k+2:
        return 0
    else:
        return (-1) ** (k + d + 1) * binomial(k + 1, d) * bernoulli(k + 1 - d) / (k + 1)

# Coefficient of x^d in Sum_{n=1}^{x-1} n^k
def psc_(k, d):
    return psc(k, d) - (1 if d == k else 0)


# Lower bound on valuation of S(r, s)
def S_min_v(r, s):
    a = [-x for x in s]
    a.sort()
    v = 0
    for i, A in enumerate(a):
        if A >= 0:
            break
        v += A * (r if i == 0 else r - 1)
    return v

# Generator for subsets of L of size k
def subsets(k, L):
    if k == 0:
        yield []
        raise StopIteration()
    for i, A in enumerate(L):
        for x in subsets(k - 1, L[i + 1:]):
            yield [A] + x

# Difference of sets
def diff(A, B):
    return [x for x in A if x not in B]

# Generator for decompositions of L into unordered pieces with sizes given by s, 
def decomps(s, L):
    if len(s) == 1:
        yield [L]
        raise StopIteration()
    for x in subsets(s[0], L):
        for y in decomps(s[1:], diff(L, x)):
            yield [x] + y

# Generator for lists of non-negative integers
def lists(k, W):
    if k == 0:
        yield []
    else:
        for i in range(W + 1):
            for L in lists(k - 1, W - i):
                yield [i] + L


# Generator for binary lists of length k
def bn(k):
    if k < 0:
        print(k)
        raise ValueError("Input needs to be a non-negative integer")
    if k == 0:
        yield []
    else:
        for L in bn(k - 1):
            yield L + [0]
            yield L + [1]

# Convert binary list to a tuple
def get_tuple(L):
    if L == []:
        return (1,)
    s = get_tuple(L[1:])
    if L[0] == 1:
        return (1,) + s
    else:
        return (s[0] + 1,) + s[1:]

# Generator for compositions of weight exactly W
def comps(W):
    if W == 0:
        yield ()
    else:
        for L in bn(W - 1):
            yield get_tuple(L)

# Generator for compositions of weight at most W
def comps_i(W):
    for i in range(W + 1):
        for L in comps(i):
            yield L

# Weight of element of mzv monomial basis
def mono_weight(i):
    return sum(a * sum(mzv_alg()[j]) for j, a in enumerate(mzv_mono()[i]))

def D_add(D,x,r):
    if x in D:
        D[x] += r
    else:
        D[x] = r
    return D

def wt(s):
    return sum(s[i]*sum(mzv_alg()[i]) for i in range(len(s)))
