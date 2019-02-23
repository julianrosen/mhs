from sage_compatibility import Int, Rat, numer, denom

def v(p, x):
    """Compute p-adic valuation of x.
    Returns 99999 if x is 0"""
    if p <= 1:
        raise ValueError('Invalid prime for computing valuation.')
        return None
    if x == 0:
        return 99999
    if isinstance(x, Int):
        x = Rat(x)
    if denom(x) % p == 0:
        return v(p, x * p) - 1
    if numer(x) % p == 0:
        return v(p, x / p) + 1
    return 0


def egcd(a, b):
    """Compute the greatest divisor of a and b."""
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    """Compute the inverse of a modulo m."""
    g, x, y = egcd(a % m, m)
    if g != 1:
        raise Exception("%i is not invertible modulo %i" % (a, m))
    else:
        return x % m

def red(x, p, r):
    """Compute the smallest fraction congruent to x modulo p^r
    whose denominator is a power of p."""
    k = v(p, x)
    if k >= r:
        return 0
    if r < 0:
        return Rat(1) / p**-r * red(x * p**-r, p, 0)
    y = x * Rat(p)**(-k)
    return (numer(y) * modinv(denom(y), p**(r-k)) * Rat(p)**k) % p**r

def MHS(n, *S):
    """Compute the multiple harmonic sum H_n(S), where S is
    allowed to have non-positive entries."""
    if len(S) == 1 and isinstance(S[0], (tuple, list)):
        S = S[0]
    if len(S)==0:
        return 1
    elif n >= 0:
        return sum(Rat(i) ** (-S[0]) * MHS(i-1, *S[1:]) for i in range(1, n + 1))
    elif len(S) == 1 and S[0] >= 1:
        raise ValueError("Illegal evaluation of harmonic sum")
    else:
        x = MHS(n,*S[1:])
        if x != 0 and S[0] >= 1 and n == -1:
            raise ValueError("Illegal evaluation of harmonic sum")
        else:
            return MHS(n+1,*S) - Rat(1,(n+1)**S[0])*x

# Like Har, but returns value modulo a prime power
def MHS_approx(p,k,n,*S):
    """Compute the multiple harmonic sum H_n(S), where S is
    allowed to have non-positive entries."""
    if len(S) == 1 and isinstance(S[0], (tuple, list)):
        S = S[0]
    if len(S)==0:
        return 1
    elif n >= 0:
        return red(sum(red(Rat(i) ** (-S[0]), p, k) * MHS_approx(p, k, i - 1, *S[1:]) for i in range(1, n + 1) if i%p != 0), p, k)
    elif len(S) == 1 and S[0] >= 1:
        raise ValueError("Illegal evaluation of harmonic sum")
    else:
        x = MHS_approx(p,k,n,*S[1:])
        if x != 0 and S[0] >= 1 and n == -1:
            raise ValueError("Illegal evaluation of harmonic sum")
        else:
            return red(MHS_approx(p,k,n+1,*S) - Rat(1,(n+1)**S[0])*x,p,k)
        
def hpQ(p,*S):
    """Compute the weighted multiple harmonic sum p^{|S|}H_{p-1}(S)."""
    if len(S) == 1 and isinstance(S[0], (tuple, list)):
        S = S[0]
    return p**sum(s for s in S)*MHS(p-1,*S)

def base(x,p):
    if x == 0:
        return []
    return [x % p] + base((x-(x%p)) / p, p)