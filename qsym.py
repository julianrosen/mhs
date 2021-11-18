import sympy as sp

from .sage_compatibility import Int, Rat, Number

# Keeps track of stuffle products that have been computed
_mul={}

class qsym:
    """This class implements the ring of quasi-symmetric functions
    over the field of rational numbers."""
    def __init__(self, s = None, err = 99999):
        self.data = {}
        self.err = err
        if isinstance(s, Number):
            self.data[()] = s
        elif isinstance(s, tuple):
            self.data[s] = 1
        elif isinstance(s, qsym):
            self.data = dict(s.data)
            self.err = s.err
    
    def clean(self):
        L = list(self.data)
        for s in L:
            if self.data[s] == 0 or sum(s) >= self.err:
                self.data.pop(s)
    
    def __add__(self, other):
        """Add two elements."""
        if isinstance(other, qsym):
            T = qsym(self)
            T.clean()
            for s in other.data:
                if s in T.data:
                    T.data[s] += other.data[s]
                elif other.data[s]!=0:
                    T.data[s] = other.data[s]
            T.err = min(self.err, other.err)
            return T
        if isinstance(other, Number):
            return self + qsym(other)
    
    def __radd__(self,other):
        """Add two elements."""
        return self + other
    
    def __sub__(self, other):
        """Subtract the second element from the first."""
        return self + -1 * other
    
    def __rsub__(self, other):
        """Subtract the first element from the second."""
        return -1 * self + other
    
    def __mul__(self,other):
        """Multiply two elements."""
        if isinstance(other, Number):
            T = qsym(self)
            for s in T.data:
                T.data[s] *= other
            T.clean()
            return T
        if isinstance(other, qsym):
            T = qsym()
            for s in self.data:
                for t in other.data:
                    c = self.data[s] * other.data[t]
                    if sum(s) + sum(t) < self.err and c != 0:
                        
                        # Checks if stuffle of s and t has already been computed, if not, computes and adds to table
                        if frozenset([(s, t)]) in _mul:
                            L = _mul[frozenset([s,t])]
                        else:
                            if len(s) == 0:
                                L = qsym(t)
                            elif len(t) == 0:
                                L = qsym(s)
                            else:
                                L = ((qsym(s[1:]) * qsym(t)).lcat(s[:1]) + (qsym(s) * qsym(t[1:])).lcat(t[:1]) + (qsym(s[1:]) * qsym(t[1:])).lcat((s[0]+t[0],)))
                            _mul[frozenset([s,t])] = L
                        T += c * L
            T.err = min(self.err, other.err)
            return T
    
    def __rmul__(self, other):
        """Multiply two elements."""
        return self * other
    
    def __pow__(self, other):
        if other == 0:
            return qsym(1)
        return self * self ** (other - 1)
    
    def lcat(self, s):
        """Concatenate s to the left of each monomial quasi-symmetric function."""
        T = qsym()
        if not isinstance(s,tuple):
            a = (s, )
        else:
            a = s
        for t in self.data:
            if sum(a) + sum(t) < self.err:
                T.data[a + t] = self.data[t]
        T.err = self.err
        return T
    
    def lcatE(self, s):
        """Concatenate s to the left of each monomial quasi-symmetric function,
        allowing for equality in nested sum"""
        T = qsym()
        if isinstance(s, Int):
            a = (s, )
        else:
            a = s
        for t in self.data:
            T += self.data[t] * qsym(a+t)
            if len(t) >= 1 and len(a) >= 1:
                T += self.data[t] * qsym(a[:-1] + (a[-1]+t[0],) + t[1:])
        T.err = self.err
        return T
