try:
    from sage.rings.integer import Integer as Int
    from sage.rings.rational import Rational as Rat
    from sage.misc.latex import latex
    from sage.functions.hypergeometric import factorial, binomial
    from sage.all import prod, bernoulli
    numer = Rat.numerator
    denom = Rat.denominator
    from sage.matrix.matrix_rational_dense import Matrix_rational_dense as Mat
    print("Sage data types loaded")
    _sage = True
except:
    Int = int
    from sympy import Rational as Rat, Matrix as Mat
    numer = lambda x:x.p
    denom = lambda x:x.q
    from sympy import latex, binomial, factorial, bernoulli, prod
    print("SymPy data types loaded")
    _sage = False
Number = (int, Int, Rat)
