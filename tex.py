from .sage_compatibility import latex

def tstr(s):
    """Converts a tuple to a string (without extra comma)"""
    if len(s) != 1:
        return str(tuple(s))
    return "(" + str(s[0]) + ")"

def error_tex(err, var='p'):
    if err > 9000:
        return None
    if err == 0:
        return "O(1)"
    elif err == 1:
        return "O(%s)"%var
    else:
        return "O(%s^{%i})"%(var,err)

    
def fix_prod(S1,S2,sign=False):
    """Returns a string that's basically the product
    S1*S2, but fixes funny 1's issue"""
    PM = ["1","-1"]
    if S1 in PM and S2 in PM:
        T = str(int(S1)*int(S2))
    elif S1 in PM:
        T = S1[:-1] + S2
    elif S2 in PM:
        T = S2[:-1] + S1
    else:
        T = S1 + S2
    if sign:
        if T[0] == "-":
            T = " - " + T[1:]
        else:
            T = " + " + T
    return T


def tex_write(D,f1,o1,f2,o2,err=99999,first=False,var='p'):
    """ 
    Input: Dictionary,
        Function from keys to outer symbol
        Function from keys to ordering of outer symbol
        Function from keys to inner symbol
        ordering of inner symbol
    Output: String of latex
    """
    K = list(D.keys())
    for x in K:
        if D[x] == 0:
            D.pop(x)
    K = list(D.keys())
    K.sort(key=lambda x:(o1(x),o2(x)))
    O1 = []
    for x in K:
        if f1(x) not in O1:
            O1.append(f1(x))
    S = ""
    for a in O1:
        Xa = [x for x in K if f1(x)==a]
        T = ""
        if len(Xa) == 1:
            w = f2(Xa[0])
            if first:
                T = fix_prod(a,w)
            else:
                T = fix_prod(w,a)
            T = fix_prod(str(latex(D[Xa[0]])),T,a!=O1[0])
        elif len(Xa) >= 2:
            for x in Xa:
                T += fix_prod(str(latex(D[x])),f2(x),x!=Xa[0])
            if a != "1":
                T = r"\left(" + T + r"\right)"
            if first:
                T = fix_prod(a,T,a!=O1[0])
            else:
                T = fix_prod(T,a,a!=O1[0])
        S += T
    if err <= 9000:
        if S == "":
            S = error_tex(err,var)
        else:
            S += " + " + error_tex(err,var)
    if S == "":
        S = "0"
    return S

def tex_write_single(D,f,o,err=99999,var='p'):
    K = list(D.keys())
    for x in K:
        if D[x] == 0:
            D.pop(x)
    K = list(D.keys())
    K.sort(key=o)
    S = ""
    for x in K:
        t = fix_prod(str(latex(D[x])),f(x))
        if S != "" and t[0] != '-':
            S += "+"
        S += t
    if err <= 9000:
        if S == "":
            S = error_tex(err,var)
        else:
            S += " + " + error_tex(err,var)
    if S == "":
        S = "0"
    return S

def monomial(L,P):
    """Returns LaTeX for a monomial with variable names in L, powers in P"""
    S = ""
    for i in range(len(L)):
        if P[i] != 0:
            S += L[i]
            if P[i] != 1:
                S += "^{%i}"%P[i]
    if S == "":
        S = "1"
    return S

def un_texify(s):
    T = ""
    T = str(s)
    T = T.replace("}{","/")
    T = T.replace(r"\zeta","zeta")
    for q in [r"\frac",r"\left",r"\right","{","}"]:
        T = T.replace(q,"")
    return T

def fix_pow(var,ex):
    if ex == 0:
        return "1"
    elif ex == 1:
        return var
    else:
        return var + "^{%i}"%ex
