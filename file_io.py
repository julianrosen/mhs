# Multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool
_max_workers = 10 # Number of cpus



from .qsym import qsym
from .global_vars import set_vars
import sympy as sp
mat = sp.Matrix
from os.path import join, isfile, isdir, realpath, dirname
from os import makedirs, getcwd
import random
import bz2
from .misc import lists, comps_i
from time import time
from .sage_compatibility import Int, Rat, Number, _sage
import datetime


MHS_DATA_DIR = join(dirname(realpath(__file__)), "mhs_data")
MZV_DATA_MINE_DIR = join(dirname(realpath(__file__)), "data_mine")


def parallel(func,max_workers=_max_workers):
    def new_func(L):
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = [executor.submit(func,q) for q in L]
            for f in as_completed(results):
                yield f.result()
    return new_func



def import_data(n,silent=False):
    max_weight = n
    filename = join(MHS_DATA_DIR, "mhs_data_%i.dat" % n)
    if not isfile(filename):
        print("Data file " + join('.', filename) + " not found.")
        print("Would you like to create this file [y/n]: ", end=' ')
        c = eval(input())
        if c != '' and c[0] in ['y', 'Y']:
            if not create_data(n):
                return
    a = open(filename, "U")
    read(a)
    mzv_alg = read(a)
    mzv_mono = read(a)
    mzv_to_mono = read(a)
    mhs_to_mono = read(a)
    mhs_basis = read(a)
    mono_to_basis = read(a)
    a.close()
    num = len(mhs_basis)
    weights = [sum(s) for s in mhs_basis]
    set_vars(max_weight, num, mhs_basis, weights, mzv_alg,
      mzv_mono, mzv_to_mono, mono_to_basis, mhs_to_mono)
    if not silent:
        print("Imported data through weight %i" % max_weight)
    return None

def tuple_to_word(s):
    """Convert composition to binary word"""
    L = []
    for i in s:
        L += (i-1)*[0] + [1]
    return L

def word_to_tuple(L):
    """Convert binary word to composition"""
    if L == []:
        return ()
    if L[-1] != 1:
        raise ValueError("Last element of word needs to be 1")
    n = L.index(1)
    return (n + 1,) + word_to_tuple(L[n + 1:])

def str_to_tuple(s):
    "Convert comma separated binary string to composition"""
    return word_to_tuple([Int(i) for i in s.split(',')])

def dual(s):
    """Convert composition to its dual"""
    L = tuple_to_word(s)
    L.reverse()
    return word_to_tuple([1 if i == 0 else 0 for i in L])


def str_to_mzv_dict(q):
    """Convert string representing mzv monomial to dictionary"""
    L = q.replace("-z","-1*z").replace("+z","+1*z").replace("+","*+").replace("-","*-").split("*")
    s = None
    r = None
    T = {}
    if L[0] == '':
        del L[0]
    if L[0][0] == 'z':
        L = ['1'] + L
    for x in L:
        if 'z' not in x:
            if s != None:
                T[s] = r
            s = ""
            if x[0] == "+":
                x = x[1:]
            r = Rat(x)
        else:
            s += '*' + x
    if s != None:
        T[s] = r
    if T == None:
        T = {}
    T[s] = r if r != None else 1
    return dict(T)

def add_line_to_dict(D, L):
    """Parses a line from MZV data mine and adds it to dictionary"""
    T = L
    if "Fill " not in T:
        return
    T = T.split("Fill ")[1]
    for w in ["Fill ", '\n', '\t', ';', ' ']:
        T = T.replace(w, '')
    T = T[T.index('(') + 1:].split('=')
    T[0] = T[0].split(')')[0]
    if T[1][0] == ' ':
        T[1] = T[1][1:]
    s = str_to_tuple(T[0])
    mzv = str_to_mzv_dict(T[1])
    D[s] = mzv
    D[dual(s)] = mzv


def get_chunk(a):
    """Reads off a chunk of a file"""
    T = ""
    for L in a:
        T += L
        if ';' in L:
            return T
    return T


def read_prc(D, filename):
    """Read prc (or bz2) file and adds to dictonary"""
    if filename[-3:] == 'bz2':
        a = bz2.BZ2File(filename, "rU")
    else:
        a = open(filename, "rU")
    L = get_chunk(a)
    while L != '':
        add_line_to_dict(D, L)
        L = get_chunk(a)
    a.close()
    
    return None


def read_hexc(D, W):
    """Read hexc.h and add to dictionary"""
    if not isfile(join(MZV_DATA_MINE_DIR, "hexc.h")):
        print("Cannot find necessary file " + join(MZV_DATA_MINE_DIR, "hexc.h"))
        print("This file can be obtained from the MZV data mine")
        print("http://www.nikhef.nl/~form/datamine/mzv/complete/complete.html")
        return False
    a = open(join(MZV_DATA_MINE_DIR, "hexc.h"), "rU")
    for L in a:
        if L == "#procedure mzv2\n":
            break
    L = get_chunk(a)
    while L != '' and "#procedure mzv%i" % (W + 1) not in L:
        add_line_to_dict(D, L)
        L = get_chunk(a)
    a.close()
    return True


def read_all(W):
    """Reads in all mzv relations through weight W"""
    D = {():{'':1}}
    if not read_hexc(D, W):
        return
    for i in range(11, W + 1):
        if isfile(join(MZV_DATA_MINE_DIR, "mzv%i.prc" % i)):
            read_prc(D, join(MZV_DATA_MINE_DIR, "mzv%i.prc" % i))
        elif isfile(join(MZV_DATA_MINE_DIR, "mzv%i.prc.bz2" % i)):
            read_prc(D, join(MZV_DATA_MINE_DIR, "mzv%i.prc.bz2" % i))
        else:
            print("Could not find file %s or %s" % (join(MZV_DATA_MINE_DIR, "mzv%i.prc" % i), join(MZV_DATA_MINE_DIR, "mzv%i.prc.bz2" % i)))
            return
    return D


def mod_z2(D):
    """Reduce modulo zeta(2)"""
    for k in D:
        T = {}
        for s in D[k]:
            if s != "*z2" and s[:4] != "*z2*" and s[:4] != "*z2^":
                T[s] = D[k][s]
        D[k] = T


def get_gens(s, B):
    """Get algebra generators from string"""
    L = s.replace('^','*').split('*')
    for i in L:
        if len(i) >= 1 and i[0] == 'z':
            t = tuple([Int(j) for j in i.split('z')[1:]])
            if t not in B:
                B.append(t)
                if t == (2,):
                    print("Added (2,)")
                    print("s: ", s)
                    print("L: ", L)


def get_gens_D(D):
    """Get algebra generators from big dictionary"""
    B = []
    for k in D:
        for q in D[k]:
            get_gens(q, B)
    def my_key(s):
        return [sum(s), len(s)] + list(s)
    B.sort(key = my_key)
    return B


def renormalize(D, W):
    """Add stuffle-regularized values (with zeta(1) = 0)"""
    for i in range(W):
        for s in list(D.keys()):
            if sum(s) >= W or (1,) + s in D:
                continue
            a = qsym((1,)) * qsym(s)
            D[(1,) + s] = {}
            for t in a.data:
                if t != (1,) + s:
                    for q in D[t]:
                        if q in D[(1,) + s]:
                            D[(1,) + s][q] += Rat(-1) / a.data[(1,) + s] * a.data[t] * D[t][q]
                        else:
                            D[(1,) + s][q] = Rat(-1) / a.data[(1,) + s] * a.data[t] * D[t][q]


def get_vs_basis(B, W):
    """Produce vector space basis from algebra basis"""
    if len(B) == 0 or W == 0:
        return [[]]
    S = get_vs_basis(B[:-1], W)
    T = []
    for i in S:
        for j in range((W - sum(i[q] * sum(B[q]) for q in range(len(B) - 1))) // sum(B[-1]) + 1):
            T.append(i + [j])
    def kk(L):
        A = list(L)
        A.reverse()
        return [sum(L[q] * sum(B[q]) for q in range(len(B)))] + A
    T.sort(key = kk)
    return T


def mzv_str_to_i(s, B, BB):
    """Take a string and returns the index of that zeta monomial"""
    a = s.split('*')[1:]
    L = len(B) * [0]
    for q in a:
        w = q.split('^')
        t = tuple([Int(d) for d in w[0].split('z')[1:]])
        if t in B:
            L[B.index(t)] += 1 if len(w) == 1 else Int(w[1])
    if L in BB:
        return BB.index(L)
    else:
        return -1


def convert(D, B, BB):
    """Convert dictionary using basis"""
    for k in D:
        L = len(BB) * [0]
        for q in D[k]:
            a = mzv_str_to_i(q, B, BB)
            if a != -1:
                L[a] += D[k][q]
        D[k] = L


def setup(n):
    """Get everything setup"""
    D = {}
    readin(D, n)
    renormalize(D, n)
    B = get_gens_D(D)
    BB = get_vs_basis(B, n)
    convert(D, B, BB)
    return D, B, BB

# Multiply two mzv's (given as lists)
def mult_mzv(L1, L2, BB):
    T = [0 for _ in BB]
    for i in range(len(L1)):
        for j in range(len(L2)):
            if list(mat(BB[i]) + mat(BB[j])) in BB:
                T[BB.index(list(mat(BB[i]) + mat(BB[j])))] += L1[i] * L2[j]
    return [[x] for x in T]

# Express an mhs in terms of p-adic MZV
def mhs_to_mzv(s, W, D, BB):
    T = mat(len(BB) * [0])
    for i in range(len(s) + 1):
        for L in lists(i, W - sum(s)):
            a = list(mat(s[:i]) + mat(L))
            a.reverse()
            T += (-1) ** sum(s[:i]) * sp.prod(sp.binomial(s[j] + L[j] - 1, L[j]) for j in range(i)) * mat(mult_mzv(D[s[i:]], D[tuple(a)], BB))
    return list(T)

# Takes a table of expressions for mhs in terms of p-adic MZV, produces a basis for mhs
def find_mhs_basis(E, W):
    L = list(E.keys())
    def kk(s):
        return [-sum(s), len(s)] + [-x for x in s]
    L.sort(key = kk)
    M = mat([E[x] for x in L]).transpose().rref()
    return L, M

# Garbage
def mhs_to_mzv_modified(s):
    return s,mhs_to_mzv(*s)

# Create the data file mhs_data_W.dat
def create_data(W,z2=False,max_workers=_max_workers,power=None):
    print("Attempting to create %s" % join(MHS_DATA_DIR, "mhs_data_%i.dat" % W))
    D = read_all(W)
    if D is None:
        return
    renormalize(D, W)
    if not z2:
        mod_z2(D)
    B = get_gens_D(D)
    BB = get_vs_basis(B, W)
    convert(D, B, BB)
    E = {}
    count = 0
    num = 2 ** W - 2
    target = 1
    print("MZV data read successfully")
    print("Computing expansion for MHS in terms of MZV (this may take some time)")
    start_time = time()
    if max_workers > 1:
        #L = [(s,W,D,BB) for s in comps_i(W)]
        #mhs_to_mzv_par = parallel(mhs_to_mzv_modified,max_workers)
        if power is None:
            if W in [10,11]:
                power = 10 ** (W-9)
            elif W >= 12:
                power = 10 ** (W-10)
        else:
            power = 10 ** power
        chunksize = max(2,2**W//(max_workers*20))
        with Pool(max_workers) as pool:
            start_time = time()
            for s in pool.imap_unordered(mhs_to_mzv_modified,((s,W,D,BB) for s in comps_i(W)),chunksize=chunksize):
        #for s in mhs_to_mzv_par(L):
                E[s[0][0]] = s[1]
                count += 1
                if W >= 10 and count * power > target * num:
                    est = int((time() - start_time) * (num - count)) // count
                    est_s = str(datetime.timedelta(seconds=est))
                    frac = target / power
                    print(f"{frac} completed (at this rate, {est_s} remaining)",flush=True)
                    target += 1
    else:
        for s in comps_i(W):
            E[s] = mhs_to_mzv(s, W, D, BB)
            count += 1
            if W in [8, 9]:
                if count * 10 > target * num:
                    print("%i0%% done" % target)
                    target += 1
            elif W >= 10:
                if count * 100 > target * num:
                    print("%i%% done" % target)
                    a.flush()
                    target += 1
    L, M = find_mhs_basis(E, W)
    make_dir("mhs_data")
    a = open(join(MHS_DATA_DIR, "mhs_data_%i.dat" % W), 'w')
    _max_weight = W
    #a.write("# Weight\n%i\n\n" % W)
    #a.write("# MHS basis\n")
    b = [L[i] for i in M[1]]
    b.reverse()
    _mhs_basis = b
    for i, x in enumerate(b):
        if i >= 1:
            #a.write(", ")
            pass
        #a.write(str(x))
    #a.write("\n\n# MHS in terms of basis\n")
    _mhs_to_basis = {}
    for i in range(len(L) - 1, -1, -1):
        b = list(M[0].col(i))
        b.reverse()
        _mhs_to_basis[L[i]] = b
        #a.write(str(L[i]) + ": " + str(b) + "\n")
    #a.write("\n# MZV algebra basis\n")
    _mzv_alg = B
    for i in range(len(B)):
        #a.write(str(B[i]))
        if i < len(B) - 1:
            #a.write(", ")
            pass
    #a.write("\n\n# MZV monomial basis\n")
    _mzv_mono = BB
    for x in BB:
        #a.write(str(x) + '\n')
        pass
    #a.write("\n# MZV in terms of monomial basis\n")
    LL = list(D.keys())
    def kk(s):
        return [sum(s), -len(s)] + list(s)
    _mzv_to_mono = D
    _mzv_to_mono_key = lambda s:[sum(s), -len(s)] + list(s)
    LL.sort(key = kk)
    for x in LL:
        #a.write(str(x) + ": " + str(D[x]) + "\n")
        pass
    #a.write("\n# MHS basis in terms of monomial basis")
    b = list(M[1])
    b.reverse()
    for i in b:
        #a.write('\n' + str(L[i]) + ": " + str(E[L[i]]))
        pass
    #a.write('\n')
    _mhs_to_mzv = E
    Q = [E[s] for s in _mhs_basis]
    a.write(f"# Weight\n{W}\n\n# MZV algebra basis\n")
    for s in _mzv_alg:
        a.write(f"{s}\n")
    a.write("\n# MZV monomial basis\n")
    for s in _mzv_mono:
        ss = tuple(s)
        a.write(f"{ss}\n")
    a.write("\n# MZV in monomial basis\n")
    for s in _mzv_to_mono:
        thing = _mzv_to_mono[s]
        a.write(str(s) + ": " + str(thing) + "\n")
    a.write("\n# MHS in MZV monomial basis\n")
    for s in _mhs_to_mzv:
        thing = _mhs_to_mzv[s]
        a.write(str(s) + ": " + str(thing) + "\n")
    a.write("\n# MHS basis\n")
    for s in _mhs_basis:
        a.write(str(s) + "\n")
    a.write("\n# MZV monomials to MHS basis\n")
    Qi = mat(Q)**(-1)
    for i in range(Qi.rows):
        x = Qi.row(i)
        a.write(str(list(x)) + "\n")
    a.write("\n")
    a.close()
    return None

# Set up a distributed version of make_data
def prep_dist(W, n):
    a = []
    make_dir("mhs_data", "dist")
    for i in range(n):
        a.append(open(join(MHS_DATA_DIR,"dist","mhs_comps_"+str(W)+"_"+str(i)),"w"))
    i = 0
    for s in comps_i(W):
        a[i].write(str(s) + '\n')
        i = (i + 1) % n
    for i in range(n):
        a[i].close()
    print("Ready for distributed computing on %i cores" % n)

# Perform distrubuted version of make_data
def dist(W, i, n):
    print("Starting computation %i of %i" % (i, n))
    D = read_all(W)
    renormalize(D, W)
    mod_z2(D)
    B = get_gens_D(D)
    BB = get_vs_basis(B, W)
    convert(D, B, BB)
    E = {}
    count = 0
    num = (2 ** W - 2) / n
    target = 1
    print("MZV data read successfully")
    print("Computing expansion for MHS in terms of MZV (this may take some time)")
    b = open(join(MHS_DATA_DIR, "dist","mhs_comps_" + str(W) + "_" + str(i)))
    a = open(join(MHS_DATA_DIR, "dist","mhs_data_" + str(W) + "_" + str(i) + ".dat"), "w")
    for L in b:
        a.write(L[:-1] + ": " + str(mhs_to_mzv(eval(L), W, D, BB)) + "\n")
        count += 1
        if count * 100 > target * num:
            print("%i%% done" % target)
            a.flush()
            target += 1
    a.close()
    b.close()
    print("Done with computation %i of %i" % (i, n))

# Create the data file mhs_data_W.dat
def finish_dist(W, n):
    D = read_all(W)
    renormalize(D, W)
    mod_z2(D)
    B = get_gens_D(D)
    BB = get_vs_basis(B, W)
    convert(D, B, BB)
    E = {}
    count = 0
    num = 2 ** W - 2
    target = 1
    for i in range(n):
        q = open(join(MHS_DATA_DIR, "dist","mhs_data_" + str(W) + "_" + str(i) + ".dat"))
        for L in q:
            E[eval(L.split(": ")[0])] = [Rat(x) for x in L.split(": ")[1][1:-2].split(",")]
        q.close()
    L, M = find_mhs_basis(E, W)
    a = open(join(MHS_DATA_DIR, "mhs_data_%i.dat" % W), 'w')
    a.write("# Weight\n%i\n\n" % W)
    a.write("# MHS basis\n")
    b = [L[i] for i in M[1]]
    b.reverse()
    for i, x in enumerate(b):
        if i >= 1:
            a.write(", ")
        a.write(str(x))
    a.write("\n\n# MHS in terms of basis\n")
    for i in range(len(L) - 1, -1, -1):
        b = list(M[0].col(i))
        b.reverse()
        a.write(str(L[i]) + ": " + str(b) + "\n")
    a.write("\n# MZV algebra basis\n")
    for i in range(len(B)):
        a.write(str(B[i]))
        if i < len(B) - 1:
            a.write(", ")
    a.write("\n\n# MZV monomial basis\n")
    for x in BB:
        a.write(str(x) + '\n')
    a.write("\n# MZV in terms of monomial basis\n")
    LL = list(D.keys())
    def kk(s):
        return [sum(s), -len(s)] + list(s)
    LL.sort(key = kk)
    for x in LL:
        a.write(str(x) + ": " + str(D[x]) + "\n")
    a.write("\n# MHS basis in terms of monomial basis")
    b = list(M[1])
    b.reverse()
    for i in b:
        a.write('\n' + str(L[i]) + ": " + str(E[L[i]]))
    a.write('\n')
    a.close()
    print("File 'mhs_data_%i.dat' created successfully!!" % W)
    return None

def make_dir(*L):
    if not isdir(join(*L)):
        makedirs(join(*L))

def is_Lyndon_word(s):
    k = 0
    while k < len(s):
        i,j = k,k+1
        while j < len(s) and s[i] <= s[j]:
            i = (s[i] == s[j]) and i+1 or k     # Python cond?yes:no syntax
            j += 1
        if k < i + 1 and j - i == len(s):
            return True
        return False

def Lyndon(W):
    for s in comps_i(W):
        if is_Lyndon_word(s):
            yield s

def lists_weighted(W, L):
    if len(L) == 0:
        yield []
    else:
        for A in lists_weighted(W, L[:-1]):
            for i in range((W - sum(A[j] * L[j] for j in range(len(L) - 1))) / L[-1] + 1):
                yield A + [i]

def tuple_to_i(s):
    return sum(tuple_to_word(s) * 2 ** i for i in range(sum(s)))

def write(fi,s,description=None,key=None):
    fi.write("# ")
    if description is not None:
        fi.write(description)
    fi.write("\n")
    if isinstance(s,dict):
        K = list(s.keys())
        if key is not None:
            K.sort(key=key)
        for x in K:
            fi.write(str(x))
            fi.write(": ")
            fi.write(str(s[x]))
            fi.write("\n")
    elif isinstance(s,(list,tuple)):
        for x in s:
            fi.write(str(x) + "\n")
    else:
        fi.write(str(s) + "\n")
    fi.write("\n")
    return None

def make_rat(s):
    T = ""
    num = False
    for x in s:
        if num and x not in [str(i) for i in range(10)]:
            T += ')'
            num = False
        T += x
        if x == '/':
            T += "Rat("
            num = True
    if num:
        T += ')'
    return T

def read(fi):
    L = fi.readline()
    if L[0] == '#':
        L = fi.readline()
    T = ""
    while L != "\n" and L != "":
        T += make_rat(L[:-1]) + ","
        L = fi.readline()
    try:
        return eval('[' + T + ']')
    except:
        return eval('{' + T + '}')
