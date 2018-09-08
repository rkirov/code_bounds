from sys import stdout as out 
from functools import partial
from operator import and_, or_
import pickle

def css_table(c, one_p=None):
    '''Return the best quantum codes using the CSS construciton 
    for all possible distances'''
    best = {'sym pure':{}, 'asym pure':{}, 'sym': {}, 'asym': {} }
    witness = {'sym pure':{}, 'asym pure':{}, 'sym': {}, 'asym': {} }
    for degC in xrange(c.MAXDEG):
        if one_p is None:
            Citer = xrange(c.m)
        else:
            Citer = [one_p] 
        for Cq in Citer:
            for degD in xrange(degC + 1, c.MAXDEG):
                if one_p is None:
                    Diter = xrange(c.m)
                else:
                    Diter = [one_p] 
                for Dq in Diter:
                    if not is_subcode(degC, Cq, degD, Dq, c):
                        continue
                    kq = k(degD, Dq, c) - k(degC, Cq, c) 
                    if kq == 0:
                        continue
                    #impure
                    x, y = d_impure(degC, Cq, degD, Dq, c), d_perp_impure(degC, Cq, degD, Dq, c)
                    if y > x:
                        x, y = y, x
                    if y > 1:
                        if best['sym'].get(kq, 0) < min(x, y):
                            best['sym'][kq] = min(x, y)
                            witness['sym'][kq] = (degC, Cq, degD, Dq)
                        if best['asym'].get((kq,x), 0) < y:
                            best['asym'][(kq, x)] = y
                            witness['asym'][(kq, x)] = (degC, Cq, degD, Dq)
                    # pure 
                    x, y = d(degD, Dq, c), d_perp(degC, Cq, c) 
                    if y > x:
                        x, y = y, x
                    if y > 1:
                        if best['sym pure'].get(kq, 0) < min(x, y):
                            best['sym pure'][kq] = min(x, y)
                            witness['sym pure'][kq] = (degC, Cq, degD, Dq)
                        if best['asym pure'].get((kq,x), 0) < y:
                            best['asym pure'][(kq, x)] = y
                            witness['asym pure'][(kq, x)] = (degC, Cq, degD, Dq)
    for key in best.keys():
        best[key], witness[key] = prune(best[key], witness[key])
    return best, witness

def tuple_compare(a, b):
    '''True if all components are as good or better, and at least one is
    strictly better'''
    return reduce(and_, [a[i] >= b[i] for i in range(len(a))]) and not a == b 

def prune(table, wtable):
    if isinstance(table.keys()[0], int):
        to_tup = lambda a, b: (a, b)
        to_key = lambda a: a[0]
    if isinstance(table.keys()[0], tuple):
        to_tup = lambda a, b: a + (b,)
        to_key = lambda a: tuple(list(a)[:-1])
    L = [to_tup(i, j) for i, j in table.iteritems()]
    to_be_removed = set([])
    for i in L:
        if i in to_be_removed:
            continue
        for j in L:
            if j in to_be_removed:
                continue
            if tuple_compare(i, j):
                to_be_removed.add(j)
    witness = {}
    best = []
    for i in L: 
        if i in to_be_removed:
            continue
        witness[i] = wtable[to_key(i)]
        best.append(i)
    return best, witness 

def seek_code(code, c, one_p = None):
    #assuming asym impure
    L = []
    for degC in xrange(c.MAXDEG):
        if one_p is None:
            Citer = xrange(c.m)
        else:
            Citer = [one_p] 
        for Cq in Citer:
            degD = degC + code[0]
            if one_p is None:
                Diter = xrange(c.m)
            else:
                Diter = [one_p] 
            for Dq in Diter:
                if not is_subcode(degC, Cq, degD, Dq, c):
                    continue
                if d_impure(degC, Cq, degD, Dq, c) == code[1] and \
                d_perp_impure(degC, Cq, degD, Dq, c) == code[2]:
                    d1 = dac(degC, Cq, c)
                    d2 = dac(degD, Dq, c)
                    if d1[1] <= d1[2] and d2[1] <= d2[2]:
                        L.append((d1,d2))
    return L

def dac(degC, Cq, c):
    Cp = degC - Cq
    return (Cp / c.m, Cp % c.m, Cq)

#############################
# Hermitian Self Dual Const #
#############################

def is_in_herm_dual(degC, Cq, c):
    # herm self-dual if (q+1)C <= K + D
    delta_deg = 2 * c.g - 2 + c.n - 2 - (c.q0 + 1) * degC
    delta_q =  (-1 - (c.q0 + 1) * Cq) % c.m 
    return delta_q <= delta_deg 

def herm_table(c, one_p=None):
    '''Return the best quantum codes using the dual hermitian construction 
    for all possible distances'''
    # no asymetric for this const.
    best = {'pure':{}, 'impure': {}}
    witness = {'pure':{}, 'impure': {}}
    for degC in xrange(c.MAXDEG):
        if one_p is None:
            Citer = xrange(c.m)
        else:
            Citer = [one_p] 
        for Cq in Citer:
            if not is_in_herm_dual(degC, Cq, c):
                continue
            kq = c.n - 2 - 2 * k(degC, Cq, c) 
            if kq <= 0:
                continue
            degD, Dq = dual(degC, Cq, c)
            x = d_impure(degC, Cq, degD, Dq, c)
            if x < 2:
                continue
            if best['impure'].get(kq, 0) < x:
                best['impure'][kq] = x 
                witness['impure'][kq] = (degC, Cq)
            # pure 
            x = d(degD, Dq, c)
            if best['pure'].get(kq, 0) < x: 
                best['pure'][kq] = x 
                witness['pure'][kq] = (degC, Cq)
    for key in best.keys():
        best[key], witness[key] = prune(best[key], witness[key])
    return best, witness

def tee(file, string):
    '''Print to file and to screen. See "man tee"'''
    out.write(string + '\n')
    file.write(string + '\n')

def better_than_anything(code, family):
    return reduce(or_, [tuple_compare(code, i) for i in family])

def tuple_diff_sum(a, b):
    return sum([a[i] - b[i] for i in range(len(a))])

def closest_in_family(code, family):
    L = [(tuple_diff_sum(code,i), i) for i in family if tuple_compare(code, i)]
    if L == []:
        return None, None
    return min(L)

def comp(b1, b2, w1, w2, c, write = False):
    ''' Symetricly compare the best codes from family b1 and b2'''
    if write:
        file = open(c.family + str(c.q) + 'css.out', 'w')
        printf = partial(tee, file) 
    else:
        printf = lambda x: out.write(x + '\n')
    good = {} 
    for code_type in b1.keys():
        good[code_type] = []
        printf(code_type)
        f1 = b1[code_type]
        f2 = b2[code_type]
        wit = w1[code_type]
        for i in range(1,3):
            printf('codes from from family%s better all of family%s.' % (i, i % 2 +1))
            for code in f1:
                dist, code2 = closest_in_family(code, f2)
                if dist:
                    printf('%s \tClosest is %13s \td: %s \nWitness: %s ' %\
                           (code, code2, dist, seek_code(code,c)))
                    good[code_type].append(code)
            f1, f2 = f2, f1
            wit = w2[code_type]
        printf('------------')
    return good

'''
b, w = css_table(curve)
b0, w0 = css_table(curve, one_p=0)
# ...
curve_comp(curve,write=True)
'''

def curve_comp(c, type='css', write=False):
    if type == 'css':
        f = css_table
    if type == 'herm':
        f = herm_table
    b, w = f(c) 
    b_small, w_small = f(c, one_p=0)
    if write:
        pickle.dump([b, w], open(c.family + str(c.q) + type + 'twoptQ', 'w')) 
        pickle.dump([b_small, w_small], open(c.family + str(c.q) + type + 'oneptQ', 'w')) 
    comp(b,b_small,w,w_small,c,write=write)
#################
# old functions #
#################


def quantp(degA, Aq, c):
    '''
    Return [k(A), d(A), d^perp(A)]
    '''
    degD = c.n - 2 
    Dq = - 1
    return (k(degA, Aq, c), d(degA, Aq, c), d_perp(degA, Aq, c))

def allquantdeg(c):
    '''
    Returns all (k, d, d_perp) pairs  
    '''
    return [[quantp(deg,i,c) for i in range(c.m)] for deg in range(len(c.DDK))]

def skew(i, j, c):
    '''
    Return [k(A), d(A), d^perp(A)]
    for inputs i = A_P and j = A_Q
    '''
    deg = j + i
    if deg not in range(len(c.DDK)):
        return '*' * 8
    else:
        return quantp(deg, j, c)

def allquant(c):
    '''
    Fill a table with [l(A), d(A), d^perp(A)]
    for all two-point divisors in RR range
    '''
    L = [[skew(i, j, c) for j in range(c.m)] for i in range(len(c.DDK))]
    fix(L)
    return L

def fix(L):
    '''
    Fix for inconsistent distantces
    Should this fix anything?!?
    '''
    m = len(L[0])
    for d in range(len(L)):
        for i in range(m):
            if d - i < 0:
                continue
            try:
                cur = L[d - i][i]
                next = L[d - i + 1][i]
                if cur[0] == next[0]:
                    L[d - i + 1][i] = (cur[0],max(cur[1],next[1]),max(cur[2],next[2]))
                next = L[d - i][i + 1]
                if cur[0] == next[0]:
                    L[d - i][i + 1] = (cur[0],max(cur[1],next[1]),max(cur[2],next[2]))
            except IndexError:
                pass

def worse(a, b):
    ''' 
    Helper function. Returns true if at least one component is worse.
    '''
    return (a[0]==b[0] and a[1]<b[1]) or (a[1]==b[1] and a[0]<b[0])

def bestq(c):
    ''' 
    Go through all pair of codes, and consturct pure CSS quantum code 
    '''
    D = [set([]) for i in range(len(c.DDK))] 
    L = allquant(c)
    for i in range(len(L)):
        for j in range(len(L[0])):
            for i2 in range(i, len(L)):
                for j2 in range(j, len(L[0])):
                    a = L[i][j]
                    b = L[i2][j2]
                    if a == '********' or b == '********':
                        continue
                    k = b[0] - a[0]
                    d1 = a[2] 
                    d2 = b[1] 
                    if d1 < d2:
                        d1, d2 = d2, d1
                    a = (d1,d2)
                    good = True
                    for b in D[k]:
                        if worse(a,b):
                            good = False
                            break
                    if good:
                        D[k].add((d1,d2))
    N = {} 
    for i in range(2*c.g + 1):
        N[i] = []
    for k, d in enumerate(D):        
        for p in d:
            defect = c.n - k - sum(p)
            N[defect].append([k] + list(p))
    return D, N

def take_until(d, L):
    i = 0
    while L[i][0] < d and i < len(L) - 1:
        i += 1
    return L[i][1]

def design(curve):
    degk = (curve.q0-2)*(curve.q0+1)
    p0 = [(curve.DDK[degk + i][0], curve.n - 2 - k(degk + i , 0, curve)) for i in range(len(curve.DDK) - degk)]
    p1 = [(curve.DDK[degk + i][1], curve.n - 2 - k(degk + i, 1, curve)) for i in range(len(curve.DDK) - degk)]
    print p0, p1
    for d in range(curve.q0 + 1, p0[-1][0] + 1, 2):
        print d, take_until(d, p0), take_until(d, p1), curve.n - 1 - curve.g - take_until(d, p1) 
