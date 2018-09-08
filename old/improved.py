'''
Functions for computing parameters of improved codes.
'''

import BAO

def rtilda(d, curve):
    '''
    Redundancy of improved codes
    '''
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    L = [0 for i in range(m)]
    for deg in range(len(BP)):
        L = [min(L[i]+BP[deg][i],L[i-1]+BQ[deg][i-1]) for i in range(m)]
    return min(L)+1

def rtilda_table(d,curve):
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    L = [0 for i in range(m)]
    N = []
    for deg in range(len(BP)):
        N.append(L)
        L = [min(L[i]+BP[deg][i],L[i-1]+BQ[deg][i-1]) for i in range(m)]
    N.append(L)
    return min(L)+1,N
    
#for testing purposes  
def rtilda_debug(d,curve):
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    ALL = []
    L = [0 for i in range(m)]
    for deg in range(len(BP)):
        ALL.append(L)
        L = [min(L[i]+BP[deg][i],L[i-1]+BQ[deg][i-1]) for i in range(m)]
    ALL.append(L)
    print 'total_min',min(L)
def rtilda(d,curve):
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!= 0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    L = [0 for i in range(m)]
    for deg in range(len(BP)):
        L = [min(L[i]+BP[deg][i],L[i-1]+BQ[deg][i-1]) for i in range(m)]
    return min(L)+1
    cq = 0
    path = ['P'] 
#0 means P, 1 means Q, defaults to last one seen
    last = 0
    for deg in range(len(ALL)-2,0,-1):
        print deg,cq,
        if(ALL[deg][cq]+BQ[deg-1][cq-1]>ALL[deg][cq]+BP[deg-1][cq]):
            path.append('P')
            last = 0
        elif(ALL[deg][cq]+BQ[deg-1][cq-1]<ALL[deg][cq]+BP[deg-1][cq]):
            path.append('Q')
            last = 1
        else:
            path.append(path[-1])
        cq = (cq-last)%m
    return path[:0:-1]
    

def rtilda_one(d,curve):
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    best = rtilda(d,curve)
    good_turns = []
    for turn in range(curve.m):
        cq = 0
        r = 1
        deg = 0
        #take turn number steps in Q, then P
        for step in range(turn):
            r+ = BQ[deg][cq]
            cq = (cq+1)%m
            deg+ = 1
        while(deg<len(BP)):
            r+ = BP[deg][cq]
            deg+=1
        if r==best:
            good_turns.append(turn)
        if r<best:
            print 'Problem with turn at %i best was %i, but got %i' % (turn,best,r)
    print good_turns
    return good_turns[-1],r

#bug!    
def rtilda_two_turn(d,curve):
    m = curve.m
    BP = [[curve.CPDK[deg][i]<d and curve.CPDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CPDK))]
    BQ = [[curve.CQDK[deg][i]<d and curve.CQDK[deg][i]!=0 for i in range(m)] for deg in range(len(curve.CQDK))]
    L_zero_t = [0 for i in range(m)]
    L_one_t = [0 for i in range(m)]
    L_two_t = [0 for i in range(m)]
    for deg in range(len(BP)):
        L_two_t = [min(L_two_t[i]+BP[deg][i],L_one_t[i]+BP[deg][i]) for i in range(m)]
        L_one_t = [min(L_one_t[i-1]+BQ[deg][i],L_zero_t[i-1]+BQ[deg][i]) for i in range(m)]
        L_zero_t = [L_zero_t[i]+BP[deg][i] for i in range(m)]
    return [min(L_zero_t),min(L_one_t),min(L_two_t)]

def r(d,curve):
    deg = 0
    best = 10*curve.MAXDEG
    for deg in range(len(curve.DDK)):
        for i in range(curve.m):        
            if(curve.DDK[deg][i]>=d and curve.l(deg,i)<best):
                best = curve.l(deg,i)
    return best

def ronepoint(d,curve):
    m = curve.m
    deg = 2*curve.g-2+d+1
    while(deg>=len(curve.DDK)-1 or curve.DDK[deg][0]>=d):
        deg -= 1
        if deg == 0:
            return 1
    return curve.l(deg+1,0)

def rtildaonepoint(d,curve):
    #add check for Q if the curve is of Klein type
    BP = [curve.CPDK[deg][0]<d and curve.CPDK[deg][0]!=0 for deg in range(len(curve.CPDK))]
    return sum(BP)+1


def bestcodesd(curve):
    #need to add klein here
    if(curve.family=='suzuki'):
        n = curve.q**2-1
    elif(curve.family=='hermitian'):
        n = curve.q0**3-1
    else:
        n = curve.n
    bestcodes = [[n,n-rtilda(d,curve),d] for d in range(curve.MAXDEG)]
    return bestcodes

def bestcodesk():
    B = bestcodesd()
    BD = []
    n = B[0][0]
    for k in range(n):
        d = 0
        for r in B:
            if(r[1]>=k and r[2]>d):
                d = r[2]
        BD.append([n,k,d])
    return BD

def improved(curve,jump=1,impr=False):
    print 'd\t1PC\t1PI\t2PC\t2PI'
    L = []
    for d in range(1,curve.MAXDEG,jump):
        v = [d,ronepoint(d,curve),rtildaonepoint(d,curve),r(d,curve),rtilda(d,curve)]
        if impr:
            v.append(min(v[1:-1])-v[-1])
        L.append(v)
        print '\t'.join([str(i) for i in v])
    return L

#only for Hermitian
def improvedBAO(curve):
    print 'd\t2PC\t2PI\t1PC\t2PI\tBAOtest\t1Pim-2Pim\t2Pc-2Pim'
    for d in range(1,curve.MAXDEG,2):
        t = (d-1)//2
        print d,'\t',r(d,curve),'\t',rtilda(d,curve),'\t',ronepoint(d,curve),'\t',rtildaonepoint(d,curve),'\t',BAO.r(t,curve.q0)==ronepoint(d,curve) and BAO.rtilda(t,curve.q0)==rtildaonepoint(d,curve),'\t',rtildaonepoint(d,curve)-rtilda(d,curve),'\t',r(d,curve)-rtilda(d,curve),'\t',rtilda_one(d,curve)[0]
