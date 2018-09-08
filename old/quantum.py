def GammaQ(degA,Aq,c):
    return degA>=c.dvalq[Aq%c.m]

def GammaP(degA,Aq,c):
    return degA>=c.dvalp[(degA-Aq)%c.m]

def conslval(c):
    lval=[[0 for i in range(c.m)]]
    for degA in range(c.MAXDEG+1):
        lval.append([lval[-1][i]+GammaP(degA,i,c) for i in range(c.m)])
    lval=lval[1:]
    return lval

def l(degA,Aq,lval,c):
    if(degA<0):
        return 0
    if(degA>2*c.g-2):
        return degA-c.g+1
    return lval[degA][Aq%c.m] 

def d(degA,Aq,c):
    if degA < 0:
        return 0
    if degA >= len(c.DDK):
        return degA-2*c.g+2
    else:
        return c.DDK[degA][Aq]

def quantp(degA,Aq,c,lval):
    #adjust depending on R = sum rational pointattrs
    degD = c.n-2 
    Dq = -1
    return (l(degA,Aq,lval,c),d(degD+2*c.g-2-degA,(Dq-Aq)%c.m,c),d(degA,Aq,c))

def allquantdeg(c):
    lval = conslval(c)
    return [[quantp(deg,i,c,lval) for i in range(c.m)] for deg in range(len(c.DDK))]

def skew(i,j,c,lval):
    deg = j+i
    if deg not in range(len(c.DDK)):
        return '*'*8
    else:
        return quantp(deg,j,c,lval)

def allquant(c):
    lval = conslval(c)
    L = [[skew(i,j,c,lval) for j in range(c.m)] for i in range(len(c.DDK))]
    fix(L)
    return L

def fix(L):
    m = len(L[0])
    for d in range(len(L)):
        for i in range(m):
            if d-i < 0:
                continue
            try:
                cur = L[d-i][i]
                next = L[d-i+1][i]
                if cur[0] == next[0]:
                    L[d-i+1][i] = (cur[0],max(cur[1],next[1]),max(cur[2],next[2]))
                next = L[d-i][i+1]
                if cur[0] == next[0]:
                    L[d-i][i+1] = (cur[0],max(cur[1],next[1]),max(cur[2],next[2]))
            except IndexError:
                pass

def worse(a,b):
    return (a[0]==b[0] and a[1]<b[1]) or (a[1]==b[1] and a[0]<b[0])

def bestq(c):
    D = [set([]) for i in range(len(c.DDK))] 
    L = allquant(c)
    for i in range(len(L)):
        for j in range(len(L[0])):
            for i2 in range(i,len(L)):
                for j2 in range(j,len(L[0])):
                    a = L[i][j]
                    b = L[i2][j2]
                    if a == '********' or b == '********':
                        continue
                    k = b[0] - a[0]
                    d1 = a[2] 
                    d2 = b[1] 
                    if d1 < d2:
                        d1,d2 = d2,d1
                    a = (d1,d2)
                    good = True
                    for b in D[k]:
                        if worse(a,b):
                            good = False
                            break
                    if good:
                        D[k].add((d1,d2))
    N = {} 
    for i in range(2*c.g+1):
        N[i] = []
    for k,d in enumerate(D):        
        for p in d:
            defect = c.n-k-sum(p)
            N[defect].append([k]+list(p))
    return D,N
