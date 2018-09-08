import families
import build
import orderbounds
import improved
import json
import gzip

def red(L,d):
    return len([i for i in L if 0<i<d])

q0 = 4

def save(curve, filename, zipped):
    jsonstring = json.dumps(curve.D)
    if zipped:
        f = gzip.open(filename + ".output.gz", "wb")
    else:
        f = open(filename + ".output", "wb")
    f.write(jsonstring)
    f.close()
    return

def makeherm(q0):
    curve = build.Curve(families.hermitian(q0))
    curve.D["MAXDEG"] = 2*curve.g
    print 'q0 =',q0,'building up to',curve.MAXDEG
    names, data = orderbounds.orderBounds(curve,'DK')
    for i,n in enumerate(names):
        curve.D[n]= data[i]
    return curve
    
q0 = 8
while q0 <= 256:
    curve = makeherm(q0)
    save(curve,'herm'+str(q0**2)+'dk',True)
    L = [curve.CPDK[i][1] for i in range(len(curve.CPDK))]
    for d in range(2,max(L)):
        if improved.rtilda(d,curve)!=red(L,d):
            print q0,d
            break
    q0 *= 2    
