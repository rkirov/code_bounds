import pylab as p

def point_plot(tuples, style, **kwds):
    ks = map(lambda t:t[0], tuples)
    ds = map(lambda t:t[1], tuples)
    p.plot(ds, ks, style, **kwds)

point_plot(best.c.get_best(), 'wo', label='Two point codes')
point_plot(best.c1.get_best(), 'kx', label='One point codes')
p.ylabel('Dimension')
p.xlabel('Distance')
#p.title('Two point vs One point codes')
p.legend()
p.show()
