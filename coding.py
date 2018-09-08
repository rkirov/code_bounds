'''Utilities for experimenting with the codes.
Assuming curve has all relevant data attached already. 
'''

class Coding():
    def __init__(self, curve):
        self.c = curve
        self.D = curve.D()
        self.K = curve.K()

    def get_k(self, C):
        '''
        Return the dimension of the code with given C and D
        k = l(D-C) - L(-C)
        '''
        return (self.D - C).l() - (-C).l()

    def constructions(self, C):
        '''Given C, return correspoding CL(D, G) and 
        COmega(D, G) that realize that code
        Gl = D - C 
        Gomega = K + C
        '''
        return {'Omega': self.K + C, 'L': self.D - C} 

    def get_d(self, C):
        return self.c.ddk(C)

    def code(self, C):
        return {'n': self.D.deg, 'k': self.get_k(C), 'd': self.get_d(C)}

    def print_code(self, C):
        print 'D = ', self.D
        print 'C = ', C
        print '[%s, %s, %s] code over F%s' % (self.D.deg, self.get_k(C),
                                              self.get_d(C), self.c.q)
        for type, div in self.constructions(C).iteritems():
            print type, '\tusing G =', div

    def design_d(self, C):
        return C.deg

    def best_codes(self, one_point=False, redun=False):
        ''' Return a list of best tuples of classic codes in the (k,d)
        format.'''
        def f(div):
            return (self.get_k(div), self.c.ddk(div))
        best = BestTuple(f=f)
        iter = self.c.two_point_iter if not one_point else self.c.one_point_iter
        for div in iter():
            best.update(div)
        best.reduce()
        return best 

    def dual(self, C):
        '''Return C' such that C(D, C) is dual to C(D, C')''' 
        return self.D - self.K - C 
         
    def is_subcode(self, C, D):
        '''Check C(_,C) \subcode C(_,D) up to equivalence.'''
        return D.less_then(C)

    def mintw_coset(self, C, D):
        '''Return the minimum wt of the words in C_omega(C) not in C_omega(D)
        Assumes D <= C'''
        assert self.is_subcode(D, C) 
        d = D - C
        d.make_std()
        temp = C.copy()
        bound = CosetBound()
        for i in range(d.P):
            bound.update(self.c.cpdk(temp))
            temp += self.c.P()
        for i in range(d.Q):
            bound.update(self.c.cqdk(temp))
            temp += self.c.Q()
        return bound.b

    def quantum(self, C1, C2, impure=True, assym=True):
        C1d = self.dual(C1) 
        C2d = self.dual(C2) 
        assert self.is_subcode(C2d, C1)
        data = {'n': self.D.deg, 'k': self.get_k(C1) - self.get_k(C2d)}
        if impure:
            data.update({ 'dx': self.mintw_coset(C2, C1d),
                    'dz': self.mintw_coset(C1, C2d)})
        else:
            data.update({ 'dx': self.get_d(C2),
                    'dz': self.get_d(C1)})
        return data

    def best_quantum(self, one_point=False):
        def q_tuple(C, C2):
            q = self.quantum(C, C2)
            return (q['k'], q['dx'], q['dz'])
        if one_point:
            get_iter = self.c.ddk.one_point_iter
        else:
            get_iter = self.c.ddk.__iter__
        best_q = BestTuple(f=q_tuple)
        #TODO: use a generic divisor iterator, since d is not used directly
        #for C, d in self.c.ddk:
        #    for C2d, d2 in self.c.ddk:
        for C, d in get_iter():
            for C2d, d2 in get_iter():
                try:
                    C2 = self.dual(C2d)
                    best_q.update(C, C2)
                except AssertionError:
                    pass
        best_q.reduce()
        return best_q

class BestTuple():
    # can you do an online version?
    def __init__(self, f=None):
        '''Object to abstract the process of turning divisors (or pairs of
        divisors) to parameters of codes and storing the best ones. The input
        function f is applied to all parameters passed to method update'''
        self.ind = [] 
        self.witness = {} 
        self.f = f

    def update(self, *arg):
        a = self.f(*arg)
        self.ind.append(a)
        self.witness[a] = arg 

    def reduce(self):
        for i in range(len(self.ind[0])):
            self.reduce_coord(i)
        self.ind.sort()

    def reduce_coord(self, i):
        d = {}
        for tup in self.ind:
            key = list(tup)
            key.pop(i)
            key = tuple(key)
            if key in d and d[key][i] >= tup[i]:
                continue
            d[key] = tup
        self.ind = d.values()
        for k in self.witness.keys():
            if k not in self.ind:
                self.witness.pop(k)

    def get_best(self):
        return self.ind

    def save(self, file):
        import pickle
        pickle.dump({'ind': self.ind, 'witness': self.witness}, file) 

    def load(self, file):
        import pickle
        data = pickle.load(file)
        self.ind = data['ind']
        self.witness = data['witness']

    def closest(self, new):
        ''' return distance to a point covered by a best tuple, and that tuple.'''
        queue = [(new, 0)]
        while True:
            cur, d = queue.pop(0)
            if self.covered(cur):
                return cur, d
            # super inefficient, exp instead of poly
            # improve by guaranteeing next doesn't give repetitions
            for new in self.next(cur):
                if new not in queue:
                    queue.append((new, d + 1))
    
    def next(self, cur):
        l = list(cur)
        out = []
        for i in range(len(l)):
            new = l[:]
            new[i] -= 1
            out.append(tuple(new))
        return out 

    def covered(self, new):
        for i in self.ind:
            if self.covered_test(new, i):
                return True
        return False

    def covered_test(self, a, b):
        return all([a[i] <= b[i] for i in range(len(a))])

def defect(best, curve):
    n = curve.D().deg
    return [(tup, n - sum(tup) + 1) for tup in best]

class CosetBound():
    def __init__(self, start=0):
        self.b = start 

    def update(self, new):
        if self.b == 0:
            self.b = new
            return
        if new > 0:
            self.b = min(new, self.b)
