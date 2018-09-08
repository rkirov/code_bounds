class BestCodes:
    def __init__(self, codes):
        self.codes = codes
        self.curve = codes.c 
        self.base_name = 'data/' + self.curve.family + str(self.curve.q)

    def compute_quantum_and_save(self): 
        q = self.codes.best_quantum()
        q1p = self.codes.best_quantum(one_point = True)
        q.save(open(self.base_name + '.q', 'w'))
        q1p.save(open(self.base_name + '.q.1p', 'w'))

    def compute_classic_and_save(self): 
        c = self.codes.best_codes()
        c1p = self.codes.best_codes(one_point = True)
        c.save(open(self.base_name + '.c', 'w'))
        c1p.save(open(self.base_name + '.c.1p', 'w'))

    def load_best_codes(self):
        import coding
        self.q = coding.BestTuple()
        self.q1 = coding.BestTuple()
        self.q.load(open(self.base_name + '.q'))
        self.q1.load(open(self.base_name + '.q.1p'))
        self.c = coding.BestTuple()
        self.c1 = coding.BestTuple()
        self.c.load(open(self.base_name + '.c'))
        self.c1.load(open(self.base_name + '.c.1p'))

def best_k_for_design_d(d, tuples):
    best = 0 
    for tuple in tuples:
        if tuple[1] >= d and tuple[0] > best:
            best = tuple[0]
    return best

def make_table(ds, best):
    table = [(d, best_k_for_design_d(d, best.c1.get_best()),
best_k_for_design_d(d, best.c.get_best())) for d in ds]
    return table

def latex_table(table):
    half_done = map(lambda x: '\t&\t'.join(map(str,x)), table)
    done = '\\\\ \n'.join(half_done)
    return done

def find_impr(A, B, opt_test=(lambda witness, dist: True)):
    ''' Finds codes from A that improve over B'''
    for code in A.get_best():
        witness, dist = B.closest(code)
        if opt_test(witness, dist):
            print code, 'closest', witness, 'distance', dist 
