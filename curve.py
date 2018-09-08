from bounds.essentials import essentials_dispatcher
from bounds.floor import floor_bounds_dispatcher
from bounds.order import order_bounds_dispatcher
import json
from itertools import imap, product
from divisor import TwoPDivisor
from divtable import DivTable 

class Curve:
    required_keys = {"q" : "field of definition",
                     "g" : "genus",
                     "m" : "smallest m such that mP ~ mQ",
                     "dvalp" : "first degree where D is a P-gap fixing iQ where 0 <= i < m",
                     "dvalq" : "first degree where D is a Q-gap fixing iP where 0 <= i < m",
                     "n" : "number of rational points",
                     "kq" : "Q value of the canonical divisor K"}
    optional_keys = {"curve_eq" : "equation of the curve",
                     "P" : "coordinates of point P",
                     "Q" : "coordinates of point Q"}
    data_keys = {}

    def __init__(self, _data={}):
        self._data = _data
        
    def __getattr__(self, name):
        if name not in self._data:
            raise TypeError("Curve doesn't have %s. Try to call build_all first" % name)
        return self._data[name]

    def __str__(self):
        return 'Curve of genus %(g)s over F%(q)s' % self._data

    def __repr__(self):
        return 'Curve of genus %(g)s over F%(q)s' % self._data

    def load_from_file(self, filename):
        ''' Load the bare minimum needed to construct all the data'''
        f = open(filename)
        input_data = json.load(f)
        missing_keys = [i for i in Curve.required_keys.keys() if i not in input_data.keys()]
        if missing_keys != []:
            for key in missing_keys:
                print("The file misses key %s which should be %s." %(key, Curve.required_keys[key]))
            raise TypeError("The input file is not properly formated. Proper JSON is required.")
        self._data = input_data
        self._data["MAXDEG"] = 2 * self.g
        return self 

    def two_point_iter(self):
        '''Iterate over 2 point divisors up to equivalence -MAXDEG < degC < MAXDEG '''
        def _(tup):
            deg, q = tup
            div = self.div(deg=deg, Q=q)
            return div
        iter = imap(_, product(xrange(-self.MAXDEG, self.MAXDEG + 1),
                               xrange(self.m)))
        return iter 

    def one_point_iter(self):
        '''Iterate over 1 point divisors up to equivalence -MAXDEG < degC < MAXDEG '''
        def _(deg):
            div = self.div(deg=deg, Q=0)
            return div
        iter = imap(_, xrange(-self.MAXDEG, self.MAXDEG + 1))
        return iter

    def build_all(self):
        ''' Build all the data using all algorithms available.'''
        offset = self._data['offset'] = {}
        essentials = ['CLP', 'FL', 'LVAL']
        floor_bound_names = ['GOP', 'BP', 'LM', 'GST', 'ABZ']
        order_bound_names = ['B0', 'B', 'ABZP', 'DP', 'DK']
        for choice in essentials + floor_bound_names + order_bound_names:
            print 'Building %s' % choice
            if choice in essentials:
                names, data = essentials_dispatcher(self, choice)                   
            elif choice in floor_bound_names:
                names, data = floor_bounds_dispatcher(self, choice)        
            elif choice in order_bound_names:
                names, data = order_bounds_dispatcher(self, choice)
            for i, name in enumerate(names):
                self._data[name] = data[i]
                #TODO to be done properly the choice of offset has the be
                #passable to all algorithms. 
                if choice != 'DK':
                    offset[name] = (0, 0)
                else:
                    K = self.K()
                    offset[name] = (-K.deg, 0)

    def wrap_table(self, name):
        offset_tup = self._data['offset'][name]
        offset = TwoPDivisor(P=offset_tup[0], Q=offset_tup[1])
        if name.lower()[0] == 'd':
            self.__dict__[name.lower()] = DivTable(
                self._data[name], self, start=offset, min=1)
        else:
            self.__dict__[name.lower()] = DivTable(
                self._data[name], self, start=offset)

    def wrap_tables(self):
        '''Wraps the tables into DivTable and CodeTable objects'''
        floor_bound_names = ['GOP', 'BP', 'LM', 'GST', 'ABZ']
        order_bound_names = ['B0', 'B', 'ABZP', 'DP', 'DK']
        for choice in floor_bound_names:
            self.wrap_table('D' + choice)
        for choice in order_bound_names: 
            self.wrap_table('D' + choice)
            self.wrap_table('CQ' + choice)
            self.wrap_table('CP' + choice)

    def serialize(self):
        '''Serialize the curve to dict object'''
        return self._data 

    def deserialize(self, dict):
        self._data = dict
        return self

    def load(self, inputfile):
        import gzip
        if inputfile[-3:] == '.gz':
            f = gzip.open(inputfile, 'rb')
        else:
            f = open(inputfile, 'r')
        self.deserialize(json.load(f))
        return self

    def save(self, filename, plain=False):
        ''' Dump the Curve dict as a .json file '''
        jsonstring = json.dumps(self.serialize())
        if not plain:
            import gzip
            f = gzip.open(filename + ".gz", "w")
        else:
            f = open(filename, "w")
        f.write(jsonstring)
        print 'Finished writing %s' % f.name 
        f.close()
        return

    def div(self,**kwds):
        return TwoPDivisor(curve=self, **kwds)

    def map_div(self, func, mindeg=None, maxdeg=None):
        maxdeg = maxdeg if maxdeg else self.MAXDEG
        mindeg = mindeg if mindeg else 0 
        L = []
        for degA in range(mindeg, maxdeg + 1):
            N = []
            for Aq in range(self.m):
                div = TwoPDivisor(Q=Aq, deg=degA, curve=self)
                N.append(func(div, self)) 
            L.append(N)
        return L

    def fill_degree_table(self, update_func, deg_neg_one, maxdeg=None,
                          mindeg=None):
        ''' Generic function that iterates by increasing degree
            from 0 to curve.MAXDEG (usually set to 2g).
            At each iteration update_func is called with -P and -Q values
            from previous iteration.
            Input: update_func, a curve.m list representing the values at -1
            degree divisors.
        '''
        maxdeg = maxdeg if maxdeg else self.MAXDEG
        mindeg = mindeg if mindeg else 0 
        L = [deg_neg_one]
        for degA in range(mindeg, maxdeg + 1):
            N = []
            for Aq in range(self.m):
                div = TwoPDivisor(Q=Aq, deg=degA, curve=self)
                minus_p_val = L[-1][Aq]
                minus_q_val = L[-1][Aq - 1] # python does the right thing with negative indices
                N.append(update_func(div, minus_p_val, minus_q_val))
            L.append(N)
        return L[1:]

    def fill_degree_table_reverse(self, update_func, deg_max, maxdeg=None,
                                  mindeg=None):
        ''' Same as fill_degree_table but in reverse degree order (i.e. going
        from DEGMAX down to 0, and updating on +P and +Q values'''
        maxdeg = maxdeg if maxdeg else self.MAXDEG
        mindeg = mindeg if mindeg else 0 
        L = [deg_max]
        for degA in range(maxdeg, mindeg - 1, -1):
            N = []
            for Aq in range(self.m):
                div = TwoPDivisor(Q=Aq, deg=degA, curve=self)
                plus_p_val = L[-1][Aq]
                plus_q_val = L[-1][(Aq + 1) % self.m] # python does the right thing with negative indices
                N.append(update_func(div, plus_p_val, plus_q_val))
            L.append(N)
        L = L[1:]
        L.reverse()
        return L

    def reduce_table(self, update_func, start, maxdeg=None, mindeg=None):
        ''' same as fill_degree_table but don't keep the table'''
        maxdeg = maxdeg if maxdeg else self.MAXDEG
        mindeg = mindeg if mindeg else 0 
        L = start 
        for degA in range(mindeg, maxdeg + 1):
            N = []
            for Aq in range(self.m):
                div = TwoPDivisor(Q=Aq, deg=degA, curve=self)
                minus_p_val = L[Aq]
                minus_q_val = L[Aq - 1] # python does the right thing with negative indices
                N.append(update_func(div, minus_p_val, minus_q_val))
            L = N
        return L

    def K(self):
        return TwoPDivisor(deg=2 * self.g - 2, Q=self.kq, curve=self)

    def P(self):
        return TwoPDivisor(deg=1, P=1, curve=self)

    def Q(self):
        return TwoPDivisor(deg=1, Q=1, curve=self)

    def D(self):
        '''Two point-divisor equivalent to sum of all rational point minus
        P,Q'''
        return TwoPDivisor(deg=self.ddeg, Q=self.dq, curve=self) 

    def magma_code(self):
        '''Code to make the curve in magma, along with points P, Q'''
        
        return '''
        P3<x,y,z>:=PolynomialRing(GF(%(q)s),3);
        f:=%(curve_eq)s;
        C:=Curve(ProjectiveSpace(GF(%(q)s),2),Homogenization(f,z));
        pl := Places(C,1);
        P := Place(C!%(P)s);
        Q := Place(C!%(Q)s);
        Exclude(~pl, P);
        Exclude(~pl, Q);
        D := &+pl;
        ''' % self._data
