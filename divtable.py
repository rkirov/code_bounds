from itertools import product, imap
class DivTable():
    '''A table which takes a divisor and outputs a number.'''
    def __init__(self, table, curve, start=None, min=0):
        '''start is the divisor corresponding to offset'''
        self.start = start if start else curve.div(deg=0, Q=0) 
        self.t = table
        self.s = start
        self.curve = curve
        self.min = min # minimum value that makes sense for the
        # table in quesiton. 1 for distances, 0 for cosets

    def __call__(self, div):
        off = div - self.start
        try:
            if off.deg < 0:
                raise IndexError
            return self.t[off.deg][off.Q % self.curve.m]
        except IndexError:
            return self.out_of_bounds(div)

    def out_of_bounds(self, div):
        ''' Return a meaningful answer of out of bound query'''
        return max(self.min, min(self.curve.D().deg, div.deg))

    def __iter__(self):
        def _(tup):
            deg, q = tup
            div = self.curve.div(deg=deg, Q=q) + self.start 
            return div, self(div)
        iter = imap(_, product(xrange(len(self.t)), xrange(self.curve.m)))
        return iter 

    def one_point_iter(self):
        def _(deg):
            div = self.curve.div(deg=deg, Q=0) + self.start 
            return div, self(div)
        iter = imap(_, xrange(len(self.t)))
        return iter 
