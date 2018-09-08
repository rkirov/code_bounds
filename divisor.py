class TwoPDivisor():
    def __init__(self, P=None, Q=None, deg=None, curve={}):
        '''
        Two point divisor class up to equivalence.
            >>> D = TwoPDivisor(P=2, Q=1) 
            >>> D = TwoPDivisor(P=2, deg=3) 
            >>> E = TwoPDivisor(3, 1) 
            >>> D += E 
            >>> D -= E + E 
            >>> print D 
            -P
        '''
        if sum(map(lambda x: x is not None, [P, Q, deg])) > 1:
            if P is None:
                P = deg - Q
            if Q is None:
                Q = deg - P
            if deg is None:
                deg = P + Q
        if P is None and Q is None and deg is not None:
            P = deg
            Q = 0
        self.P = P if P else 0
        self.Q = Q if Q else 0
        self.deg = deg if deg else P + Q 
        self.curve = curve
    
    def deg(self):
        return self.deg

    def dict_rep(self):
        return {'P': self.P, 'Q': self.Q}

    def __str__(self):
        l = []
        d = self.dict_rep()
        points = d.keys()
        points.sort() 
        points = map(lambda p: str(d[p]) + p, points)
        def sign(old, new):
            if new[0] == '0':
                return old
            elif new[0] == '-':
                return old + ' - ' + new[1:]
            else:
                return old + ' + ' + new
        return reduce(sign, points)

    def __add__(self, other):
        return TwoPDivisor(P=self.P + other.P, Q=self.Q + other.Q, curve=self.curve)

    def __sub__(self, other):
        return TwoPDivisor(P=self.P - other.P, Q=self.Q - other.Q, curve=self.curve)

    def __rmul__(self, other):
        return TwoPDivisor(P=other * self.P, Q=other * self.Q, curve=self.curve)

    def make_std(self):
        self.Q = self.Q % self.curve.m
        self.P = self.deg - self.Q

    def __iadd__(self, other):
        self.P += other.P
        self.Q += other.Q
        self.deg += other.deg
        return self 

    def __isub__(self, other):
        self.P -= other.P
        self.Q -= other.Q
        self.deg -= other.deg
        return self 

    def __getstate__(self):
        return {'P': self.P, 'Q': self.Q}

    def __setstate__(self, state):
        self.P = state['P']
        self.Q = state['Q']

    def __neg__(self):
        return TwoPDivisor(P=-self.P, Q=-self.Q, curve=self.curve)

    def less_then(self, other):
        '''Up to equality.'''
        d = other - self
        d.make_std()
        if d.P == 0 and d.Q == 0:
            return False 
        if d.P >= 0: 
            return True 
        else:
            return False

    def __eq__(self, other):
        d = self - other
        d.make_std()
        return d.P == 0 and d.Q == 0

    def copy(self):
        return TwoPDivisor(P=self.P, Q=self.Q, curve=self.curve)

    def l(self):
        '''
        l(D) - dimension of the Riemann-Roch space of D.
        '''
        if self.deg < 0:
            return 0
        if self.deg > 2 * self.curve.g - 2:
            return self.deg - self.curve.g + 1
        return self.curve.LVAL[self.deg][self.Q % self.curve.m] 

    def floor(self):
        '''
        TODO: reference maharaj-etal.
        '''
        if self.deg < 0:
            raise ValueError('Floor of negative divisor is -inf, which is not implemented yet')
        if self.deg > 2 * self.curve.g - 2:
            tup = (self.deg, self.Q)
        else:
            tup = self.curve.FL[self.deg][self.Q % self.curve.m] 
        return TwoPDivisor(deg=tup[0], Q=tup[1], curve=self.curve)

    def is_P_nongap(self):
        return self.deg >= self.curve.dvalp[self.P % self.curve.m]

    def is_Q_nongap(self):
        return self.deg >= self.curve.dvalq[self.Q % self.curve.m]

    def is_nongap(self, point):
        if point == 'P':
            return self.is_P_nongap()
        elif point == 'Q':
            return self.is_Q_nongap()

    def to_tuple(self):
        return (self.deg, self.Q)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
