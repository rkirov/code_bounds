from curve import Curve

def suzuki(q0):
    '''
    TODO: add reference
    '''
    data = {}
    data['q0'] = q0
    data['family'] = 'suzuki'
    data['q'] = q = 2 * q0 ** 2
    data['m'] = m = q + 2 * q0+1
    data['g'] = g = q0 * (q - 1)
    #generating values for d(i) where i is taken mod m
    dlist=[]
    for a in range(-q0,q0 + 1):
        for b in range(-q0 + abs(a),q0-abs(a) + 1):
            dlist.append([(a * (q0 + 1) + b * q0-q0 * (q0 + 1))%m,
                          (q0-a) * (q-1)])
    dlist.sort()
    data['dvalp'] = data['dvalq'] = [i[1] for i in dlist]
    data['kq'] = 0
    data['n'] = q + 1 + 2 * g * q0
    data['curve_eq_latex'] = 'y^%(q)s + y = x^%(q0)s (x^%(q)s + x)' % data 
    data['curve_eq_code'] = 'y^%(q)s + y - x^%(q0)s * (x^%(q)s + x)' % data 
    data['P'] = '[0 , 0 , 1]'
    data['Q'] = '[1 , 0 , 0]'
    data['ddeg'] = q ** 2 - 1
    data['dq'] = -1 
    c = Curve()
    c._data = data
    return c 
    
def hermitian(q0):
    '''
    TODO: add reference to Stichtenoth paper
    '''
    data = {}
    data['q0'] = q0
    data['family'] = 'hermitian'
    data['q'] = q = q0 ** 2
    data['m'] = m = q0 + 1
    data['g'] = g = q0 * (q0 - 1) / 2
    #generating values for d(i) where i is taken mod m
    dlist=[]
    for d in range(m):
        dlist.append([(-d)%m,d*(q0-1)])
    dlist.sort()
    data['dvalp'] = data['dvalq'] = [i[1] for i in dlist]
    data['kq'] = 0
    data['n'] = q + 1 + 2*g*q0
    data['curve_eq_latex'] = 'y^{q + 1} = x^q + x' 
    data['curve_eq'] = 'y^(%(q0)s + 1) - x^%(q0)s - x' % data 
    data['P'] = '[0 , 0 , 1]'
    data['Q'] = '[1 , 0 , 0]'
    data['ddeg'] = q0 ** 3 - 1
    data['dq'] = -1 
    c = Curve()
    c._data = data
    return c 
