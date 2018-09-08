# TODO: make 'if point' better
from functools import partial
from divtable import DivTable

def B0coset(*arg):
    return Bcoset_gen(*arg, B0=True)

def Bcoset(point, C, curve, B0=False):
    start = [0 for i in range(curve.m)]
    if point == 'P':
        def update(div, minus_p, minus_q):
            delta_p = div.is_P_nongap() and not (div - C).is_P_nongap()
            return minus_p + delta_p
    elif point == 'Q':
        def update(div, minus_p, minus_q):
            delta_q = div.is_Q_nongap() and not (div - C).is_Q_nongap()
            return minus_q + delta_q
    L = curve.reduce_table(update, start, mindeg=min(0, C.deg),
                           maxdeg=curve.MAXDEG + max(0, C.deg))
    return max(L) if not B0 else L[0]

def ABZPcoset(point, C, curve):
    start = [(0,0,0) for i in range(curve.m)]
    if point == 'P':
        def update(div, minus_p, minus_q):
            delta_p = div.is_P_nongap() and not (div - C).is_P_nongap()
            return (minus_p[0] + delta_p, max(minus_p[1], minus_p[2]) + delta_p,
                   max(minus_q[2], minus_q[0]))
    elif point == 'Q':
        def update(div, minus_p, minus_q):
            delta_q = div.is_Q_nongap() and not (div - C).is_Q_nongap()
            return (minus_q[0] + delta_q, max(minus_q[1], minus_q[2]) + delta_q,
                   max(minus_p[2], minus_p[0]))
    L = curve.reduce_table(update, start, mindeg=min(0, C.deg),
                           maxdeg=curve.MAXDEG + max(0, C.deg))
    return max(map(max, L))

def DPcoset(point, C, curve):
    start = [0 for i in range(curve.m)]
    if point == 'P':
        def update(div, minus_p, minus_q):
            delta_p = div.is_P_nongap() and not (div - C).is_P_nongap()
            delta_q = 0
            return max(minus_p + delta_p, minus_q + delta_q)
    elif point == 'Q':
        def update(div, minus_p, minus_q):
            delta_p = 0
            delta_q = div.is_Q_nongap() and not (div - C).is_Q_nongap()
            return max(minus_p + delta_p, minus_q + delta_q)
    L = curve.reduce_table(update, start, mindeg=min(0, C.deg),
                           maxdeg=curve.MAXDEG + max(0, C.deg))
    return max(L)

def DKsubset(C, curve):
    start = [0 for i in range(curve.m)]
    def update(div, minus_p, minus_q):
        delta_p = div.is_P_nongap() and not (div - C).is_P_nongap()
        delta_q = div.is_Q_nongap() and not (div - C).is_Q_nongap()
        return max(minus_p + delta_p, minus_q + delta_q)
    L = curve.reduce_table(update, start, mindeg=min(0, C.deg),
                           maxdeg=curve.MAXDEG + max(0, C.deg))
        #L = [max(L[i]+(curve.dvalp[degAm-i] <= degA and degAC < curve.dvalp[bunchM-i]),L[i-1]+(curve.dvalq[i] <= degA and degAC < curve.dvalq[(i-cq)%m])) for i in range(m)]
    return max(L)

def order_bounds_dispatcher(curve, choice):
    if choice == 'DK':
        return DKOrderBound(curve)
    elif choice == 'B0':
        coset = Bcoset
    elif choice == 'B':
        coset = Bcoset
    elif choice == 'ABZP':
        coset = ABZPcoset
    elif choice == 'DP':
        coset = DPcoset
    else:
        raise TypeError("Choose between 'B0','B','ABZP','DP'.")
    return coset_bound(curve, choice, coset)

def min_star(old, new):
    if new <= 0:
        return old
    else:
        return min(new, old)

def dist_from_cosets(coset_p_table, coset_q_table, curve, mindeg=0):
    def update(div, plus_p_value, plus_q_value):
        return max(min_star(plus_p_value, coset_p_table(div)),
                   min_star(plus_q_value, coset_q_table(div)))
    D = [curve.MAXDEG for i in range(curve.m)]
    return curve.fill_degree_table_reverse(update, D, mindeg=mindeg)

def cosets_from_subsets(subset_table, curve, mindeg=0):
    def update_P(div, plus_p_value, plus_q_value):
        return min(plus_q_value, subset_table(div))
    def update_Q(div, plus_p_value, plus_q_value):
        return min(plus_p_value, subset_table(div))
    D = [curve.MAXDEG for i in range(curve.m)]
    return (curve.fill_degree_table_reverse(update_P, D, mindeg=mindeg),
           curve.fill_degree_table_reverse(update_Q, D, mindeg=mindeg))

def coset_bound(curve, choice, coset):
    coset_P = partial(coset, 'P')
    CP = curve.map_div(coset_P)
    CP_table = DivTable(CP, curve) 
    coset_Q = partial(coset, 'Q')
    CQ = curve.map_div(coset_Q)
    CQ_table = DivTable(CQ, curve) 
    D = dist_from_cosets(CP_table, CQ_table, curve)
    return ['CP' + choice, 'CQ' + choice, 'D' + choice], [CP, CQ, D]

def DKOrderBound(curve):
    ''' Unlike other order bounds DK is computer for negative C too.'''
    m = curve.m
    min_deg = -curve.K().deg
    #FIXME this assumes kq=0 
    minusK = curve.div(deg=min_deg, Q=0)
    DKsubsets = curve.map_div(DKsubset)
    #DKN contains the negative subsets i.e. DKN[a][b] has degC = -a and Cq = b
    DKNsubsets = [[DKsubsets[deg][-i] - deg for i in range(m)] 
                  for deg in range(1, -min_deg + 1)]
    DKNsubsets.reverse()
    DKsubsets = DKNsubsets + DKsubsets
    DKsubsets_table = DivTable(DKsubsets, curve, start=minusK)
    Pcosets, Qcosets = cosets_from_subsets(DKsubsets_table, curve, mindeg=min_deg) 
    Pcosets_table = DivTable(Pcosets, curve, start=minusK)
    Qcosets_table = DivTable(Qcosets, curve, start=minusK)
    D = dist_from_cosets(Pcosets_table, Qcosets_table, curve, mindeg=min_deg)
    return ['SDK', 'CPDK', 'CQDK', 'DDK'], [DKsubsets, Pcosets, Qcosets, D]
