def LMbound(C, curve):
    m = curve.m
    K = curve.K()
    bestc = 0
    for degA in range(0, curve.MAXDEG + 1):
        for aq in range(m):
            A = curve.div(deg=degA, Q=aq)
            floor = A.floor()
            flq = (aq - floor.Q) % m
            flp = (A.deg - floor.deg - flq) % m
            B = K + C - A
            degZ = 0
            rq = min(flq, curve.CLQ[min(B.deg, curve.MAXDEG)][B.Q])
            for dq in range(1 + rq):
                if 0 <= B.deg + dq  <=  K.deg:
                    dp = min(flp, curve.CLP[B.deg + dq][(dq + B.Q)%m])
                else:
                    dp = 0
                degZ = dq + dp
                if degZ > bestc:
                    bestc = degZ       
    return C.deg + bestc    

def GSTbound(C, curve):
    bestB = 0
    for degA in range(curve.MAXDEG + 1):
        for aq in range(curve.m):
            A = curve.div(deg=degA, Q=aq)
            floor = A.floor() 
            newval = floor.l() - (floor - C).l() - (A.l() - (A - C).l())
            if newval > bestB:
                bestB = newval
    return C.deg + bestB

def build_best_table(C, curve):
    '''Table of best f(A') = l(A') - l(A' - C), where A' <= A, for each A.'''
    def update(div, minus_p_val, minus_q_val):
        return max([minus_p_val, minus_q_val, div.l() - (div - C).l() ]) 
    return curve.fill_degree_table(update, [0 for _ in range(curve.m)],
                                   maxdeg=curve.MAXDEG + C.deg)

def ABZbound(C, curve):
    m = curve.m
    K = curve.K()
    Best = build_best_table(C, curve)
    def f(A, curve):
        B = K + C - A
        B.make_std()
        return Best[A.deg][A.Q] + Best[B.deg][B.Q]
    val = curve.map_div(f, maxdeg=K.deg + C.deg)
    return max(map(max, val))

def floor_bounds_dispatcher(curve, choice):
    m = curve.m
    if choice == 'ABZ':
        func = ABZbound
    elif choice == 'LM':
        func = LMbound
    elif choice == 'GST':
        func = GSTbound
    elif choice == 'GOP':
        def func(div, curve):
            return div.deg
    elif choice == 'BP':
        def func(div, curve):
            add = (curve.FL[div.deg][div.Q][0] < div.deg 
                   or curve.FL[div.deg][div.Q] == (0, 1))
            return div.deg + add 
    else:
        print "choose between 'DP', 'LM', 'GST', 'GOP', 'BP'"
        return
    data = curve.map_div(func)
    return ['D' + choice], [data]
