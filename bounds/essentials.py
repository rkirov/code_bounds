def l_values(c):
    minus_one_deg = [0 for _ in range(c.m)]
    def update(div, minus_p_val, minus_q_val):
        return minus_p_val + div.is_P_nongap() 
    lval = c.fill_degree_table(update, minus_one_deg)
    return ['LVAL'], [lval]

def build_floor_table(c):
    #if l(D) = 0, floor is assigned to be [0,1], XXX: why not [0,0] (2011)
    minus_one_deg = [(0,1) for _ in range(c.m)]
    def update(div, minus_p_val, minus_q_val):
        jumpP = div.is_P_nongap() 
        jumpQ = div.is_Q_nongap() 
        if jumpP and jumpQ:
            return div.to_tuple()
        elif jumpQ:
            return minus_p_val
        else:
            # this includes two cases
            return minus_q_val
    floor_table = c.fill_degree_table(update, minus_one_deg)
    return ['FL'], [floor_table]

#the table should be read T[deg][a] where a is the multiplicity of the Q
def build1DCeilingTable(c):
    '''entry for A is max k s.t. l(A) = l(A+kP) and l(A+kQ) '''
    max_deg = [0 for _ in range(c.m)]
    def update(div, plus_p_val, plus_q_val):
        div_plus_P = div + c.div(P=1,deg=1) 
        return (plus_p_val + 1) if not div_plus_P.is_P_nongap() else 0 
    CLP = c.fill_degree_table_reverse(update, max_deg)
    def update(div, plus_p_val, plus_q_val):
        div_plus_Q = div + c.div(Q=1,deg=1) 
        return (plus_q_val + 1) if not div_plus_Q.is_Q_nongap() else 0 
    CLQ = c.fill_degree_table_reverse(update, max_deg)
    return ['CLP','CLQ'], [CLP,CLQ]

def essentials_dispatcher(curve, choice):
    if choice == 'CLP' or choice == 'CLQ':
        return build1DCeilingTable(curve)
    elif choice == 'FL':
        return build_floor_table(curve)
    elif choice == 'LVAL':
        return l_values(curve)
