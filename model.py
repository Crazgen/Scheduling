import pulp as lp

Clerk_m_s_h = 'Clerk minimum shifting hour'
Clerk_m_s = 'Clerk multiple shift'
Clerk_w_t_b = 'Clerk work time balance'
Clerk_w_d_b = 'Clerk work day balance'


def solve_single_floor(data):
    # constraints in data: data['constraint name'] = (is_involved, must_meet, priority, optional_data_dicts)
    hour_num, day_num, clerk_num = data['working hours'], data['wording days'], data['clerks']
    time_slots_num = hour_num * day_num
    clerks, time_slots, hours, days = (range(1, clerk_num+1), range(1, time_slots_num+1), range(1, hour_num+1),
                                       range(1, day_num+1))
    clerk_time = [(n, t) for n in clerks for t in time_slots]
    clerk_day = [(n, d) for n in clerks for d in days]

    prob = lp.LpProblem("Scheduling", lp.LpMinimize)
    obj = lp.LpAffineExpression()
    # This part for clerk constraints
    x_c = lp.LpVariable.dicts('Time on duty', clerk_time, cat='Binary')
    z_c = lp.LpVariable.dicts('Day on duty', clerk_day, cat='Binary')
    kesi_c = lp.LpVariable.dicts('Switch on duty', clerk_time, cat='Binary')
    # first is the relations between different variables
    for n, d in clerk_day:
        prob += z_c[(n, d)] <= lp.lpSum(x_c[(n, hour_num*(d-1)+h)] for h in hours)
        prob += z_c[(n, d)] >= lp.lpSum(x_c[(n, hour_num*(d-1)+h)] for h in hours)/(hour_num*1.0)
        prob += kesi_c[(n, hour_num*(d-1)+1)] == x_c[(n, hour_num*(d-1)+1)]
        for h in range(2, hour_num+1):
            prob += kesi_c[(n, hour_num*(d-1)+h)] >= x_c[(n, hour_num*(d-1)+h)] - x_c[(n, hour_num*(d-1)+h-1)]
    # constraint for the minimum shifting hour
    if data[Clerk_m_s_h][0]:
        if data[Clerk_m_s_h][1]:
            for n, d in clerk_day:
                prob += (data[Clerk_m_s_h][3]['minimum hour'] * z_c[(n, d)] <=
                         lp.lpSum(x_c[(n, hour_num*(d-1)+h)] for h in hours))
        else:
            cmsh = lp.LpVariable.dicts('Number of exceeding hours', clerk_day, 0)
            obj += lp.lpSum(cmsh[n_d] for n_d in clerk_day) * (10**data[Clerk_m_s_h][2])
            for n, d in clerk_day:
                prob += (data[Clerk_m_s_h][3]['minimum hour'] * z_c[(n, d)] <=
                         lp.lpSum(x_c[(n, hour_num * (d - 1) + h)] for h in hours) + cmsh[(n, d)])
    # constraint for no cross shifty
    if data[Clerk_m_s][0]:
        if data[Clerk_m_s][1]:
            for n, d in clerk_day:
                prob += lp.lpSum(kesi_c[(n, hour_num*(d-1)+h)] for h in hours) <= 1
        else:
            cms = lp.LpVariable.dicts('Number of cross shifting', clerk_day, 0)
            obj += lp.lpSum(cms[n_d] for n_d in clerk_day) * (10**data[Clerk_m_s][2])
            for n, d in clerk_day:
                prob += lp.lpSum(kesi_c[(n, hour_num * (d - 1) + h)] for h in hours) <= 1 + cms[(n, d)]
    # constraint for work time balance
    if data[Clerk_w_t_b][0]:
        clerk_wtb = data[Clerk_w_t_b][3]['Clerks included']
        dif_wtb = data[Clerk_w_t_b][3]['Difference bound']
        clerk_inter_wtb = [(i, j) for i in clerk_wtb for j in clerk_wtb if i != j]
        if data[Clerk_w_t_b][1]:
            for i, j in clerk_inter_wtb:
                prob += lp.lpSum(x_c[(i, t)] for t in time_slots) - lp.lpSum(x_c[(j, t)] for t in time_slots) <= dif_wtb
        else:
            cwtb = lp.LpVariable.dicts('Number of exceeding work hours', clerk_inter_wtb, 0)
            obj += lp.lpSum(cwtb[i_j] for i_j in clerk_inter_wtb) * (10**data[Clerk_w_t_b][2])
            for i, j in clerk_inter_wtb:
                prob += (lp.lpSum(x_c[(i, t)] for t in time_slots) - lp.lpSum(x_c[(j, t)] for t in time_slots) <=
                         dif_wtb + cwtb[(i, j)])
    # constraint for work day balance
    if data[Clerk_w_d_b][0]:
        clerk_wdb = data[Clerk_w_d_b][3]['Clerks included']
        dif_wdb = data[Clerk_w_d_b][3]['Difference bound']
        clerk_inter_wdb = [(i, j) for i in clerk_wdb for j in clerk_wdb if i != j]
        if data[Clerk_w_d_b][1]:
            for i, j in clerk_inter_wdb:
                prob += lp.lpSum(z_c[(i, d)] for d in days) - lp.lpSum(z_c[(j, d)] for d in days) <= dif_wdb
        else:
            cwdb = lp.LpVariable.dicts('Number of exceeding work days', clerk_inter_wdb, 0)
            obj += lp.lpSum(cwdb[i_j] for i_j in clerk_inter_wdb) * (10**data[Clerk_w_d_b][2])
            for i, j in clerk_inter_wdb:
                prob += (lp.lpSum(z_c[(i, d)] for d in days) - lp.lpSum(z_c[(j, d)] for d in days) <=
                         dif_wdb + cwdb[(i, j)])

    # add object
    prob += obj
    return None

