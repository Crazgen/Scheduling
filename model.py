# encoding=utf-8
import pulp as lp

Clerk_m_s_h = 'Clerk minimum shifting hour'
Clerk_m_s = 'Clerk multiple shift'
Clerk_w_t_b = 'Clerk work time balance'
Clerk_w_d_b = 'Clerk work day balance'
Clerk_m_o_d = 'Clerk minimum off days'
Clerk_m_w_d = 'Clerk most whole day'
Clerk_m_r = 'minimum clerks'
Clerk_m_t = 'Clerk minimum work time'
Clerk_o_t = 'Clerk most over time'
Clerk_i_d = 'Clerks included'
Clerk_c_o = 'Clerk continuous off days'
Clerk_w_b = 'Clerk weekly workday balance'
Clerk_e_l_s = 'Clerk early late switch'


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
        clerk_wtb = data[Clerk_w_t_b][3][Clerk_i_d]
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
        clerk_wdb = data[Clerk_w_d_b][3][Clerk_i_d]
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
    # constraint for minimum off days in the period
    if data[Clerk_m_o_d][0]:
        clerk_mod = data[Clerk_m_o_d][3][Clerk_i_d]
        mod = data[Clerk_m_o_d][3]['Min off day']
        if data[Clerk_m_o_d][1]:
            for n in clerk_mod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) <= day_num - mod
        else:
            cmod = lp.LpVariable.dicts('Number of exceeding days', clerk_mod, 0)
            obj += lp.lpSum(cmod[n] for n in clerk_mod) * (10**data[Clerk_m_o_d][2])
            for n in clerk_mod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) <= day_num - mod + cmod[n]
    # constraint for 最多连续上全班天数
    if data[Clerk_m_w_d][0]:
        clerk_mwd = data[Clerk_m_w_d][3][Clerk_i_d]
        mwd = data[Clerk_m_w_d][3]['Most whole day'] * hour_num
        mwd_time_slots = range(1, time_slots_num - mwd + 2)
        clerk_time_mwd = [(n, t) for n in clerk_mwd for t in mwd_time_slots]
        if data[Clerk_m_w_d][1]:
            for n, t in clerk_time_mwd:
                prob += lp.lpSum(x_c[(n, t + k - 1)] for k in range(1, mwd+1)) <= mwd - 1
        else:
            cmwd = lp.LpVariable.dicts('Number of exceeding whole days', clerk_time_mwd, 0)
            obj += lp.lpSum(cmwd[n_t] for n_t in clerk_time_mwd) * (10**data[Clerk_m_w_d][2])
            for n, t in clerk_time_mwd:
                prob += lp.lpSum(x_c[(n, t + k - 1)] for k in range(1, mwd + 1)) <= mwd - 1 + cmwd[(n, t)]
    # constraint for minimum clerks required, combined with traffic..
    if data[Clerk_m_r][0]:
        mini_clerks = data[Clerk_m_r][3]['minimum clerk num']
        if data[Clerk_m_r][1]:
            for t in time_slots:
                prob += lp.lpSum(x_c[(n, t)] for n in clerks) >= mini_clerks[t]
        else:
            cmr = lp.LpVariable.dicts('Number of short of clerks', time_slots, 0)
            obj += lp.lpSum(cmr[t] for t in time_slots) * (10**data[Clerk_m_r][2])
            for t in time_slots:
                prob += lp.lpSum(x_c[(n, t)] for n in clerks) + cmr[t] >= mini_clerks[t]
    # constraint for standard work time and overtime
    if data[Clerk_m_t][0]:
        clerks_mt = data[Clerk_m_t][3][Clerk_i_d]
        mt = data[Clerk_m_t][3]['minimum work time']
        if data[Clerk_m_t][1]:
            for n in clerks_mt:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) >= mt
        else:
            cmt = lp.LpVariable.dicts('Number of less time', clerks_mt, 0)
            obj += lp.lpSum(cmt[n] for n in clerks_mt) * (10**data[Clerk_m_t][2])
            for n in clerks_mt:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) + cmt[n] >= mt
    if data[Clerk_o_t][0]:
        clerk_ot = data[Clerk_o_t][3][Clerk_i_d]
        ot = data[Clerk_o_t][3]['maximum work time']
        if data[Clerk_o_t][1]:
            for n in clerk_ot:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) <= ot
        else:
            cot = lp.LpVariable.dicts('Number of overtime', clerk_ot, 0)
            obj += lp.lpSum(cot[n] for n in clerk_ot) * (10**data[Clerk_o_t][2])
            for n in clerk_ot:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) <= ot + cot[n]
    # constraint for most continuous off days
    if data[Clerk_c_o][0]:
        clerk_co = data[Clerk_c_o][3][Clerk_i_d]
        cod = data[Clerk_c_o][3]['maximum off days']
        cod_days = range(1, day_num-cod+1)
        cod_clerk_days = [(n, d) for n in clerk_co for d in cod_days]
        if data[Clerk_c_o][1]:
            for n, d in cod_clerk_days:
                prob += lp.lpSum(z_c[(n, d+i-1)] for i in range(1, cod+2)) >= 1
        else:
            codd = lp.LpVariable.dicts('Number of cod', cod_clerk_days, 0)
            obj += lp.lpSum(codd[n_d] for n_d in cod_clerk_days) * (10**data[Clerk_c_o][2])
            for n, d in cod_clerk_days:
                prob += lp.lpSum(z_c[(n, d + i - 1)] for i in range(1, cod + 2)) + codd[(n, d)] >= 1
    # constraint for work day balance in different weeks
    if data[Clerk_w_b][0]:
        clerk_wb = data[Clerk_w_b][3][Clerk_i_d]
        wb_s_d = data[Clerk_w_b][3]['Week start days']
        wb_dif_b = data[Clerk_w_b][3]['Weekly work day differences']
        wb_clerk_d = [(n, d1, d2) for n in clerk_wb for d1 in wb_s_d for d2 in wb_s_d if d1 != d2]
        if data[Clerk_w_b][1]:
            for n, d1, d2 in wb_clerk_d:
                prob += (lp.lpSum(z_c[(n, d1+i)] for i in range(7)) -
                         lp.lpSum(z_c[(n, d2+i)] for i in range(7)) <= wb_dif_b)
        else:
            cwb = lp.LpVariable.dicts('Number of work day balance', wb_clerk_d, 0)
            obj += lp.lpSum(cwb[ndd] for ndd in wb_clerk_d) * (10**data[Clerk_w_b][2])
            for n, d1, d2 in wb_clerk_d:
                prob += (lp.lpSum(z_c[(n, d1 + i)] for i in range(7)) -
                         lp.lpSum(z_c[(n, d2 + i)] for i in range(7)) <= wb_dif_b + cwb[(n, d1, d2)])
    # early late switch
    if data[Clerk_e_l_s][0]:
        clerk_els = data[Clerk_e_l_s][3][Clerk_i_d]
        els_exception = data[Clerk_e_l_s][3]['Switch permission']
        els_clerk_day = [(n, d) for n in clerk_els for d in range(1, day_num) if d not in els_exception]
        if data[Clerk_e_l_s][1]:
            for n in clerk_els:
                for d in range(1, day_num):
                    if d in els_exception:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 2
                    else:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 1
        else:
            cels = lp.LpVariable.dicts('Number of unfavorable switch', els_clerk_day, 0)
            obj += lp.lpSum(cels[n_d] for n_d in els_clerk_day) * (10**data[Clerk_e_l_s][2])
            for n in clerk_els:
                for d in range(1, day_num):
                    if d in els_exception:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 2
                    else:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 1 + cels[(n, d)]
    # add object
    prob += obj
    prob.solve()
    return None
