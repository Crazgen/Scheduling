# encoding=utf-8
import pulp as lp
import xlwt
from datetime import datetime, timedelta, time
Clerk_i_d = 'Clerks included'

Clerk_m_s_h = 'Clerk minimum shifting hour'
Clerk_m_s = 'Clerk multiple shift'
Clerk_w_t_b = 'Clerk work time balance'
Clerk_w_d_b = 'Clerk work day balance'
Clerk_m_o_d = 'Clerk minimum off days'
Clerk_m_w_d = 'Clerk most whole day'
Clerk_m_r = 'minimum clerks'
Clerk_m_t = 'Clerk minimum work time'
Clerk_o_t = 'Clerk most over time'
Clerk_c_o = 'Clerk continuous off days'
Clerk_w_b = 'Clerk weekly workday balance'
Clerk_e_l_s = 'Clerk early late switch'
CONSTRAINTS = [Clerk_m_s_h, Clerk_m_s, Clerk_w_t_b, Clerk_w_d_b, Clerk_m_o_d, Clerk_m_w_d, Clerk_m_r, Clerk_m_t,
               Clerk_o_t, Clerk_c_o, Clerk_w_b, Clerk_e_l_s]


# constraints in data: data['constraint name'] = (is_involved, must_meet, priority, optional_data_dicts)
def solve_single_floor(data):
    hour_num, day_num, clerk_num = data['working hours'], data['working days'], data['clerks']
    time_slots_num = hour_num * day_num
    clerks, time_slots, hours, days = (range(1, clerk_num+1), range(1, time_slots_num+1), range(1, hour_num+1),
                                       range(1, day_num+1))
    clerk_time = [(n, t) for n in clerks for t in time_slots]
    clerk_day = [(n, d) for n in clerks for d in days]

    prob = lp.LpProblem("Scheduling", lp.LpMinimize)
    # This part for clerk constraints
    x_c = lp.LpVariable.dicts('Time on duty', clerk_time, cat='Binary')
    z_c = lp.LpVariable.dicts('Day on duty', clerk_day, cat='Binary')
    kesi_c = lp.LpVariable.dicts('Switch on duty', clerk_time, cat='Binary')
    cons_relax = {c_str: None for c_str in CONSTRAINTS}
    obj = lp.LpAffineExpression()
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
            cons_relax[Clerk_m_s_h] = lp.LpVariable.dicts('Number of exceeding hours', clerk_day, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_s_h][n_d] for n_d in clerk_day) * (10 ** data[Clerk_m_s_h][2])
            for n, d in clerk_day:
                prob += (data[Clerk_m_s_h][3]['minimum hour'] * z_c[(n, d)] <=
                         lp.lpSum(x_c[(n, hour_num * (d - 1) + h)] for h in hours) + cons_relax[Clerk_m_s_h][(n, d)])
    # constraint for no cross shifty
    if data[Clerk_m_s][0]:
        if data[Clerk_m_s][1]:
            for n, d in clerk_day:
                prob += lp.lpSum(kesi_c[(n, hour_num*(d-1)+h)] for h in hours) <= 1
        else:
            cons_relax[Clerk_m_s] = lp.LpVariable.dicts('Number of cross shifting', clerk_day, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_s][n_d] for n_d in clerk_day) * (10 ** data[Clerk_m_s][2])
            for n, d in clerk_day:
                prob += (lp.lpSum(kesi_c[(n, hour_num * (d - 1) + h)] for h in hours) <=
                         1 + cons_relax[Clerk_m_s][(n, d)])
    # constraint for work time balance
    if data[Clerk_w_t_b][0]:
        clerk_wtb = data[Clerk_w_t_b][3][Clerk_i_d]
        dif_wtb = data[Clerk_w_t_b][3]['Difference bound']
        clerk_inter_wtb = [(i, j) for i in clerk_wtb for j in clerk_wtb if i != j]
        if data[Clerk_w_t_b][1]:
            for i, j in clerk_inter_wtb:
                prob += lp.lpSum(x_c[(i, t)] for t in time_slots) - lp.lpSum(x_c[(j, t)] for t in time_slots) <= dif_wtb
        else:
            cons_relax[Clerk_w_t_b] = lp.LpVariable.dicts('Number of exceeding work hours', clerk_inter_wtb, 0)
            obj += lp.lpSum(cons_relax[Clerk_w_t_b][i_j] for i_j in clerk_inter_wtb) * (10 ** data[Clerk_w_t_b][2])
            for i, j in clerk_inter_wtb:
                prob += (lp.lpSum(x_c[(i, t)] for t in time_slots) - lp.lpSum(x_c[(j, t)] for t in time_slots) <=
                         dif_wtb + cons_relax[Clerk_w_t_b][(i, j)])
    # constraint for work day balance
    if data[Clerk_w_d_b][0]:
        clerk_wdb = data[Clerk_w_d_b][3][Clerk_i_d]
        dif_wdb = data[Clerk_w_d_b][3]['Difference bound']
        clerk_inter_wdb = [(i, j) for i in clerk_wdb for j in clerk_wdb if i != j]
        if data[Clerk_w_d_b][1]:
            for i, j in clerk_inter_wdb:
                prob += lp.lpSum(z_c[(i, d)] for d in days) - lp.lpSum(z_c[(j, d)] for d in days) <= dif_wdb
        else:
            cons_relax[Clerk_w_d_b] = lp.LpVariable.dicts('Number of exceeding work days', clerk_inter_wdb, 0)
            obj += lp.lpSum(cons_relax[Clerk_w_d_b][i_j] for i_j in clerk_inter_wdb) * (10 ** data[Clerk_w_d_b][2])
            for i, j in clerk_inter_wdb:
                prob += (lp.lpSum(z_c[(i, d)] for d in days) - lp.lpSum(z_c[(j, d)] for d in days) <=
                         dif_wdb + cons_relax[Clerk_w_d_b][(i, j)])
    # constraint for minimum off days in the period
    if data[Clerk_m_o_d][0]:
        clerk_mod = data[Clerk_m_o_d][3][Clerk_i_d]
        mod = data[Clerk_m_o_d][3]['Min off day']
        if data[Clerk_m_o_d][1]:
            for n in clerk_mod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) <= day_num - mod
        else:
            cons_relax[Clerk_m_o_d] = lp.LpVariable.dicts('Number of exceeding days', clerk_mod, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_o_d][n] for n in clerk_mod) * (10 ** data[Clerk_m_o_d][2])
            for n in clerk_mod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) <= day_num - mod + cons_relax[Clerk_m_o_d][n]
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
            cons_relax[Clerk_m_w_d] = lp.LpVariable.dicts('Number of exceeding whole days', clerk_time_mwd, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_w_d][n_t] for n_t in clerk_time_mwd) * (10 ** data[Clerk_m_w_d][2])
            for n, t in clerk_time_mwd:
                prob += (lp.lpSum(x_c[(n, t + k - 1)] for k in range(1, mwd + 1)) <=
                         mwd - 1 + cons_relax[Clerk_m_w_d][(n, t)])
    # constraint for minimum clerks required, combined with traffic..
    if data[Clerk_m_r][0]:
        mini_clerks = data[Clerk_m_r][3]['minimum clerk num']
        if data[Clerk_m_r][1]:
            for t in time_slots:
                prob += lp.lpSum(x_c[(n, t)] for n in clerks) >= mini_clerks[t]
        else:
            cons_relax[Clerk_m_r] = lp.LpVariable.dicts('Number of short of clerks', time_slots, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_r][t] for t in time_slots) * (10 ** data[Clerk_m_r][2])
            for t in time_slots:
                prob += lp.lpSum(x_c[(n, t)] for n in clerks) + cons_relax[Clerk_m_r][t] >= mini_clerks[t]
    # constraint for standard work time and overtime
    if data[Clerk_m_t][0]:
        clerks_mt = data[Clerk_m_t][3][Clerk_i_d]
        mt = data[Clerk_m_t][3]['minimum work time']
        if data[Clerk_m_t][1]:
            for n in clerks_mt:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) >= mt
        else:
            cons_relax[Clerk_m_t] = lp.LpVariable.dicts('Number of less time', clerks_mt, 0)
            obj += lp.lpSum(cons_relax[Clerk_m_t][n] for n in clerks_mt) * (10 ** data[Clerk_m_t][2])
            for n in clerks_mt:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) + cons_relax[Clerk_m_t][n] >= mt
    if data[Clerk_o_t][0]:
        clerk_ot = data[Clerk_o_t][3][Clerk_i_d]
        ot = data[Clerk_o_t][3]['maximum work time']
        if data[Clerk_o_t][1]:
            for n in clerk_ot:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) <= ot
        else:
            cons_relax[Clerk_o_t] = lp.LpVariable.dicts('Number of overtime', clerk_ot, 0)
            obj += lp.lpSum(cons_relax[Clerk_o_t][n] for n in clerk_ot) * (10 ** data[Clerk_o_t][2])
            for n in clerk_ot:
                prob += lp.lpSum(x_c[(n, t)] for t in time_slots) <= ot + cons_relax[Clerk_o_t][n]
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
            cons_relax[Clerk_c_o] = lp.LpVariable.dicts('Number of cod', cod_clerk_days, 0)
            obj += lp.lpSum(cons_relax[Clerk_c_o][n_d] for n_d in cod_clerk_days) * (10 ** data[Clerk_c_o][2])
            for n, d in cod_clerk_days:
                prob += lp.lpSum(z_c[(n, d + i - 1)] for i in range(1, cod + 2)) + cons_relax[Clerk_c_o][(n, d)] >= 1
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
            cons_relax[Clerk_w_b] = lp.LpVariable.dicts('Number of work day balance', wb_clerk_d, 0)
            obj += lp.lpSum(cons_relax[Clerk_w_b][ndd] for ndd in wb_clerk_d) * (10 ** data[Clerk_w_b][2])
            for n, d1, d2 in wb_clerk_d:
                prob += (lp.lpSum(z_c[(n, d1 + i)] for i in range(7)) -
                         lp.lpSum(z_c[(n, d2 + i)] for i in range(7)) <= wb_dif_b + cons_relax[Clerk_w_b][(n, d1, d2)])
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
            cons_relax[Clerk_e_l_s] = lp.LpVariable.dicts('Number of unfavorable switch', els_clerk_day, 0)
            obj += lp.lpSum(cons_relax[Clerk_e_l_s][n_d] for n_d in els_clerk_day) * (10 ** data[Clerk_e_l_s][2])
            for n in clerk_els:
                for d in range(1, day_num):
                    if d in els_exception:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 2
                    else:
                        prob += (x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <=
                                 1 + cons_relax[Clerk_e_l_s][(n, d)])
    # add object
    prob += obj
    prob.solve(lp.PULP_CBC_CMD(msg=1, maxSeconds=data['Time limit']))
    print(lp.LpStatus[prob.status])
    if 'Not Solved' != lp.LpStatus[prob.status] != 'Infeasible':
        obj_res = lp.value(prob.objective)
        return obj_res, cons_relax, x_c, z_c, kesi_c, lp.LpStatus[prob.status]
    else:
        return None, cons_relax, None, None, None, lp.LpStatus[prob.status]


def output_excel(filename, data, result):
    wb = xlwt.Workbook()
    analysis_s = wb.add_sheet(u'结果分析')
    if not (result['x_c'] is None):
        # write the scheduling result
        schedule_s = wb.add_sheet(u'排班')
        on_st = xlwt.easyxf('pattern: pattern solid, pattern_fore_colour 17;' +
                            ' borders: left THIN, right THIN, top THIN, bottom THIN')
        off_st = xlwt.easyxf('pattern: pattern solid, pattern_fore_colour 73;' +
                             ' borders: left THIN, right THIN, top THIN, bottom THIN')
        date_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN;' +
                              ' align: wrap on, vert center, horiz center',
                              num_format_str=u'YYYY年M月D日 AAAA')
        hour_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN',
                              num_format_str='h:mm')
        reg_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN')
        hour_num, day_num, clerk_num = data['working hours'], data['working days'], data['clerks']
        day_start_row, date_delta, hour_delta = 0, timedelta(days=1), timedelta(hours=1)
        start_day, start_hour = data['Start date'], data['Start hour']
        clerk_workdays, clerk_work_time = [0*n for n in range(clerk_num+1)], [0*n for n in range(clerk_num+1)]
        for d in range(1, day_num+1):
            schedule_s.write_merge(day_start_row, day_start_row, 1, hour_num, start_day + timedelta(days=d-1), date_st)
            day_start_row += 1
            for h in range(1, hour_num+1):
                schedule_s.write(day_start_row, h,
                                 (datetime.combine(datetime.today(), start_hour) + timedelta(hours=h-1)).time(),
                                 hour_st)
            day_start_row += 1
            for n in range(1, clerk_num+1):
                clerk_workdays[n] += result['z_c'][(n, d)].varValue
                schedule_s.write(day_start_row, 0, data['Clerk names'][n], reg_st)
                for h in range(1, hour_num+1):
                    if result['x_c'][(n, hour_num*(d-1)+h)].varValue == 0:
                        schedule_s.write(day_start_row, h, 0, off_st)
                    else:
                        clerk_work_time[n] += 1
                        schedule_s.write(day_start_row, h, 1, on_st)
                day_start_row += 1
            day_start_row += 1
        # write the result analysis
        print(clerk_workdays)
        print(clerk_work_time)
    else:
        pass
    wb.save(filename)


def test(timelim):
    data = dict()
    data['working hours'], data['working days'], data['clerks'] = 12, 31, 4
    data['Time limit'] = timelim
    data['Start date'] = datetime(2016, 3, 14)
    data['Start hour'] = time(10)
    data['Clerk names'] = [None, u'工', u'了', u'要', u'在']
    clerks = range(1, 5)
    data[Clerk_m_s_h] = (True, False, 1, {'minimum hour': 6})
    data[Clerk_m_s] = (True, False, 1, None)
    data[Clerk_w_t_b] = (True, False, 1, {Clerk_i_d: clerks, 'Difference bound': 3})
    data[Clerk_w_d_b] = (True, False, 1, {Clerk_i_d: clerks, 'Difference bound': 2})
    data[Clerk_m_o_d] = (True, False, 1, {Clerk_i_d: clerks, 'Min off day': 31*2.0/7})
    data[Clerk_m_w_d] = (True, False, 1, {Clerk_i_d: clerks, 'Most whole day': 3})
    data[Clerk_m_r] = (True, False, 1, {'minimum clerk num': [None, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2,
                                                              2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
                                                              1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2,
                                                              2, 2,
                                                              2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2,
                                                              2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
                                                              1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2,
                                                              2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2,
                                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
                                                              3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2,
                                                              2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                                                              1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2,
                                                              2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                              2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
                                                              2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1,
                                                              2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
                                                              2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2,
                                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                                                              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2,
                                                              2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1,
                                                              1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1]})
    data[Clerk_m_t] = (True, False, 1, {Clerk_i_d: clerks, 'minimum work time': 167})
    data[Clerk_o_t] = (True, False, 0, {Clerk_i_d: clerks, 'maximum work time': 167+9})
    data[Clerk_c_o] = (True, False, 1, {Clerk_i_d: clerks, 'maximum off days': 3})
    data[Clerk_w_b] = (True, False, 1, {Clerk_i_d: clerks, 'Weekly work day differences': 1, 'Week start days': [1, 8,
                                                                                                                 15,
                                                                                                                 22]})
    tttjust = [None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    p_days = set()
    for d in range(1, 32):
        for h in range(1, 13):
            t = (d-1) * 12 + h
            if tttjust[t] > 0:
                p_days.add(d)
    data[Clerk_e_l_s] = (True, False, 0, {Clerk_i_d: clerks, 'Switch permission': p_days})
    result = dict()
    (result['objective'], result['relax'], result['x_c'],
     result['z_c'], result['kesi_c'], result['status']) = solve_single_floor(data)
    output_excel('test_schedule.xls', data, result)
    return result
