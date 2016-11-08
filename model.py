# encoding=utf-8
import pulp as lp
import xlwt
import xlrd
import os
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
Clerk_ma_o_d = 'Clerk maximum off day'
CONSTRAINTS = [Clerk_m_s_h, Clerk_m_s, Clerk_w_t_b, Clerk_w_d_b, Clerk_m_o_d, Clerk_m_w_d, Clerk_m_r, Clerk_m_t,
               Clerk_o_t, Clerk_c_o, Clerk_w_b, Clerk_e_l_s, Clerk_ma_o_d]
CONSTRAINTS_NAME = [u'每班次最少工时', u'无串班', u'员工间工时平衡', u'员工间工作天数平衡', u'排班周期最少休息天数',
                    u'最多连续全班数', u'各时段最少员工数', u'排班周期内员工最少工时', u'排班周期内员工最多工时',
                    u'最多连续休息天数', u'同一员工每周工作天数平衡', u'消除晚、早班连上', u'最多休息天数']
CONSTRAINTS_IND = {cons_sss: '' for cons_sss in CONSTRAINTS}


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
            CONSTRAINTS_IND[Clerk_m_s_h] = u'(员工, 天数): '
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
            CONSTRAINTS_IND[Clerk_m_s] = u'(员工, 天数): '
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
            CONSTRAINTS_IND[Clerk_w_t_b] = u'(员工, 员工): '
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
            CONSTRAINTS_IND[Clerk_w_d_b] = u'(员工, 员工): '
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
            CONSTRAINTS_IND[Clerk_m_o_d] = u'员工: '
            obj += lp.lpSum(cons_relax[Clerk_m_o_d][n] for n in clerk_mod) * (10 ** data[Clerk_m_o_d][2])
            for n in clerk_mod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) <= day_num - mod + cons_relax[Clerk_m_o_d][n]
    # constraint for the maximum off days
    if data[Clerk_ma_o_d][0]:
        clerk_maod = data[Clerk_ma_o_d][3][Clerk_i_d]
        maod = data[Clerk_ma_o_d][3]['Max off day']
        if data[Clerk_ma_o_d][1]:
            for n in clerk_maod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) >= day_num - maod
        else:
            cons_relax[Clerk_ma_o_d] = lp.LpVariable.dicts('Exceeding off days', clerk_maod, 0)
            CONSTRAINTS_IND[Clerk_ma_o_d] = u'员工：'
            obj += lp.lpSum(cons_relax[Clerk_ma_o_d][n] for n in clerk_maod) * (10 ** data[Clerk_ma_o_d][2])
            for n in clerk_maod:
                prob += lp.lpSum(z_c[(n, d)] for d in days) >= day_num - maod - cons_relax[Clerk_ma_o_d][n]
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
            CONSTRAINTS_IND[Clerk_m_w_d] = u'(员工, 时间): '
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
            CONSTRAINTS_IND[Clerk_m_r] = u'时间: '
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
            CONSTRAINTS_IND[Clerk_m_t] = u'员工: '
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
            CONSTRAINTS_IND[Clerk_o_t] = u'员工: '
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
            CONSTRAINTS_IND[Clerk_c_o] = u'(员工, 天数): '
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
            CONSTRAINTS_IND[Clerk_w_b] = u'(员工, 天数, 天数): '
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
                    if d not in els_exception:
                        prob += x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <= 1
        else:
            cons_relax[Clerk_e_l_s] = lp.LpVariable.dicts('Number of unfavorable switch', els_clerk_day, 0)
            CONSTRAINTS_IND[Clerk_e_l_s] = u'(员工, 天数): '
            obj += lp.lpSum(cons_relax[Clerk_e_l_s][n_d] for n_d in els_clerk_day) * (10 ** data[Clerk_e_l_s][2])
            for n in clerk_els:
                for d in range(1, day_num):
                    if d not in els_exception:
                        prob += (x_c[(n, hour_num * d)] + x_c[(n, hour_num * d + 1)] <=
                                 1 + cons_relax[Clerk_e_l_s][(n, d)])
    # add object
    prob += obj
    # prob.solve(lp.PULP_CBC_CMD(msg=1, maxSeconds=data['Time limit']))
    prob.solve(lp.COIN_CMD(msg=1, maxSeconds=data['Time limit'], path=os.getcwdu() + '\\cbc.exe'))
    print(lp.LpStatus[prob.status])
    if lp.LpStatus[prob.status] == 'Infeasible':
        return None, cons_relax, None, None, None, lp.LpStatus[prob.status], prob.constraints
    else:
        obj_res = lp.value(prob.objective)
        return obj_res, cons_relax, x_c, z_c, kesi_c, lp.LpStatus[prob.status], prob.constraints


def output_excel(filename, data, result):
    is_feasible = check_feasible(result)
    wb = xlwt.Workbook()
    analysis_s = wb.add_sheet(u'结果分析')
    analysis_s.write(0, 0, u'结果状态：')
    off_st = xlwt.easyxf('pattern: pattern solid, pattern_fore_colour 73;' +
                         ' borders: left THIN, right THIN, top THIN, bottom THIN')
    if result['status'] != 'Optimal':
        analysis_s.write(0, 1, u'未找到最优状态，可尝试增加求解时间得到更好的排班。')
    if (not (result['x_c'] is None)) and is_feasible:
        # write the scheduling result
        schedule_s = wb.add_sheet(u'排班')
        on_st = xlwt.easyxf('pattern: pattern solid, pattern_fore_colour 17;' +
                            ' borders: left THIN, right THIN, top THIN, bottom THIN')
        date_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN;' +
                              ' align: wrap on, vert center, horiz center',
                              num_format_str=u'YYYY年M月D日 AAAA')
        hour_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN',
                              num_format_str='h:mm')
        reg_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN')
        reg_m_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN;' +
                               ' align: wrap on, vert center, horiz center')
        hour_num, day_num, clerk_num = data['working hours'], data['working days'], data['clerks']
        day_start_row, date_delta, hour_delta = 0, timedelta(days=1), timedelta(hours=1)
        start_day, start_hour = data['Start date'], data['Start hour']
        clerk_workdays, clerk_work_time = [0*n for n in range(clerk_num+1)], [0*n for n in range(clerk_num+1)]
        clerk_num_at_shop = {(d, h): 0 for d in range(1, day_num+1) for h in range(1, hour_num+1)}
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
                        clerk_num_at_shop[(d, h)] += 1
                        schedule_s.write(day_start_row, h, 1, on_st)
                day_start_row += 1
            day_start_row += 1
        # calculate the constraints
        constraint_re_num = {con_str: [con_name, 0, CONSTRAINTS_IND[con_str]] for con_str, con_name
                             in zip(CONSTRAINTS, CONSTRAINTS_NAME)}
        for con_str in CONSTRAINTS:
            cons_var = result['relax'][con_str]
            if not (cons_var is None):
                for the_cons_k in cons_var.keys():
                    the_cons = cons_var[the_cons_k].varValue
                    if abs(the_cons) > 0.01:
                        constraint_re_num[con_str][2] += str(the_cons_k) + ', '
                    constraint_re_num[con_str][1] += the_cons
        unmet_cons = len(CONSTRAINTS) - [constraint_re_num[con_str][1] for con_str in CONSTRAINTS].count(0)
        # write the header
        analysis_s.write(1, 0, u'得到排班方案，违反约束类别数为' + str(unmet_cons))
        work_time_rc, cons_rc, clerk_num_rc = (3, 0), (3, 5), (6 + clerk_num, 0)
        analysis_s.write_merge(work_time_rc[0], work_time_rc[0], work_time_rc[1], work_time_rc[1] + 3, u'员工工时总结',
                               reg_m_st)
        analysis_s.write(work_time_rc[0] + 1, work_time_rc[1] + 1, u'总工时', reg_st)
        analysis_s.write(work_time_rc[0] + 1, work_time_rc[1] + 2, u'工作天数', reg_st)
        analysis_s.write(work_time_rc[0] + 1, work_time_rc[1] + 3, u'休息天数', reg_st)
        analysis_s.write_merge(cons_rc[0], cons_rc[0], cons_rc[1], cons_rc[1] + 5, u'排班未满足约束', reg_m_st)
        analysis_s.write_merge(cons_rc[0] + 1, cons_rc[0] + 1, cons_rc[1], cons_rc[1] + 2, u'未满足约束', reg_m_st)
        analysis_s.write(cons_rc[0] + 1, cons_rc[1] + 3, u'优先级', reg_st)
        analysis_s.write(cons_rc[0] + 1, cons_rc[1] + 4, u'违反量', reg_st)
        analysis_s.write(cons_rc[0] + 1, cons_rc[1] + 5, u'详情', reg_st)
        analysis_s.write_merge(clerk_num_rc[0], clerk_num_rc[0], clerk_num_rc[1], clerk_num_rc[1] + 3,
                               u'不同时段上班员工数', reg_m_st)
        analysis_s.write(clerk_num_rc[0] + 1, clerk_num_rc[1] + 0, u'日期', reg_st)
        analysis_s.write(clerk_num_rc[0] + 1, clerk_num_rc[1] + 1, u'星期', reg_st)
        analysis_s.write(clerk_num_rc[0] + 1, clerk_num_rc[1] + 2, u'时间', reg_st)
        analysis_s.write(clerk_num_rc[0] + 1, clerk_num_rc[1] + 3, u'员工数', reg_st)
        # write work time summary
        for n in range(1, clerk_num+1):
            analysis_s.write(work_time_rc[0] + 1 + n, work_time_rc[1], data['Clerk names'][n], reg_st)
            analysis_s.write(work_time_rc[0] + 1 + n, work_time_rc[1] + 1, clerk_work_time[n], reg_st)
            analysis_s.write(work_time_rc[0] + 1 + n, work_time_rc[1] + 2, clerk_workdays[n], reg_st)
            analysis_s.write(work_time_rc[0] + 1 + n, work_time_rc[1] + 3, day_num - clerk_workdays[n], reg_st)
        # write the constraint summary
        cons_row = cons_rc[0] + 2
        for con_str in CONSTRAINTS:
            if constraint_re_num[con_str][1] > 0:
                analysis_s.write_merge(cons_row, cons_row, cons_rc[1], cons_rc[1] + 2, constraint_re_num[con_str][0],
                                       reg_m_st)
                analysis_s.write(cons_row, cons_rc[1] + 3, data[con_str][2], reg_st)
                analysis_s.write(cons_row, cons_rc[1] + 4, constraint_re_num[con_str][1], reg_st)
                analysis_s.write(cons_row, cons_rc[1] + 5, constraint_re_num[con_str][2], reg_st)
                cons_row += 1
        # write clerk num summary
        c_n_row = clerk_num_rc[0] + 2
        date_st1 = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN;', num_format_str=u'M月D日')
        week_st = xlwt.easyxf('borders: left THIN, right THIN, top THIN, bottom THIN;', num_format_str=u'AAAA')
        for d in range(1, day_num+1):
            for h in range(1, hour_num+1):
                analysis_s.write(c_n_row, clerk_num_rc[1], start_day + timedelta(days=d-1), date_st1)
                analysis_s.write(c_n_row, clerk_num_rc[1] + 1, start_day + timedelta(days=d - 1), week_st)
                analysis_s.write(c_n_row, clerk_num_rc[1] + 2,
                                 (datetime.combine(datetime.today(), start_hour) + timedelta(hours=h-1)).time(),
                                 hour_st)
                analysis_s.write(c_n_row, clerk_num_rc[1] + 3, clerk_num_at_shop[(d, h)], reg_st)
                c_n_row += 1
    elif is_feasible:
        analysis_s.write(1, 0, u'未找到可行排班，请尝试调整以下某些约束的优先级：')
        start_row = 2
        temp_l = zip(CONSTRAINTS_NAME, CONSTRAINTS)
        temp_l.reverse()
        for cons_name, cons_str in temp_l:
            if data[cons_str][0] and data[cons_str][1]:
                analysis_s.write_merge(start_row, start_row, 0, 4, str(start_row-1) + '. ' + cons_name, off_st)
                start_row += 1
    else:
        analysis_s.write(1, 0, u'未找到可行排班，请增加求解时间。')
    try_str, try_num, ori_f_name = '', 0, filename
    while try_num < 100000:
        try:
            wb.save(filename + '.xls')
        except IOError, e:
            try_num += 1
            if e.errno == 13:
                filename = ori_f_name + str(try_num)
            else:
                break
        else:
            os.startfile(os.path.abspath(filename + '.xls'))
            break


def input_excel(filename):
    filename = os.path.abspath(filename)
    try:
        w_f = xlrd.open_workbook(filename)
    except IOError:
        print(u'无法找到参数文件。请确认文件在正确位置。')
        os.system('pause')
        return None
    else:
        basic = w_f.sheet_by_name(u'基本信息')
        data = dict()
        data['working days'], data['clerks'] = int(basic.cell_value(2, 1)), int(basic.cell_value(1, 1))
        data['working hours'], data['Time limit'] = int(basic.cell_value(5, 1)), int(basic.cell_value(6, 1))
        data['Start date'] = datetime(*xlrd.xldate_as_tuple(basic.cell_value(3, 1), w_f.datemode)[:])
        data['Start hour'] = time(*xlrd.xldate_as_tuple(basic.cell_value(4, 1), w_f.datemode)[3:])
        flags = w_f.sheet_by_name('FLAGS')
        flag_dict = {str(f_str): fl == 1 for f_str, fl in zip(flags.col_values(0), flags.col_values(1))}
        if flag_dict['Clerk names']:
            c_list_sheet = w_f.sheet_by_name(u'员工列表')
            data['Clerk names'] = [None]
            for n in range(1, data['clerks']+1):
                try:
                    data['Clerk names'].append(c_list_sheet.cell_value(n, 0))
                except IndexError:
                    data['Clerk names'].append(u'员工' + str(n))
        else:
            data['Clerk names'] = [None]
            for n in range(1, data['clerks'] + 1):
                data['Clerk names'].append(u'员工' + str(n))
        clerks = range(1, data['clerks'] + 1)
        cons_sheet = w_f.sheet_by_name(u'约束管理')
        priority_dict = {str(f_str): pri for f_str, pri in zip(flags.col_values(0)[1:13],
                                                               cons_sheet.col_values(1)[1:13])}
        priority_dict[str(flags.col_values(0)[14])] = cons_sheet.col_values(1)[14]
        tight_dict = dict()
        for con_str in CONSTRAINTS:
            if priority_dict[con_str] == u'必须满足':
                tight_dict[con_str] = True
            else:
                tight_dict[con_str] = False
        cons_option_dict = dict()
        cons_option_dict[Clerk_m_s_h] = {'minimum hour': cons_sheet.cell_value(1, 4)}
        cons_option_dict[Clerk_m_s] = None
        cons_option_dict[Clerk_w_t_b] = {Clerk_i_d: clerks, 'Difference bound': cons_sheet.cell_value(3, 4)}
        cons_option_dict[Clerk_w_d_b] = {Clerk_i_d: clerks, 'Difference bound': cons_sheet.cell_value(4, 4)}
        cons_option_dict[Clerk_m_o_d] = {Clerk_i_d: clerks, 'Min off day': cons_sheet.cell_value(5, 4)}
        cons_option_dict[Clerk_m_w_d] = {Clerk_i_d: clerks, 'Most whole day': int(cons_sheet.cell_value(6, 4)) + 1}
        cons_option_dict[Clerk_m_t] = {Clerk_i_d: clerks, 'minimum work time': cons_sheet.cell_value(8, 4)}
        cons_option_dict[Clerk_o_t] = {Clerk_i_d: clerks, 'maximum work time': cons_sheet.cell_value(9, 4)}
        cons_option_dict[Clerk_c_o] = {Clerk_i_d: clerks, 'maximum off days': int(cons_sheet.cell_value(10, 4))}
        cons_option_dict[Clerk_ma_o_d] = {Clerk_i_d: clerks, 'Max off day': cons_sheet.cell_value(14, 4)}
        week_s_d, try_sd = [], 1
        while try_sd + 6 <= data['working days']:
            week_s_d.append(try_sd)
            try_sd += 7
        cons_option_dict[Clerk_w_b] = {Clerk_i_d: clerks, 'Weekly work day differences': cons_sheet.cell_value(11, 4),
                                       'Week start days': week_s_d}
        els_sheet = w_f.sheet_by_name('els')
        p_days, d_st = set(), 1
        for exc_str in els_sheet.col_values(1)[2:]:
            if d_st > data['working days']:
                break
            if exc_str == u'是':
                p_days.add(d_st)
            d_st += 1
        cons_option_dict[Clerk_e_l_s] = {Clerk_i_d: clerks, 'Switch permission': p_days}
        minic_sheet = w_f.sheet_by_name('minic')
        mini_c_num = list()
        mini_c_num.append(None)
        for d in range(1, data['working days'] + 1):
            for h in range(1, data['working hours'] + 1):
                if minic_sheet.cell_type(d + 1, h) != 2:
                    mini_c_num.append(0)
                else:
                    mini_c_num.append(minic_sheet.cell_value(d + 1, h))
        if flag_dict['Consider traffic']:
            service_level = cons_sheet.cell_value(13, 4)
            traff_sheet = w_f.sheet_by_name('traffic')
            for d in range(1, data['working days'] + 1):
                for h in range(1, data['working hours'] + 1):
                    if traff_sheet.cell_type(d + 1, h) == 2:
                        t = (d-1)*data['working hours'] + h
                        trf = traff_sheet.cell_value(d + 1, h)*1.0/service_level
                        if mini_c_num[t] < trf:
                            mini_c_num[t] = trf
        cons_option_dict[Clerk_m_r] = {'minimum clerk num': mini_c_num}
        for con_str in CONSTRAINTS:
            data[con_str] = (flag_dict[con_str], tight_dict[con_str], priority_dict[con_str], cons_option_dict[con_str])
        return data


def check_feasible(result, eps=0.01):
    for c in result['constraints'].values():
        if not c.valid(eps):
            return False
    for x in result['x_c'].values():
        xv = x.varValue
        if abs(xv) > eps < abs(xv-1.0):
            return False
    for x in result['z_c'].values():
        xv = x.varValue
        if abs(xv) > eps < abs(xv - 1.0):
            return False
    for x in result['kesi_c'].values():
        xv = x.varValue
        if abs(xv) > eps < abs(xv - 1.0):
            return False
    return True


def test(timelim):
    data = dict()
    data['working hours'], data['working days'], data['clerks'] = 12, 31, 4
    data['Time limit'] = timelim
    data['Start date'] = datetime(2016, 3, 14)
    data['Start hour'] = time(10)
    data['Clerk names'] = [None, u'工', u'了', u'要', u'在']
    clerks = range(1, 5)
    data[Clerk_m_s_h] = (True, True, 1, {'minimum hour': 6})
    data[Clerk_m_s] = (True, True, 1, None)
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
    data[Clerk_m_t] = (True, True, 1, {Clerk_i_d: clerks, 'minimum work time': 167})
    data[Clerk_o_t] = (True, True, 0, {Clerk_i_d: clerks, 'maximum work time': 167+5})
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
    output_excel('../Schedule result/test_schedule', data, result)
    return result


if __name__ == '__main__':
    data = input_excel(u'../排班工具.xlsm')
    result = dict()
    (result['objective'], result['relax'], result['x_c'],
     result['z_c'], result['kesi_c'], result['status'], result['constraints']) = solve_single_floor(data)
    output_excel('../Schedule result/test_schedule', data, result)
