import pulp as lp


def solve_single_floor():
    x = lp.LpVariable('duty', range(10), cat='Binary')
    return x

