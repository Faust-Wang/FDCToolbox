from .utils import fix, sign


# Ref [1] page 637
def interp_1(data_1d, a):

    s = .2 * a
    k = fix(s)

    if k <= -2:
        k = -1

    if k >= 9:
        k = 8

    da = s - k
    l = k + fix(1.1 * sign(da))
    k = k + 3
    l = l + 3

    d = data_1d[k-1] + abs(da) * (data_1d[l-1] - data_1d[k-1])

    return d


# Ref [1] page 638
def interp_2(data_2d, a, b):

    s = .2 * a
    k = fix(s)
    if k <= -2:
        k = -1

    if k >= 9:
        k = 8

    da = s - k
    l = k + fix(1.1 * sign(da))
    s = b / 12
    m = fix(s)
    if m <= -2:
        m = -1

    if m >= 2:
        m = 1

    de = s - m
    n = m + fix(1.1 * sign(de))
    k = k + 3
    l = l + 3
    m = m + 3
    n = n + 3
    t = data_2d[k-1, m-1]
    u = data_2d[k-1, n-1]
    v = t + abs(da) * (data_2d[l-1, m-1] - t)
    w = u + abs(da) * (data_2d[l-1, n-1] - u)
    d = v + (w - v) * abs(de)

    return d


# Ref [1] page 639
def interp_3(data_2d, a, b):

    s = .2 * a
    k = fix(s)

    if k <= -2:
        k = -1

    if k >= 9:
        k = 8

    da = s - k
    l = k + fix(1.1 * sign(da))
    s = .2 * abs(b)
    m = fix(s)
    if m == 0:
        m = 1

    if m >= 6:
        m = 5

    db = s - m
    n = m + fix(1.1 * sign(db))
    l = l + 3
    k = k + 3
    m = m + 1
    n = n + 1
    t = data_2d[k-1, m-1]
    u = data_2d[k-1, n-1]
    v = t + abs(da) * (data_2d[l-1, m-1] - t)
    w = u + abs(da) * (data_2d[l-1, n-1] - u)
    d = v + (w - v) * abs(db)

    return d * sign(b)


# Ref [1] page 640
def interp_4(data_2d, a, b):

    s = .2 * a
    k = fix(s)
    if k <= -2:
        k = -1

    if k >= 9:
        k = 8

    da = s - k
    l = k + fix(1.1 * sign(da))
    s = .1 * b
    m = fix(s)
    if m <= -3:
        m = -2

    if m >= 3:
        m = 2

    db = s - m
    n = m + fix(1.1 * sign(db))
    l = l + 3
    k = k + 3
    m = m + 4
    n = n + 4
    t = data_2d[k-1, m-1]
    u = data_2d[k-1, n-1]
    v = t + abs(da) * (data_2d[l-1, m-1] - t)
    w = u + abs(da) * (data_2d[l-1, n-1] - u)

    d =  v + (w - v) * abs(db)

    return d

