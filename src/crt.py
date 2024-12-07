"""."""

import functools
import math


# https://www3.cuc.ac.jp/~miyata/mathprg/nt/xeuclid.html
def xgcd(a, b):
    r0, r1 = a, b
    x0, x1 = 1, 0
    y0, y1 = 0, 1

    while r1 != 0:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1

    return (x0, y0)

def pack_ms(ms):
    k=0
    while 2**k < len(ms):
        k += 1
    r = 2**k
    return ms + ([1] * (r-len(ms))), r, k

def alg10_3(ms):
    r = len(ms)
    k = int(math.log2(r))
    M = [[1]*r for _ in range(k+1)]

    for i in range(r):
        M[0][i] = ms[i]

    for i in range(1, k+1):
        for j in range(0, 2**(k-i)):
            M[i][j] = M[i-1][2*j] * M[i-1][2*j+1]

    return M

def alg10_14(ms, f):
    if len(ms) == 1:
        return [f]
    else:
        M = alg10_3(ms)
        f0 = f % M[-2][0]
        f1 = f % M[-2][1]
    return alg10_14(ms[:len(ms)//2], f0) + alg10_14(ms[len(ms)//2:], f1)

def alg10_14x(i_ms_s, i_ms_e, f, M, i_f0, j_f0):
    if i_ms_e - i_ms_s == 1:
        return [f]
    else:
        f0 = f % M[i_f0][j_f0]
        f1 = f % M[i_f0][j_f0+1]
        return alg10_14x(i_ms_s, i_ms_s+(i_ms_e-i_ms_s)//2, f0, M, i_f0-1, 2*j_f0) \
             + alg10_14x(i_ms_s+(i_ms_e-i_ms_s)//2, i_ms_e, f1, M, i_f0-1, 2*j_f0+2)

def alg10_18(ms):
    m = functools.reduce(lambda a, b: a*b, ms, 1)
    ms2 = [mi**2 for mi in ms]

    # 1. call alg10_16 (alg10_16はalg10_3とalg10_14を併せたもの)
    M = alg10_3(ms2)
    m_rem_mi2 = alg10_14x(0, len(ms2), m, M, -2, 0)

    # 2. m/mi rem mi
    m_rem_mi = [m_rem2//mi for mi, m_rem2 in zip(ms, m_rem_mi2)]

    # 3. extended euclidean algo
    si = []
    for mi, mi_rem in zip(ms, m_rem_mi):
        a, b = xgcd(mi, mi_rem)
        if 0 < b:
            si.append(b)
        else:
            si.append(b+mi)
    return si

def alg10_20(cs, i_s, i_e, M, i_r0, j_r0):
    if i_e - i_s == 1:
        return cs[i_s]
    else:
        return M[i_r0][j_r0+1] * alg10_20(cs, i_s, i_s+(i_e-i_s)//2, M, i_r0-1, 2*j_r0) \
             + M[i_r0][j_r0]  *  alg10_20(cs, i_s+(i_e-i_s)//2, i_e, M, i_r0-1, 2*j_r0+2)

def alg10_22(ms, vs):
    assert(len(ms) == len(vs))

    M = alg10_3(ms)
    ss = alg10_18(ms)
    vs_ss = [(vi*si) % mi for vi, si, mi in zip(vs, ss, ms)]
    return alg10_20(vs_ss, 0, len(ms), M, -2, 0)

def crt(ms, vs):
    return alg10_22(pack_ms(ms)[0], pack_ms(vs)[0])
