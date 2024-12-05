"""."""

import math
import sympy
import functools
from fractions import Fraction

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

def get_primitive_root_(p):
    g = 3
    for i in range(3, p):
        # 原始根が奇数である判定を追加
        if i%2 != 1:
            continue
        elements = set([(i**j) % p for j in range(1, p)])
        if (len(elements) == p-1):
            g = i
            break
    return g

def get_primitive_root(p):
    for i in sympy.primerange(3, p+1):
        x = i
        for j in range(2, p):
            x *= i
            x %= p
            if x == 1:
                break
        if j == p-1:
            return i

def calc_mod(p, k):
    # 原始根(生成元)を求める
    g = get_primitive_root(p)

    r = pow(g, k-1, p)
    u = ((g - 1) // 2) % p
    S = 0; X = 1; Y = r
    for i in range(1, (p//2)+1):
        q = math.floor(g*X/p)
        S = (S + (u - q)*Y) % p
        X = (g * X) % p
        Y = (r * Y) % p

    # 1 - g**k の逆元を求める
    inv_1_gk = pow(1 - pow(g, k, p), -1, p)

    return ((2 * k * S) * inv_1_gk) % p

def rational(k):
    Y = max(37, math.ceil(k+0.5*math.log2(k)))
    primes = list(sympy.primerange(Y+1))
    Dk = 1
    for p in primes:
        if k % (p-1) == 0:
            Dk *= p
    beta = math.ceil((k+0.5)*math.log2(k)-4.094*k+2.470+math.log2(Dk))
    p=3; M_prime = 1
    while M_prime < 2**(beta+1):
        p = sympy.nextprime(p)
        if (k % (p-1)) != 0:
           M_prime = p * M_prime
    X = p

    rp = {}
    for p in sympy.primerange(X+1):
        if k % (p-1) != 0:
            rp[p] = calc_mod(p, k)
    # print(rp)
    R = alg10_22(pack_ms(list(rp.keys()))[0], pack_ms(list(rp.values()))[0])
    M = 1
    for p in sympy.primerange(X+1):
        if k % (p-1) != 0:
            M *= p
    N_prime = (Dk * R) % M

    if (k % 4) == 2:
        Nk = N_prime
    else:
        Nk = N_prime - M

    return Fraction(Nk, Dk)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int)

    args = parser.parse_args()
    print(rational(args.n))
