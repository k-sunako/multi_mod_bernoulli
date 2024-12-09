"""."""

import math
from fractions import Fraction

import sympy

from src.crt import crt


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

def get_primitive_root_old(p):
    for i in sympy.primerange(3, p+1):
        x = i
        for j in range(2, p):
            x *= i
            x %= p
            if x == 1:
                break
        if j == p-1:
            return i


def get_primitive_root(m):
    if m == 2: return 1
    if m == 167772161: return 3
    if m == 469762049: return 3
    if m == 754974721: return 11
    if m == 998244353: return 3

    # m-1の素因数抽出
    divs = [2]
    x = (m - 1) // 2
    while x % 2 == 0: x //= 2
    i = 3
    while i * i <= x:
        if x % i == 0:
            divs.append(i)
            while x % i == 0: x //= i
        i += 2
    if x > 1: divs.append(x)

    # 全ての素因数で1と合同でない最小のgを探す
    g = 2
    while True:
        if all(pow(g, (m - 1) // div, m) != 1 for div in divs): return g
        g += 1

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

    R = crt(list(rp.keys()), list(rp.values()))
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
