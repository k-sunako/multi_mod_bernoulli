"""."""

from src.bernoulli import calc_mod
from concurrent.futures import ThreadPoolExecutor
import sympy


def test_threads():
    X = 1000
    k = 10
    rp = {}
    for p in sympy.primerange(X+1):
        if k % (p-1) != 0:
            rp[p] = calc_mod(p, k)

    executor = ThreadPoolExecutor(max_workers=4)
    futures = []
    rp_thread = {}
    for p in sympy.primerange(X+1):
        if k % (p-1) != 0:
            future = executor.submit(calc_mod, p, k)
            futures.append(future)
            rp_thread[p] = 1
    for p, future in zip(rp.keys(), futures):
        rp_thread[p] = future.result()

    for p in rp.keys():
        assert rp[p] == rp_thread[p]
