from src.bernoulli import rational


def test_do():
    b100 = rational(100)
    assert (b100.numerator % 263) == 0
    assert (b100.numerator % 379) == 0
    assert (b100.numerator % 28717943) == 0
    assert (b100.numerator % 65677171692755556482181133) == 0

    b200 = rational(200)
    assert (b200.numerator % 389) == 0
    assert (b200.numerator % 691) == 0

    b300 = rational(300)
    assert (b300.numerator % 353) == 0

    b1000 = rational(1000)
    assert (b1000.numerator % 19919) == 0

    b2000 = rational(2000)
    assert (b2000.numerator % 1217) == 0
    assert (b2000.numerator % 306953) == 0

    b3000 = rational(3000)
    assert (b3000.numerator % 547) == 0
    assert (b3000.numerator % 40609) == 0
