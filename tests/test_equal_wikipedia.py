import pytest

from fractions import Fraction
from src.main import rational

def test_n():
    assert rational(2) == Fraction(1, 6)
    assert rational(4) == Fraction(-1, 30)
    assert rational(6) == Fraction(1, 42)
    assert rational(8) == Fraction(-1, 30)
    assert rational(10) == Fraction(5, 66)
    assert rational(12) == Fraction(-691, 2730)
    assert rational(14) == Fraction(7, 6)
    assert rational(16) == Fraction(-3617, 510)
    assert rational(18) == Fraction(43867, 798)
    assert rational(20) == Fraction(-174611, 330)
