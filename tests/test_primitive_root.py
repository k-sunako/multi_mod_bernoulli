from src.bernoulli import get_primitive_root_01, get_primitive_root_02

def test_do():
    assert get_primitive_root_01(13) == get_primitive_root_02(13)
