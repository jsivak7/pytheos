from pytheos.utils import sum


def test_sum() -> None:
    assert sum([1, 2, 3]) == 6


def test_sum_2() -> None:
    assert sum((1, 2)) == 3
