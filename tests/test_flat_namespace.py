import typing

from ALifeStdDev import ALifeStdDev as asd


def test_access():
    assert hasattr(asd, "load_phylogeny_to_networkx")
    assert isinstance(asd.load_phylogeny_to_networkx, typing.Callable)
