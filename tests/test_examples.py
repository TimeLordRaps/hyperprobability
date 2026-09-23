import runpy
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

EXPECTED = {
    "coin_and_clock.py": [
        "hyperprobability from 'tock'  1",
        "hyperprobability from 'T'  ω",
        "ordinal probability         1/ω",
        "Kac bridge                  1/2 <= 1/2, tight: True",
        "Kac bridge                  1/ω <= 1/2, tight: False",
        "with probability 1/2, 1/4, 1/8, 1/16, 1/32, ...",
    ],
    "strange_loop.py": [
        "hyperprobability              ∞",
        "hyperprobability              ω",
        "hyperprobability              ω + 3",
        "hyperprobability              ω^2",
        "ordinal probability           1/(ω + 3)",
        "first reached at stage ω^2    1",
    ],
}


def test_every_example_is_checked():
    assert sorted(path.name for path in EXAMPLES.glob("*.py")) == sorted(EXPECTED)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_example(name, capsys):
    runpy.run_path(str(EXAMPLES / name), run_name="__main__")
    output = capsys.readouterr().out
    for line in EXPECTED[name]:
        assert line in output
