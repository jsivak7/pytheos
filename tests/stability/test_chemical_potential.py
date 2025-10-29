from pytheos.stability import chemical_potential


def test_calc_overlap():

    fake_ranges = [
        (-1, 1),
        (-5, 5),
        (-0.5, 2),
    ]

    overlap_dict = chemical_potential.calc_overlap(ranges=fake_ranges)

    assert overlap_dict["overlap"] == 1.5
    assert overlap_dict["lower bound"] == -0.5
    assert overlap_dict["upper bound"] == 1
