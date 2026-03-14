"""
Regression tests for simplify.py aggregation functions.
"""

import pytest

from osm_chordify.osm.simplify import bool_all, mean_maxspeed, most_restrictive_access


class TestBoolAll:
    """
    Regression (L1): bool_all used to return None for empty input.
    Python's all([]) returns True (vacuous truth); bool_all should match.
    """

    def test_all_true(self):
        assert bool_all([True, True, True]) is True

    def test_any_false(self):
        assert bool_all([True, False, True]) is False

    def test_single_false(self):
        assert bool_all([False]) is False

    def test_single_true(self):
        assert bool_all([True]) is True

    # --- regression: empty list must return True, not None ------------------

    def test_empty_list_returns_true(self):
        """
        Regression: bool_all([]) returned None instead of True.
        Vacuous truth: no False value means allowed (True).
        """
        result = bool_all([])
        assert result is True, (
            f"bool_all([]) returned {result!r}, expected True (vacuous truth)"
        )

    def test_return_type_is_bool_not_none(self):
        assert isinstance(bool_all([True]), bool)
        assert isinstance(bool_all([False]), bool)
        assert isinstance(bool_all([]), bool)


class TestMeanMaxspeed:
    def test_single_value(self):
        assert mean_maxspeed(["30 mph"]) == "30 mph"

    def test_averages_two_values(self):
        assert mean_maxspeed(["20 mph", "40 mph"]) == "30 mph"

    def test_empty_returns_none(self):
        assert mean_maxspeed([]) is None

    def test_all_none_returns_none(self):
        assert mean_maxspeed([None, None]) is None

    def test_skips_none_elements(self):
        assert mean_maxspeed([None, "60 mph"]) == "60 mph"


class TestMostRestrictiveAccess:
    def test_no_is_most_restrictive(self):
        assert most_restrictive_access(["yes", "no", "permissive"]) == "no"

    def test_private_more_restrictive_than_permissive(self):
        assert most_restrictive_access(["private", "permissive"]) == "private"

    def test_yes_when_only_option(self):
        assert most_restrictive_access(["yes"]) == "yes"

    def test_empty_returns_none(self):
        assert most_restrictive_access([]) is None

    def test_all_none_returns_none(self):
        assert most_restrictive_access([None]) is None