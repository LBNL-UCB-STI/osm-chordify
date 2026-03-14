"""
Regression tests for tags.py standardization functions.

Each test class maps to one function and covers the specific bug that was
fixed.  Comments mark the exact regression being guarded.
"""

import pytest

from osm_chordify.osm.tags import (
    standardize_access,
    standardize_hgv,
    standardize_maxspeed,
    standardize_motor_vehicle,
    standardize_oneway,
    standardize_weight,
)


# ---------------------------------------------------------------------------
# standardize_oneway
# ---------------------------------------------------------------------------


class TestStandardizeOneway:
    """
    BEAM/MATSim require "-1" for reverse-direction oneway edges, not "reverse".
    Regression: the function used to return "reverse" for OSM value "-1".
    """

    def test_forward_yes(self):
        assert standardize_oneway("yes") == "yes"

    def test_forward_true(self):
        assert standardize_oneway("true") == "yes"

    def test_forward_1(self):
        assert standardize_oneway("1") == "yes"

    def test_forward_bool_true(self):
        assert standardize_oneway(True) == "yes"

    # --- regression: must be "-1", not "reverse" ----------------------------

    def test_reverse_minus_one_string(self):
        """OSM "-1" must produce the string "-1", not "reverse"."""
        assert standardize_oneway("-1") == "-1"

    def test_reverse_reverse_string(self):
        """OSM "reverse" must also produce "-1"."""
        assert standardize_oneway("reverse") == "-1"

    # --- bidirectional ------------------------------------------------------

    def test_no(self):
        assert standardize_oneway("no") == "no"

    def test_false(self):
        assert standardize_oneway("false") == "no"

    def test_none(self):
        assert standardize_oneway(None) == "no"

    def test_empty_string(self):
        assert standardize_oneway("") == "no"

    def test_unknown_value(self):
        assert standardize_oneway("alternating") == "no"

    # --- semicolon-separated ------------------------------------------------

    def test_semicolon_all_yes(self):
        assert standardize_oneway("yes;yes") == "yes"

    def test_semicolon_all_reverse(self):
        assert standardize_oneway("-1;reverse") == "-1"

    def test_semicolon_mixed(self):
        assert standardize_oneway("yes;no") == "no"


# ---------------------------------------------------------------------------
# standardize_maxspeed
# ---------------------------------------------------------------------------


class TestStandardizeMaxspeed:
    """
    Regression: list-valued maxspeed (from simplify_graph merging) was not
    handled — the function would crash or return None instead of averaging.
    """

    # --- scalar inputs ------------------------------------------------------

    def test_numeric_kph_default_unit(self):
        """Bare number assumed kph, converted to mph."""
        result = standardize_maxspeed("50")
        assert result == "31 mph"

    def test_explicit_kph(self):
        assert standardize_maxspeed("80 km/h") == "50 mph"

    def test_explicit_mph(self):
        assert standardize_maxspeed("25 mph") == "25 mph"

    def test_signals_uses_default(self):
        assert standardize_maxspeed("signals", default_kph=50) == "31 mph"

    def test_none_returns_none_without_default(self):
        assert standardize_maxspeed(None) is None

    def test_none_returns_default_when_given(self):
        assert standardize_maxspeed(None, default_kph=50) == "31 mph"

    # --- regression: list input must average, not crash ---------------------

    def test_list_two_elements_averages(self):
        """["50", "100"] in kph → each rounded (31, 62 mph) → average 46 mph.
        Note: averaging happens on already-rounded integer mph values, so
        round((31+62)/2) = round(46.5) = 46 (Python banker's rounding)."""
        result = standardize_maxspeed(["50", "100"])
        assert result == "46 mph"

    def test_list_single_element(self):
        result = standardize_maxspeed(["60"])
        assert result == "37 mph"

    def test_list_with_none_element(self):
        """None inside list falls back to default (31 mph); "60" kph → 37 mph.
        Average of (31, 37) = 34 mph."""
        result = standardize_maxspeed([None, "60"], default_kph=50)
        assert result == "34 mph"

    def test_list_all_none_falls_back_to_default(self):
        result = standardize_maxspeed([None], default_kph=50)
        assert result == "31 mph"

    def test_list_mixed_units(self):
        """["30 mph", "60 km/h"] → average of 30 and 37 → 34 mph."""
        result = standardize_maxspeed(["30 mph", "60 km/h"])
        assert result == "34 mph"


# ---------------------------------------------------------------------------
# standardize_weight
# ---------------------------------------------------------------------------


class TestStandardizeWeight:
    """
    Regression: OSM convention is that a bare numeric weight (no suffix) means
    metric tons.  The old code returned the raw float as-is, so "5" was treated
    as 5 kg rather than 5000 kg.
    """

    # --- regression: bare numeric = metric tons ----------------------------

    def test_bare_integer_treated_as_metric_tons_to_kg(self):
        """'5' must be interpreted as 5 metric tons = 5000 kg."""
        assert standardize_weight("5", "kg") == pytest.approx(5000.0)

    def test_bare_integer_to_lbs(self):
        """5 metric tons × 2204.62 lbs/ton ≈ 11023.1 lbs."""
        assert standardize_weight("5", "lbs") == pytest.approx(5 * 2204.62, rel=1e-3)

    def test_bare_decimal_treated_as_metric_tons(self):
        assert standardize_weight("3.5", "kg") == pytest.approx(3500.0)

    # --- explicit unit strings ----------------------------------------------

    def test_explicit_tons_to_kg(self):
        assert standardize_weight("5t", "kg") == pytest.approx(5000.0)

    def test_explicit_kg_passthrough(self):
        assert standardize_weight("500kg", "kg") == pytest.approx(500.0)

    def test_explicit_lbs_to_kg(self):
        assert standardize_weight("1000lbs", "kg") == pytest.approx(453.592, rel=1e-3)

    def test_explicit_lbs_to_lbs(self):
        # lbs→kg→lbs is a two-step conversion; accumulated floating point error
        # means the result is ~999.998, not exactly 1000.  Loosen tolerance.
        assert standardize_weight("1000lbs", "lbs") == pytest.approx(1000.0, rel=1e-3)

    def test_stone_to_kg(self):
        """1 stone = 6.35029 kg."""
        assert standardize_weight("1stone", "kg") == pytest.approx(6.35029, rel=1e-4)

    def test_stone_to_lbs(self):
        """1 stone = 14 lbs."""
        assert standardize_weight("1stone", "lbs") == pytest.approx(14.0, rel=1e-3)

    def test_kg_to_stone(self):
        """
        Regression (M1): 'stone' was missing from from_kg dict so kg→stone
        returned the raw kg value instead of converting.
        6.35029 kg should be ~1 stone.
        """
        assert standardize_weight("6.35029kg", "stone") == pytest.approx(1.0, rel=1e-3)

    def test_none_returns_none(self):
        assert standardize_weight(None, "kg") is None

    def test_non_numeric_returns_none(self):
        assert standardize_weight("heavy", "kg") is None


# ---------------------------------------------------------------------------
# standardize_motor_vehicle
# ---------------------------------------------------------------------------


class TestStandardizeMotorVehicle:
    """
    Regression: "private" was not in the restrictive set, so
    motor_vehicle=private was treated as "yes" (allowed) instead of "no".
    """

    def test_yes(self):
        assert standardize_motor_vehicle("yes") == "yes"

    def test_no(self):
        assert standardize_motor_vehicle("no") == "no"

    def test_false(self):
        assert standardize_motor_vehicle("false") == "no"

    # --- regression ---------------------------------------------------------

    def test_private_is_restrictive(self):
        """motor_vehicle=private must map to "no"."""
        assert standardize_motor_vehicle("private") == "no"

    # --- edge cases ---------------------------------------------------------

    def test_none_defaults_to_yes(self):
        assert standardize_motor_vehicle(None) == "yes"

    def test_permissive_is_yes(self):
        assert standardize_motor_vehicle("permissive") == "yes"

    def test_semicolon_with_private(self):
        assert standardize_motor_vehicle("yes;private") == "no"

    def test_semicolon_all_allowed(self):
        assert standardize_motor_vehicle("yes;permissive") == "yes"


# ---------------------------------------------------------------------------
# standardize_access
# ---------------------------------------------------------------------------


class TestStandardizeAccess:
    """
    Regression: same as motor_vehicle — "private" was missing from the
    restrictive set.
    """

    def test_yes(self):
        assert standardize_access("yes") == "yes"

    def test_no(self):
        assert standardize_access("no") == "no"

    # --- regression ---------------------------------------------------------

    def test_private_is_restrictive(self):
        """access=private must map to "no"."""
        assert standardize_access("private") == "no"

    # --- edge cases ---------------------------------------------------------

    def test_none_defaults_to_yes(self):
        assert standardize_access(None) == "yes"

    def test_permissive_is_yes(self):
        assert standardize_access("permissive") == "yes"

    def test_customers_is_yes(self):
        """Non-private restricted-but-allowed tags stay "yes"."""
        assert standardize_access("customers") == "yes"

    def test_semicolon_with_private(self):
        assert standardize_access("no;private") == "no"


# ---------------------------------------------------------------------------
# standardize_hgv
# ---------------------------------------------------------------------------


class TestStandardizeHgv:
    """
    Regression: the old guard was ``if not value`` which treated boolean False
    as missing and returned True (allowed), flipping ferry edges that had
    hgv=False back to allowed.
    """

    def test_yes_string(self):
        assert standardize_hgv("yes") is True

    def test_no_string(self):
        assert standardize_hgv("no") is False

    def test_bool_true(self):
        assert standardize_hgv(True) is True

    # --- regression ---------------------------------------------------------

    def test_bool_false_is_not_treated_as_missing(self):
        """
        boolean False must return False (prohibited), NOT True (allowed).
        The old ``if not value`` guard incorrectly treated False as missing.
        """
        assert standardize_hgv(False) is False

    # --- edge cases ---------------------------------------------------------

    def test_none_defaults_to_allowed(self):
        assert standardize_hgv(None) is True

    def test_empty_string_defaults_to_allowed(self):
        assert standardize_hgv("") is True

    def test_zero_string_is_restricted(self):
        assert standardize_hgv("0") is False

    def test_list_all_yes(self):
        assert standardize_hgv(["yes", "yes"]) is True

    def test_list_with_no(self):
        assert standardize_hgv(["yes", "no"]) is False

    def test_empty_list_defaults_to_allowed(self):
        assert standardize_hgv([]) is True

    def test_semicolon_with_no(self):
        assert standardize_hgv("yes;no") is False