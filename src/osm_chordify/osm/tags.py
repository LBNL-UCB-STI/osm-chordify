"""Tag standardization functions for OSM data."""

import re

import pandas as pd


def parse_other_tags(other_tags):
    """Parse the 'other_tags' column from OSM PBF file to extract key-value pairs.

    Parameters
    ----------
    other_tags : str or None
        Raw other_tags string in ``"key"=>"value"`` format.

    Returns
    -------
    dict
        Parsed key-value pairs, or empty dict for null/empty input.
    """
    if not other_tags or pd.isna(other_tags):
        return {}
    pattern = r'"([^"]+)"=>"([^"]+)"'
    matches = re.findall(pattern, other_tags)
    return {key: value for key, value in matches}


def extract_tag_as_float(tags_dict, key):
    """Extract a numeric value from a parsed tags dictionary.

    Parameters
    ----------
    tags_dict : dict
        Dictionary of tag key-value pairs (from :func:`parse_other_tags`).
    key : str
        Tag key to look up (e.g. ``"length"``).

    Returns
    -------
    float or None
        The numeric value, or ``None`` if the key is missing or not numeric.
    """
    value_str = tags_dict.get(key, None)
    if value_str is None:
        return None
    try:
        return float(value_str)
    except (ValueError, TypeError):
        return None


def convert_weight(value: float, from_unit: str, to_unit: str) -> float:
    """Convert weight between different units."""
    # Conversion factors
    conversions = {
        "lbs_to_kg": 0.453592,
        "kg_to_lbs": 2.20462,
        "tons_to_kg": 1000,
        "kg_to_tons": 0.001
    }

    if from_unit == to_unit:
        return value

    conversion_key = f"{from_unit}_to_{to_unit}"
    if conversion_key in conversions:
        return value * conversions[conversion_key]

    # Handle two-step conversions if needed
    if from_unit == "lbs" and to_unit == "tons":
        return value * conversions["lbs_to_kg"] * conversions["kg_to_tons"]
    if from_unit == "tons" and to_unit == "lbs":
        return value * conversions["tons_to_kg"] * conversions["kg_to_lbs"]

    raise ValueError(f"Unsupported conversion from {from_unit} to {to_unit}")


def standardize_weight(weight_str: str, target_unit: str) -> float:
    """Convert weight string to numeric value in target unit."""
    if pd.isna(weight_str):
        return None

    # Handle numeric-only strings (assume they're in target unit)
    if str(weight_str).replace('.', '').isdigit():
        return float(weight_str)

    # Enhanced pattern to match more formats:
    # - Numbers with optional decimal point
    # - Various unit formats: tons, ton, t, kg, lbs, lb, st (stone)
    match = re.match(r'(\d+\.?\d*)\s*(tons?|t|kg|lbs?|st|stone)', str(weight_str).lower())
    if not match:
        # Try again with a simpler pattern in case the unit is missing or unusual
        number_match = re.match(r'(\d+\.?\d*)', str(weight_str))
        if number_match:
            # If we can extract just a number, assume it's in the target unit
            return float(number_match.group(1))
        return None

    value, unit = match.groups()
    value = float(value)

    # Standardize unit names
    unit_mapping = {
        't': 'tons',
        'ton': 'tons',
        'lb': 'lbs',
        'kg': 'kg',
        'st': 'stone',
        'stone': 'stone'
    }
    unit = unit_mapping.get(unit, unit)

    # Convert to standard unit first (kg), then to target unit
    # Conversion factors to kg
    to_kg = {
        'lbs': 0.453592,
        'kg': 1.0,
        'tons': 1000.0,
        'stone': 6.35029  # 1 stone = 14 lbs = 6.35029 kg
    }

    # Convert to kg first
    weight_in_kg = value * to_kg.get(unit, 1.0)

    # Then convert from kg to target unit
    from_kg = {
        'lbs': 2.20462,
        'kg': 1.0,
        'tons': 0.001
    }

    # If target unit isn't recognized, default to kg
    conversion_factor = from_kg.get(target_unit, 1.0)

    return weight_in_kg * conversion_factor


def standardize_oneway(value):
    """
    Standardize oneway tag to "yes", "reverse", or "no" strings to match MATSim expectations.
    - "yes" for forward direction oneway
    - "reverse" for backward direction oneway
    - "no" for bidirectional
    """
    # Return None if the value is None or empty
    if value is None or value == '':
        return "no"

    # Values that explicitly mean "yes" (forward oneway)
    valid_yes = {'yes', 'true', '1', True, 1}

    # Values that explicitly mean reverse oneway
    valid_reverse = {'-1', 'reverse'}

    # Values that explicitly mean "no"
    valid_no = {'no', 'false', '0', False, 0}

    # Handle strings
    if isinstance(value, str):
        value = value.lower().strip()
        # Handle semicolon-separated values
        if ';' in value:
            parts = [part.strip().lower() for part in value.split(';')]
            if all(part in valid_yes for part in parts):
                return "yes"
            elif all(part in valid_reverse for part in parts):
                return "reverse"
            else:
                return "no"

        # Handle single string
        if value in valid_yes:
            return "yes"
        elif value in valid_reverse:
            return "reverse"
        elif value in valid_no:
            return "no"
        else:
            # If we can't interpret it, MATSim logs a warning and ignores it
            return "no"

    # Handle boolean and numeric
    if isinstance(value, (bool, int)):
        if value in valid_yes:
            return "yes"
        else:
            return "no"

    # Handle lists (if that's a use case)
    if isinstance(value, list):
        if all(str(v).lower().strip() in valid_yes for v in value if v):
            return "yes"
        elif all(str(v).lower().strip() in valid_reverse for v in value if v):
            return "reverse"
        else:
            return "no"

    # Default case
    return "no"


def standardize_motor_vehicle(value):
    """
    Standardize motor_vehicle tag to "yes" or "no" strings, focusing on a defined set of restrictive values.

    Parameters
    ----------
    value : any
        The motor_vehicle tag value

    Returns
    -------
    str
        "no" if motor vehicles are restricted (no, false, 0, private)
        "yes" otherwise
    """
    # Define restrictive values
    restrictive_values = {"no", "false", "0"}

    # If value is None or empty, assume motor vehicles are allowed
    if value is None or pd.isna(value) or (isinstance(value, str) and not value.strip()):
        return "yes"

    # Convert to string and lowercase for consistent processing
    if not isinstance(value, str):
        value = str(value)

    value = value.lower().strip()

    # Handle special cases with multiple values (separated by semicolons or vertical bars)
    if ';' in value or '|' in value:
        # Split by either semicolon or vertical bar
        parts = re.split(r'[;|]+', value)
        parts = [p.strip() for p in parts if p.strip()]

        # If any part is in the restrictive values, the overall value is "no"
        if any(p in restrictive_values for p in parts):
            return "no"
        else:
            return "yes"

    # Check if the value is in the restrictive set
    if value in restrictive_values:
        return "no"

    # All other values indicate some form of access
    return "yes"


def standardize_maxspeed(value, default_kph=None):
    """
    Standardize maxspeed values and return them in the format "25 mph".

    Parameters
    ----------
    value : any
        The maxspeed tag value
    default_kph : int, optional
        Default speed in kph to use if the value can't be parsed

    Returns
    -------
    str or None
        Speed in format "XX mph", or None if the value can't be parsed and no default is provided
    """
    if value is None or pd.isna(value) or (isinstance(value, str) and not value.strip()):
        if default_kph is not None:
            return f"{round(default_kph / 1.60934)} mph"  # Convert kph to mph
        return None

    # Convert to string for processing
    if not isinstance(value, str):
        value = str(value)

    value = value.lower().strip()

    # Handle special cases
    if value == "signals" or value == "none" or value == "variable":
        if default_kph is not None:
            return f"{round(default_kph / 1.60934)} mph"  # Convert kph to mph
        return None

    # Try to extract numeric value and unit
    match = re.match(r'^(\d+(?:\.\d+)?)\s*(mph|kmh|km/h|kph)?$', value)
    if match:
        speed_val = float(match.group(1))
        unit = match.group(2) if match.group(2) else "kph"  # Default to kph if no unit

        # Convert to mph if necessary
        if unit in ["kmh", "km/h", "kph"]:
            speed_mph = round(speed_val / 1.60934)  # Convert kph to mph
        else:
            # Already in mph
            speed_mph = round(speed_val)

        return f"{speed_mph} mph"

    # If we can't parse the value and have a default
    if default_kph is not None:
        return f"{round(default_kph / 1.60934)} mph"  # Convert kph to mph

    # If we can't parse the value and don't have a default
    return None


def standardize_access(value):
    """
    Standardize access tag to "yes" or "no" strings, focusing on a defined set of restrictive values.

    Parameters
    ----------
    value : any
        The access tag value

    Returns
    -------
    str
        "no" if access is restricted (no, private, forestry, permit, etc.)
        "yes" otherwise
    """
    # Define restrictive values - values that indicate restricted access
    restrictive_values = {"no", "false", "0"}

    # If value is None or empty, assume access is allowed
    if value is None or pd.isna(value) or (isinstance(value, str) and not value.strip()):
        return "yes"

    # Convert to string and lowercase for consistent processing
    if not isinstance(value, str):
        value = str(value)

    value = value.lower().strip()

    # Handle special cases with multiple values (separated by semicolons or vertical bars)
    if ';' in value or '|' in value:
        # Split by either semicolon or vertical bar
        parts = re.split(r'[;|]+', value)
        parts = [p.strip() for p in parts if p.strip()]

        # If any part is in the restrictive values, the overall value is "no"
        if any(p in restrictive_values for p in parts):
            return "no"
        else:
            return "yes"

    # Check if the value is in the restrictive set
    if value in restrictive_values:
        return "no"

    # All other values (yes, permissive, etc.) indicate general access
    return "yes"


def standardize_hgv(value):
    """
    Standardize HGV access values to boolean (True/False)
    Returns True if HGVs are allowed, False if they are not
    """
    if not value:
        return True  # Default to allowed if no value

    # Values that indicate HGV prohibition
    restrictive_values = {"no", "false", "0"}

    # Handle boolean inputs
    if isinstance(value, bool):
        return value

    # Handle semicolon-separated string values
    if isinstance(value, str) and ';' in value:
        # If any part is "no", the whole is restricted
        for part in value.split(';'):
            if part.strip().lower() in restrictive_values:
                return False
        return True

    # Handle list case
    if isinstance(value, list):
        if not value:
            return True
        # If any value is "no", the whole is restricted
        for v in value:
            if str(v).strip().lower() in restrictive_values:
                return False
        return True

    # Handle single string value
    if isinstance(value, str):
        return str(value).strip().lower() not in restrictive_values

    # For any other case, convert to string and check
    return str(value).strip().lower() not in restrictive_values
