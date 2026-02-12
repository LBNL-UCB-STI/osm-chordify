"""Custom aggregation functions for graph simplification."""

import math
import re
from statistics import mean, median

import pandas as pd


def most_restrictive_bool_str(values):
    """Returns "no" if any value is "no", otherwise returns "yes" (or True/False logic)."""
    # Filters out None/NaN/empty strings, converts booleans/numbers to strings for safety
    valid_values = [str(v).strip().lower() for v in values if pd.notna(v) and str(v).strip()]
    if not valid_values:
        return None

    # If any value is a clear restriction, enforce restriction.
    return "no" if "no" in valid_values or "false" in valid_values or "0" in valid_values else "yes"


def min_numeric_or_string(values):
    """
    Finds the minimum numeric value from a list, ignoring non-numeric strings.
    If no numeric value is found, returns the first non-NaN string.
    Used for 'maxweight' where the minimum constraint is the most restrictive.
    """
    numeric_values = []
    first_string = None

    for v in values:
        if pd.isna(v):
            continue

        try:
            # Try to convert to float (handles strings like "1000")
            numeric_v = float(v)
            if not math.isnan(numeric_v):
                numeric_values.append(numeric_v)
        except (ValueError, TypeError):
            # If conversion fails, it's a string (e.g., "30 tons" or "5000 kg").
            if first_string is None and isinstance(v, str):
                first_string = v
            continue

    if numeric_values:
        return min(numeric_values)

    # Fallback: if no valid numeric value, return the first encountered string (e.g., "5000 kg")
    return first_string if first_string is not None else None


def first_valid_value(values):
    """
    Returns the first non-NaN, non-empty value encountered.
    Used for dimensions like 'maxheight', 'maxwidth' where a single value is sufficient
    and aggregation is complex.
    """
    for v in values:
        if pd.notna(v) and str(v).strip():
            return v
    return None


def median_lanes(values):
    """
    Calculate median after converting string values to numbers.
    Handles:
    - Lists of values
    - Semicolon-separated values
    - Mixed numeric types
    """
    # Initialize empty list for numeric values
    numeric_values = []

    # Handle case where values is already a single value, not an iterable
    if isinstance(values, (int, float)):
        return int(values)
    elif isinstance(values, str):
        values = [values]

    # Process each value in the iterable
    for v in values:
        # Skip None values
        if v is None:
            continue

        # Handle different types
        if isinstance(v, (int, float)):
            numeric_values.append(int(v))
            continue

        if not isinstance(v, str):
            v = str(v)

        # Split by semicolon to handle multiple values
        parts = v.split(';')
        for part in parts:
            part = part.strip()
            try:
                numeric_values.append(int(part))
            except (ValueError, TypeError):
                # Skip non-numeric parts
                continue

    if not numeric_values:
        return None
    return int(median(numeric_values))


def most_restrictive_access(values):
    """
    Returns the most restrictive access value from a list based on a predefined priority order.

    Parameters
    ----------
    values : list
        List of access values

    Returns
    -------
    str
        The most restrictive access value, or None if no valid values
    """
    # Define a priority order for access restrictions (from most to least restrictive)
    priority = {
        "no": 1,  # Most restrictive
        "private": 2,
        "permit": 3,
        "destination": 4,
        "delivery": 5,
        "customers": 6,
        "forestry": 7,
        "agricultural": 8,
        "discouraged": 9,
        "permissive": 10,
        "yes": 11  # Least restrictive
    }

    # Default priority for unknown values - place between "discouraged" and "permissive"
    default_priority = 9.5

    if not values:
        return None

    # Process each value and find the most restrictive
    most_restrictive = None
    highest_priority = float('inf')  # Lower number = higher priority

    for value in values:
        if value is None or value == "nan" or pd.isna(value):
            continue

        if isinstance(value, str):
            value = value.strip().lower()
            if not value or value == "nan":
                continue

        # Get priority for this value
        value_priority = priority.get(value, default_priority)

        # Update most restrictive if this has higher priority (lower number)
        if value_priority < highest_priority:
            most_restrictive = value
            highest_priority = value_priority

    return most_restrictive


def bool_all(values):
    """
    Returns False if any value is False, otherwise returns True.
    Expects only boolean values (True or False).

    Parameters
    ----------
    values : list
        List of boolean values

    Returns
    -------
    bool
        False if any value is False, True otherwise
    """
    if not values:
        return None

    # If any value is False, return False
    return all(values)


def mean_maxspeed(speed_values):
    """
    Calculate the mean speed from a list of speed values in the format "XX mph".

    Parameters
    ----------
    speed_values : list
        List of speed values in format "XX mph"

    Returns
    -------
    str
        Mean speed in format "XX mph", or None if no valid speeds found
    """
    if not speed_values:
        return None

    # Extract numeric values
    speeds_mph = []

    for value in speed_values:
        if not value or pd.isna(value):
            continue

        # Convert to string if needed
        if not isinstance(value, str):
            value = str(value)

        # Extract the numeric part
        match = re.match(r'^(\d+(?:\.\d+)?)\s*mph$', value.lower().strip())
        if match:
            speeds_mph.append(float(match.group(1)))

    # Calculate mean if we have valid values
    if speeds_mph:
        mean_speed = mean(speeds_mph)
        return f"{round(mean_speed)} mph"

    return None


def yes_no_all(values):
    """
    Returns "no" if any value is "no", otherwise returns "yes".
    Expects string values ("yes" or "no").

    Parameters
    ----------
    values : list
        List of string values ("yes" or "no")

    Returns
    -------
    str
        "no" if any value is "no", "yes" otherwise
    """
    if not values:
        return None

    # If any value is "no", return "no"
    return "no" if "no" in values else "yes"
