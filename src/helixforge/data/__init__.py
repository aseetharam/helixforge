"""Package data files for HelixForge.

This package contains data files used by HelixForge, including:
- pwm_plant.json: Plant splice site position weight matrices
"""

from importlib import resources
from pathlib import Path


def get_data_path(filename: str) -> Path:
    """Get the path to a data file.

    Args:
        filename: Name of the data file.

    Returns:
        Path to the data file.
    """
    return resources.files(__name__).joinpath(filename)


def load_pwm_plant() -> dict:
    """Load the plant PWM data.

    Returns:
        Dictionary with donor and acceptor PWM data.
    """
    import json

    data_path = get_data_path("pwm_plant.json")
    with resources.as_file(data_path) as path:
        with open(path) as f:
            return json.load(f)


__all__ = ["get_data_path", "load_pwm_plant"]
