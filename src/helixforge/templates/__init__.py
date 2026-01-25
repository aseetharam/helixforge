"""Template management for HelixForge.

This module provides access to templates for generating job scripts
and wrapper scripts for parallel execution.

Available templates:
    - chunk_wrapper.sh: Generic wrapper script for chunk processing
    - slurm/*.j2: SLURM-specific templates (for advanced users)

Example:
    >>> from helixforge.templates import get_wrapper_template
    >>> wrapper = get_wrapper_template()
    >>> print(wrapper[:100])  # View beginning of template
"""

from pathlib import Path
from typing import Any

try:
    from jinja2 import Environment, FileSystemLoader, Template
    HAS_JINJA2 = True
except ImportError:
    HAS_JINJA2 = False
    Environment = None
    FileSystemLoader = None
    Template = None

# Template directory
TEMPLATE_DIR = Path(__file__).parent


def get_template_dir() -> Path:
    """Get the templates directory path.

    Returns:
        Path to the templates directory.
    """
    return TEMPLATE_DIR


def get_environment() -> "Environment":
    """Get Jinja2 environment for template rendering.

    Returns:
        Configured Jinja2 Environment.

    Raises:
        ImportError: If jinja2 is not installed.
    """
    if not HAS_JINJA2:
        raise ImportError(
            "jinja2 is required for template rendering. "
            "Install with: pip install jinja2"
        )

    return Environment(
        loader=FileSystemLoader(str(TEMPLATE_DIR)),
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )


def get_template(name: str) -> "Template":
    """Get a template by name.

    Args:
        name: Template name relative to templates directory.
              e.g., "slurm/submit.sh.j2"

    Returns:
        Jinja2 Template object.

    Raises:
        ImportError: If jinja2 is not installed.
        FileNotFoundError: If template doesn't exist.
    """
    env = get_environment()
    return env.get_template(name)


def render_template(name: str, **kwargs: Any) -> str:
    """Render a template with given variables.

    Args:
        name: Template name relative to templates directory.
        **kwargs: Template variables.

    Returns:
        Rendered template string.

    Raises:
        ImportError: If jinja2 is not installed.
        FileNotFoundError: If template doesn't exist.
    """
    template = get_template(name)
    return template.render(**kwargs)


def list_templates(subdir: str | None = None) -> list[str]:
    """List available templates.

    Args:
        subdir: Optional subdirectory to list (e.g., "slurm").

    Returns:
        List of template names.
    """
    base = TEMPLATE_DIR
    if subdir:
        base = base / subdir

    templates = []
    if base.exists():
        for path in base.rglob("*.j2"):
            rel_path = path.relative_to(TEMPLATE_DIR)
            templates.append(str(rel_path))

    return sorted(templates)


def get_wrapper_template() -> str:
    """Get the chunk wrapper script template.

    Returns:
        Content of chunk_wrapper.sh template.

    Example:
        >>> wrapper = get_wrapper_template()
        >>> # Customize and save
        >>> Path("my_wrapper.sh").write_text(wrapper)
    """
    wrapper_path = TEMPLATE_DIR / "chunk_wrapper.sh"
    return wrapper_path.read_text()


def copy_wrapper_template(output_path: Path | str) -> Path:
    """Copy wrapper template to specified location.

    Args:
        output_path: Where to copy the template.

    Returns:
        Path to the copied template.
    """
    import shutil

    output_path = Path(output_path)
    wrapper_path = TEMPLATE_DIR / "chunk_wrapper.sh"
    shutil.copy(wrapper_path, output_path)
    return output_path


__all__ = [
    "get_template_dir",
    "get_environment",
    "get_template",
    "render_template",
    "list_templates",
    "get_wrapper_template",
    "copy_wrapper_template",
    "TEMPLATE_DIR",
]
