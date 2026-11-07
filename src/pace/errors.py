"""Errors safe to present at the command line."""


class PaceError(ValueError):
    """An invalid scientific input, unsupported contract, or unavailable asset."""
