"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class NetCleaveExecutionError(Exception):
    """Failed to run NetCleave locally: missing installation, failed/
    timed-out subprocess, or the output .xlsx was not generated / does not
    match the expected format."""
