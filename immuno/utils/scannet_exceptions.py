"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class ScanNetExecutionError(Exception):
    """Failed to run ScanNet locally: missing installation (either
    runtime), failed/timed-out subprocess, or the output CSV was not
    generated / does not match the expected format."""
