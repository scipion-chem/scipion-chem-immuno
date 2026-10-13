"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class DiscoTopeExecutionError(Exception):
    """Failed to run DiscoTope-3.0 locally: missing installation, failed/
    timed-out subprocess, or the output CSV was not generated / does not
    match the expected format."""
