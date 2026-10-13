"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class SignalPExecutionError(Exception):
    """Failed to run SignalP-6.0 locally: missing installation, failed/
    timed-out subprocess, or the output was not generated / does not
    match the expected format."""
