"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class EpiDopeExecutionError(Exception):
    """Failed to run EpiDope locally: missing installation, failed/timed-out
    subprocess, or the per-accession score CSV was not generated / does not
    match the expected format."""
