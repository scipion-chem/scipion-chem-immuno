"""Exception hierarchy for this plugin: never let a raw
FileNotFoundError/CalledProcessError escape to the Scipion GUI without an
actionable message.
"""


class TMbedExecutionError(Exception):
    """Failed to run TMbed locally: missing installation, failed subprocess,
    or the output prediction file was not generated."""


class TMbedParseError(Exception):
    """The TMbed output file does not match the expected 3-line
    (header/sequence/prediction) format for the requested --out-format."""
