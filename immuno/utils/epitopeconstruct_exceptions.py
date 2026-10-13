"""Exception hierarchy for this plugin."""


class EpitopeConstructError(Exception):
    """Failed to assemble a multi-epitope construct: no candidates
    available in any class, or an internal invariant was violated."""
