"""Mass-error windows behind the mass-match sub-score.

Standard library only, and separate from ``predict`` on purpose. ``predict``
pulls the inference stack, so anything that merely needs to *describe* the
scoring windows — the exported parameter table, its tests, CI — cannot import
from there. Splitting the two numbers out keeps the description and the scorer
reading from one definition instead of restating each other.
"""

from __future__ import annotations

# Error at or below which the mass sub-score is a flat 1.0.
MASS_SCORE_PLATEAU_DA: float = 0.001

# Error at or above which the mass sub-score is 0.0. Between the two the score
# falls linearly, so it lands on 0.80 the moment it leaves the plateau: no
# candidate can score between 0.80 and 1.00. Worth knowing before reading any
# cutoff placed in that interval as if it were reachable.
MASS_SCORE_TOLERANCE_DA: float = 0.005
