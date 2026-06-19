"""Background dispatcher for async discovery work items.

Drains queued ``UnitOfWork`` jobs, runs miner discovery with caching, and
delegates persistence and notification to companion modules.
"""

