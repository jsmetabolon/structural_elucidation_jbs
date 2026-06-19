"""Asynchronous discovery API package for the NCE REST service.

Queues compound-discovery jobs by curation set and library id, then hands
work to the background dispatcher for miner execution and persistence.
"""

