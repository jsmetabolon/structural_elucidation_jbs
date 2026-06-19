"""Flask REST API package for Metabolon compound discovery (Sally service).

Vendored HTTP layer that queues and serves discovery jobs against
curation sets and HD4 libraries. Not invoked directly by spec2mol
inference, but shares the miner/elucidation stack under ``vendor/nce/``.
"""
