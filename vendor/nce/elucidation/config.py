"""Filesystem paths for NCE discovery and delivery output.

``BASE_DIR`` and ``DELIVERY_DIR`` default to Metabolon network shares and
can be overridden via environment variables for local or EC2 runs.
"""

import os

BASE_DIR = os.getenv("BASE_DIR", "L:\\New Compound Discovery")
DELIVERY_DIR = os.getenv("DELIVERY_DIR", "L:\\New Compound Discovery")
REFERENCE_DIR = os.getenv("REFERENCE_DIR", "L:\\New Compound Discovery\\references\\production")
