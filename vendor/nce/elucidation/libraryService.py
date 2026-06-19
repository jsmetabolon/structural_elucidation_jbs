"""LibraryService facade selecting the active MLIMS data provider.

Re-exports MlimsLibraryService as LibraryService for live Oracle queries during
spec2mol v21 inference. Alternate cached or testing backends are commented out.
"""

from elucidation.library_services.mlimsLibraryService import MlimsLibraryService as LibraryService
