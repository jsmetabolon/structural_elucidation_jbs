"""The manuscript's identification-confidence tier claims must track the audit.

Nothing checked them before: ``audit_statistics`` covers recovery, the
shown/withheld split and pathway agreement, and Table S2 is regenerated from the
same artifact, so the tier sentence in Results and the Figure 4 caption beside it
were the only tier claims in the package that could drift silently.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


V = load_tool("verify_paper_package")

PRESENTATION = json.loads(
    (ROOT / "docs" / "paper" / "data" / "presentation" / "audit_presentation_summary.json").read_text()
)
MANUSCRIPT = (ROOT / "docs" / "paper" / "manuscript.md").read_text()


def presentation(**tiers):
    return {"tier_distribution": tiers}


# --- the phrase builder -------------------------------------------------------


def test_phrase_renders_tiers_in_the_order_results_uses():
    phrase = V.tier_distribution_phrase(presentation(Confirmed=44, Moderate=1416, Low=1344, Tentative=1093))
    assert phrase == "Confirmed 44 / High 0 / Moderate 1,416 / Low 1,344 / Tentative 1,093"


def test_an_absent_tier_reads_as_zero_rather_than_raising():
    """``High`` is missing from the artifact whenever nothing reaches it."""
    assert "High 0" in V.tier_distribution_phrase(presentation(Confirmed=1))


def test_phrase_uses_thousands_separators_as_the_prose_does():
    assert "Moderate 1,416" in V.tier_distribution_phrase(presentation(Moderate=1416))


# --- the guard ----------------------------------------------------------------


def test_the_committed_manuscript_matches_the_committed_audit():
    assert V.tier_errors(MANUSCRIPT, PRESENTATION) == []


def test_a_drifted_count_is_caught():
    moved = json.loads(json.dumps(PRESENTATION))
    moved["tier_distribution"]["Moderate"] += 1
    problems = V.tier_errors(MANUSCRIPT, moved)
    assert problems and "tier distribution drifted" in problems[0]


def test_a_filled_high_tier_is_caught_by_the_count_check():
    moved = json.loads(json.dumps(PRESENTATION))
    moved["tier_distribution"]["High"] = 12
    assert any("tier distribution drifted" in p for p in V.tier_errors(MANUSCRIPT, moved))


def test_a_filled_high_tier_also_flags_the_caption_that_argues_from_emptiness():
    """The caption is prose, so renumbering the sentence above it is not enough."""
    moved = json.loads(json.dumps(PRESENTATION))
    moved["tier_distribution"]["High"] = 12
    problems = V.tier_errors(MANUSCRIPT, moved)
    assert any("argues from that emptiness" in p for p in problems), problems
    assert any("leaves the strongest confidence tier empty" in p for p in problems)


# --- artifacts CI cannot regenerate -------------------------------------------

PAPER = ROOT / "docs" / "paper"


def test_the_committed_downstream_artifacts_agree_with_the_audit():
    """baseline_claims and claim_evidence are bundle-gated, so nothing else checks them."""
    assert V.downstream_tier_errors(PAPER, PRESENTATION) == []


def test_a_drifted_downstream_claim_is_caught():
    moved = json.loads(json.dumps(PRESENTATION))
    moved["tier_distribution"]["Moderate"] += 7
    problems = V.downstream_tier_errors(PAPER, moved)
    assert any("baseline_claims.json" in p for p in problems), problems
    assert any("claim_evidence.csv" in p for p in problems), problems


def test_a_stale_figure_is_caught():
    """Figure 4 bakes its counts in as text, and readers take numbers from it."""
    moved = json.loads(json.dumps(PRESENTATION))
    moved["tier_distribution"]["Tentative"] += 3
    assert any("fig4" in p for p in V.downstream_tier_errors(PAPER, moved))


def test_the_caption_guard_stays_quiet_while_high_is_empty():
    for claim in V.EMPTY_HIGH_CLAIMS:
        assert claim in MANUSCRIPT, f"caption wording changed; the guard now watches nothing: {claim}"
    assert not any("argues from that emptiness" in p for p in V.tier_errors(MANUSCRIPT, PRESENTATION))
