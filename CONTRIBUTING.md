# Contributing to Open Exoplanet Hunt

An optional [preregistered TLS detector evaluation](discovery/docs/tls-independent-controls.md)
provides independent synthetic transit detection and adversarial controls.
Its manual hosted workflow keeps TLS dependencies separate from ordinary
verification. It does not scan astronomical targets, change the ten gates,
establish completeness or authorize new candidate claims.

Human astronomers, citizen scientists, developers, and autonomous AI coding agents (`Claude Code`, `Codex`, `Antigravity`, `Cursor`, etc.) are invited to contribute!

---

## Quick Start

```bash
git clone https://github.com/Th0rgal/open-exoplanet-hunt.git
cd open-exoplanet-hunt
uv venv .venv
uv pip install --python .venv/bin/python astropy numpy scipy matplotlib requests pytest

# 1. Verify existing catalog index & candidate certificates
make verify

# 2. Run the 10 positive & negative adversarial falsification tests
make test

# 3. Scan a new batch of unreported multi-sector TESS K/M dwarfs
.venv/bin/python -B scripts/hunt_sectors.py --max-targets 10 --min-sectors 3 --page 7
```

---

## How to Submit a Pull Request

### Option 1: Submit a Target Scan Batch or New Candidate Dossier
1. Check [`registry/scanned_targets.jsonl`](registry/scanned_targets.jsonl) and [`catalog/known_tics.json.gz`](catalog/known_tics.json.gz) so you do not duplicate already-scanned or already-known stars (`hunt_sectors.py` checks both automatically).
2. Run `scripts/hunt_sectors.py` (varying `--page`, `--tmag-min`, `--tmag-max`, `--teff-min`, `--teff-max`, or `--dist-max`) or vet a specific target with:
   ```bash
   .venv/bin/python -B scripts/vet_candidate.py --tic <TIC_ID>
   ```
3. Run `make verify && make test` to update `registry/leaderboard.json` and confirm all checks pass.
4. Open a Pull Request! Both negative scan logs (`registry/scanned_targets.jsonl` entries showing which stars were scanned and falsified) and new passing/watchlist dossiers (`candidates/TIC-<ID>/`) are valuable scientific contributions.

### Option 2: Improve the Falsification Pipeline or Index New External Catalogs
1. If adding a new vetting check to [`scripts/vet_candidate.py`](scripts/vet_candidate.py), you **must** add a synthetic negative-control unit test in [`tests/test_vetting_controls.py`](tests/test_vetting_controls.py) demonstrating that the check rejects the intended false-positive scenario without breaking the positive control.
2. If indexing new external papers, CTOI batches, or ML catalogs, add them to [`catalog/external_discoveries.json`](catalog/external_discoveries.json) and run `make sync-catalogs && make verify`.

---

## Automated & Maintainer PR Review Rubric

Every Pull Request is reviewed (by humans and autonomous review agents) against five criteria:
1. **Deterministic Verification**: `python3 -B scripts/verify_all.py` and `python3 -B -m unittest discover -s tests -v` pass with zero errors.
2. **Catalog De-duplication**: No candidate marked `unreported_candidate` or `candidate_pass_all_10` collides with an existing NASA confirmed planet, official `TOI`, or ExoFOP `CTOI`.
3. **Raw MAST FITS Replay**: For any new candidate dossier in `candidates/TIC-<ID>/`, the reviewer independently downloads the raw MAST FITS files, checks the SHA-256 hashes, and replays `scripts/vet_candidate.py --tic <ID>`.
4. **Adversarial Scrutiny**: We actively try to *disprove* every candidate (checking nearby Gaia DR3 stars within $60''$, momentum dump timestamps, eclipsing binary secondaries at eccentric phases, and stellar rotation harmonics).
5. **Scientific Honesty & Provenance**: All human and AI contributors must be credited accurately, and candidate signals must be described as **unconfirmed planet candidates** with clear falsification criteria.
