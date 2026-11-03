"""Shared project paths for the Cornus DQD/SDH primer pipeline.

Every step script imports this instead of hard-coding file names, so the scripts run
from any working directory and the repository layout has a single source of truth.

    data/raw/         as received (curated cDNA, GenBank records, alignment tables)
    data/processed/   derived tables (UTR/CDS extracts, k-mer scans, coordinate tables)
    data/results/     deliverables (final primers, backup pair) + specificity/ evidence
    data/external/    large inputs, NOT committed: transcriptome FASTA/GFF, BLAST DBs
                      (see data/external/README.md for what to place there)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, 'data')

RAW = os.path.join(DATA, 'raw')
PROC = os.path.join(DATA, 'processed')
RESULTS = os.path.join(DATA, 'results')
SPEC = os.path.join(RESULTS, 'specificity')
EXT = os.path.join(DATA, 'external')
BLASTDB = os.path.join(EXT, 'blastdb')


def _j(base, name):
    return os.path.join(base, name)


def raw(name):
    return _j(RAW, name)


def align(name):
    return _j(os.path.join(RAW, 'alignment'), name)


def proc(name):
    return _j(PROC, name)


def results(name):
    return _j(RESULTS, name)


def spec(name):
    """Evidence file. Creating the folder is safe; scripts write their reports here."""
    os.makedirs(SPEC, exist_ok=True)
    return _j(SPEC, name)


def ext(name):
    return _j(EXT, name)


def blastdb(name):
    os.makedirs(BLASTDB, exist_ok=True)
    return _j(BLASTDB, name)


def require(path, hint=None):
    """Fail early with a message that says where the file should come from."""
    if not os.path.exists(path):
        sys.exit('FATAL: missing input %s%s' % (path, '\n       ' + hint if hint else ''))
    return path


def genbank_records(prefix='PZ5186', folder=RAW):
    return sorted(_j(folder, f) for f in os.listdir(folder)
                  if f.startswith(prefix) and f.endswith('.gb'))
