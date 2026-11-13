"""Shared project paths for family-level primer design."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUTDIR = os.path.join(ROOT, 'data', 'housekeeping', 'primers')
os.makedirs(OUTDIR, exist_ok=True)


RAW = os.path.join(OUTDIR, 'raw')
PROC = os.path.join(OUTDIR, 'processed')
RESULTS = os.path.join(OUTDIR, 'results')
SPEC = os.path.join(OUTDIR, 'specificity')
BLASTDB = os.path.join(OUTDIR, 'blastdb')

for directory in [RAW, PROC, RESULTS, SPEC, BLASTDB]:
    os.makedirs(directory, exist_ok=True)

print("ROOT:", ROOT)
print("RAW:", RAW)
print("PROC:", PROC)
print("RESULTS:", RESULTS)
print("OUTDIR:", OUTDIR)

# Transcriptome (symlink in external/)
TRANSCRIPTOME = os.path.join(ROOT, 'data','external', 'transcriptome.fasta')
# Filtered family candidate FASTA files live one level up
FAMILY_CANDIDATES = os.path.join(ROOT, 'data', 'housekeeping', 'family_candidates_filtered')

def _j(base, name):
    return os.path.join(base, name)

def raw(name):
    return _j(RAW, name)

def proc(name):
    os.makedirs(PROC, exist_ok=True)
    return _j(PROC, name)

def results(name):
    os.makedirs(RESULTS, exist_ok=True)
    return _j(RESULTS, name)

def spec(name):
    os.makedirs(SPEC, exist_ok=True)
    return _j(SPEC, name)

def ext(name):
    return _j(EXT, name)

def blastdb(name):
    os.makedirs(BLASTDB, exist_ok=True)
    return _j(BLASTDB, name)

def require(path, hint=None):
    if not os.path.exists(path):
        sys.exit('FATAL: missing input %s%s' % (path, '\n       ' + hint if hint else ''))
    return path
