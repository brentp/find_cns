"""
load stuff from cns pipeline into proofing database
"""
import sys
import sqlite3
from flatfeature import Bed
import os.path as op
sys.path.insert(0, op.dirname(__file__))

def type_lookup(dup_info, accn):
    """given the dup_info, accn values from the cns pipeline,
    return the of gene.
    """
    dup_lookup_type = {
        'P': 'gene.parent',
        '': 'gene',
        'I': 'gene.interruptor',
    }
    type = dup_lookup_type.get(dup_info, 'gene.localdup')
    if type != 'gene': return type
    if accn.endswith('_protein'): return 'gene.protein'
    if accn.endswith('_rna'): return 'gene.rna'
    tokens = accn.split("_")
    if len(tokens) < 4: return 'gene'
    if tokens[-1].isdigit() and tokens[-2].isdigit():
        return 'gene.coanno'
    return type


def _load_bed(bed, org, comparison):
    seen = {}
    for i, accn in enumerate(bed):
        ftype = 'gene'
        if (org, accn['seqid'], accn['start'], accn['end'], ftype) in seen: continue
        seen[(org, accn['seqid'], accn['start'], accn['end'], ftype)] = True

        yield (org, accn['seqid'], accn['accn'], int(accn["start"]), int(accn["end"]), \
                    accn["strand"], ftype, comparison)

        for start, end in accn['locs']:
            if (org, accn['seqid'], start, end, 'CDS') in seen: continue
            seen[(org, accn['seqid'], start, end, 'CDS')] = True
            yield (org, accn['seqid'], accn['accn'], start, end, \
                    accn["strand"], "CDS", comparison)

def load_bed(bed, org, comparison, db):
    bed = Bed(bed)
    cur = db.cursor()
    cur.execute("PRAGMA syncronous='OFF'")
    cur.execute("PRAGMA locking_mode=EXCLUSIVE")
    cur.executemany("INSERT INTO simple "
                   "(organism, seqid, name, start, end, strand, type, group_name) "
                   "values (?, ?, ?, ?, ?, ?, ?, ?)",
                _load_bed(bed, org, comparison))
    db.commit()

def update_dups(fdups, org, group, db):
    def gen_dupcmds():
        for line in open(fdups):
            dups = line.rstrip().split()
            yield ("gene.parent", org, dups[0], group)
            for dup in dups[1:]:
                yield ("dup", org, dup, group)
    cur = db.cursor()
    cur.execute("PRAGMA syncronous='OFF'")
    cur.execute("PRAGMA locking_mode=EXCLUSIVE")
    cur.executemany("UPDATE simple set type = ? WHERE type IN ('gene', 'pseudo_gene') AND organism = ? AND name = ? AND group_name = ?",
                       gen_dupcmds())
    print >>sys.stderr, "dups updated"
    db.commit()

def _add_cns(cnsfile, org, qors, strand, group):
    # if strand is none, use the qors strand.
    fh = open(cnsfile)
    header = fh.readline().strip().split(",")
    seen = {}
    type = "CNS_" + group
    for i, line in enumerate(l.strip() for l in fh):
        d = dict(zip(header, line.split(",")))
        if not i % 5000:
            sys.stdout.write(".")
            sys.stdout.flush()
        for h in ('start', 'stop'):
            d[qors + h] = int(d[qors + h])

        name = "cns-%s-%s" % (d['qaccn'], d['saccn'])
        tup = (org, d[qors + "chr"], d[qors + "start"], d[qors + "stop"])
        if tup in seen: continue
        seen[tup] = None
        yield (org, d[qors + "chr"], name, d[qors + "start"], d[qors + "stop"], \
                strand or d[qors + "strand"], type, group)

def add_cns(db, cnsfile, org, qors, strand, comparison):
    cur = db.cursor()
    cur.execute("PRAGMA syncronous='OFF'")
    cur.execute("PRAGMA locking_mode=EXCLUSIVE")
    cur.executemany("INSERT INTO simple "
                   "(organism, seqid, name, start, end, strand, type, group_name)"
                   "values (?, ?, ?, ?, ?, ?, ?, ?)",
                _add_cns(cnsfile, org, qors, strand, comparison))
    db.commit()



def create(path):
    if not op.exists(path):
        db = sqlite3.connect(path)
        #group is the comparison e.g. brachy_v1_sorghum_v1.4
        db.executescript("""
CREATE TABLE simple (
   id INTEGER PRIMARY KEY,
   organism TEXT,
   seqid TEXT,
   name TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   type TEXT,
   group_name TEXT,
   UNIQUE(organism, seqid, strand, start, end, type, group_name)
);
    """)
        db.commit()
    else:
        db = sqlite3.connect(path)
    return db


if __name__ == "__main__":

    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option('--db', dest='db', help='path to sqlitedb to create or update')
    p.add_option('--prefix', dest='prefix', help='organism to add stuff to. e.g. rice_v6')
    p.add_option('--comparison', dest='comparison', help='comparison. e.g: rice_v6_sorghum_v1.4')
    p.add_option('--qors', dest='qors', help='cns are coming from query or subject?')
    p.add_option('--assigned-cns', dest='assigned_cns', help="path to cnss assigned to genes")
    p.add_option('--cns-strand', dest='cns_strand', help='render cnss on the plus or minus strand. used to separate out'
                    " e.g. sorghum_rice cnss from sorghum_brachy_cnss")

    opts, _ = p.parse_args()
    if not opts.db:
        print >>sys.stderr, "must specify database"
        sys.exit(p.print_help())
    db = create(opts.db)

    if not opts.prefix and opts.comparison: sys.exit(p.print_help())

    org = op.basename(opts.prefix.strip())

    for suffix in (".bed", ".localdups"):
        assert op.exists(opts.prefix + suffix), opts.prefix + suffix
    dups = opts.prefix + ".localdups"
    bed = opts.prefix + ".bed"


    load_bed(bed, org, opts.comparison, db)
    print "%s loaded" % (bed, )


    update_dups(dups, org, opts.comparison, db)

    if not opts.assigned_cns: sys.exit(0)
    add_cns(db, opts.assigned_cns, org, opts.qors, opts.cns_strand, opts.comparison)
