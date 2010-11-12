import web
from web.contrib.template import render_mako
import urllib
import sqlite3
import numpy as np
import sys
import os
import os.path as op
from pyfasta import Fasta, complement
from commands import getoutput
from itertools import combinations

path = op.dirname(op.abspath(op.join(__file__, "..")))
sys.path.insert(0, "/opt/src/quota-alignment/scripts")
TMPDIR='/tmp/blasts'
URL = 'icns'
if not op.exists(TMPDIR):
    os.makedirs(TMPDIR)

from bed_utils import BlastLine

fastas = {
    'rice_v6': Fasta(path + '/data/rice_v6_sorghum_v1.4/rice_v6.fasta'),
    'sorghum_v1.4': Fasta(path + '/data/rice_v6_sorghum_v1.4/sorghum_v1.4.fasta'),
    'brachy_v1': Fasta(path + '/data/brachy_v1_sorghum_v1.4/brachy_v1.fasta'),
}


BL2SEQ = "/usr/bin/bl2seq -p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %(mask)s" \
           " -Y 812045000 -d 26195 -e 2.11 -i %(iseq_file)s -j %(jseq_file)s  " \
           "| grep -v '#' | grep -v 'WARNING' | grep -v 'ERROR' "


render = render_mako(directories=[op.join(op.dirname(__file__), 'templates')])

class Loc(object):
    __slots__ = ('org', 'start', 'end', 'seqid', 'rc')
    def __init__(self, org, seqid, start, end):
        self.org = org
        self.start = int(start)
        self.end = int(end)
        self.seqid = seqid
        self.rc = start > end

    @classmethod
    def from_input(self, winput):
        for loc in winput(locs=[]).locs:
            r = loc.split("..")
            start, end = map(int, r[2:4])
            fields = ('org', 'seqid', 'start', 'end', 'rc')
            yield Loc(r[0], r[1], start, end)

    def __str__(self):
        return "..".join(map(str, (self.org, self.seqid, self.start, self.end)))

    __repr__ = __str__

def get_annos(loc, db):
    gstart, gend = (loc.end, loc.start) if loc.rc else (loc.start, loc.end)
    for i, row in enumerate(db.execute("SELECT * from simple where organism = ? "
                          "AND seqid = ? AND start < ? AND end > ?",
                          (loc.org, loc.seqid, gend, gstart))):
        start, end, strand = row['start'], row['end'], row['strand']
        if loc.rc:
            strand = "+" if strand == '-' else '-'
            start, end = end, start
        yield (str(i) +":" + str(loc), str(loc), start, end, row['type'].replace(' ', '-'), strand, row['name'])

def write_fasta(loc, tmpdir=TMPDIR):
    fasta = fastas[loc.org]
    seq = fasta[loc.seqid]
    fname = "%s/%s_%i_%i.fasta" % (tmpdir, loc.seqid, loc.start, loc.end)
    if op.exists(fname):
        return fname
    start, end = sorted((loc.start, loc.end))
    seq = seq[start - 1: end]
    if loc.rc: seq = complement(seq)[::-1]

    fh = open(fname, "w")
    print >>fh, ">%s\n%s" % (fh.name, seq)
    return fh.name

def url_params():
    params = ('gc',)
    p = []
    for param in params:
        val = web.input(**{param:False})[param]
        if not val: continue
        p.append('%s=%s' % (param, val))
    return "&".join(p)


class index(object):
    def GET(self, dbname):
        dbname = dbname or "thaliana_v9"
        web.header('content-type', 'text/html')
        wi = web.input(locs=[])
        anno_url = ("/%s/%s/annos/?locs=" % (URL, dbname)) + "&locs=".join(wi.locs)
        anno_url += "&" + url_params()
        return render.index(anno_url=anno_url,
                            locs=list(Loc.from_input(web.input)))

class annos(object):
    def GET(self, dbname):
        web.header('Content-type', 'text/plain')
        gc = web.input(gc='f').gc[0] not in ('f', 'F', '0')
        locs = list(Loc.from_input(web.input))
        dbpath = op.abspath(op.join(path, "data/db/%s.db") % dbname)
        web.debug(dbpath)
        db = sqlite3.connect(dbpath)
        db.row_factory = sqlite3.Row
        for i, loc in enumerate(locs):
            yield ",".join(map(str, (loc, loc, loc.start, loc.end, "track"))) + "\n"
            for a in get_annos(loc, db):
                yield ",".join(map(str, a)) + "\n"

        for loc_pair in combinations(locs, 2):
            for qanno, sanno in run_blast(loc_pair):
                yield ",".join(map(str, qanno)) + "\n"
                yield ",".join(map(str, sanno)) + "\n"


def run_blast(loc_pair, idx=[0]):
    loci, locj = loc_pair
    iseq_file = write_fasta(loci)
    jseq_file = write_fasta(locj)
    mask = web.input(mask='F').mask

    # the start is the end if it's flipped.
    starti = loci.start - 1
    startj = locj.start - 1
    flipi = -1 if loci.rc else 1
    flipj = -1 if locj.rc else 1

    results = getoutput(BL2SEQ % locals())
    #web.debug(BL2SEQ % locals())
    # assure unique-ness of id's.
    i = idx[0]
    for line in results.split("\n"):
        b = BlastLine(line)
        strand = '+' if b.sstart < b.sstop else '-'
        qanno = (i, str(loci), starti + flipi * b.qstart, starti + flipi * b.qstop, 'HSP', strand)
        i += 1
        sanno = (i, str(locj), startj + flipj * b.sstart, startj + flipj * b.sstop, 'HSP', strand)
        i += 1
        yield qanno, sanno
    idx[0] = i

class query:
    def GET(self, name):
        xml = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "zome_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

    <Dataset name = "phytozome" interface = "default" >
        <Filter name = "gene_name_filter" value = "%s"/>
        <Attribute name = "gene_description" />
        <Attribute name = "go_id" />
        <Attribute name = "go_desc" />
    </Dataset>
</Query>
""".replace("\n", "") % name
        url = "http://www.phytozome.net/biomart/martservice"
        data = urllib.urlencode({'query': xml})
        val = urllib.urlopen(url + "?" +  data).read()
        try:
            desc, gid, godesc = val.split("\t")[:3]
            if gid:
                gid = "GO:%06d %s" % (int(gid), godesc)
            return "%s<br/>%s" % (gid, desc)
        except Exception, e:
            return "molecular function:" + str(e) + "<br/>" + val

urls = (
    '/query/(.+)', 'query',
    '/(.+)/annos/?', 'annos',
    #'/(.+)/plot/?', 'plot',
    '/(.+)?/?', 'index',
)

application = web.application(urls, globals()).wsgifunc()
