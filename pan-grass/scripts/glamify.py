GLAM="/usr/local/src/meme/meme_4.4.0/src/glam2"
"""
Main alphabets: p = proteins, n = nucleotides
Main options (default settings):
-h: show all options and their default settings
-o: output directory; will not clobber existing files
-O: output directory (glam2_out); allow clobbering
-r: number of alignment runs (10)
-n: end each run after this many iterations without improvement (10000)
-2: examine both strands - forward and reverse complement
-z: minimum number of sequences in the alignment (2)
-a: minimum number of aligned columns (2)
-b: maximum number of aligned columns (50)
-w: initial number of aligned columns (20)
-Q: run quietly
"""

from itertools import groupby
import collections
import sys
import subprocess
fh = open(sys.argv[1])

def sh(cmd):
    p = subprocess.Popen(cmd, shell=True)
    return p.communicate()

def pairs(iterable):
    "s -> (s0,s1), (s2,s3), ..."
    for i in range(0, len(iterable), 2):
        yield iterable[i], iterable[i + 1]


def combs(*args):
    """
    >>> combs([1], [2])
    [[1, 2]]

    >>> combs([1], [2, 3])
    [[1, 2], [1, 3]]

    >>> combs(['a'], ['b', 'c'], ['d'])
    [['a', 'b', 'd'], ['a', 'c', 'd']]

    >>> combs(['a'], ['b', 'c'], [1, 2])
    [['a', 'b', 1], ['a', 'b', 2], ['a', 'c', 1], ['a', 'c', 2]]
    """
    results = []
    if len(args) == 1:
        #print >>sys.stderr, "SHORT args[0]", tuple(args[0])
        return tuple(args[0])

    for item in args[0]:
        for other in combs(*args[1:]):
            if isinstance(other, list):
                results.append([item] + other)
            else:
                results.append([item] + [other])
    return results





orgs = ('brachy_v1', 'sorghum_v1.4', 'rice_v6')

glamin = "t.glam.fa"
for group in (list(g) for b, g in groupby(fh, lambda line: bool(line.strip())) if b):
    ########################3
    if len(group) < 3: continue
    fh = open(glamin, "w")
    print >>fh, "".join(group)
    fh.close()
    sh("%s -z 3 n %s 2>/dev/null" % (GLAM, glamin))
    res = [r.strip() for r in open("glam2_out/glam2.re").readlines()]
    gs = collections.defaultdict(list)
    for header, seq in pairs(group):
        header = header[1:].strip() # remove ">"
        org = header.split()[0]
        gs[org].append((header, seq.strip()))


    args = []
    for org in orgs:
        args.append(gs[org])

    for (bheader, bseq), (sheader, sseq), (rheader, rseq) in combs(*args): 
        line = []
        for header, seq in (bheader, bseq), (sheader, sseq), (rheader, rseq):
            line.extend(header.split() + [seq.upper()])
        line.append("\t".join(res))
        print "\t".join(line).strip()


    """
    lines = []
    for org in orgs:
        org_set = []
        for header, seq in gs[org]:
            org_set.extend(header.split() + [seq])
            n += 1

        org_set = "\t".join(org_set)
        line.append(org_set)
    if n < 4: continue
    line.append("\t".join(res))
    print "\t".join(line).strip()
    """


