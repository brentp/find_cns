import sys
import os
import os.path as op
import numpy as np
import commands
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed

NCPU = 8
from processing import Pool
pool = Pool(NCPU)


EXPON = 0.90

def get_feats_in_space(locs, ichr, bpmin, bpmax, bed):
    """ locs == [start, stop]
    bpmin is the lower extent of the window, bpmax ... 
    """
    assert bpmin < bpmax, (locs, ichr, bpmin, bpmax)
    feats = bed.get_features_in_region(str(ichr), bpmin, bpmax)
    feats = [f for f in feats if not (f['start'] == locs[0] and f['end'] == locs[1])]
    if len(feats) != 0:
        assert feats[0]['seqid'] == str(ichr)
    return [(f['start'], f['end'], f['accn']) for f in feats]
    

def parse_blast(blast_str, orient, qfeat, sfeat, qbed, sbed, pad):
    blast = []
    slope = orient

    qgene = [qfeat['start'], qfeat['end']]
    sgene = [sfeat['start'], sfeat['end']]
    qcds = qfeat['locs']
    scds = sfeat['locs']


    sgene = sgene[::slope]
    center = sum(qgene)/2., sum(sgene)/2.

    EXP = EXPON
    if abs(abs(qgene[1] - qgene[0]) - abs(sgene[1] - sgene[0])) > 3000:
        EXP = 0.94


    #intercept = (sgene[0] + sgene[1])/2.  - slope * (qgene[0] + qgene[1])/2.
    intercept = center[1] - slope * center[0]
    rngx = qgene[1] - qgene[0]
    rngy = abs(sgene[1] - sgene[0])

    x = np.linspace(qgene[0] - pad, qgene[1] + pad, 50)
    y = slope * x + intercept


    xb = x + -slope * rngx/3. + -slope * np.abs(x - center[0])**EXP
    yb = y + rngy/3. + np.abs(y - center[1])**EXP

    xy = x + slope * rngx/3. + slope * np.abs(x - center[0])**EXP
    yy = y - rngy/3. - np.abs(y - center[1])**EXP

    if slope == 1:
        xall = np.hstack((xy[::-1], xb[::slope], xy[-1]))
        yall = np.hstack((yy[::-1],yb, yy[-1]))
    if slope == -1:
        xall = np.hstack((xy, xb[::-1], xy[0]))
        yall = np.hstack((yy,yb[::-1], yy[0]))

    feats_nearby = {}
    feats_nearby['q'] = get_feats_in_space(qgene, qfeat['seqid'], int(x.min()), int(x.max()), qbed)
    feats_nearby['s'] = get_feats_in_space(sgene, sfeat['seqid'], int(y.min()), int(y.max()), sbed)



    genespace_poly = Polygon(zip(xall, yall))
    
    for sub in ('q', 's'):
        if len(feats_nearby[sub]) !=0:
            feats_nearby[sub] = MultiLineString([[(0, c0),(0, c1)] for c0, c1, fname in feats_nearby[sub]])
        else:
            feats_nearby[sub] = None

    cnss = set([])

    qgene_poly = LineString([(0.0, qgene[0]), (0.0, qgene[1])])
    sgene_poly = LineString([(0.0, sgene[0]), (0.0, sgene[1])])
    intronic_removed = 0

    for line in blast_str.split("\n"):
        line = line.split("\t")
        locs = map(int, line[6:10])
        locs.extend(map(float, line[10:]))

        xx = locs[:2]
        yy = locs[2:4]
        
        # get rid of stuff on the wrong strand
        if slope == 1 and locs[2] > locs[3]: continue
        if slope == -1 and locs[2] < locs[3]: continue

        # to be saved. a hit must either be in an intron in both
        # genes, or in neither.

        ##########################################################
        # DEAL WITH INTRONIC cnss in the gene of interest.
        ##########################################################
        xls = LineString([(0, locs[0]), (0, locs[1])])
        yls = LineString([(0, locs[2]), (0, locs[3])])

        locs = tuple(locs) # make it hashable.
        if qgene_poly.intersects(xls) and sgene_poly.intersects(yls): 
            cnss.update((locs,))
            continue

        # has to be both or neither.
        if qgene_poly.intersects(xls) or sgene_poly.intersects(yls): 
            intronic_removed += 1
            continue

        ##########################################################

        ###############################################################
        # for all other genes, if it's in an intron, we dont keep it.
        ###############################################################
        intronic = False 
        # get rid of stuff that overlaps another gene:
        for sub, (start, stop) in (('q', locs[:2]), ('s', locs[2:4])):
            feats = feats_nearby[sub]
            if feats is None: continue
            # the current hsp is overlapping another gene. we dont want that...
            if feats.contains(Point(0, start)) or feats.contains(Point(0, stop)):
                intronic = True
                break

        if intronic: continue

        ##########################################################

        # this is the bowtie.
        if not genespace_poly.contains(LineString(zip(xx, yy))): continue
        cnss.update((locs,))

    # cant cross with < 2 cnss.
    # get rid of the eval, bitscore stuff.
    if len(cnss) < 2: return [l[:4] for l in cnss]

    cnss = list(cnss)
    # need to flip to negative so the overlapping stuff still works.
    if orient == -1:
        cnss = list(cnss)
        for i, cns in enumerate(cnss):
            cns = list(cns)
            cns[2] *= - 1
            cns[3] *= - 1
            cnss[i] = tuple(cns)
        sgene[0] *= -1
        sgene[1] *= -1

    cnss = [l[:4] for l in remove_crossing_cnss(cnss, qgene, sgene)]
    if orient == -1:
        cnss = [(c[0], c[1], -c[2], -c[3]) for c in cnss]
    return cnss


def remove_overlapping_cnss(cnss):
    """for cases when there is nearly the same cns, but with 1
    basepair shfit up/down. that create many cnss stacked on top
    of each other. this reduces those down to one."""
    qcnss = [LineString([(0, cns[0]), (0, cns[1])]) for i, cns in enumerate(cnss)]
    scnss = [LineString([(0, cns[2]), (0, cns[3])]) for i, cns in enumerate(cnss)]

    remove = []
    for zcnss in (qcnss, scnss):
        for i, csi in enumerate(zcnss[:-1]):
            for _j, csj in enumerate(zcnss[i + 1:]):
                j = i + _j + 1 # cause enumerate starts at 0
                if csi.overlaps(csj):
                    if cnss[i][-2] < cnss[j][-2] or cnss[i][-1] > cnss[j][-2] or csi.y < csj.y:
                        remove.append(j)
                    else:
                        remove.append(i)
    remove = frozenset(remove)
    return [cns for i, cns in enumerate(cnss) if not i in remove]


def remove_crossing_cnss(cnss, qgene, sgene):
    diff = qgene[0] - sgene[0] # adjust subject so it's in same range as query
    cns_shapes = [LineString([((c[0] + c[1])/2., 0 ), ((c[2] + c[3])/2. + diff, 1000)]) for c in cnss]

    overlapping = len(cnss)
    cnss = remove_overlapping_cnss(cnss)
    overlapping -= len(cnss)
    cns_shapes = [LineString([((c[0] + c[1])/2., 0 ), ((c[2] + c[3])/2. + diff, 1000)]) for c in cnss]


    # and save a reference to the orginal cnss as that's the data we want.
    for cs, cns in zip(cns_shapes, cnss):
        cs.cns = cns
        # hold the number of times an hsp crosses any other.
        cs.cross_list = set([])
        # mark for removal.
        cs.do_remove = False


    for csi in cns_shapes:
        for csj in cns_shapes: 
            if csi == csj: continue
            if csi.crosses(csj):
                csi.cross_list.update([csj])
                csj.cross_list.update([csi])

    ###################################################################### 
    # first remove anything that cross more than 5 other cnss. 
    ###################################################################### 
    nremoved = 0
    any_removed = True
    while any_removed:
        # need this outer loop to refresh the sorting.
        cns_shapes = sorted(cns_shapes, reverse=True, cmp=lambda a, b: cmp(len(a.cross_list), len(b.cross_list)))[:]
        any_removed = False
        for i, cs in enumerate(cns_shapes):
            if len(cs.cross_list) > 3:
                # remove this from all other lists as it's a bad guy.
                for crossed in cs.cross_list:
                    crossed.cross_list.difference_update(set([cs]))
                cs.do_remove = True
                any_removed = True
                nremoved += 1
                del cns_shapes[i]
                break


    ###################################################################### 
    # then remove crosses one-by-one, keeping the < evalue, > bitscore.
    ###################################################################### 
    for csi in cns_shapes:
        if csi.do_remove or len(csi.cross_list) == 0: continue
        for csj in cns_shapes:
            if csj.do_remove or len(csj.cross_list) == 0: continue
            if csi.do_remove or len(csi.cross_list) == 0: continue
            if csi.crosses(csj):

                # access the assocated cns.
                # evalue: less is better       bitscore: more is better
                if csi.cns[-2] < csj.cns[-2] or csi.cns[-1] > csj.cns[-1]:
                    csj.do_remove = 1
                    map(lambda crossed: crossed.cross_list.difference_update(set([csj])), csj.cross_list)

                else:
                    csi.do_remove = 1
                    map(lambda crossed: crossed.cross_list.difference_update(set([csi])), csi.cross_list)
                    break

    for c in cns_shapes:
        if not c.do_remove: continue
        nremoved += 1
    return [c.cns for c in cns_shapes if not c.do_remove]
        

def get_pair(pair_file, fmt, seen={}):
    """ read a line and make sure it's unique handles
    dag, cluster, and pair formats."""
    fh = open(pair_file)
    for line in open(pair_file):
        if line[0] == "#": seen.clear(); continue
        line = line.split("\t")
        if fmt == 'dag':
            assert len(line) > 5, line
            key = tuple(line[:-1])
            pair = line[1], line[5]
        elif fmt in ('cluster', 'qa', 'raw'):
            assert len(line) == 5, line
            key = tuple(line[:-1])
            pair = line[1], line[3]
        elif fmt == 'pairs':
            if len(line) == 1:
                line = line.split(",")
            assert len(line) >= 2, "dont know how to handle %s" % line
            pair = line[0], line[1]
            key = tuple(line)

        if fmt in ('qa', 'raw'):
            pair = int(pair[0]), int(pair[1])

        if key in seen:
            continue
        else:
            seen[tuple(key)] = True
            if isinstance(pair[0], (int, long)) and \
                           isinstance(pair[1], (int, long)):
                   yield int(pair[0]), int(pair[1])

            assert len(pair) == 2, (pair, line)
            yield pair


def get_masked_fastas(bed):
    """
    create the masked fasta files per chromosome. needed to run bl2seq.
    """
    f = bed.fasta.fasta_name
    fname = op.splitext(op.basename(f))[0]
    d = op.dirname(f) + "/%s_split" % fname
    try: os.mkdir(d)
    except OSError: pass

    fastas = {}
    for seqid, seq in bed.mask_cds():
        f = d + "/%s.fasta" % seqid
        fastas[seqid] = f
        if op.exists(f): continue
        fh = open(f, "wb")
        print >>fh, ">%s" % seqid
        print >>fh, seq.tostring()
        fh.close()
    return fastas

def main(qbed, sbed, pairs_file, pad, pair_fmt):
    """main runner for finding cnss"""

     
    bl2seq = "/usr/bin/bl2seq -p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F F " \
           " -Y 812045000 -d 26195 -e 2.11 -i %(qfasta)s -j %(sfasta)s \
              -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
            | grep -v 'WARNING' | grep -v 'ERROR' "

    fcnss = sys.stdout
    print >> fcnss, "#qseqid,qaccn,sseqid,saccn,[qstart,qend,sstart,send...]"

    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed)

    pairs = [True]
    _get_pair_gen = get_pair(pairs_file, pair_fmt)
    # need this for parallization stuff.
    def get_pair_gen():
        try: return _get_pair_gen.next()
        except StopIteration: return None

    while any(pairs):
        pairs = [get_pair_gen() for i in range(NCPU)]

        # this helps in parallelizing.
        def get_cmd(pair):
            if pair is None: return None
            if pair_fmt == 'dag':
                qfeat, sfeat = qbed.d[pair[0]], sbed.d[pair[1]]
            elif pair_fmt in ('cluster', 'raw', 'qa'):
                # qseqid, qidx, sseqid, sidx
                qfeat, sfeat = qbed[pair[0]], sbed[pair[1]]
            elif pair_fmt == 'pair':
                qq, ss = pair
                if isinstance(qq, int) and isinstance(ss, int):
                    qfeat, sfeat = qbed[qq], sbed[ss]
                else:
                    qfeat, sfeat = qbed.d[qq], sbed.d[ss]

            qfasta = qfastas[qfeat['seqid']]
            sfasta = sfastas[sfeat['seqid']]

            qstart, qstop = max(qfeat['start'] - pad, 1), qfeat['end'] + pad
            sstart, sstop = max(sfeat['start'] - pad, 1), sfeat['end'] + pad

            assert qstop - qstart > 2 * pad or qstart == 1, (qstop, qstart)
            assert sstop - sstart > 2 * pad or sstart == 1, (sstop, sstart)

            cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta, qstart=qstart,
                                sstart=sstart, qstop=qstop, sstop=sstop)
            return cmd, qfeat, sfeat

        cmds = [c for c in map(get_cmd, [l for l in pairs if l]) if c]
        results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))
        #results = (r for r in map(commands.getoutput, [c[0] for c in cmds]))

        for res, (cmd, qfeat, sfeat) in zip(results, cmds):
            if not res.strip(): continue
            print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn']),
            orient = qfeat['strand'] == sfeat['strand'] and 1 or -1

            cnss = parse_blast(res, orient, qfeat, sfeat, qbed, sbed, pad)
            print >>sys.stderr, "(%i)" % len(cnss)
            if len(cnss) == 0: continue

            qname, sname = qfeat['accn'], sfeat['accn']
            print >> fcnss, "%s,%s,%s,%s,%s" % (qfeat['seqid'], qname, sfeat['seqid'], sname,
                             ",".join(map(lambda l: ",".join(map(str,l)),cnss)))

    return None 

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-q", dest="qfasta", help="path to genomic query fasta")
    parser.add_option("--qbed", dest="qbed", help="query bed file")
    parser.add_option("-s", dest="sfasta", help="path to genomic subject fasta")
    parser.add_option("--sbed", dest="sbed", help="subject bed file")
    parser.add_option("-p", dest="pairs", help="the pairs file. output from dagchainer")
    choices = ("dag", "cluster", "pair", 'qa', 'raw')
    parser.add_option("--pair_fmt", dest="pair_fmt", default='raw',
                      help="format of the pairs, one of: %s" % str(choices),
                      choices=choices)
    parser.add_option("--pad", dest="pad", type='int', default=12000,
                      help="how far from the end of each gene to look for cnss")
    (options, _) = parser.parse_args()


    if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
        sys.exit(parser.print_help())
    
    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()

    main(qbed, sbed, options.pairs, options.pad, options.pair_fmt)
