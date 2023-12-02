# NOT ONLY 0 and 16 in $2 !!!

import sys
from collections import defaultdict
import re
import bisect
from microdict import mdict
import io
import argparse

cigar_pattern = re.compile(r'(\d+)([MIDNS])')
DNA_BASES = ('A', 'C', 'G', 'T')

umi_to_code = dict(zip([''.join([base_1, base_2, base_3, base_4])
                        for base_1 in DNA_BASES for base_2 in DNA_BASES for base_3 in DNA_BASES for base_4 in
                        DNA_BASES],
                       range(256)))
ind = iter(range(16))
for a in DNA_BASES:
    for b in DNA_BASES:
        umi_to_code[''.join([a, b])] = next(ind)


def hash_umi(umi: str) -> int:
    """
    :param umi: String with UMI; it's assumed that the actual UMI has length 12, and that the first
    5 characters in variable umi is the samtools prefix ('UB:Z:'). Therefore, the variable umi is expected
    to have a length of 5+12 = 17
    :return: int representing the hash of UMI
    """
    return umi_to_code[umi[5:9]] * 1000000 + umi_to_code[umi[9:13]] * 1000 + umi_to_code[umi[13:17]]


def cigar2pos(start_cigar):
    pos, cigar = start_cigar.split("_")
    pos = int(pos)
    if cigar == "93M":
        return ([pos], [pos + 92])
    starts = []
    ends = []
    cigar_ops = cigar_pattern.findall(cigar)
    if cigar_ops[0][1] == 'S':
        cigar_ops = cigar_ops[1:]
    for length, op in cigar_ops:
        length = int(length)
        if op == 'M':
            starts.append(pos)
            ends.append(pos + length - 1)
        elif op == 'I':
            continue
        pos += length
    return (starts, ends)


def find_intersect(starts, ends, th):
    n = len(starts)
    starts.append(ends[-1] + 1)

    result_starts = []
    result_ends = []
    i = 0
    j = 0
    overlap_count = 0

    while i <= n and j < n:
        if starts[i] <= ends[j]:
            overlap_count += 1
            if overlap_count == th:
                result_starts.append(starts[i])
            i += 1
        else:
            overlap_count -= 1
            if overlap_count == th - 1:
                result_ends.append(ends[j])
                if i == n:
                    break
            j += 1

    if result_starts:
        return result_starts, result_ends
    else:
        return


def intersect(startss_endss, th):
    starts = []
    ends = []
    for starts_ends in startss_endss:
        for start in starts_ends[0]:
            bisect.insort(starts, start)
        for end in starts_ends[1]:
            bisect.insort(ends, end)
    return find_intersect(starts, ends, th)


def proc(chrom, dup_min):
    global out
    for bc in out:
        for umi in list(out[bc].keys()):  # changed size bc of +out[bc][chrom]
            if len(out[bc][umi]) > 1:
                bad2bcumis(bc, umi)
                continue
            strand = list(out[bc][umi].keys())[0]
            if out[bc][umi][strand][0] < dup_min:
                bad2bcumis(bc, umi)
                continue
            starts_ends_dups = list()
            start_cigars = out[bc][umi][strand][1]
            start_cigars = start_cigars.split(
                ',')  # [1::] #ONLY IF out defaultdict str(), ','.joni(out, first start_cigar) => ",.." => empty first, but we switched to normal dict!
            for start_cigar in start_cigars:
                starts_ends_dup = cigar2pos(start_cigar)
                starts_ends_dups.append(starts_ends_dup)
            starts_ends_umi = intersect(starts_ends_dups, th=dup_min)  # -1
            if starts_ends_umi:
                if chrom in out[bc]:
                    out[bc][chrom][umi] = (strand, starts_ends_umi)
                else:
                    out[bc][chrom] = {umi: (strand, starts_ends_umi)}
            else:
                bad2bcumis(bc, umi)

        if chrom not in out[bc]:
            continue

        if len(out[bc][chrom]) < cov_min:
            for umi in out[bc][chrom]:
                bad2bcumis(bc, umi)
            del out[bc][chrom]
            continue

        starts_ends_umis = defaultdict(lambda: list())  # [[],[]] #
        for umi in out[bc][chrom]:
            strand, starts_ends_umi = out[bc][chrom][umi]
            starts_ends_umis[strand].append(starts_ends_umi)
        starts_ends_chrbc = {}
        for strand in starts_ends_umis:
            starts_ends_chrbc[strand] = intersect(starts_ends_umis[strand], th=cov_min)  # -1

            if not starts_ends_chrbc[strand]:
                for umi in list(out[bc][chrom].keys()):
                    if out[bc][chrom][umi][0] == strand:
                        bad2bcumis(bc, umi)
                        del out[bc][chrom][umi]
                        continue
        for umi in list(out[bc][chrom].keys()):
            strand, starts_ends_umi = out[bc][chrom][umi]
            coords = intersect([starts_ends_chrbc[strand], starts_ends_umi], th=2)
            if coords:
                out[bc][chrom][umi] = (strand, coords)
            else:
                bad2bcumis(bc, umi)
                del out[bc][chrom][umi]

        if len(out[bc][chrom]) < cov_min:
            for umi in out[bc][chrom]:
                bad2bcumis(bc, umi)
            del out[bc][chrom]
            continue

        for umi in out[bc][chrom]:
            strand, coords = out[bc][chrom][umi]
            value = ' '.join([chrom, strand, str(coords[0]), str(coords[1]), ";"])  # True #
            good2bcumis(bc, umi, value)

        del out[bc][chrom]


def add2out(bc, umi, strand, start_cigar):
    global out
    if bc in out:
        if umi in out[bc]:
            if strand in out[bc][umi]:
                out[bc][umi][strand][0] += 1
                out[bc][umi][strand][1] = ','.join([out[bc][umi][strand][1], start_cigar])
            else:
                bad2bcumis(bc, umi)
        else:
            out[bc][umi] = {strand: [1, start_cigar]}
    else:
        out[bc] = {umi: {strand: [1, start_cigar]}}


def bad2bcumis(bc, umi):
    global bcumis
    bcumis[bc][umi] = 0


def good2bcumis(bc, umi, value):
    global bcumis, keeper, keeper_ind
    bcumis[bc][umi] = keeper_ind
    keeper.write(value)
    keeper_ind += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filt_bcs')
    parser.add_argument('out_file')
    parser.add_argument('-dup_min', type=int, default=5)
    parser.add_argument('-cov_min', type=int, default=5)
    args = parser.parse_args()

    cov_min = args.cov_min

    out = {}

    bcumis = {}
    keeper = io.StringIO()
    keeper.write(';')
    keeper_ind = 1
    with open(args.filt_bcs) as file:
        for line in file.readlines():
            barcode = line.strip()
            bcumis[''.join(['CB:Z:', barcode])] = mdict.create("i32:i32")

    chromosome_old = ''

    for line in sys.stdin:
        line = line.strip().split()
        if line[-2] not in bcumis or line[2] == "chrM":
            continue

        barcode = line[-2]
        umi = line[-1]
        hashed_umi = hash_umi(umi)

        if hashed_umi in bcumis[barcode]:
            bad2bcumis(barcode, hashed_umi)
            continue

        strand = line[1][0]
        if (barcode in out) and (hashed_umi in out[barcode]) and (strand not in out[barcode][hashed_umi]):
            bad2bcumis(barcode, hashed_umi)
            continue

        chrom = line[2]
        start_cigar = '_'.join([line[3], line[5]])

        if chrom != chromosome_old:
            if chromosome_old:
                proc(chromosome_old, args.dup_min)
                out = {}
            chromosome_old = chrom
        add2out(barcode, hashed_umi, strand, start_cigar)

    proc(chrom, args.dup_min)

    keeper.seek(0)
    result = keeper.read().split(';')
    with open(args.out_file, 'w') as f:
        for barcode in bcumis:
            for hashed_umi in bcumis[barcode]:
                if bcumis[barcode][hashed_umi]:
                    print(barcode, hashed_umi, result[bcumis[barcode][hashed_umi]], file=f)
    exit()
