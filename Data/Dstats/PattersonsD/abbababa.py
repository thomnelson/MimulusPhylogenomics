#!/usr/bin/env python
#
# PARSE A VCF AND KEEP A SUBSET OF REGIONS/INDIVIDUALS
#

import sys
newpath = '/home/thom_nelson/modules/'
if newpath not in sys.path:
    sys.path.append(newpath)
from thomtools import readVCFregion
import argparse
import subprocess

# parse command line options
parser = argparse.ArgumentParser(description='Enumerate ABBA/BABA sites from a series of tabix-indexed VCFs and a specified tree of ((P1,P2),P3),P0')

parser.add_argument('-c','--chroms', required=True,nargs='+', help='list of chromosomes to process.')
parser.add_argument('-p','--prefix', required=True, help='prefix to VCF files (e.g. prefix.[chromosomeID].suffix).')
parser.add_argument('-s','--suffix', required=True, help='suffix to VCF files (e.g. prefix.[chromosomeID].suffix).')
parser.add_argument('-0','--P0', required=True, help='Sample representing P0, the outgroup with which to polarize SNPs.')
parser.add_argument('-1','--P1', required=True, help='Sample representing P1.')
parser.add_argument('-2','--P2', required=True, help='Sample representing P2.')
parser.add_argument('-3','--P3', required=True, help='Sample representing P3.')
# parser.add_argument('-o','--out', required=False, default="filtered.vcf.gz", help='compressed file to write output.')
args=parser.parse_args()

chroms= args.chroms
prefix= args.prefix
suffix= args.suffix
P1    = args.P1
P2    = args.P2
P3    = args.P3
P0    = args.P0
taxa  = [P0,P1,P2,P3]
###
### EXTRACT VCF HEADER FOR LATER EXPORT AND TO INDEX COLS2KEEP IF SAMPLE SUBSET SPECIFIED
###

cmd_line = ("tabix","-H",prefix+chroms[0]+suffix)
header   = subprocess.Popen(cmd_line, stdout = subprocess.PIPE)
header   = header.communicate()[0].decode("utf-8")
p0    = 0
p1    = 0
p2    = 0
p3    = 0
for line in header.split('\n'):
    if "#CHROM" not in line:
        continue
    # parse field names
    line = line.split('\t')
    ncols = len(line)
    ### GET GT INDEX FOR EACH SAMPLE
    for i in range(9,ncols):
        test = line[i]
        if test == P0:
            p0 = i - 9
        if test == P1:
            p1 = i - 9
        if test == P2:
            p2 = i - 9
        if test == P3:
            p3 = i - 9
### NEED A STDERR MESSAGE FOR POP TREE IN NEWICK
sys.stderr.write("Analyzing genotype information for samples:\n")
for i in taxa:
    sys.stderr.write("  %s\n"%(i))
sys.stderr.write("Assuming population tree:\n  (((%s,%s),%s),%s)\n"%(P1,P2,P3,P0))


###
nentries = 0
sys.stderr.write("\n")
ABBA = 0
BABA = 0
BBAA = 0
### SAVE CHR/POS/CLASSIFICATION TO A LIST FOR EXPORT
positional = []
for chrom in chroms:
    # Extract entries from VCF
    chrentries = 0
    vcf = prefix+chrom+suffix
    vcfSNPs = readVCFregion(vcf, chrom, '0', '50000000', True, True)
    nSNPs = len(vcfSNPs)
    chrABBA = 0
    chrBABA = 0
    chrBBAA = 0
    for SNP in vcfSNPs:
        CHROM = SNP['CHROM']
        POS   = SNP['POS']
        SNPid = CHROM+POS
        ID    = SNP['ID']
        REF   = SNP['REF']
        ALT   = SNP['ALT']
        QUAL  = SNP['QUAL']
        FILTER= SNP['FILTER']
        INFO  = SNP['INFO']
        FORMAT= SNP['FORMAT']
        GTs   = SNP['GTs']
        ### CHECK IF HETS OR MISSING, PASS FOR NOW
        if '0/1' in GTs or './.' in GTs:
            continue
        A = GTs[p0]
        nentries +=1
        chrentries += 1
        siteBBAA = 0
        siteABBA = 0
        siteBABA = 0
        ### CHECK IF SITE MATCHES POP TREE
        if GTs[p1] != A and GTs[p2] != A and GTs[p3] == A:
            chrBBAA += 1
            BBAA    += 1
            siteBBAA+= 1
        ### CHECK IF AN ABBA SITE
        if GTs[p1] == A and GTs[p2] != A and GTs[p3] != A:
            chrABBA += 1
            ABBA    += 1
            siteABBA+= 1
        ### CHECK IF A BABA SITE
        if GTs[p1] != A and GTs[p2] == A and GTs[p3] !=A:
            chrBABA += 1
            BABA    += 1
            siteBABA+= 1
        if nentries % 100 == 0:
            sys.stderr.write("Chromosome %s: %s SNPs analyzed...\r"%(chrom,chrentries))
        positional.append([CHROM,POS,siteBBAA,siteABBA,siteBABA])
    patsD = (chrABBA - chrBABA) / (chrABBA+chrBABA)
    patsD = round(patsD, len(str(chrentries)))
    sys.stderr.write("\n  BBAA sites: %s\n"%(chrBBAA))
    sys.stderr.write("  ABBA sites: %s\n"%(chrABBA))
    sys.stderr.write("  BABA sites: %s\n"%(chrBABA))
    sys.stderr.write("  Patterson's D: %s\n\n"%(patsD))
sys.stderr.write("GRAND TOTAL:")
patsD = (ABBA - BABA) / (ABBA+BABA)
patsD = round(patsD, len(str(nentries)))
sys.stderr.write("\n  BBAA sites: %s\n"%(BBAA))
sys.stderr.write("  ABBA sites: %s\n"%(ABBA))
sys.stderr.write("  BABA sites: %s\n"%(BABA))
sys.stderr.write("  Patterson's D: %s\n"%(patsD))

### WRITE POSITIONAL DATA TO STDOUT

sys.stderr.write("\nWriting positional ABBA-BABA classification to stdout...\n")
sys.stdout.write("CHR\tPOS\tBBAA\tABBA\tBABA\n")
for site in positional:
    sys.stdout.write("%s\t%s\t%s\t%s\t%s\n"%(site[0],site[1],site[2],site[3],site[4]))
