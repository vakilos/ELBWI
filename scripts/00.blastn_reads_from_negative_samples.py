from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re,sys,os, argparse
import gzip
pd.set_option('display.max_rows',None)

parser=argparse.ArgumentParser(description="This script will merge reads (forward and reverse orientation) for read 1 and read 2 for each sample. Run in an python environment with the latest version of blast+ installed")
parser.add_argument("-i", dest="inDir", help="Directory with all fastqs from SRA. (use full path)", required=True)
parser.add_argument('-o', dest="outDir", help="Directory to export files. (use full path)", required=True, default='')
args=parser.parse_args()

inDir = args.inDir
outDir = args.outDir
if outDir.endswith("/") is False:
    outDir += "/"
files = sorted([x for x in os.listdir(inDir) if x.startswith('neg')])
os.chdir(inDir)
fasta = []

os.chdir(inDir)
if os.path.isfile(outDir+"negative_control_reads.fasta") is False:
    for file in files:
        name = re.match(re.compile(r'(neg\d).+'),file).group(1)
        idx = 0
        with gzip.open(file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                idx += 1
                newRec = SeqRecord(seq=record.seq,id=name+"_"+str(idx), description="")
                fasta.append(newRec)

    with open(outDir+"negative_control_reads.fasta", "w") as out_handle:
        SeqIO.write(fasta, out_handle, "fasta")

os.chdir(outDir)
if os.path.isfile(outDir+"neg_control_hits.csv") is False:
    from Bio.Blast import NCBIWWW
    from Bio.Blast.Applications import NcbiblastnCommandline as blastn
    from Bio.Blast import NCBIXML
    blastn_cm = blastn(query='negative_control_reads.fasta',db='16S_ribosomal_RNA', evalue=0.001,outfmt=5,out='nc_reads_blast.xml')
    stout,sterr = blastn_cm()
    result_handle = open('nc_reads_blast.xml')
    blast_records = NCBIXML.parse(result_handle)
    hits = pd.DataFrame()
    for record in blast_records:
        if len(record.alignments) > 0:
            best_alignment=record.alignments[0]
            print(record.query)
        else:
            print(record.query+" has NO ALIGNMENT")
        for hsp in best_alignment.hsps:
            #if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence:", best_alignment.title)
            print("length:", best_alignment.length)
            print("e value:", hsp.expect)
            print("score:", hsp.score)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
            print('\n')
            hit = re.match(re.compile(r'gi\|\d+\|ref\|[^\|]+\|\s([^\s]+)\s[^\s].+16S ribosomal RNA, (partial|complete) sequence'),best_alignment.title).group(1)
            tf = pd.DataFrame({'sample':[record.query.split("_")[0]], 'hit':[hit]} )
            hits = hits.append(tf)

    hits.to_csv("neg_control_hits.csv",index=False)

#tax = pd.read_csv("qiime2_dada2/taxonomy_table.tsv", sep='\t')
hf = pd.read_csv("neg_control_hits.csv")
#hits = pd.Series([re.match(re.compile(r'([^\s]+).*'),x).group(1) for x in hf.hit])
for n,g in hf.groupby("sample"):
    print(n)
    print(g.hit.value_counts())
    print('\n')

print(hf.hit.value_counts())
