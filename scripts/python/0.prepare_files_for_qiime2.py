import sys, argparse, re, subprocess, os
import pandas as pd
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

parser=argparse.ArgumentParser(description="This script will merge reads (forward and reverse orientation) for read 1 and read 2 for each sample.")
parser.add_argument("-i", dest="inDir", help="Directory with all fastqs from SRA. (use full path)", required=True)
parser.add_argument('-o', dest="outDir", help="Export directory for the concatenated fastq files. (use full path)", required=True, default='')
args=parser.parse_args()

"""This script will merge reads (forward and reverse orientation) for read 1 and read 2 for each sample."""
pd.set_option('display.expand_frame_repr', False)
df = pd.DataFrame()


in_dir = args.inDir
out_dir  = args.outDir

os.chdir(in_dir)

def merge_samples(in_dir, out_dir):
    ### will merge forward and reverse orientation reads, for the first and the second read respectively
    if os.path.isdir(out_dir) is False:
        os.mkdir(out_dir)
    def extract_samples(in_dir):
        os.chdir(in_dir)
        samples = [re.match(re.compile(r'(\w+)\.\w\.\d\.fastq\.gz'), file).group(1) for file in os.listdir(in_dir) if
                   file.endswith('fastq.gz')]

        return list(set(samples))

    if out_dir[-1] !="/":
        out_dir +='/'
    import glob
    for sample in sorted(extract_samples(in_dir)):
        os.chdir(in_dir)
        R1 = sorted(glob.glob(sample+'*.1.*'))
        R2 = sorted(glob.glob(sample + '*.2.*'))

        if len(R1) == 2:
            subprocess.run('cat {} {} > '.format(*R1) +out_dir+sample+'_R1.fastq.gz', shell=True)
            print(('Concatenating '+color.BOLD+color.BLUE+'{}'+color.END+' and '+color.BOLD+color.BLUE+'{}'+color.END+' --> '+color.BOLD+color.BLUE+'{}'+color.END).format(*R1,sample+'_R1.fastq.gz'))
        else:
            subprocess.run('cat {} > '.format(*R1) +out_dir+sample+'_R1.fastq.gz', shell=True)
            print("Expected 2 files, but {} found!!".format(len(R1)))
            print(('Concatenating '+color.BOLD+color.RED+'{}'+color.END+' --> '+color.BOLD+color.RED+'{}'+color.END).format(*R1,sample+'_R1.fastq.gz'))

        if len(R2) == 2:
            subprocess.run('cat {} {} > '.format(*R2) +out_dir+sample+'_R2.fastq.gz', shell=True)
            print(('Concatenating '+color.BOLD+color.BLUE+'{}'+color.END+' and '+color.BOLD+color.BLUE+'{}'+color.END+' --> '+color.BOLD+color.BLUE+'{}'+color.END).format(*R2,sample+'_R2.fastq.gz'))

        else:
            subprocess.run('cat {} > '.format(*R2) +out_dir+sample+'_R2.fastq.gz', shell=True)
            print("Expected 2 files, but {} found!!".format(len(R2)))
            print(('Concatenating '+color.BOLD+color.RED+'{}'+color.END+' --> '+color.BOLD+color.RED+'{}'+color.END).format(*R2,sample+'_R2.fastq.gz'))

merge_samples(in_dir,out_dir)
