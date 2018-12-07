#!/usr/bin/env python

import sys
import os
import argparse
import datetime
import platform
import subprocess
import shutil
from Bio import SeqIO
from natsort import natsorted

__version__ = "0.0.2"
#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='redmask.py',
    description = '''Wraper for Red - repeat identification and masking for genome annotation''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--genome', required=True, help='genome assembly FASTA format')
parser.add_argument('-o', '--output', required=True, help='Output basename')
parser.add_argument('-m', '--min', default=3, type=int, help='Minimum number of observed k-mers')
parser.add_argument('--training', default=1000, type=int, help='Min length for training')
parser.add_argument('-l', '--word_len', type=int, help='word length (kmer length)')
parser.add_argument('-t', '--threshold', type=int, help='threshold of low adjusted scores of non-repeats')
parser.add_argument('-g', '--gaussian', type=int, help='Gaussian smoothing width')
parser.add_argument('-c', '--markov_order', type=int, help='Order of background markov chain')
parser.add_argument('--debug', action='store_true', help='Keep intermediate files')
parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
args=parser.parse_args()

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None
    
#via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def group(L):
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

#via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def list2groups(L):
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

def maskingstats2bed(input, counter):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    masked = []
    maskedSize = 0
    bedfilename = input.replace('.msk', '.bed')
    with open(input, 'rU') as infile:
        for header, Seq in SimpleFastaParser(infile):
            if ' ' in header:
                ID = header.split(' ')[0]
            else:
                ID = header
            for i,c in enumerate(Seq):
                if c.islower():
                    masked.append(i) #0 based
                    maskedSize += 1
    if maskedSize > 0: #not softmasked, return False
        with open(bedfilename, 'w') as bedout:
            repeats = list(list2groups(masked))
            for item in repeats:
                if len(item) == 2:
                    bedout.write('{:}\t{:}\t{:}\tRepeat_{:}\n'.format(ID, item[0], item[1], counter))
                    counter += 1
    return maskedSize, counter

def n50(input):
    lengths = []
    with open(input, 'rU') as infile:
        for rec in SeqIO.parse(infile, 'fasta'):
            lengths.append(len(rec.seq))
    lengths.sort()
    nlist = []
    for x in lengths:
        nlist += [x]*x
    if len(nlist) % 2 == 0:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos-1]) / 2)
    else:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
    return N50, len(lengths), sum(lengths)

def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return
        
def softwrap(string, every=80):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

sys.stdout.write('[{:}] Running Python v{:} \n'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), platform.python_version()))
dependencies = ['Red']
for x in dependencies:
    if not which_path(x):
        print('{:} is not properly installed, install and re-run script'.format(x))
        sys.exit(1)

pid = os.getpid()
inputDir = 'redmask_contigs_'+str(pid)
trainDir = 'redmask_train_'+str(pid)
outputDir = 'redmask_output_'+str(pid)
logfile = 'redmask_'+str(pid)+'.log'
if os.path.isfile(logfile):
    os.remove(logfile)
os.makedirs(inputDir)
os.makedirs(outputDir)
os.makedirs(trainDir)

Scaffolds = {}
with open(logfile, 'w') as log:
    calcN50, numContigs, lenContigs = n50(args.genome)
    sys.stdout.write('[{:}] Loading assembly:\n    Contigs: {:,}\n    Length:  {:,} bp\n    nN50:    {:,} bp\n'.format(
    					datetime.datetime.now().strftime('%b %d %I:%M %p'), numContigs, lenContigs, calcN50))
    sys.stdout.write('[{:}] Splitting genome assembly into training set (contigs > {:} bp)\n'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), args.training))
    with open(args.genome, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if not rec.id in Scaffolds:
                Scaffolds[rec.id] = len(rec.seq)
            if len(rec.seq) >= args.training:
                with open(os.path.join(trainDir,rec.id+'.fa'), 'w') as output:
                    SeqIO.write(rec, output, 'fasta')
            else:
                with open(os.path.join(inputDir,rec.id+'.fa'), 'w') as output:
                    SeqIO.write(rec, output, 'fasta')
    sys.stdout.write('[{:}] Finding repeats with Red (REpeat Detector)\n'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    cmd = ['Red', '-gnm', trainDir, '-dir', inputDir, '-min', str(args.min), '-msk', outputDir]
    if args.gaussian:
        cmd = cmd + ['-gau', str(args.gaussian)]
    if args.threshold:
        cmd = cmd + ['-thr', str(args.threshold)]
    if args.word_len:
        cmd = cmd + ['-len', str(args.word_len)]
    if args.markov_order:
        cmd = cmd + ['-ord', str(args.markov_order)]
    subprocess.call(cmd, stdout=log, stderr=log)
    
    sys.stdout.write('[{:}] Collecting results from Red\n'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    maskedOut = args.output+'.softmasked.fa'
    maskedfiles = []
    for file in os.listdir(outputDir):
        if file.endswith('.msk'):
            maskedfiles.append(os.path.join(outputDir, file))
    
    sys.stdout.write('[{:}] Summarizing results and converting to BED format\n'.format(datetime.datetime.now().strftime('%b %d %I:%M %p')))
    maskedBED = args.output+'.repeats.bed'
    maskedFA = args.output+'.repeats.fasta'
    maskedSize = 0
    num = 1
    with open(maskedOut, 'w') as outfile:
        for fname in natsorted(maskedfiles):
            masksize, num = maskingstats2bed(fname, num)
            maskedSize += masksize
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    with open(maskedBED, 'w') as bedout:
        for file in natsorted(os.listdir(outputDir)):
            if file.endswith('.bed'):
                with open(os.path.join(outputDir, file), 'rU') as infile:
                    bedout.write(infile.read())
    #generate masked fasta sequences
    SeqRecords = SeqIO.to_dict(SeqIO.parse(maskedOut, 'fasta'))
    with open(maskedFA, 'w') as outfile:
    	with open(maskedBED, 'rU') as infile:
    		for line in infile:
    			if line.startswith('#') or line.startswith('\n'):
    				continue
    			line = line.strip()
    			contig, start, end, name = line.split('\t')
    			Seq = str(SeqRecords[contig][int(start):int(end)].seq)
    			outfile.write('>{:}\n{:}\n'.format(name, softwrap(Seq)))

    #get masked repeat stats
    GenomeLength = sum(Scaffolds.values())
    percentMask = maskedSize / float(GenomeLength)
    sys.stdout.write('\nMasked genome: {:}\nRepeat BED file: {:}\nRepeat FASTA file: {:}\nnum scaffolds: {:,}\nassembly size: {:,} bp\nmasked repeats: {:,} bp ({:.2f}%)\n\n'.format(
    				 os.path.abspath(maskedOut), os.path.abspath(maskedBED), os.path.abspath(maskedFA), 
    				 len(Scaffolds), GenomeLength, maskedSize, percentMask*100))

#clean up
if not args.debug:
    SafeRemove(inputDir)
    SafeRemove(trainDir)
    SafeRemove(outputDir)