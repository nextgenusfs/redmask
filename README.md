# redmask
Genome assembly soft-masking using Red (Repeat Detector)

Requires Red (REpeat Detector) software which you can get [http://toolsmith.ens.utulsa.edu](http://toolsmith.ens.utulsa.edu).  Install one of the pre-built binaries for your system or build from source.  The `Red` executable needs to be in your $PATH. 

You can then run the Red mediated genome masking:
```
redmask.py -i genome.fa -o mygenome
```
This will run the masking generating the following files:

mygenome.softmasked.fa -- soft-masked genome
mygenome.repeats.bed -- BED file of repeats
mygenome.repeats.fasta -- FASTA file of masked sequences (putative repeats)


