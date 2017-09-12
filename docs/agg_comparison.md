We want the output from initial versions of this tool to be on parity with agg. This just a basic sanity check. We expect the output to (rapidly) evolve such that it becomes inconsistent with agg's crude representations. At that point, we need to start evaluate verses GATK/Sention joint-called output.


agg commands:

```
ls data/*tiny.vcf.gz > gvcfs.txt
python ~/workspace/agg/make_chunk.py  gvcfs.txt -ref data/tiny.ref.fa -o chunk
~/workspace/agg/agg genotype chunk.bcf -Ob -o agg.output.bcf
```

gvcfmerge commands:

```
./gvcfmerge -f data/tiny.ref.fa -l gvcfs.txt -Ob -o gvcfmerge.output.bcf
```

check variant set is exactly the same:

```
bcftools index agg.output.bcf 
bcftools index gvcfmerge.output.bcf 
bcftools stats agg.output.bcf gvcfmerge.output.bcf | grep SN
```

output: 

```
# SN, Summary numbers:
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	3
SN	1	number of samples:	3
SN	0	number of records:	0
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	0
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0
SN	0	number of multiallelic SNP sites:	0
SN	1	number of records:	0
SN	1	number of no-ALTs:	0
SN	1	number of SNPs:	0
SN	1	number of MNPs:	0
SN	1	number of indels:	0
SN	1	number of others:	0
SN	1	number of multiallelic sites:	0
SN	1	number of multiallelic SNP sites:	0
SN	2	number of records:	495
SN	2	number of no-ALTs:	0
SN	2	number of SNPs:	398
SN	2	number of MNPs:	0
SN	2	number of indels:	89
SN	2	number of others:	8
SN	2	number of multiallelic sites:	0
SN	2	number of multiallelic SNP sites:	0

```

check consistency of genotypes:

```
$ vcftools --bcf agg.output.bcf --diff-bcf gvcfmerge.output.bcf  --diff-discordance-matrix

VCFtools - 0.1.15
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--bcf agg.output.bcf
	--diff-bcf gvcfmerge.output.bcf
	--diff-discordance-matrix

Using zlib version: 1.2.8
After filtering, kept 3 out of 3 Individuals
Using zlib version: 1.2.8
Outputting Discordance Matrix
	For bi-allelic loci, called in both files, with matching alleles only...
Non-matching REF. Skipping all such sites.
Non-matching ALT. Skipping all such sites.
Found 485 sites common to both files.
Found 0 sites only in main file.
Found 0 sites only in second file.
After filtering, kept 495 out of a possible 495 Sites
Run Time = 0.00 seconds
joconnell@ubuntu:~/workspace/gvcfmerge/debug$ cat out.diff.discordance_matrix 
-	N_0/0_file1	N_0/1_file1	N_1/1_file1	N_./._file1
N_0/0_file2	317	0	0	0
N_0/1_file2	0	1043	0	0
N_1/1_file2	0	0	71	0
N_./._file2	0	0	0	0
```