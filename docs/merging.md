# GVCFGenotyper merging procedure

## Merging variants

This problem is a straightforward [K-way Merge](https://en.wikipedia.org/wiki/K-Way_Merge_Algorithms), although does require variants to be pre-processed before merging to ensure consistent representation across samples.  We have implemented the "Direct K-Way Merge" described on wikipedia, there are more efficient algorithms but it is unlikely that this is a bottleneck at the moment. 

## Canonicalising single sample input

We "canonicalise" variants from our input GVCFs prior to merging them (this is handled on-the-fly by the GVCFReader class). To simplify merging, we also ensure that one VCF row represents one alternative allele eg. a row with two alternate alleles is split into two rows on-the-fly. 

The canonicalisation steps are as follows:

1. Decomposition of multi-nucleotide polymorphisms (MNP) 
2. Variant normalisation (left-shifting and trimming) 
3. Multi-allelic variant splitting
4. Removal of duplicated variants (very rare)

These steps are applied in the `Normaliser::unarise` function. Note that all these processes are performed in memory, there are no intermediate files involved, we only use VCF representations below to describe the process.

### Step 1:

MNP decomposition is straightforward, any mismatching base between the `REF` and (possibly multiple) `ALT` becomes a new SNP and FORMAT/INFO fields are directly copied. 
 
Input:

```
7	150648198	.	ATA	GTG	49	PASS	.	GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL	0/1:14:4:39:6,10:1,3:5,7:PASS:88,0,11

```

Output:

```
7	150648198	.	A	G	49	PASS	OLD_CLUMPED=7:150648198:ATA/GTG	GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL	0/1:14:4:39:6,10:1,3:5,7:PASS:88,0,11
7	150648200	.	A	G	49	PASS	OLD_CLUMPED=7:150648198:ATA/GTG	GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL	0/1:14:4:39:6,10:1,3:5,7:PASS:88,0,11
```

This is handled in the `Normaliser::mnp_decompose` function.

### Step 2.

 We apply the normalisation routine described in the [vt paper](https://academic.oup.com/bioinformatics/article/31/13/2202/196142) to all indels.  Note that in that all Strelka bi-allelic variants (one `ALT`) appear normalised so this step is mostly redundant. However a small number of multi-allelics still appear unnormalised, this seems to be occurring at overlapping variants. 

For example, take the following variant:

```
chr1	197095592	.	AAAAAAG	AAA,A	738	LowGQX	.	GT:GQ:GQX:DPI:AD:FT:PL	1/2:151:0:29:1,15,14:1,7,8:0,8,6:LowGQX:816,279,187,308,0,231
``` 

or represented as haplotypes:

```
REF:  AAAAAAGGGAGGGGGACTT
ALT1: AAA----GGAGGGGGACTT
ALT2: A------GGAGGGGGACTT
```

the first alternate allele can have two bases trimmed and be moved to position `197095594`. The problem with this shifting is that the variants still overlap, meaning that values of `FORMAT/AD` and `FORMAT/PL` strictly should involve the overlapping allele but there is no good way to do this in VCF. Our pragmatic solution is to collapse the values into the reference call:

```
chr1	95593	.	AAAAAAG	A	738	.	.	GT:GQ:GQX:DP:DPF:AD:PL	0/1:151:0:30:.:16,14:187,0,231
chr1	95595	.	AAAAG	A	738	.	.	GT:GQ:GQX:DP:DPF:AD:PL	1/0:151:0:30:.:15,15:231,0,187
```
The above result is achieved by splitting an *N* alternate allele VCF row into *N* bi-allelic rows, normalising each row independently, and then collapsing alleles with the same normalised position into a multi-allelic row. In practice this often results in the output being identical to the input, but occasionally results into multiple new rows as in the above example.   

Finally, we would prefer to introduce the symbolic deletion allele `*`  here:

```
chr1	95593	.	AAAAAAG	A,*	738	.	.	GT:GQ:GQX:DP:DPF:AD	2/1:151:0:30:.:1,14,15
chr1	95595	.	AAAAG	A,*	738	.	.	GT:GQ:GQX:DP:DPF:AD	1/2:151:0:30:.:1,15,14
```
 
and plan to do this at a later date (https://jira.illumina.com/browse/GG-2).


### Step 3.

Our merging procedure only ever looks at the first alternate allele. To handle a row with *N>1* alternate alleles, we create *N* duplicate rows with each allele shuffled to the front. 

For example:

```
chr1:62430033:C:T:A	GT=1/2	AD=0,16,12	PL=370,214,166,242,0,206
```

becomes

```
chr1:62430033:C:T:A	GT=1/2	AD=0,16,12	PL=370,214,166,242,0,206
chr1:62430033:C:A:T	GT=2/1	AD=0,12,16	PL=370,242,206,214,0,166

```

Whilst there is some inefficiency in duplicating rows like this, very few variants in a single sample GVCF  are multi-allelc ( 61864/5042595=1.2%  of variants in NA12888) so the cost is negligble. The (debatable) advantage of this approach is that it simplifies the variant comparison when merging. I might revisit this logic at a later date but it is working for now!

Steps 2 and 3 are both implemented in `Normaliser::multi_split`.

### Step 4.

We still occasional see duplicated variants with different representations, this appears to be an artefact of forced genotyping (usually one of the variants has a clinvar annotation).  Sometimes the genotypes disagree between duplicated variants creating a conflict, in these cases we prefer the variant with an alternate genotype call over one with a homozygous reference call. If both samples are homref or alternate, we take the first copy of the variant that appears.

For example:

```
chr20	6065729	.	CTT	TTC	0	PASS	.	GT:GQ:GQX:AD:PL	0/0:71:71:25,0:0,74,530
chr20	6065729	.	C	T	177	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:210:30:45:1:25,20:211,0,265
chr20	6065731	.	T	C	187	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:220:30:47:0:26,21:222,0,272
```

After decomposition the two SNPs are duplicated with conflicting genotype calls:

```
chr20	6065729	.	C	T	0	PASS	.	GT:GQ:GQX:AD:PL	0/0:71:71:25,0:0,74,530
chr20	6065729	.	C	T	177	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:210:30:45:1:25,20:211,0,265
chr20	6065731	.	T	C	0	PASS	.	GT:GQ:GQX:AD:PL	0/0:71:71:25,0:0,74,530
chr20	6065731	.	T	C	187	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:220:30:47:0:26,21:222,0,272
```

we drop the homref calls in the output:

```
chr20	6065729	.	C	T	177	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:210:30:45:1:25,20:211,0,265
chr20	6065731	.	T	C	187	PASS	.	GT:GQ:GQX:DP:DPF:AD:PL	0|1:220:30:47:0:26,21:222,0,272
```
