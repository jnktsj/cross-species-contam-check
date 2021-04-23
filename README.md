## Cross-species contamination check

### settings

The repo contains python scripts that check cross-species contamination
by mapping alt supporting reads to 57 mammalian genomes.  Before running
any python script, make sure to download the 57 mammalian genomes from
the UCSC Genome Database. `species_list.txt` contains the direct links to
download the genome 2bit files.

```shell
$ head species_list.txt

human	hg19	https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
mouse	mm39	https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit
alpaca	vicPac2	https://hgdownload.soe.ucsc.edu/goldenPath/vicPac2/bigZips/vicPac2.2bit
armadillo	dasNov3	https://hgdownload.soe.ucsc.edu/goldenPath/dasNov3/bigZips/dasNov3.2bit
baboon		papAnu4	https://hgdownload.soe.ucsc.edu/goldenPath/papAnu4/bigZips/papAnu4.2bit
bison		bisBis1	https://hgdownload.soe.ucsc.edu/goldenPath/bisBis1/bigZips/bisBis1.2bit
bonobo		panPan3	https://hgdownload.soe.ucsc.edu/goldenPath/panPan3/bigZips/panPan3.2bit
bushbaby	otoGar3	https://hgdownload.soe.ucsc.edu/goldenPath/otoGar3/bigZips/otoGar3.2bit
cat		felCat9	https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/felCat9.2bit
chimpanzee	panTro6	https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.2bit
...
```

### analysis workflow

The analysis workflow is built as per-chromosome base.  `extract_alt_reads.py` requires a MAF file that contains mutations.

```shell
# extract alt supporting reads in ${chrom}
python src/extract_alt_reads.py ${chrom} mutations.maf sample.bam | sort -u > ${chrom}_alt_reads.txt

# create fasta from the text file above
awk '{print ">"$2":"$3"\n"$4}' ${chrom}_alt_reads.txt > ${chrom}_alt_reads.fa

# run blat on the chr21 alt supporting reads against mammalian genomes
for i in `ls -1 twoBit/*.2bit`;
do
  qsub -w e -j y -b y -V -cwd -r y -l h_vmem=8G -l h_rt=7:00:00 blat $i ${chrom}_alt_reads.fa `basename $i .2bit`.psl
done

# download blat score calculation script and then compute blat scores
wget http://genome-source.soe.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl

# strip psl headers
mkdir tmp && mv *.psl tmp
for i in `ls -1 tmp/*.psl`;
do
  grep ^[0-9] $i > `basename $i`
done

# calculate BLAT score
for i in `ls -1 psl_headerless/*.psl`;
do
  perl pslScore.pl $i > `basename $i .psl`.score.txt
done

# clean up some blat files...
rm -r tmp psl_headerless && mkdir blatOut && mv *.score.txt blatOut

# create a matrix
python parse_blat_score.py ${chrom}_alt_reads.fa `ls -1 blatOut/*.score.txt` > ${chrom}_alt_blat_scores.txt
```

How to visualize the results can be found in `plot.R`.

