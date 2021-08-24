**Supplementary Notes**
# 1. Commands used for downloading resource data.
## (1) Download the nstd162 and nstd137 callsets.
```
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd162.GRCh37.variant_call.vcf.gz
gzip -d nstd162.GRCh37.variant_call.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd137.GRCh37.variant_call.vcf.gz
gzip -d nstd137.GRCh37.variant_call.vcf.gz
```
## (2) Download hg19 reference genome.
```
curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > human_hs37d5.fasta.gz
gunzip human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' human_hs37d5.fasta
```

# 2. Commands used for simulating diverse settings of long-read.
## 2.1 Integration of the SVs in CHM1 human sample from the nstd162 into the reference genome using VISOR.
### (1) Extraction of deletion.bed.
```
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DEL | awk -F '\t' '{print $1"\t"$2}' > del_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DEL | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\tdeletion\tNone\t0"}' > del_col2.txt
paste del_col1.txt del_col2.txt | uniq > del.bed
```
### (2) Extraction of insertion.bed.
```
grep EXPERIMENT=2 nstd137.GRCh37.variant_call.vcf | grep INS |awk -F '\t' '{print $1"\t"$2-1"\t"$2}' > ins_col1.txt
grep EXPERIMENT=2 nstd137.GRCh37.variant_call.vcf | grep INS |awk -F '\t' '{print $8}' | awk -F ';SEQ=' '{print "insertion\t" $2 "\t0"}' >ins_col2.txt
paste ins_col1.txt ins_col2.txt | uniq > del.bed
```
### (3) Extraction of duplication.bed.
```
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DUP | awk -F '\t' '{print $1"\t"$2}' > dup_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep DUP | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\ttandem duplication\t2\t0"}' > dup_col2.txt
paste dup_col1.txt dup_col2.txt | uniq > dup.bed
```
### (4) Extraction of inversion.bed.
```
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep INV | awk -F '\t' '{print $1"\t"$2}' > inv_col1.txt
grep CHM1 nstd162.GRCh37.variant_call.vcf | grep INV | awk -F '\t' '{print $8}' | awk -F ';END=' '{print $2}' | awk -F ';' '{print $1"\tinversion\tNone\t0"}' > inv_col2.txt
paste inv_col1.txt inv_col2.txt | uniq > inv.bed
```
### (5) Creation of simulated donor genome.
```
cat del.bed ins.bed dup.bed inv.bed | sort -k 1,1 -k 2,2n > chm1.bed
git clone https://github.com/SQLiu-youyou/SV_evaluation.git
VISOR HACk -g human_hs37d5.fasta -bed SHORtS.LASeR.bed -o chm1
```
## 2.2 Generation of the synthetic diplontic long-reads within various sequencing attributes via VISOR.
### (1) Generate long-reads under various sequencing coverages.
```
VISOR LASeR -g human_hs37d5.fasta -s chm1/ -bed SHORtS.LASeR.bed -o [COVERAGE]_20k_90 -c [COVERAGE] -a 0.9 -l 20000 --readstype PB --threads 8 --noaddtag 
where [COVERAGE] represents 3×, 5×, 10×, 20×, 30×, 40×, and 50×, respectively.
```
### (2) Generate long-reads under various mean read lengths.
```
VISOR LASeR -g human_hs37d5.fasta -s chm1/ -bed SHORtS.LASeR.bed -o 50×_[READLENGTH]_90 -c 50 -a 0.9 -l  [READLENGTH] --readstype PB --threads 8 --noaddtag 
where [READLENGTH] represents 2.5k, 5k, 7.5k, 10k, 15k, 20k, 50k, 100k and 500k, respectively.
```
### (3) Generate long-reads under various sequencing error rates.
```
VISOR LASeR -g human_hs37d5.fasta -s chm1/ -bed SHORtS.LASeR.bed -o 50×_20k_[ACCURATE] -c 50 -a [ACCURATE] -l 20000 --readstype PB --threads 8 --noaddtag 
where [ACCURATE] represents 0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.99, and 0.998 respectively.
When [ACCURATE] represents 0.99 and 0.998 -ccs parameter should be settled
```

# 3. Commands used for SV calling with each tool.
## (1) cuteSV
```
cuteSV sim.srt.bam human_hs37d5.fasta cutesv.vcf ./ -s [sup_read] -l 30 –genotype -mi 0
where [sup_read] are 1, 2, 3, 4, 5, 5, and 5 for 3×, 5×, 10×, 20×, 30×, 40×, and 50× coverage datasets, respectively.
```
## (2) Sniffles
```
sniffles -m sim.srt.bam -v sniffles.vcf -s [sup_read] -l 30 –genotype
where [sup_read] are 1, 2, 3, 4, 5, 5, and 5 for 3×, 5×, 10×, 20×, 30×, 40×, and 50× coverage datasets, respectively.
```
## (3) SVIM
```
svim alignment svim sim.srt.bam human_hs37d5.fasta --min_sv_size 30 --minimum_score 0 --minimum_depth 1
grep -v “SUPPORT=1;\|……\|SUPPORT=[sup_read-1];” svim/variant.vcf > ./svim.vcf
where [sup_read] are 1, 2, 3, 4, 5, 5, and 5 for 3×, 5×, 10×, 20×, 30×, 40×, and 50× coverage datasets, respectively.
```
## (4) PBSV
```
pbsv discover sim.srt.bam pbsv.svsig.gz -s pbsv && pbsv call human_hs37d5.fasta pbsv.svsig.gz pbsv.vcf
```
## (5) NanoSV
```
NanoSV sim.srt.bam -s samtools -o nanosv.vcf
```
## (6) NanoVar
```
nanovar sim.srt.bam human_hs37d5.fasta nanovar -x pacbio-clr
```
## (7) Ensemble calling
```
ls cutesv.vcf sniffles.vcf svim.vcf pbsv.vcf > sample_files
./SURVIVOR merge sample_files 1000 2 1 0 0 30 merge.vcf
```

# 4. Commands used for benchmarking.
## (1) Generation of various base calling for Truvari
```
git clone https://github.com/SQLiu-youyou/SV_evaluation.git
bash covert_vcf.sh
grep '[sv_type]\|#' chm1.vcf > [sv_type].chm1.vcf
where [sv_type] are DEL, INS, DUP and INV
bgzip -c [all_type].chm1.vcf > [all_type].chm1.vcf.gz
tabix [all_type].chm1.vcf.gz
where [all_type] are TOTAL, DEL, INS, DUP and INV
```
## (2) Generation of comparison calling for Truvari
```
Sniffles.vcf should be sorted advanced, so does NanoSV.vcf
grep '#' sniffles.vcf > head
grep -v '#' sniffles.vcf > body
sort -k 1,1 -k 2,2n body > Body
cat head Body > sniffles.vcf
rm head Body body
```
```
grep '[sv_type]\|#' [caller].vcf > [sv_type].[caller].vcf
where [sv_type] are DEL, INS, DUP and INV, [caller] are cuteSV, Sniffles, SVIM, PBSV, NanoSV and NanoVar for coverage datasets and apart from NanoSV and NanoVar for read length and error rate datasets
bgzip -c [all_type].[caller].vcf > [all_type].[caller].vcf.gz
tabix [all_type].[caller].vcf.gz
where [all_type] are TOTAL, DEL, INS, DUP and INV, [caller] are cuteSV, Sniffles, SVIM, PBSV, NanoSV and NanoVar
```
## (3) Benchmarking via Truvari
```
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o all -p 0 --sizemax 100000000
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o 0_99 -p 0 -s 0 --sizemax 99
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o 100_499 -p 0 -s 100 --sizemax 499
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o 500_999 -p 0 -s 500 --sizemax 999
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o 1k_5k -p 0 -s 1000 --sizemax 5000
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o 5k_10k -p 0 -s 5000 --sizemax 10000
truvari bench -b [all_type].chm1.vcf.gz -c [all_type].[caller].vcf.gz -o over_10k -p 0 -s 10000 --sizemax 100000000
where [all_type] are TOTAL, DEL, INS, DUP and INV, [caller] are cuteSV, Sniffles, SVIM, PBSV, NanoSV and NanoVar
```


