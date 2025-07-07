# Software

* [JCVI](https://github.com/tanghaibao/jcvi)
```
pip install git+https://github.com/tanghaibao/jcvi.git
conda install -c conda-forge numpy scipy matplotlib pandas networkx
conda install -c bioconda diamond samtools
conda install -c conda-forge imagemagick
conda install -c conda-forge tectonic
jcvi --version
```
* [FastOMA](https://github.com/DessimozLab/FastOMA)
* [OrthoFinder](https://github.com/davidemms/OrthoFinder)
```
conda install orthofinder -c bioconda
```


# Data sources

* [GDR](https://www.rosaceae.org/)
* [ensembl-compara](https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-60/tsv/ensembl-compara/homologies/)
* [PLAZA](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_dicots/download/download)
* [Grapedia](https://grapedia.org/files-download/)
* [OrthoDB](https://www.orthodb.org/)

# Genomes, proteomes, and General Feature Format files

Reference genomes/proteomes in **bold**

* _Arabidopsis thaliana_: thale cress (ath)
  - [**ARAPORT11**](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts
  - [TAIR10](https://phytozome-next.jgi.doe.gov/)
  - [C24](https://www.orthodb.org/)
* _Malus domestica_: apple (mdo)
  - [**Golden delicious GDDH13 v1.1**](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts
  - [Honeycrisp Genome v1.1.a1, haplotype 1](https://www.rosaceae.org/species/malus/malus_x_domestica)
* _Prunus amygdalus_, _Prunus dulcis_: almond (pdul)
  - [**Texas v3.0, F1**](https://www.rosaceae.org/organism/24336)
* _Prunus armeniaca_: apricot (parm)
  - [**Marouch n14 v1.0**](https://www.rosaceae.org/organism/24338)
* _Prunus avium_: wild cherry (pavi)
  - [**Tieton v2.0**](https://www.rosaceae.org/organism/24337)
* _Prunus cerasifera_: cherry plum (pcer)
  - [**Montmorency' v1.0.a2**](https://www.rosaceae.org/organism/24337)
* _Prunus persica_: peach (ppe)
  - [**Lovell v2.0.a1**](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts
* _Prunus sibirica_: sibearian apricot (psib)
  - [**F106 v1.0**](https://www.rosaceae.org/organism/26133)
* _Pyrus communis_: pear (pcox)
  - [**d'Anjou v2.3.a1, hap 1**](https://www.rosaceae.org/organism/24590)
  - [Pyrus ussuriensis x Pyrus communis](https://www.orthodb.org/)
* _Solanum lycopersicum_: tomato (sly)
  - [**Heinz 1706 ITAG 4.0**](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts 
  - [ITAG 5.0](https://phytozome-next.jgi.doe.gov/)
  - [ITAG 3.0](https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/annotation/)
* _Solanum tuberosum_: potato (stu)
  - [**UniTato**](https://unitato.nib.si/)
  - [PGSC 4.03](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts
  - [PGSC DM v6.1](https://spuddb.uga.edu/dm_v6_1_download.shtml)
* _Vitis vinifera_: grapevine (vvi)
  - [**T2T (v5)**](https://grapedia.org/files-download/)
  - [12Xv2](https://bioinformatics.psb.ugent.be/plaza.dev/_dev_instances/feedback/download/download) primary/selected transcripts
  - [v4](https://grapedia.org/files-download/)
  - [12Xv0](https://grapedia.org/files-download/)
  - [12X](https://grapedia.org/files-download/)
  - [8x](https://grapedia.org/files-download/)



# Approaches

## JCVI

between species

* keeping one isoform per gene
```
python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only X.gff3.gz -o X.bed;
```
or
```
python -m jcvi.formats.gff bed --type=mRNA --key=ID --primary_only X.gff3 -o X.bed;
```
depending on the .gff

* CDS fasta

```
python -m jcvi.formats.fasta format X.cds.fa.gz X.cds
```

* pairwise synteny search
```
python -m jcvi.compara.catalog ortholog plant1 plant2 --no_strip_names
```

* pairwise synteny visualization
```
python -m jcvi.graphics.dotplot plant1.plant2.anchors
python -m jcvi.compara.synteny depth --histogram plant1.plant2.anchors
```
* macrosynteny visualization
```
python -m jcvi.compara.synteny screen --minspan=30 --simple plant1.plant2.anchors plant1.plant2.anchors.new
python -m jcvi.graphics.karyotype seqids layout
```

## FastOMA

between species

```
nextflow run dessimozlab/FastOMA -profile conda \
         --input_folder ./input  \
         --omamer_db ./input/Viridiplantae.h5 \
         --output_folder ./output \
         --report | tee run.log
```

## OrthoFinder

within genomes

```
#!/bin/bash

# Number of threads to use
THREADS=28

# File containing input/output directory pairs
PAIR_FILE="dir_pairs.txt"

while read -r input_dir output_dir; do
  echo "Running OrthoFinder for input: $input_dir, output: $output_dir"
  
  orthofinder -f "$input_dir" \
              -M msa \
              -S diamond_ultra_sens \
              -T iqtree \
              -t "$THREADS" \
              -n orthofinder \
              -o "$output_dir"
  
  # Check exit status
  if [ $? -ne 0 ]; then
    echo "OrthoFinder failed for $input_dir"
  else
    echo "OrthoFinder completed for $input_dir"
  fi

done < "$PAIR_FILE"

```

