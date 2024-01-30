# Usage

## STEP 0: Select RNA-seq samples

<div align="justify"> We select all RNA-seq samples from Sequence Read Archive database (SRA) and from a particular species. To this end, we use the <a href="https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/">E-utilities</a>. <br /><br /></div>

```
cd 01-Sample_processing_and_selection/cme/
sbatch 00-Samples.sh
```

<div align="justify"> <br />As a result, we obtain a list of RNA-seq sample accessions and a metadata table with all the information by RNA-seq sample accession. </div>

## STEP 1: Download RNA-seq samples

## STEP 2: Quality control RNA-seq samples

## STEP 3: Select strand-specific RNA-seq samples

## STEP 4: Mapping

## STEP 5: Assembly

## STEP 6: Transcript annotation

## STEP 7: LncRNA prediction

## STEP 8: Quantification

## STEP 9: Get random intergenic regions

## STEP 10: Molecular properties comparison

## STEP 11: Genomic distribution

## STEP 12: Comparative genomics

## STEP 13: Tissue-specificity analysis

## STEP 14: Differential expression analysis


<div align="justify"> <br />All RNA-seq samples from Sequence Read Archive database and from a particular species are collected (SW: SRA Toolkit). Next, we remove adapters and filter the reads by quality (SW: Fastp). Then, we deduce whether the library is strand-specific or not, and select the strand-specific samples (SW: Salmon).<br /><br />Once the data have been preprocessed, we map the clean data to the reference genome (SW: Hisat2) and assemble the transcriptome using a genome-guided assembly approach (SW: Stringtie2). Finally, we classify the assembled transcripts by their genomic position relative to the protein-coding genes (SW: Gffcompare). As a result, a class code is assigned to each transcript.</div><br />

- <div align="justify"> <b>LncRNA prediction</b><br /><br />We select transcripts that have any of the following five class codes: “u” (intergenic), “x” (antisense), “i” (intronic) and, “o” (sense) or “e” (sense). Next, transcripts are filtered by length (longer than 200 nucleotides) and expression level (more than 0.3 FPKM).<br /><br />Then, we assess the coding potential of each transcript using three alignment-free computational tools (SW: CPC2, FEELnc and CPAT) and two protein databases (Swissprot and Pfam). In the following step, we classify the transcripts into three confidence-levels (High-, medium- and low-confidence) according to the results of the previous step. After that, we annotate transcripts using different ncRNAs databases and those annotated as miRNA precursors (miRBAse and PmiREN) or housekeeping ncRNAs (RNAcentral) are discarded. We also annotate transcripts using databases of known potential lncRNAs in order to provide additional information (CANTATAdb, PLncDB, GreeNC). The program used to align the transcripts to the different ncRNAs databases is blastn. In the case of miRNA precursors, the MIReNA program is also used to validate them. </b><br /><br />Finally, we discard redundant transcripts (SW: CGAT) and create the database that will contain all potential lncRNAs. In addition, the different classes of potential lncRNAs will be renamed from intergenic, antisense, intronic and sense to lincRNA, NAT-lncRNA, int-lncRNA and SOT-lncRNA, respectively.</div><br />

- <div align="justify"> <b>Downstream analysis</b> </div><br />Now, there are some downstream analyses that we can perform:<br /><br />
 
    + Molecular properties comparison.<br />
    + Comparative genomics.</b>
    + Tissue-specificity analysis.<br />
    + Differential expression analysis (development and environment).<br />

<br />
<br />


<div align="justify"> This pipeline has been configured to be executed in a HPC cluster environment using the open-source workload manager Slurm. So, it was executed in <a href="https://garnatxadoc.uv.es/">Garnatxa HPC Cluster</a> located at Data Center of the Institute for Integrative Systems Biology (<a href="https://www.uv.es/institute-integrative-systems-biology-i2sysbio/en/institute-integrative-systems-biology-i-sysbio.html">I2SysBio</a>). All the scripts that compose the pipeline are stored in <a href="Scripts">Scripts</a>. As the scripts for each species are the same, only the scripts for cucumis melo (cme) have been uploaded. </div>

<br />

<div align="justify">Initially, this pipeline was designed for cucurbit species, but it could be adapted to other species. </div>

<br />

<div align="justify"> For more information you can read the paper <b>"Identification, characterization and transcriptional analysis of long non-coding RNAs in Cucurbits"</b>. </div>

<br />

## Software

- [E-utilities](https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools) v.2.11.2
- [FastQC](https://github.com/s-andrews/FastQC) v.0.11.9
- [Multiqc](https://github.com/MultiQC/MultiQC) v.1.11
- [Fastp](https://github.com/OpenGene/fastp) v.0.23.2
- [RSEM](https://github.com/deweylab/RSEM) v.1.3.3
- [Salmon](https://github.com/COMBINE-lab/salmon) v.1.6.0
- [HISAT2](https://github.com/DaehwanKimLab/hisat2) v2.2.1
- [StringTie2](https://github.com/gpertea/stringtie) v2.2.0
- [GffCompare](https://github.com/gpertea/gffcompare) v.0.12.6
- [CPC2](https://github.com/gao-lab/CPC2_standalone) v.1.0.1
- [FEELnc](https://github.com/tderrien/FEELnc) v.0.2
- [CPAT](https://cpat.readthedocs.io/en/latest/) v.3.0.2
- [DIAMOND](https://github.com/bbuchfink/diamond) v.2.0.14
- [Transdecoder](https://github.com/TransDecoder/TransDecoder) v.5.5.0
- [HMMER](https://github.com/EddyRivasLab/hmmer) v.2.0.14
- [NCBI-BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) v.2.13.0+
- [MIReNA](https://www.lcqb.upmc.fr/mirena/index.html) v.2.0
- [CGAT](https://cgat.readthedocs.io/en/latest/cgat.html) v.1.0
- [Bedtools](https://bedtools.readthedocs.io/en/latest/) v.2.27.1
- [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) v.2.0.3
- [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) v.4.1.3-p1
- [OrthoFinder](https://github.com/davidemms/OrthoFinder) v.2.5.4
- [MEME](https://meme-suite.org/meme/) v.5.5.1
- [Tspex](https://apcamargo.github.io/tspex/) (tissue-specificity calculator tool)

<br />

**R packages:**

<div align="justify"> circlize, colorspace, ComplexHeatmap, data.table, DESeq2, dplyr, GenomicFeatures, ggbreak, ggExtra, ggplot2, ggpubr, ggradar, ggridges, ggtext, ggvenn, grid, htmltools, limma, palmerpenguins, pheatmap, plotly, png, pRoloc, ragg, RColorBrewer, rtracklayer, scales, stringr, tibble, tidyr, tximport and UpSetR. </div>

<br />

## Databases

- [Swissprot](https://www.uniprot.org/help/downloads)
- [Pfam](https://www.ebi.ac.uk/interpro/download/pfam/)
- [RNAcentral](https://rnacentral.org/)
- [miRBase](https://mirbase.org/)
- [PmiREN](https://www.pmiren.com/)
- [CANTATAdb](http://cantata.amu.edu.pl/)
- [PLncDB](https://www.tobaccodb.org/plncdb/)
- [GreeNC](http://greenc.sequentiabiotech.com/wiki2/Main_Page)

<br />

## Authors

* **Pascual Villalba-Bermell** - *Initial work* - [pasviber](https://github.com/pasviber) (pascual.villalba@csic.es) at [ncRNA-lab](https://github.com/ncRNA-lab).

<br />



