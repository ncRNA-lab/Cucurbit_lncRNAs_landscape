# Cucurbit lncRNAs landscape

 <div align="justify"> Long non-coding RNAs (lncRNAs) constitute a fascinating class of regulatory RNAs, widely distributed in eukaryotes. In plants, they exhibit features such as tissue-specific expression, spatiotemporal regulation and responsiveness to environmental changes, suggesting their involvement in specific biological processes. Although an increasing number of studies support the regulatory role of lncRNAs in model plants, our knowledge about these transcripts in relevant crops is limited. </div>

<br />

<div align="justify"> In the paper <b>"Identification, characterization and transcriptional analysis of long non-coding RNAs in Cucurbits"</b>, we identified and analyze potential lncRNAs in nine representative species of the family <em>Cucurbitaceae</em> comprising a dataset of more than 1,000 RNA-seq studies. </div>

<br />
<br />

<div align="center"> <sub><b>Table 1: Information on the number of RNA-seq samples and projects used to predict and analyze lncRNAs in the nine selected cucurbit species.</b></sub> </div>

<br />

| <sub>Species</sub>                                  | <sub>Genome Size<br />(bases; Mb)</sub> | <sub>Download<br />(Samples)</sub> | <sub>QC<br />(Samples)</sub> | <sub>Strand Specific<br />(Samples)</sub> | <sub>Download<br />(Projects)</sub> | <sub>QC<br />(Projects)</sub> | <sub>Strand Specific<br />(Projects)</sub> | <sub>Final Data Size <br />(bytes; Gb)</sub> |
|-----------------------------------------------------|----------------------------------------------|-------------------------------|-------------------------|------------------------------------------------|--------------------------------|--------------------------|-------------------------------------------------|---------------------------------------------|
| <sub>*Cucumis sativus* (csa)</sub>                  | <sub>226.21</sub>                            | <sub>1388</sub>               | <sub>1334</sub>         | <sub>360</sub>                                 | <sub>129</sub>                 | <sub>127</sub>           | <sub>35</sub>                                   | <sub>1167.43</sub>                          |
| <sub>*Cucumis melo* (cme)</sub>                     | <sub>357.86</sub>                            | <sub>802</sub>                | <sub>776</sub>          | <sub>383</sub>                                 | <sub>45</sub>                  | <sub>44</sub>            | <sub>16</sub>                                   | <sub>820.87</sub>                           |
| <sub>*Citrullus lanatus* (cla)</sub>                | <sub>365.45</sub>                            | <sub>717</sub>                | <sub>711</sub>          | <sub>231</sub>                                 | <sub>55</sub>                  | <sub>54</sub>            | <sub>17</sub>                                   | <sub>663.61</sub>                           |
| <sub>*Lagenaria siceraria* (lsi)</sub>              | <sub>313.81</sub>                            | <sub>92</sub>                 | <sub>92</sub>           | <sub>9</sub>                                   | <sub>7</sub>                   | <sub>7</sub>             | <sub>3</sub>                                    | <sub>27.07</sub>                            |
| <sub>*Cucurbita moschata* (cmo)</sub>               | <sub>273.42</sub>                            | <sub>127</sub>                | <sub>126</sub>          | <sub>39</sub>                                  | <sub>17</sub>                  | <sub>16</sub>            | <sub>6</sub>                                    | <sub>102.73</sub>                           |
| <sub>*Cucurbita argyrosperma* (car)</sub>           | <sub>231.58</sub>                            | <sub>10</sub>                 | <sub>10</sub>           | <sub>9</sub>                                   | <sub>2</sub>                   | <sub>2</sub>             | <sub>2</sub>                                    | <sub>30.36</sub>                            |
| <sub>*Cucurbita pepo* (cpe)</sub>                   | <sub>263.38</sub>                            | <sub>143</sub>                | <sub>142</sub>          | <sub>50</sub>                                  | <sub>13</sub>                  | <sub>13</sub>            | <sub>7</sub>                                    | <sub>112.61</sub>                           |
| <sub>*Cucurbita maxima* (cma)</sub>                 | <sub>279.69</sub>                            | <sub>50</sub>                 | <sub>50</sub>           | <sub>27</sub>                                  | <sub>10</sub>                  | <sub>10</sub>            | <sub>4</sub>                                    | <sub>43.36</sub>                            |
| <sub>*Momordica charantia* (mch)</sub>              | <sub>294.01</sub>                            | <sub>74</sub>                 | <sub>74</sub>           | <sub>8</sub>                                   | <sub>5</sub>                   | <sub>5</sub>             | <sub>2</sub>                                    | <sub>27.73</sub>                            |

<br />
<br />

<div align="justify">All the lncRNAs identified for each of the species can be found in the <a href="Results">Results</a> folder in fasta, tsv and gtf format so that they can be used by the scientific community. </div>

<br />
<br />


## Pipeline

To identify and analyze potential lncRNAs, we used a custom pipeline consisting of three parts:

- <div align="justify"> <b>Data preprocessing and assembly</b><br /><br />All RNA-seq samples from Sequence Read Archive (SRA) database and from a particular species are collected (SRA Toolkit). Next, we remove adapters and filter the reads by quality (Fastp). Then, we deduce whether the library is strand-specific or not, and select the strand-specific samples (Salmon).<br /><br />Once the data have been preprocessed, we map the clean data to the reference genome (Hisat2) and assemble the transcriptome using a genome-guided assembly approach (Stringtie2). Finally, we classify the assembled transcripts by their genomic position relative to the protein-coding genes (Gffcompare). As a result, a class code is assigned to each transcript.</div><br />

- <div align="justify"> <b>LncRNA prediction</b><br /><br />We select transcripts that have any of the following five class codes: “u” (intergenic), “x” (antisense), “i” (intronic) and, “o” (sense) or “e” (sense). Next, transcripts are filtered by length (longer than 200 nucleotides) and expression level (more than 0.3 FPKM).<br /><br />Then, we assess the coding potential of each transcript using three alignment-free computational tools (CPC2, FEELnc and CPAT) and two protein databases (Swissprot and Pfam). In the following step, we classify the transcripts into three confidence-levels (High-, medium- and low-confidence) according to the results of the previous step. Those transcripts that don’t meet any of the present scenarios are not classified. After that, we annotate transcripts using different ncRNAs databases and those annotated as miRNA precursors or housekeeping ncRNAs are discarded. The program used to align the transcripts to the different ncRNAs databases is blastn. In the case of miRNA precursors, the MIReNA program is also used to validate them. We also annotate transcripts using databases of known potential lncRNAs in order to provide additional information.</b><br /><br />Finally, we discard redundant transcripts (CGAT) and create the database that will contain all potential lncRNAs. In addition, the different classes of potential lncRNAs will be renamed from intergenic, antisense, intronic and sense to lincRNA, NAT-lncRNA, int-lncRNA and SOT-lncRNA, respectively.</div><br />

- <div align="justify"> <b>Downstream analysis</b> </div><br />Now, there are some downstream analyses that we can perform:<br /><br />
 
    + Molecular properties comparison (LncRNAs vs Protein-coding genes).
    + Comparative genomics (Sequence, position and motif level).
    + Tissue-specificity analysis.
    + Differential expression (development and environment).

<br />
<br />

<div align="justify"> <sub><b>Figure 1: Graphic representation of the bioinformatic workflow used for the prediction, classification and analysis of lncRNAs. (A) Data recovering, preprocessing and transcriptome assembly. (B) Prediction and categorization of lncRNAs: intergenic (lincRNAs), natural antisense (NAT-lncRNAs), intronic (int-lncRNAs) and sense overlapping (SOT-lncRNAs) lncRNAs. (C) Downstream analysis to compare lncRNA features, conservation at three levels (sequence, genomic position and motifs) and differential expression (related to tissue, development and environment).</b></sub> </div>

<br />

![Image Alt text](Figure_1.png)

<br />
<br />

<div align="justify">The <a href="Scripts">scripts</a> stored in this github repository. </div>

<br />

<div align="justify"> For more information you can read the paper <b>"Identification, characterization and transcriptional analysis of long non-coding RNAs in Cucurbits"</b>. </div>

<br />

nf-core/sarek is a workflow designed to detect variants on whole genome or targeted sequencing data. Initially designed for Human, and Mouse, it can work on any species with a reference genome. Sarek can also handle tumour / normal pairs and could include additional relapses.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The Nextflow DSL2 implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from nf-core/modules in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the nf-core website.

## Softwares

- SRA Toolkit
- Fastp
- Salmon
- Hisat2
- Stringtie2
- Gffcompare


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
