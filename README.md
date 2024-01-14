# Cucurbit lncRNAs landscape

Long non-coding RNAs (lncRNAs) constitute a fascinating class of regulatory RNAs, widely distributed in eukaryotes. In plants, they exhibit features such as tissue-specific expression, spatiotemporal regulation and responsiveness to environmental changes, suggesting their involvement in specific biological processes. Although an increasing number of studies support the regulatory role of lncRNAs in model plants, our knowledge about these transcripts in relevant crops is limited. 

This repository includes the scripts used to identify the lncRNAs shown in the paper "Identification, characterization and transcriptional analysis of long non-coding RNAs in Cucurbits".


| Species                 | ID  | Genome Size (bases; Mb) | Download | Trimming and QC | Strand Specific | Download | Trimming and QC | Strand Specific | Final Data Size (Gb) |
|-------------------------|-----|------------------|---------|-----------------|----------|------------------|-----------------|----------|----------------------|
| *Cucumis sativus*       | csa | 226.21           | 1388    | 1334            | 360      | 129              | 127             | 35       | 1167.43              |
| *Cucumis melo*          | cme | 357.86           | 802     | 776             | 383      | 45               | 44              | 16       | 820.87               |
| *Citrullus lanatus*     | cla | 365.45           | 717     | 711             | 231      | 55               | 54              | 17       | 663.61               |
| *Lagenaria siceraria*   | lsi | 313.81           | 92      | 92              | 9        | 7                | 7               | 3        | 27.07                |
| *Cucurbita moschata*    | cmo | 273.42           | 127     | 126             | 39       | 17               | 16              | 6        | 102.73               |
| *Cucurbita argyrosperma*| car | 231.58           | 10      | 10              | 9        | 2                | 2               | 2        | 30.36                |
| *Cucurbita pepo*        | cpe | 263.38           | 143     | 142             | 50       | 13               | 13              | 7        | 112.61               |
| *Cucurbita maxima*      | cma | 279.69           | 50      | 50              | 27       | 10               | 10              | 4        | 43.36                |
| *Momordica charantia*   | mch | 294.01           | 74      | 74              | 8        | 5                | 5               | 2        | 27.73                |



In this study we employ a custom pipeline on a dataset of over 1,000 RNA-seq studies across nine representative species of the family Cucurbitaceae to predict 91,209 non-redundant lncRNAs. LncRNAs were predicted according to three confidence levels and classified into intergenic, natural antisense, intronic, and sense overlapping. Predicted lncRNAs have lower expression levels compared to protein-coding genes but a more specific behavior when considering plant tissues, developmental stages, and the response to environmental changes, emphasizing their potential roles in regulating various aspects of plant biology. Additionally, the evolutionary analysis indicates higher positional conservation than sequence conservation, which may be linked to the presence of conserved modular motifs within syntenic lncRNAs. In short, this research provides a comprehensive map of lncRNAs in the agriculturally relevant Cucurbitaceae family, offering a valuable resource for future investigations in crop improvement.



nf-core/sarek is a workflow designed to detect variants on whole genome or targeted sequencing data. Initially designed for Human, and Mouse, it can work on any species with a reference genome. Sarek can also handle tumour / normal pairs and could include additional relapses.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The Nextflow DSL2 implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from nf-core/modules in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the nf-core website.
<br />
<br />
<br />
<br />
<br />
<br />

![Image Alt text](Figure_1.png)

## Table of contents

  1. Sample processing
  2. X
  3. gg
  4. ff

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

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
