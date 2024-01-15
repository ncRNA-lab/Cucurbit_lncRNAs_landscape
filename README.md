# Cucurbit lncRNAs landscape

 <div align="justify"> Long non-coding RNAs (lncRNAs) constitute a fascinating class of regulatory RNAs, widely distributed in eukaryotes. In plants, they exhibit features such as tissue-specific expression, spatiotemporal regulation and responsiveness to environmental changes, suggesting their involvement in specific biological processes. Although an increasing number of studies support the regulatory role of lncRNAs in model plants, our knowledge about these transcripts in relevant crops is limited. </div>

<br />

<div align="justify"> This repository includes the scripts used in the paper "Identification, characterization and transcriptional analysis of long non-coding RNAs in Cucurbits" to identify and analyze the lncRNAs for the nine species present in the table below. </div>

<br />
<br />

| <sub>Species</sub>                            | <sub>ID</sub>  | <sub>Genome<br />Size<br />(bases; Mb)</sub> | <sub>Download (Samples)</sub> | <sub>QC (Samples)</sub> | <sub>Strand<br />Specific<br />(Samples)</sub> | <sub>Download (Projects)</sub> | <sub>QC (Projects)</sub> | <sub>Strand<br />Specific<br />(Projects)</sub> | <sub>Final<br />Data Size (bytes; Gb)</sub> |
|-----------------------------------------------|----------------|------------------|---------|-----------------|----------|------------------|-----------------|----------|----------------------|
| <sub>*Cucumis sativus*</sub>                  | <sub>csa</sub> | <sub>226.21</sub>           | <sub>1388</sub>    | <sub>1334</sub>            | <sub>360</sub>      | <sub>129</sub>              | <sub>127</sub>             | <sub>35</sub>       | <sub>1167.43</sub>              |
| <sub>*Cucumis melo*</sub>                     | <sub>cme</sub> | <sub>357.86</sub>           | <sub>802</sub>     | <sub>776</sub>             | <sub>383</sub>      | <sub>45</sub>               | <sub>44</sub>              | <sub>16</sub>       | <sub>820.87</sub>               |
| <sub>*Citrullus lanatus*</sub>                | <sub>cla</sub> | <sub>365.45</sub>           | <sub>717</sub>     | <sub>711</sub>             | <sub>231</sub>      | <sub>55</sub>               | <sub>54</sub>              | <sub>17</sub>       | <sub>663.61</sub>               |
| <sub>*Lagenaria siceraria*</sub>              | <sub>lsi</sub> | <sub>313.81</sub>           | <sub>92</sub>      | <sub>92</sub>              | <sub>9</sub>        | <sub>7</sub>                | <sub>7</sub>               | <sub>3        | <sub>27.07</sub>                |
| <sub>*Cucurbita moschata*</sub>               | <sub>cmo</sub> | <sub>273.42</sub>           | <sub>127</sub>     | <sub>126</sub>             | <sub>39</sub>       | <sub>17</sub>               | <sub>16</sub>              | <sub>6        | <sub>102.73</sub>               |
| <sub>*Cucurbita argyrosperma*</sub>           | <sub>car</sub> | <sub>231.58</sub>           | <sub>10</sub>      | <sub>10</sub>              | <sub>9</sub>        | <sub>2</sub>                | <sub>2</sub>               | <sub>2        | <sub>30.36</sub>                |
| <sub>*Cucurbita pepo*</sub>                   | <sub>cpe</sub> | <sub>263.38</sub>           | <sub>143</sub>     | <sub>142</sub>             | <sub>50</sub>       | <sub>13</sub>               | <sub>13</sub>              | <sub>7        | <sub>112.61</sub>               |
| <sub>*Cucurbita maxima*</sub>                 | <sub>cma</sub> | <sub>279.69</sub>           | <sub>50</sub>     | <sub>50</sub>              | <sub>27</sub>       | <sub>10</sub>               | <sub>10</sub>              | <sub>4        | <sub>43.36</sub>                |
| <sub>*Momordica charantia*</sub>              | <sub>mch</sub> | <sub>294.01</sub>           | <sub>74</sub>      | <sub>74</sub>              | <sub>8</sub>        | <sub>5</sub>                | <sub>5</sub>               | <sub>2        | <sub>27.73</sub>                |

<br />
<br />








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
