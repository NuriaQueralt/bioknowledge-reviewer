# BioKnowledge reviewer library
This is a Python library to create structured reviews integrating knowledge and data from different types. We built a library for a dynamic, interactive and evolving network construction and hypothesis generation process. It is designed to build the review network based on the research question hypothesis. In Figure workflow we show the workflow of network-based review and hypothesis generation process. 

##### Prerequisites
Python 3. We provide a requirements.txt file to set a virtual environment to run the library for the creation of structured reviews around the NGLY1 Deficiency. The library + the environment runs without problems in an Ubuntu 18.04 distribution.


Neo4j needs Java 8. It won't work with superior java versions.

### Input / Output
#### Input
Edges to build the structured review.

Data used for Graph v3.2:
* curation/data/v20180118
* transcriptomics/ngly1-fly-chow-2018
* regulation/msigdb/data
* regulation/tftargets
* ontologies/mondo.owl


http://edamontology.org/data_2080

#### Output
Structured review network.

http://edamontology.org/data_2600


### Library architecture
The library is currently under testing.  

#### Edges
##### Curation
Edges from curation.

##### Monarch
Edges from the Monarch Knowledge Graph.

##### Transcriptomics
Edges from RNA-seq graph.

##### Regulation
Edges from regulation graph.

#### Graph
Build the review network.

#### Neo4jlib
Store the network in the Neo4j graph database.

#### Hypothesis
Compute review-based explanations.

#### Summary
Summarize extracted explanations.

#### Utils
Utils.py

#### Ontologies
mondo_class.py

### Usage
This sections showcase examples of use.

To run the library to reproduce the generation of the NGLY1 Deficiency Knowledge Graph v3.2, the user can use either the jupyter notebook or the python script provided in this repository. To run the jupyter notebook, the user should have installed the Jupyter framework and the ipython 3 kernel.


#### 1. Build a review network
Build the network by compiling edges.


bzip2 -d ontologies/mondo.owl.bz2
#### License
GPL v3.0

#### DOI
[![DOI](https://zenodo.org/badge/132827298.svg)](https://zenodo.org/badge/latestdoi/132827298)
