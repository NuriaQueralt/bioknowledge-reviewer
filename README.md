# BioKnowledge reviewer library
This is a Python library to create structured reviews integrating knowledge and data from different types. We built a library for a dynamic, interactive and evolving network construction and hypothesis generation process. It is designed to build the review network based on the research question hypothesis. In Figure workflow we show the workflow of network-based review and hypothesis generation process. 

##### Prerequisites
Python 3

### Input / Output
#### Input
Edges to build the structured review.

http://edamontology.org/data_2080

#### Output
Structured review network.

http://edamontology.org/data_2600

### Library architecture
The library is currently in development. 

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


### Usage
This sections showcase examples of use.

#### 1. Build a review network
Build the network by compiling edges.

#### License
GPL v3.0

#### DOI
[![DOI](https://zenodo.org/badge/132827298.svg)](https://zenodo.org/badge/latestdoi/132827298)
