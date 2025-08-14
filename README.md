# Mosaic Patterns in the Cortex

**Code for the paper _“Evidence from spatial transcriptomics for the mosaic hypothesis and pure cell types in the cortex”_**

## Overview
This repository contains the code to test the statistical validity of the mosaic hypothesis introduced in the paper ([BioRxiv link](https://www.biorxiv.org/content/10.1101/2024.08.09.607193v2)). While the code provided is intended to reproduce figures and analysis in the paper, the method can be adapted to other spatial transcriptomic dataset with cell typing and 2D cell coordinates. Please see instruction under **Usage**.


## Getting Started
### Prerequisites
- Python (>=3.8)
- R (>=3.5.0)
- spatstat (>=3.4-0)

We suggest you export the repository to Code Ocean for a fully interactive run. Alternatively, we suggest you download the repository (git clone https://github.com/yiliuw/mosaics.git) so you can access it locally with Anaconda. You can also click the "download" button, rather than using git.


### Usage
We provide data and several Jupyter Notebooks in `code/` folder. The notebooks are arranged in the order of main figures in the paper. We recommend reading the following steps along with the main paper. If you choose to download the code to your own computer, please also check the `environment/` folder for environment setups. 
1. **F1 - Data and mechanism** introduces the exemplar dataset (*illustration.ipynb*) and key mathematical concepts (*mechanism.ipynb*). While the data is publicly available at [Data link](https://alleninstitute.github.io/abc_atlas_access/intro.html), sample of data used for demonstration is included in this folder (for details, see *sagittal.ipynb*).
2. **F2 - Hypothesis testing** includes key algorithms for cell type filtering (*filter.ipynb*) and hypothesis testing(*hypothesis-x.ipynb*). The outputs for excitatory and inhibitory cell types are separated. This folder is key if you are interested in testing the hypothesis for your own spatial transcriptomic data.
3. **F3 - Segregation effects** introduces the concept of segregation index (*index.ipynb*) and discusses the differences of excitatory and inhibitory cell types at the subclass level (*distribution.ipynb*).
4. **F4 - ECC** introduces the effective cluster count and the construction of the tree hierarchy diagram as we merge cell types. 

As an example of adapting the method to another spatial transcriptomic dataset, see notebooks in the `code/Supplementary dataset/` folder.

---


## Citation
If you use this code, please cite:

> Wang, Y., Koch, C., & Sümbül, U. (2025). *Evidence from spatial transcriptomics for the mosaic hypothesis and pure cell types in the cortex.* Cell reports & [BioRxiv](https://www.biorxiv.org/content/10.1101/2024.08.09.607193v2).

---


