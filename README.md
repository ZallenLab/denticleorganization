# Scaling of cytoskeletal organization with cell size in <em>Drosophila</em>

This repository contains the code used in [Spencer, Schaumberg and Zallen, Molecular Biology of the Cell, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28404752) to investigate patterns of subcellular denticle organization in the <em>Drosophila</em> embryo. It is divided into two parts: 

## 1. Calculate denticle positions and separation distances
The first part is the MATLAB pipline for evaluating denticle organization in <em>in vivo</em> samples. This pipeline takes in the output of the [ImageJ](https://imagej.nih.gov/ij/) [CellCounter plugin](https://imagej.nih.gov/ij/plugins/cell-counter.html), calculates relative positions, and outputs a number of csv files, which can be used in other programs for further analysis [denticle spacing readme](/analyze_denticle_spacing_MATLAB/README_DenticleDist.md). 

## 2. Test models for the denticle organization pattern
The second part is the statistical model that tests key features of the denticle organization pattern by comparing simulated spacing distributions (generated within) to the measured <em>in vivo</em> data from the first part [statistical modeling readme](/statistical_modeling_PYTHON/README_StatModel.md).

