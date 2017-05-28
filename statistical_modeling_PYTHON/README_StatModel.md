# A statistical model testing the parameters of structure organization

### Overview 
To run via command line: 
python StatisticalModel_DenticleOrganization.py \<name of input file\> \<number of iterations\>

where the input file is a csv file (\*\_CellbyCell.csv), structured as described below.


### References 
For more details, see [Spencer, Schaumberg and Zallen, Molecular Biology of the Cell, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28404752). If this code has been helpful to you, please cite our paper.


## Usage
Written & tested primarily in python 3.5 (Anaconda).
Runs in python 2.7, but has not been extensively tested.


###### Required libraries 
- numpy
- pandas
- scipy (scipy.stats)
- matplotlib
- seaborn
(see requirements.txt)


To run via command line: 
python StatisticalModel_DenticleOrganization.py \<name of input file\> \<number of iterations\>


###### Input file 
Takes the output \*\_CellbyCell.csv from the MATLAB script set "Analyze Denticle Organization.m".


format: 

column | label | description
------ | ----- | -----------
01 | embryoID | unique # for each embryo, so that the script finds individual embryos, and so I can trace problems; embryo pattern is 0xxx
02 | row | aka cell column; [1:6]        
03 | belt | 'band' number, corresponding to the abdominal segment number; I focused on denticle belts 3-7, as these are most similarly strucured.    
04 | cell | from left to right (embryo tail-up) (or top to bottom, embryo with head pointing left) , the # of the cell within that specific belt/row (e.g. if there are 6 cells in that row/belt, there will be 6 rows of data, with this column [1:6] 
05 | dentincell | # of denticles in that cell
06 | Dvlen | length of the cell     
07 | dentEdgeL | distance between the left/top-most denticle and the left/top-most cell DV edge
08 | dentEdgeR | ditto for right/bottom-most
09 | AdjL | for a given cell, the distance from the left/top-most denticle in that cell to the closest denticle in the adjacent cell moving left/up
10 | AdjR | ditto but for right/down
11 | Intra2 | distance between denticle 1 and denticle 2 (from left to right)
12 | Intra3 | distance between denticle 2 and denticle 3
13 | Intra4 | distance between denticle 3 and denticle 4
.. | ...    | etc


### Options 

* choose number of repeats/iterations 
* choose to use the 'absolute' cell DV length or the summed total of denticle-edge and denticle-denticle distances (as in the equation below)
    * To swich between these optisons, comment/uncomment the sections of code under the headers 
        > # using the 'absolute' DV length (dist between edge markers)
        or 
        > using the summed DVlength (sum of dent-edge, dent-dent ... dent-dent, dent-edge); sum to get the additive, rather than absolute, DV length
        in StatisticalModel_DenticleOrganization*.py



#### Additional Information

Basic methodology:  
The equation L = α\*D + D(N-1) + α\*D describes the general case for determining the distance between denticles (D) given the known cell length (L) and the number of denticles (N) in the cell, where α is the spacing ratio (the ratio of the average denticle-to-cell-edge distance to the average denticle-to-denticle distance).



### Authors
Written by Alison Spencer in Jennifer Zallen's lab at the Sloan Kettering Institute. 

