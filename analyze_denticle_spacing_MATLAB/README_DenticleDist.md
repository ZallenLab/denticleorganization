# Overview

This set of matlab scripts takes manually identified points (as a series of (x,y) coordinates with a 'type' label, as output by the ImageJ CellCounter plugin), calculates the absolute separation between each pair, then performs a variety of additional calculations and rearrangements that are output to a series of csv files. 


### References 
For more details, see [Spencer, Schaumberg and Zallen, Molecular Biology of the Cell, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28404752). If this code has been helpful to you, please cite our paper.



# Requirements 
- Matlab R2013
*This pipeline has also been run successfully on versions as late as R2016a, but not guaranteed to work as written on versions other than R2013*


### INPUT: 
.txt files containing (x,y) coordinates of points of interest
These were manually identified in [ImageJ](https://imagej.nih.gov/ij/index.html) using the [CellCounter](https://imagej.nih.gov/ij/plugins/cell-counter.html) plugin 
(if these are saved from ImageJ as '.xls' files, just replace '.xls' with '.txt' in Automator)

**These files MUST be named as follows:**

\<ID \#\>\_\<stage\>\_\<genotype\>\_\<\*optional\*\>\_b\<belt\#\>r\<row\#\>.txt

where
* \<ID #\> is a 4-digit identifier unique to each embryo
* \<stage\> is a 2-digit identifier
* \<genotype\>\_\<\*optional\*\> is where the code looks for the string of characters that you give it initially
* \<belt\#\> and \<row\#\> are one-digit integers (belt[3,8], row[0,7] typically)

This can be any length, but all embryos in a set should include at least some of the same strings


the following are all acceptable names
- 0151\_15\_yw\_b3r2.txt
- 0151\_15\_yw\_****\_b3r2.txt
- 0151\_15\_yellowwhite\_isoA\_b3r2.txt

but these are not 
- 151\_15\_yw\_b3r2.txt
- 0151\_15\_yw\_32.txt
- 151\_15\_yw\_b32.txt


**The integrity of the first 7 characters and the last 4 characters are essential for the program to work correctly!**


### OUTPUT:
- \<genotype\>\_stats.txt:  descriptive statistics of the dataset - embryos analyzed, cell & denticle numbers
- \<name\>.csv:             Data file 


__\*\_CellbyCell.csv  *(version w/ headers)* and \*\_cellbyCellStack.csv  *(version w/o headers)*__
this is all the data you can really want from the raw (x,y,label) txt files
contains the data from every cell analyzed, one row per cell

column | label | description
------ | ----- | -----------
01 | embryoID | unique # for each embryo, so that the script finds individual embryos, and so I can trace problems; embryo pattern is 0xxx
02 | row | aka cell column; [1:6]        
03 | belt | 'band' number, corresponding to the abdominal segment number; [3:8]        
04 | cell | from left to right (embryo tail-up) (or top to bottom, embryo with head pointing left) , the # of the cell within that specific belt/row (e.g. if there are 6 cells in that row/belt, there will be 6 rows of data, with this column [1:6] 
05 | dentincell | # of denticles in that cell
06 | Dvlen | length of the cell     
07 | dentEdgeL | distance between the left/top-most denticle and the left/top-most cell DV edge
08 | dentEdgeR | ditto for right/bottom-most
09 | AdjL | for a given cell, the distance from the left/top-most denticle in that cell to the closest denticle in the adjacent cell moving left/up
10 | AdjR | ditto but for right/down
11 | Intra2 | distance between denticle 1 and denticle 2 (from left to right)
12 | Intra3 |     
13 | Intra4 |     
14 | Intra5 |
15 | Intra6 |
16 | Intra7 |
17 | Intra8 |
18 | Intra9 |

__\*\_IntrabyCell.csv *(version w/ headers)* and \*\_intrabyCellStack.csv *(version w/o headers)*__
contains the data for cells with n≥2 denticles, with n lines for each cell (where n is the # of denticles that cell posesses)

column | label | description
------ | ----- | -----------
01 | embryoID | unique # for each embryo, so that the script finds individual embryos, and so I can trace problems; embryo pattern is 0xxx
02 | row | aka cell column; [1:6]        
03 | belt | 'band' number, corresponding to the abdominal segment number; [3:8]        
04 | cell | from left to right (embryo tail-up) (or top to bottom, embryo with head pointing left) , the # of the cell within that specific belt/row (e.g | if there are 6 cells in that row/belt, there will be 6 rows of data, with this column [1:6] 
05 | dentincell | # of denticles in that cell
06 | Dvlen | length of the cell     
07 | Intra | distance between the left/top-most denticle and the left/top-most cell DV edge


__\*\_CellOrder.csv *(version w/ headers)* and \*\_cellOrderStack.csv *(version w/o headers)*__
*contains the data for cells with n≥2 denticles, with n lines for each cell (where n is the # of denticles that cell posesses)

column | label | description
------ | ----- | -----------
01 | embryoID | unique # for each embryo, so that the script finds individual embryos, and so I can trace problems; embryo pattern is 0xxx
02 | row | aka cell column; [1:6]        
03 | belt | 'band' number, corresponding to the abdominal segment number; [3:8]        
04 | cell | from left to right (embryo tail-up) (or top to bottom, embryo with head pointing left) , the # of the cell within that specific belt/row (e.g if there are 6 cells in that row/belt, there will be 6 rows of data, with this column [1:6] 
05 | dentincell | # of denticles in that cell
06 | Dvlen | length of the cell     
07 | dentEdgeL | distance between the left/top-most denticle and the left/top-most cell DV edge
08 | AdjL | for a given cell, the distance from the left/top-most denticle in that cell to the closest denticle in the adjacent cell moving left/up
10 | Intra2 | distance between denticle 1 and denticle 2 (from left to right)
11 | Intra3 |    
12 | Intra4 |   
13 | Intra5 |
14 | Intra6 |
15 | Intra7 |    
16 | Intra8 |    
17 | Intra9 |
18 | dentEdgeR | distance between the right/bottom-most denticle and the left/top-most cell DV edge
19 | AdjR | for a given cell, the distance from the right/bottom-most denticle in that cell to the closest denticle in the adjacent cell moving right/down 




#### FILES & DEPENDENCIES:
    
*primary scripts*   
* __AnalyzeDenticleOrganization.m__  this is the main script that calls all the others; creates the primary .csv files that are used by all other functions
* __DenticleOrganization\_ImportData.m__ parses \*.csv files for a single genotype, as specified by the user in _AnalyzeDenticleOrganization.m_, creates statistics files, calls _DataDivider.m_ & _DenticleCalculations.m_
* __DataDivider.m__ imports \*.csv from _AnalyzeDenticleOrganization.m_ to output derivative calculations
* __DenticleCalculations.m__ called by _DenticleOrganization\_ImportData.m_ to output derivitive calculations

*generic functions*
* AddHeaders.m
* ColumnShifter.m
* EmbryoAvg.m
* GimmeFractions.m
* GroupAvg.m
* SpiffyName.m



## USAGE & DIRECTIONS

Things that the user will need to change (make sure you save the file before running it):

* In AnalyzeDenticleOrganization.m:
    - `(line 33)`       set the file path to the folder where these scripts and your data files reside
    - `(lines 53-55)`   add scale at which your images were taken to the list & comment out the others


Input files and .m files need to be in the same place (the dir set above) 



### Running

Open AnalyzeDenticleOrganization.m

Run (F5 is the shortcut when you are in the Editor window)
Make sure you are looking at the 'Command Window' tab in Matlab so you can respond to user prompts

\>Print-out of the current set scale, and the date & time of the run

1. user prompt: What is the genotype & stage?
    * This is the string that will be used to name your OUTPUT files. you can make it as long & verbose as you like. 
    * The majority of the time, you can just input a short string here (that will be in all your files), then hit [Enter] for all other prompts (use the default values for each)
    
2. user prompt: What is the genotype in the filenames?
    * This is the string that the program will look for in all your data files. 
    * IF the string you put in above is in all your desired input files, this is unnecssary. just hit [Enter]
    * IF you put in additional info above, or the string is NOT found in all your data files, enter the search string here  
    * for example: 
        - if all your text files are named with \*\*\*\*\_\*\*\_yw\_b\*r\*, and you enter at (1) 'yw', then at (2) you can hit [enter]
        - but if you enter 1. yellowwhite_stage15, then at prompt (2) you will have to enter yw so that it can find your files 

3. user prompt: Which calculations? 1=[pool]  2=[pool+indiv]  3=[indiv]?
    Calculations will look at.... 
    *    1=[pool]  = all data in a big group
    *    2=[pool+indiv]  = all data in a big group AND each embryo alone
    *    3=[indiv]  =  each embryo alone

4. user prompt: Make auto-gen figures? (DIC, etc)  [0 = no];  1 = yes 
[default = 0/no]
    * This can be ignored in this version (just hit enter/0)


### Authors
Written by Alison Spencer in Jennifer Zallen's lab at the Sloan Kettering Institute. 