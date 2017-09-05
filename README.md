# esiCancer
Darlan Conterno Minussi, 
Bernardo Henz,
Mariana dos Santos Oliveira,
Eduardo C. Filippi-Chiela,
Manuel M. Oliveira,
Guido Lenz

### What is esiCancer?
esiCancer is a cell-autonomous software tool to model early evolution of cancer. This program uses real genomic information of tumors to construct a biallelic genome and applies stochastic mutations to this genome. Mutations can produce several alterations to the probability of proliferation, death, mutation rate, and maximal divisions of the cells containing mutation in one allele of oncogenes or two alleles of tumor suppressor genes (TSGs).

esiCancer has a user friendly interface and produces .cvs outputs that can be easily analyzed by computer professionals or amateurs alike.




### Compiling the code
Source code is written in C++, using Qt library for interface design. For compiling the code you will need to: (1) having a compiler installed (we suggest GCC or MinGW, but any compiler should work fine); (2) having the QT library, which can be download on [QT website](https://www.qt.io/download/). We strongly recommend the use of QtCreator (installed with QT), as we provide the .pro project file. Any problem compiling the code, please let me know.


### Start simulating clonal evolution
We provide two files with standard configurations. Put the *config.xml* file in the same folder as your runnable, run the program and load the *esiTable.csv*. By clicking in *Automatic Runs*, the program will simulate several runs (accordingly to the values on the spinners), and generating output files. Now, you can modify your parameters to test different hypothesis, checking the outcome in no time. 

### Parameters
The parameters are described next:

|   | **Parameter** | **Action** |
| --- | --- | --- |
|Essential parameters | Seed | Controls the pseudorandom number generator |
| | Number of cells | Starting number of cells |
| | Mutations per division | Number of mutations to be inserted into each cell at each division |
| | Probability of proliferation | Default probability of proliferation for all cells |
| | Probability of death | Default probability of proliferation for all cells |
| | Maximum divisions | Number of divisions that a cell can go through |
|Automatic runs | Parameter to iterate | Parameter that will be automatically iterated |
| | Minimum value | Starting iterated parameter value |
| | Maximum value | Ending iterated parameter value |
| | Increment | Defines the increment to be used from minimum to maximum |
| | Output Files | Defines the format of output files |
|Stop Conditions | Number of generations | Number of maximum generations allowed per run |
| | Number of Cells | Maximum number of cells |
| | Mutated Cells | Maximum number of mutated cells |
|Manual Mutations | Generation | Generation that mutation rule is added |
| | Percentage | Percent of cells that is added mutation rule |
| | Mutation Name | Gene that is added mutation rule (according to mutation table genes) |
| | Allels Changed | Defines the mutation rule format (mono or bi-allelic) |

### Output Files
When simulating different runs on *Automatic Runs*, our program generates output files informing important info about the runs. Each output file presents different data about the runs:

| **Output** | **Data presented** |
| --- | --- |
| ancestralResults.csv | Returns the number of affected cells for each initial cell per generation, as the total number of affected cell per generation |
| GeneMutationResults.csv | Returns the absolute number of cells mutated for each gene per generation |
| numberOfCellsWithNMutations.csv | Returns the absolute number of cells with N mutations by generation, the maximum number of mutations per cell is not fixed |
| parameters.txt | Returns the input parameters (seed, Number of Cells, Proliferation Rate, Death Rate, Maximum Division, Mutations per Division) |
| automatics\_ancestralResults.csv | Returns the last row of ancestralResults.csv for each seed in one unique file |
| automatics\_GeneMutationResults.csv | Returns the last row of GeneMutationResults.csv for each seed in one unique file |
| automatics\_numberOfCellsWithNMutations.csv | Returns the last row of numberOfCellsWithNMutations.csv for each seed in one unique file |


### Running into issues?
Contact Bernardo Henz <bernardohenz@gmail.com>
