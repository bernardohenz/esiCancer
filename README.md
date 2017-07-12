# esiCancer
Darlan Conterno Minussi, 
Bernardo Henz,
Mariana dos Santos Oliveira,
Eduardo C. Filippi-Chiela,
Manuel M. Oliveira,
Guido Lenz

### What is esiCancer?
Knowledge about cancer biology and its consequent evolution is mainly inferred from tissue sampled at very limited time points, thus only representing a glimpse of the complete tumor evolutionary path. In this work, we developed a software to simulate the initial steps of cancer evolution. Our model randomly applies mutations to a real oncogenome of a subset of diploid cells, thus altering their probabilities of proliferation and death. Accumulation of mutations over hundreds of generations leads to a gradual change in cell fitness and establishment of a fast growing population. esiCancer reproduces essential stochastic characteristics of the evolution process and provides a tool to predict the impact on fitness of oncogenes and tumor suppressor genes. It also predicts the impact that mutation rate, fitness levels, and mono or bi-allelic changes in specific genes have on Kaplan-Meier-like curves. Thus, esiCancer offers a unique and realistic model for the initial tumor evolution process. 


### Compiling the code
Source code is written in C++, using Qt library for interface design. We strongly recommend the use of QtCreator, as we provide the .pro project file. We suggest using GCC or MinGW as compiler, but any other should work just fine. Any problem compiling the code, please let me know.


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