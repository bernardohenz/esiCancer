# esiCancer
Darlan Conterno Minussi, 
Bernardo Henz,
Mariana dos Santos Oliveira,
Eduardo C. Filippi-Chiela,
Manuel M. Oliveira,
Guido Lenz

### What is esiCancer?

esiCancer is a cell-autonomous software tool to model evolution of cancer. This program uses real genomic information of tumors to construct a biallelic genome and applies stochastic events, such as mutations, indels and other alterations to this genome. Each event can produce several alterations to the probability of division, death, event rate per division, and/or maximal divisions, which impact the population of cells over time.

esiCancer has a user friendly interface and produces .cvs outputs that can be easily analyzed by computer professionals or amateurs alike.


### Compiling the code
Source code is written in C++, using Qt library for interface design. For compiling the code you will need to: (1) have a compiler installed (we suggest GCC or MinGW, but any compiler should work fine); (2) have the QT library, which can be downloaded from the [QT website](https://www.qt.io/download/). We strongly recommend the use of QtCreator (installed with QT), as we provide the .pro project file. Any problem compiling the code, please let me know.


### Start simulating clonal evolution
In the esiCancer interface, set the parameters and load the esiTable.csv. By clicking in Automatic Runs, the program will simulate several runs (according to the values on the spinners), generating output files. By modifying the parameters, you can test
different hypothesis, checking the outcome quickly.

### Parameters
The parameters are described next:

|   | **Parameter** | **Action** |
| --- | --- | --- |
|Simulation parameters | Seed | Controls the pseudorandom number generator |
| | Number of cells | Starting number of cells |
|Standard Values | Probability of division | Initial probability of division for all cells |
| | Probability of death | Initial probability of death for all cells |
| | Maximum divisions | Number of divisions that a cell can go through |
| | Events per division | Number of events to be inserted into each cell at each division |
|Automatic runs | Parameter to iterate | Parameter that will be automatically iterated |
| | Minimum value | Starting iterated parameter value |
| | Maximum value | Ending iterated parameter value |
| | Increment | Defines the increment to be used from minimum to maximum |
| | Output Files | Defines the format of output files |
|Stop Conditions | Number of generations | Number of maximum generations allowed per run |
| | Number of Cells | Maximum number of cells |
| | Mutated Cells | Maximum number of mutated cells |
|Tumor Environment | Maximal Tumor Growth Rate | Maximal population growth rate allowed per generation |
|Manual Events | Generation | Generation that event rule is added |
| | Percentage | Percent of cells that is added event rule |
| | Gene-event Name | Gene-event that is added event rule (according to esiTable) |
| | Allels Changed | Number of alleles affected by event rule(mono or bi-allelic) |

### esiTable Parameters
The esiTable parameters are described next:

|   | **Parameter Meaning** | **Action** |  
| --- | --- | --- |
|genome size | Genome Size | Sum of the oncogenome and the normal genome |
|GENES | Gene Name | Name of a given gene |
|EVENTS| Evet Name | Name of a given event |
|PROBEVENT| Probability of Event | Probability of a specific event to occur (Prob = PROBEVENT / genome size) |
|DOMINANCE| Dominance | Determines the impact of a specific event if it occurs in one allele. |
|DIVFUNC| Function for division rate | Determines if the event adds or multiplies a specific value to its division rate. |
|DIVRATE| Division Rate | Determines the impact of a specific event on division rate. |
|DEFUNC| Function for death rate | Determines if the event adds or multiplies a specific value to its death rate. |
|DEATHRATE| Death Rate | Determines the impact of a specific event on death rate. |
|MUTFUNC| Function for mutation rate | Determines if the event adds or multiplies a specific value to its event per division rate. |
|MUT| Mutation Rate | Determines the impact of a specific event on the number of events per division rate. |
|MAXDIVFUNC| Function for maximum division rate | Determines if the event adds or multiplies a specific value to its maximum number of division rate. |
|MAXDIVRATE| Maximum Division Rate | Determines the impact of a specific event on the maximum number of divisions rate. |
|MEFUNC| Function for microenvironment   | Determines if the event adds or multiplies a specific value to MTGR. |
|MICROENVIROMENT| Microenvironment  | Determines the impact of a specific event on MTGR. |

### interactionTable  Parameters
The interactionTable  parameters are described next:

|   | **Parameter Meaning** | **Action** |  
| --- | --- | --- |
|Gene_Before | Gene Name | Defines the gene name affected first. |
|Event_Before | Event Name | Defines the event name affected first. |
|Gene_After| Gene Name | Defines the gene name affected after the first. |
|Event_After| Event Name | Defines the event name affected after the first. |
|DIVMOD| Division rate modifier | Determines the impact of interaction on division rate. |
|DEMOD| Death rate modifier | Determines the impact of interaction on death rate. |
|MAXDIVMOD| Maximum division rate modifier | Determines the impact of interaction on maximum division rate. |
|Link| Gene-Event Link | Determines if an event affects more than one gene (1 = affect; 0 = does not affect). |


### Output Files
When simulating different runs on *Automatic Runs*, our program generates output files informing important info about the runs. Each output file presents different data about the runs:

| **Output** | **Data presented** |
| --- | --- |
| \_seed\_ancestralResults.csv | Returns the number of descendants cells for each initial cell per generation. |
| \_seed\_eventsMutationResults.csv | Returns the number of cells affected on each event in one or two alleles in each generation. |
| \_seed\_genesMutationResults.csv | Returns the number of cells affected on each gene in one or two alleles in each generation. |
| \_seed\_parameters.txt | Returns the input parameters (seed, Number of Cells, Proliferation Rate, Death Rate, Maximum Division, Events per Division). |
| \_automatics\_ancestralResults.csv | Returns the number of descendants cells for each initial cell in the final population for each seed. Also, returns the total number of divisions for each seed. |
| \_automatics\_eventsMutationResults.csv | Returns the number of cells affected on each event in one or two alleles in the final population for each seed.|
| \_automatics\_genesMutationResults.csv | Returns the number of cells affected on each gene in one or two alleles in the final population for each seed.|
| \_automatics\_sequenceEachCell.csv | Returns the single cell genes in the two alleles as two independent lists in a binary form (1 = affected; 0 = not affected). |

### Documentation
We haven't created a documentation for all classes and fucntions of the code, but a friendly documentation can be found [here](./esiCancerDocumentation.pdf).

### Running into issues?
Contact Bernardo Henz <bernardohenz@gmail.com>

Other informations can be obtained from [www.ufrgs.br/labsinal/esiCancer](https://www.ufrgs.br/labsinal/esiCancer).
