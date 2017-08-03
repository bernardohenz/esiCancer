#ifndef CELLSYSTEM_H
#define CELLSYSTEM_H

#include <vector>
#include <queue>
#include <random>
#include <QString>
#include "mutationhistogram.h"
#include "mutationcounter.h"
#include "ancestorhistogram.h"
#include "mutationmanual.h"

class Cell;
class MutationTable;



class CellSystem
{
private:
    //Overal management
    unsigned int generation;

    //Cell Collection
    float stdProliferation,stdDeath;
    unsigned int startingNumberOfCells,stdMutationsPerDivision,stdTelomeres;
    unsigned int maxMutationsPerDivision;
    float maxProliferation, maxDeath;
    std::vector<Cell *> cellVec;
    std::queue<unsigned int> freePositions;
    unsigned int lastPositionCell;
    unsigned int generationLastPositionCell;
    unsigned int numberOfCells;
    long counterId;
    unsigned int numberOfDivisions;

    //Probability control
    float counterCellId;
    int seed;
    std::mt19937 randomGenerator;
    std::default_random_engine probabilityGenerator;
    std::uniform_real_distribution<double> probabilityDistribution;
    QString debugRandomString;
    long genome_size;


    MutationTable *myMutationTable;
    //Manual Mutations
    std::vector<MutationManual *> rulesMutationManual;


    //Log Information
    std::vector<unsigned int> populationHistory;
    //Mutations
    std::vector<unsigned int> historyNumberOfAffectedCells;
    std::vector<MutationHistogram> historyMutationHistogram;
    std::vector<MutationCounter> historyMutationCounters;
    std::vector<AncestorHistogram> historyAncestorsCounters;

    //Stop Conditions
    int stopGenerations;
    int stopNumberOfCells;
    int stopMutatedCells;



    //Ancestor
    unsigned int stdNumberOfAncestors;
public:
    CellSystem();

    //Cell Management
    Cell *addNewCell(Cell * father);
    void updateOncsNowTumorSuppressor(Cell *tcell);
    void informNaturalKilled(long id);

    //Simulation
    void reset();
    float getRandomProb();
    long getRandomInteger();
    bool process();
    void startSimulation();

    //Load Params
    void loadMutationTableFile(QString filename);
    void loadDefaultValues(float tprobProl,float tprobDeath,int tTelomeres,int mutPerDiv);
    void loadSimulationParameters(int tseed, int tnumberOfCells);
    void loadStopConditions(int tgenerations,int tCells,int tMutatedCells);

    //Manual Mutations
    void addMutationRule(MutationManual *tMut);
    void updateMutationRule(int id, unsigned int generation, int percentage, QString mutationname, int mutationId, manualMutationAllelsChanged tallelschanged);
    void removeMutationRule(int index);
    void removeAllRules();
    MutationManual *getMutationManualFromIndex(int id) const;
    void processRuleOfManualMutations(unsigned int id);

    //Interface
    QStringList getMutationListNames() const;

    //Export
    void exportInfo(QString tname);
    void exportLastLine(QString tname);
    void exportInputParams(QString tname);
    void exportMutationCounters(QString tname);
    void exportMutationHistograms(QString tname);
    void exportAncestralResults(QString tname);
    void exportMutationCountersLastLine(QString tname);
    void exportMutationHistogramsLastLine(QString tname);
    void exportAncestralResultsLastLine(QString tname);
    void exportMutsEachCellLastLine(QString tname);

    void show();
    int getStopGenerations() const;
    void setStopGenerations(int value);
    unsigned int getGeneration() const;
    void setGeneration(unsigned int value);
    unsigned int getMaxMutationsPerDivision() const;
    void setMaxMutationsPerDivision(unsigned int value);
    void setMaxProliferation(float value);
    void setMaxDeath(float value);
    void incrementNumberOfDivisions();
};

#endif // CELLSYSTEM_H
