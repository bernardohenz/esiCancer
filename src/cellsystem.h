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
    unsigned long genome_size;
    //Tumor Environment
    float stdGrowthRate;

    MutationTable *myMutationTable;
    //Manual Mutations
    std::vector<MutationManual *> rulesMutationManual;


    //Log Information
    std::vector<unsigned int> populationHistory;
    //Mutations
    std::vector<unsigned int> historyNumberOfAffectedCells;
    std::vector<MutationHistogram> historyMutationGeneHistogram;

    std::vector<MutationCounter> historyMutationGeneCounters;
    std::vector<MutationCounter> historyFirstTapeGeneCounters;
    std::vector<MutationCounter> historySecondTapeGeneCounters;

    std::vector<MutationCounter> historyMutationEventCounters;
    std::vector<MutationCounter> historyFirstTapeEventCounters;
    std::vector<MutationCounter> historySecondTapeEventCounters;
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
    void addMutationEventToNewCell(int tmpMutationEventId,int tmpMutationGeneId, std::vector<unsigned int> &newFirstTapeMutationGenes, std::vector<unsigned int> &newSecondTapeMutationGenes,
                                   std::vector<unsigned int> &newFirstTapeMutationEvents, std::vector<unsigned int> &newSecondTapeMutationEvents,
                                   std::vector<unsigned int> &newMutations, std::vector<unsigned int> &newMutationEvents,int tTapeChoosen,
                                   float &proliferationChild, float &deathChild, float &tmpMutationsPerRepChild,
                                   float &telomeresChild, float &tmpMicroEnvironment);
    Cell *addNewCell(Cell * father);

    //Simulation
    void reset();
    float getRandomProb();
    unsigned long getRandomInteger();
    long getRandomIntegerInRange(int max);
    bool process();
    void startSimulation();
    int killAtRandom(int numberOfCellsToDie);

    //Load Params
    void loadMutationTableFile(QString filename);
    void loadSinergyTableFile(QString filename);
    void loadDefaultValues(float tprobProl,float tprobDeath,int tTelomeres,int mutPerDiv);
    void loadSimulationParameters(int tseed, int tnumberOfCells);
    void loadStopConditions(int tgenerations,int tCells,int tMutatedCells);

    //Manual Mutations
    void addMutationRule(MutationManual *tMut);
    void updateMutationRule(int id, unsigned int generation, int percentage, QString mutationFullname, int mutationId, manualMutationAllelsChanged tallelschanged);
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
    void exportMutationGeneCounters(QString tname);
    void exportMutationEventCounters(QString tname);
    void exportMutationHistograms(QString tname);
    void exportAncestralResults(QString tname);
    void exportMutationEventCountersLastLine(QString tname);
    void exportMutationGeneCountersLastLine(QString tname);
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
    void setStdGrowthRate(float trate);
};

#endif // CELLSYSTEM_H
