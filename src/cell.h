#ifndef CELL_H
#define CELL_H

#include <vector>
#include <string>

class CellSystem;
class MutationTable;

class Cell
{
private:
    CellSystem *mySystem;

    bool alive;
    float mProbDeath;
    float mProbReproduction;
    float mMutationsPerRep;

    long id;
    unsigned int ancestor;
    unsigned int numberOfReproductions;
    std::vector<long> parentsInOrder;

    unsigned short int lifeTime;
    unsigned int mTelomeres;

    //List of mutations in order
    std::vector<unsigned int> firstTapeMutationGenes;
    std::vector<unsigned int> secondTapeMutationGenes;

    std::vector<unsigned int> firstTapeMutationEvents;
    std::vector<unsigned int> secondTapeMutationEvents;

    std::vector<unsigned int> mutationsInOrder; //only fullMutations
    std::vector<unsigned int> mutationEventsInOrder;

    void setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres);
    void setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres,int maxIdMutation);
public:
    Cell(long tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres,int maxIdMutation);
    Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, int maxIdMutation);
    Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, std::vector<unsigned int> tFirstTapeMutationGenes, std::vector<unsigned int> tSecondTapeMutationGenes, std::vector<long> tparentsOrder);
    Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder,std::vector<unsigned int> tmutationEventsInOrder,std::vector<unsigned int> tFirstTapeMutationGenes, std::vector<unsigned int> tSecondTapeMutationGenes, std::vector<unsigned int> tFirstTapeMutationEvents, std::vector<unsigned int> tSecondTapeMutationEvents, std::vector<long> tparentsOrder);

    void registerCellSystem(CellSystem *tsystem);

    //Getters & Setters
    float getProbDeath() const;
    void setProbDeath(float probDeath);
    float getProbReproduction() const;
    void setProbReproduction(float probReproduction);
    std::vector<unsigned int> getMutationsInOrder() const;
    void setMutationsInOrder(const std::vector<unsigned int> &value);
    std::vector<unsigned int> getMutationEventsInOrder() const;
    void setMutationEventsInOrder(const std::vector<unsigned int> &value);

    std::vector<unsigned int> getFirstTapeMutationGenes() const;
    void setFirstTapeMutationGenes(const std::vector<unsigned int> &value);
    std::vector<unsigned int> getSecondTapeMutationGenes() const;
    void setSecondTapeMutationGenes(const std::vector<unsigned int> &value);

    std::vector<unsigned int> getFirstTapeMutationEvents() const;
    void setFirstTapeMutationEvents(const std::vector<unsigned int> &value);
    std::vector<unsigned int> getSecondTapeMutationEvents() const;
    void setSecondTapeMutationEvents(const std::vector<unsigned int> &value);

    void addInFirstTapeMutationEvents(unsigned int tmutevent);
    void addInSecondTapeMutationEvents(unsigned int tmutevent);

    float getMutationsPerRep() const;
    void setMutationsPerRep(float mutationsPerRep);
    unsigned int getTelomeres() const;
    void setTelomeres(unsigned int telomeres);
    std::vector<long> getParentsInOrder() const;
    void setParentsInOrder(const std::vector<long> &value);
    long getId() const;
    void setId(long value);
    unsigned int getAncestor() const;
    void setAncestor(unsigned int value);
    unsigned int getNumberOfReproductions() const;
    void setNumberOfReproductions(unsigned int value);
    void die();

    void addMutation(unsigned int mutid);

    bool isAlive() const;
    bool reproduce();
    bool hasThisMut(int tMutId) const;
    bool hasThisMutGeneOnTape(int tMutId,int tTape) const;
    bool hasThisMutEventOnTape(int tMutId,int tTape) const;




};

#endif // CELL_H
