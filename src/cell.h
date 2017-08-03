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
    bool hasSuppressorMut;
    float mProbDeath;
    float mProbReproduction;
    float mMutationsPerRep;
    float sinergy;

    long id;
    unsigned int ancestor;
    unsigned int numberOfReproductions;
    std::vector<long> parentsInOrder;

    unsigned short int lifeTime;
    unsigned int mTelomeres;

    //Mutations, size = number of mutations
    //Returning to ints
    std::vector<unsigned int> firstTapeMutations;
    std::vector<unsigned int> secondTapeMutations;

    std::vector<unsigned int> mutationsInOrder;

    void setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres);
    void setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres,int maxIdMutation);
public:
    Cell(long tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres,int maxIdMutation);
    Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, int maxIdMutation);
    Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, std::vector<unsigned int> tFirstTapeMutations, std::vector<unsigned int> tSecondTapeMutations, std::vector<long> tparentsOrder);

    void registerCellSystem(CellSystem *tsystem);

    //Getters & Setters
    float getProbDeath() const;
    void setProbDeath(float probDeath);
    float getProbReproduction() const;
    void setProbReproduction(float probReproduction);
    std::vector<unsigned int> getMutationsInOrder() const;
    void setMutationsInOrder(const std::vector<unsigned int> &value);
    std::vector<unsigned int> getFirstTapeMutations() const;
    void setFirstTapeMutations(const std::vector<unsigned int> &value);
    std::vector<unsigned int> getSecondTapeMutations() const;
    void setSecondTapeMutations(const std::vector<unsigned int> &value);
    float getMutationsPerRep() const;
    void setMutationsPerRep(float mutationsPerRep);
    unsigned int getTelomeres() const;
    void setTelomeres(unsigned int telomeres);
    std::vector<long> getParentsInOrder() const;
    void setParentsInOrder(const std::vector<long> &value);
    bool getHasSuppressorMut() const;
    void setHasSuppressorMut(bool value);
    long getId() const;
    void setId(long value);
    float getSinergy() const;
    void setSinergy(float value);
    unsigned int getAncestor() const;
    void setAncestor(unsigned int value);
    unsigned int getNumberOfReproductions() const;
    void setNumberOfReproductions(unsigned int value);

    void addMutation(unsigned int mutid);


    bool isAlive() const;
    bool reproduce();
    bool hasThisMut(int tMutId) const;


};

#endif // CELL_H
