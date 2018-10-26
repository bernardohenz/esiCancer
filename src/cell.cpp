#include "cell.h"
#include "cellsystem.h"
#include <QDebug>
#include "mutationtable.h"


float Cell::getMutationsPerRep() const
{
    return mMutationsPerRep;
}

void Cell::setMutationsPerRep(float mutationsPerRep)
{
    mMutationsPerRep = mutationsPerRep;
}

unsigned int Cell::getTelomeres() const
{
    return mTelomeres;
}

void Cell::setTelomeres(unsigned int telomeres)
{
    mTelomeres = telomeres;
}

std::vector<long> Cell::getParentsInOrder() const
{
    return parentsInOrder;
}

void Cell::setParentsInOrder(const std::vector<long> &value)
{
    parentsInOrder = value;
}


bool Cell::isAlive() const
{
    return alive;
}

long Cell::getId() const
{
    return id;
}

void Cell::setId(long value)
{
    id = value;
}

unsigned int Cell::getAncestor() const
{
    return ancestor;
}

void Cell::setAncestor(unsigned int value)
{
    ancestor = value;
}

unsigned int Cell::getNumberOfReproductions() const
{
    return numberOfReproductions;
}

void Cell::setNumberOfReproductions(unsigned int value)
{
    numberOfReproductions = value;
}

void Cell::die()
{
    alive = false;
}

void Cell::addMutation(unsigned int mutid)
{
    mutationsInOrder.push_back(mutid);
}


void Cell::setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres)
{
    id = tid;
    mProbReproduction = tprobOfRep;
    mProbDeath = tprobOfDeath;
    mMutationsPerRep = tmutationPerRep;
    mTelomeres = ttelomeres;
    alive=true;
    numberOfReproductions=0;
    lifeTime=0;
    ancestor=id;
}

void Cell::setupParams(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, int maxIdMutation)
{
    setupParams(tid,tprobOfRep,tprobOfDeath,tmutationPerRep,ttelomeres);
    firstTapeMutationGenes = std::vector<unsigned int>();//(maxIdMutation,false);
    secondTapeMutationGenes = std::vector<unsigned int>();//(maxIdMutation,false);
    firstTapeMutationEvents = std::vector<unsigned int>();//(maxIdMutation,false);
    secondTapeMutationEvents = std::vector<unsigned int>();//(maxIdMutation,false);
    mutationsInOrder = std::vector<unsigned int>();
    parentsInOrder = std::vector<long>();
}

Cell::Cell(long tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, int maxIdMutation)
{
    setupParams(tid,tprobOfRep,tprobOfDeath,tmutationPerRep,ttelomeres,maxIdMutation);
}

Cell::Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder,int maxIdMutation)
{
    setupParams(tid,tprobOfRep,tprobOfDeath,tmutationPerRep,ttelomeres,maxIdMutation);
    //initialize mutations
    mutationsInOrder =tmutationsInOrder;

}

Cell::Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, std::vector<unsigned int> tFirstTapeMutationGenes, std::vector<unsigned int> tSecondTapeMutationGenes, std::vector<long> tparentsOrder)
{
    setupParams(tid,tprobOfRep,tprobOfDeath,tmutationPerRep,ttelomeres);
    //initialize mutations
    mutationsInOrder =tmutationsInOrder;
    firstTapeMutationGenes = tFirstTapeMutationGenes;
    secondTapeMutationGenes = tSecondTapeMutationGenes;
    parentsInOrder = tparentsOrder;
}

Cell::Cell(float tid, float tprobOfRep, float tprobOfDeath, float tmutationPerRep, unsigned int ttelomeres, std::vector<unsigned int> tmutationsInOrder, std::vector<unsigned int> tmutationEventsInOrder, std::vector<unsigned int> tFirstTapeMutationGenes, std::vector<unsigned int> tSecondTapeMutationGenes, std::vector<unsigned int> tFirstTapeMutationEvents, std::vector<unsigned int> tSecondTapeMutationEvents, std::vector<long> tparentsOrder)
{
    setupParams(tid,tprobOfRep,tprobOfDeath,tmutationPerRep,ttelomeres);
    //initialize mutations
    mutationsInOrder =tmutationsInOrder;
    firstTapeMutationGenes = tFirstTapeMutationGenes;
    secondTapeMutationGenes = tSecondTapeMutationGenes;

    mutationEventsInOrder = tmutationEventsInOrder;
    firstTapeMutationEvents = tFirstTapeMutationEvents;
    secondTapeMutationEvents = tSecondTapeMutationEvents;
    parentsInOrder = tparentsOrder;

}

void Cell::registerCellSystem(CellSystem *tsystem)
{
    mySystem = tsystem;
    /*
    if(ancestor==mySystem->getDominantGeneAncestor())
        mTag = QString::number(ancestor+1).toStdString();
    else
        mTag = "";
    */
}


bool Cell::reproduce()
{
    lifeTime++;
    if(mySystem==NULL || !alive)
        return false;

    float tmpProbReproduction = mProbReproduction;
    float tmpProbDeath = mProbDeath;
    if (mySystem->getRandomProb() < tmpProbReproduction){
        if (mTelomeres>0){
            mTelomeres--;
            //New system
            mySystem->addNewCell(this);
            mySystem->addNewCell(this);
            //numberOfReproductions++;
            alive = false;
            mySystem->incrementNumberOfDivisions();
        }
    }
    if(mySystem->getRandomProb() < tmpProbDeath){
        //mySystem->informNaturalKilled(id);
        alive = false;
    }
    return alive;
}

bool Cell::hasThisMut(int tMutId) const
{
    bool hasIt=false;
    for (unsigned int i=0; i<mutationsInOrder.size(); i++)
        if (mutationsInOrder[i]==(unsigned int)tMutId){
            hasIt=true;
            break;
        }
    return hasIt;
}

bool Cell::hasThisMutGeneOnTape(int tMutId, int tTape) const
{
    bool hasIt=false;
    if (tTape==0){
        for (unsigned int i=0; i<firstTapeMutationGenes.size(); i++)
            if (firstTapeMutationGenes[i]==(unsigned int) tMutId){
                hasIt=true;
                break;
            }
    } else if(tTape==1){
        for (unsigned int i=0; i<secondTapeMutationGenes.size(); i++)
            if (secondTapeMutationGenes[i]==(unsigned int) tMutId){
                hasIt=true;
                break;
            }
    }
    return hasIt;
}

bool Cell::hasThisMutEventOnTape(int tMutId, int tTape) const
{
    bool hasIt=false;
    if (tTape==0){
        for (unsigned int i=0; i<firstTapeMutationEvents.size(); i++)
            if (firstTapeMutationEvents[i]==(unsigned int) tMutId){
                hasIt=true;
                break;
            }
    } else if(tTape==1){
        for (unsigned int i=0; i<secondTapeMutationEvents.size(); i++)
            if (secondTapeMutationEvents[i]==(unsigned int) tMutId){
                hasIt=true;
                break;
            }
    }
    return hasIt;
}

float Cell::getProbDeath() const
{
    return mProbDeath;
}

void Cell::setProbDeath(float probDeath)
{
    if(probDeath>=0)
        mProbDeath = probDeath;
}

float Cell::getProbReproduction() const
{
    return mProbReproduction;
}

void Cell::setProbReproduction(float probReproduction)
{
    mProbReproduction = probReproduction;
}


std::vector<unsigned int> Cell::getMutationsInOrder() const
{
    return mutationsInOrder;
}

void Cell::setMutationsInOrder(const std::vector<unsigned int> &value)
{
    mutationsInOrder = value;
}

std::vector<unsigned int> Cell::getMutationEventsInOrder() const
{
    return mutationEventsInOrder;
}

void Cell::setMutationEventsInOrder(const std::vector<unsigned int> &value)
{
    mutationEventsInOrder = value;
}


std::vector<unsigned int> Cell::getFirstTapeMutationGenes() const
{
    return firstTapeMutationGenes;
}

void Cell::setFirstTapeMutationGenes(const std::vector<unsigned int> &value)
{
    firstTapeMutationGenes = value;
}

std::vector<unsigned int> Cell::getSecondTapeMutationGenes() const
{
    return secondTapeMutationGenes;
}

void Cell::setSecondTapeMutationGenes(const std::vector<unsigned int> &value)
{
    secondTapeMutationGenes = value;
}


std::vector<unsigned int> Cell::getFirstTapeMutationEvents() const
{
    return firstTapeMutationEvents;
}

void Cell::setFirstTapeMutationEvents(const std::vector<unsigned int> &value)
{
    firstTapeMutationEvents = value;
}

std::vector<unsigned int> Cell::getSecondTapeMutationEvents() const
{
    return secondTapeMutationEvents;
}

void Cell::setSecondTapeMutationEvents(const std::vector<unsigned int> &value)
{
    secondTapeMutationEvents = value;
}

void Cell::addInFirstTapeMutationEvents(unsigned int tmutevent)
{
    //if (tmutevent<firstTapeMutationEvents.size())
    if(!Contains(firstTapeMutationEvents,tmutevent))
        firstTapeMutationEvents.push_back(tmutevent);
}

void Cell::addInSecondTapeMutationEvents(unsigned int tmutevent)
{
    //if (tmutevent<secondTapeMutationEvents.size())
    if(!Contains(secondTapeMutationEvents,tmutevent))
        secondTapeMutationEvents.push_back(tmutevent);
}
