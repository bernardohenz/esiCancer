#include "cellsystem.h"
#include "cell.h"
#include "mutationtable.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include "mutationmanual.h"
#include <QFileInfo>

bool fileExists(QString path) {
    QFileInfo check_file(path);
    // check if file exists and if yes: Is it really a file and no directory?
    if (check_file.exists() && check_file.isFile()) {
        return true;
    } else {
        return false;
    }
}

int CellSystem::getStopGenerations() const
{
    return stopGenerations;
}

void CellSystem::setStopGenerations(int value)
{
    stopGenerations = value;
}

unsigned int CellSystem::getGeneration() const
{
    return generation;
}

void CellSystem::setGeneration(unsigned int value)
{
    generation = value;
}

unsigned int CellSystem::getMaxMutationsPerDivision() const
{
    return maxMutationsPerDivision;
}

void CellSystem::setMaxMutationsPerDivision(unsigned int value)
{
    maxMutationsPerDivision = value;
}

void CellSystem::setMaxProliferation(float value)
{
    maxProliferation = value;
}

void CellSystem::setMaxDeath(float value)
{
    maxDeath = value;
}

void CellSystem::incrementNumberOfDivisions()
{
    numberOfDivisions++;
}

void CellSystem::setStdGrowthRate(float trate)
{
    stdGrowthRate = trate;
}

CellSystem::CellSystem()
{
    probabilityDistribution = std::uniform_real_distribution<double> (0.0,1.0);

    myMutationTable = new MutationTable();
    startingNumberOfCells = 1000;
    reset();
    //default values
    seed=13;
    genome_size = 4568;
    stdMutationsPerDivision = 20;
    stdProliferation = 0.01;
    stdDeath = 0.01;
    stdTelomeres = 30;
    stopGenerations = stopNumberOfCells = stopMutatedCells = -1;
    maxMutationsPerDivision = 2000;
    maxProliferation = maxDeath = 1;
    stdGrowthRate = -1;
}


void CellSystem::reset()
{
    generation=1;
    generationLastPositionCell =0;
    counterId=0;
    numberOfDivisions=0;

    //cell informations
    for(unsigned int i=0;i<cellVec.size();i++)
        if(cellVec[i] != NULL)
            delete cellVec[i];
    cellVec.clear();
    while(!freePositions.empty())
        freePositions.pop();

    //random number generators
    randomGenerator = std::mt19937((unsigned int)seed);
    probabilityGenerator = std::default_random_engine((unsigned int)seed);
    probabilityDistribution.reset();

    //log containers
    populationHistory.clear();


    //mutations
    historyNumberOfAffectedCells.clear();
    historyNumberOfAffectedCells.push_back(0);
    historyMutationGeneHistogram.clear();
    historyMutationGeneHistogram.push_back(MutationHistogram(myMutationTable->getMaxIdMutationGene()));

    historyMutationGeneCounters.clear();
    historyMutationGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));
    historyFirstTapeGeneCounters.clear();
    historyFirstTapeGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));
    historySecondTapeGeneCounters.clear();
    historySecondTapeGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));

    historyMutationEventCounters.clear();
    historyMutationEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));
    historyFirstTapeEventCounters.clear();
    historyFirstTapeEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));
    historySecondTapeEventCounters.clear();
    historySecondTapeEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));

    stdNumberOfAncestors = startingNumberOfCells;
    historyAncestorsCounters.clear();
    historyAncestorsCounters.push_back(AncestorHistogram(stdNumberOfAncestors));
}

void CellSystem::addMutationEventToNewCell(int tmpMutationEventId,int tmpMutationGeneId, std::vector<unsigned int> &newFirstTapeMutationGenes, std::vector<unsigned int> &newSecondTapeMutationGenes,
                               std::vector<unsigned int> &newFirstTapeMutationEvents, std::vector<unsigned int> &newSecondTapeMutationEvents,
                               std::vector<unsigned int> &newMutations, std::vector<unsigned int> &newMutationEvents,int tTapeChoosen,
                               float &proliferationChild, float &deathChild, float &tmpMutationsPerRepChild,
                               float &telomeresChild, float &tmpMicroEnvironment){
    if(tmpMutationEventId>-1){
        int activationLevel = 0; //0=noActivation 1=partialActivation 2=fullActivation  3=fullActivationFullDominance
        int tapeChosen;
        if (tTapeChoosen>0)
            tapeChosen=tTapeChoosen;
        else {
            if (getRandomProb()>0.5)
                tapeChosen = 1;
            else
                tapeChosen = 2;
        }
        if ((tapeChosen==1) || (tapeChosen==3)){
            activationLevel = myMutationTable->checkIfAddMutationByReference(tmpMutationGeneId,newFirstTapeMutationGenes,newSecondTapeMutationGenes,1);
            if(!Contains(newFirstTapeMutationGenes,(unsigned int)tmpMutationGeneId))
                newFirstTapeMutationGenes.push_back(tmpMutationGeneId);
            if (activationLevel>0) {
                if(!Contains(newFirstTapeMutationEvents,(unsigned int)tmpMutationEventId))
                    newFirstTapeMutationEvents.push_back(tmpMutationEventId);
                if(activationLevel==2){
                    if(!Contains(newSecondTapeMutationEvents,(unsigned int)tmpMutationEventId)) //accounting for different events of the same gene
                        activationLevel=3;
                }
            }

            //newFirstTapeMutations[tmpMutationId] = true;
        }
        if ((tapeChosen==2)||(tapeChosen==3)) {
            activationLevel += myMutationTable->checkIfAddMutationByReference(tmpMutationGeneId,newFirstTapeMutationGenes,newSecondTapeMutationGenes,2);
            if(!Contains(newSecondTapeMutationGenes,(unsigned int)tmpMutationGeneId))
                newSecondTapeMutationGenes.push_back(tmpMutationGeneId);
            if (activationLevel>0) {
                if(!Contains(newSecondTapeMutationEvents,(unsigned int)tmpMutationEventId))
                    newSecondTapeMutationEvents.push_back(tmpMutationEventId);
                if(activationLevel==2){
                    if(!Contains(newFirstTapeMutationEvents,(unsigned int)tmpMutationEventId)) //accounting for different events of the same gene
                        activationLevel=3;
                }
            }
            //newSecondTapeMutations[tmpMutationId] = true;
        }
        if(Contains(newFirstTapeMutationEvents,(unsigned int)tmpMutationEventId) && Contains(newSecondTapeMutationEvents,(unsigned int)tmpMutationEventId) && !Contains(newMutationEvents,(unsigned int)tmpMutationEventId))
            newMutationEvents.push_back(tmpMutationEventId);

        if (activationLevel>0){
            MutationEvent *activatedMutation = myMutationTable->getMutationEvent(tmpMutationEventId);
            proliferationChild = activatedMutation->computeNewProlRate(proliferationChild,activatedMutation->computeProlSinergyModifier(newFirstTapeMutationEvents,newSecondTapeMutationEvents),activationLevel);
            deathChild = activatedMutation->computeNewDeathRate(deathChild,activatedMutation->computeDeathSinergyModifier(newFirstTapeMutationEvents,newSecondTapeMutationEvents),activationLevel);
            tmpMutationsPerRepChild = activatedMutation->computeNewMoreMutRate(tmpMutationsPerRepChild,1);
            telomeresChild = activatedMutation->computeNewTelomeresRate(telomeresChild,activatedMutation->computeTelSinergyModifier(newFirstTapeMutationEvents,newSecondTapeMutationEvents),activationLevel);
            tmpMicroEnvironment = activatedMutation->computeNewMicroEnvironmentModifier(tmpMicroEnvironment,activationLevel);

            if (activationLevel>=2){ //mutations on both tapes
                if (!Contains(newMutations,(unsigned int)tmpMutationGeneId))
                    newMutations.push_back(tmpMutationGeneId);
            }
            for (int i=0;i<activatedMutation->listOfChainReactionlMutationEvents.size(); i++){
                unsigned int cur_eventId = activatedMutation->listOfChainReactionlMutationEvents[i];
                unsigned int cur_geneId = myMutationTable->getGeneOfEvent(cur_eventId);
                addMutationEventToNewCell(cur_eventId,cur_geneId,newFirstTapeMutationGenes,newSecondTapeMutationGenes,
                                 newFirstTapeMutationEvents,newSecondTapeMutationEvents,
                                 newMutations,newMutationEvents,tapeChosen,proliferationChild,deathChild,tmpMutationsPerRepChild,
                                 telomeresChild,tmpMicroEnvironment);
            }

        }

    }
}


Cell *CellSystem::addNewCell(Cell *father){
    unsigned int tmpMutatedPosition;
    int tmpMutationEventId,tmpMutationGeneId;
    //copy from father
    std::vector<unsigned int> newMutations = father->getMutationsInOrder();
    std::vector<unsigned int> newMutationEvents = father->getMutationEventsInOrder();
    std::vector<unsigned int> newFirstTapeMutationGenes = father->getFirstTapeMutationGenes();
    std::vector<unsigned int> newSecondTapeMutationGenes = father->getSecondTapeMutationGenes();
    std::vector<unsigned int> newFirstTapeMutationEvents = father->getFirstTapeMutationEvents();
    std::vector<unsigned int> newSecondTapeMutationEvents = father->getSecondTapeMutationEvents();
    std::vector<long> newParentsOrder = father->getParentsInOrder();

    float proliferationChild = father->getProbReproduction();
    float deathChild = father->getProbDeath();
    float telomeresChild = father->getTelomeres();
    float tmpMutationsPerRepChild = father->getMutationsPerRep();
    float tmpMutationsPerRepFather = father->getMutationsPerRep();
    float tmpMicroEnvironment = 1;
//    if (generation>520){
//        qDebug()<<"Adding New Cell:  mutPerRep: "<<tmpMutationsPerRepFather<<"    mutations: "<<newMutations.size()<<"    Muts:";
//        for (int i=0;i<newMutations.size();i++)
//            qDebug()<<"   Mut "<<newMutations[i];
//    }
    for(unsigned int i=0; i<tmpMutationsPerRepFather; i++){
        //get a random basis position
        tmpMutatedPosition = getRandomInteger();

        //get the mutation id corresponding to the mutated base
        tmpMutationEventId = myMutationTable->getNumberofMutationEvent(tmpMutatedPosition);
        tmpMutationGeneId = myMutationTable->getNumberOfMutationGene(tmpMutatedPosition);
        addMutationEventToNewCell(tmpMutationEventId,tmpMutationGeneId,newFirstTapeMutationGenes,newSecondTapeMutationGenes,
                         newFirstTapeMutationEvents,newSecondTapeMutationEvents,
                         newMutations,newMutationEvents,-1,proliferationChild,deathChild,tmpMutationsPerRepChild,
                         telomeresChild,tmpMicroEnvironment);
    }
    if (stdGrowthRate>0){
        stdGrowthRate *= tmpMicroEnvironment;
    }
    newParentsOrder.push_back(father->getId());
    tmpMutationsPerRepChild = std::min((unsigned int)tmpMutationsPerRepChild,maxMutationsPerDivision);
    proliferationChild = std::min(proliferationChild,maxProliferation);
    deathChild = std::min(deathChild,maxDeath);
    //qDebug()<<"Creating Cell";
    Cell *newCell = new Cell(counterCellId,proliferationChild, deathChild, tmpMutationsPerRepChild,
                                 telomeresChild, newMutations,newMutationEvents, newFirstTapeMutationGenes, newSecondTapeMutationGenes,newFirstTapeMutationEvents, newSecondTapeMutationEvents, newParentsOrder);
    newCell->setAncestor(father->getAncestor());

    newCell->registerCellSystem(this);
    //qDebug()<<"Cell created and registered";
    counterCellId++;
    if (freePositions.empty()){
        if (lastPositionCell == cellVec.size()){
            cellVec.resize(cellVec.size()+300,NULL);
        }
        cellVec[lastPositionCell] = newCell;
        lastPositionCell++;
    } else{
        cellVec[freePositions.front()] = newCell;
        freePositions.pop();
    }
    numberOfCells++;
    //qDebug()<<"Memory handled";
    return newCell;
}


float CellSystem::getRandomProb()
{
    return (float)probabilityDistribution(probabilityGenerator);
    //float a = (float)probabilityDistribution(probabilityGenerator);
    //return a;
}

unsigned long CellSystem::getRandomInteger()
{
    unsigned long randomNumber = ((randomGenerator)());
    //qDebug()<<randomNumber%genome_size;
    return randomNumber%genome_size;
    //long randomNumber = (randomGenerator)()%genome_size;
    //return randomNumber;
}

long CellSystem::getRandomIntegerInRange(int max)
{
    return (randomGenerator)()%max;
}

bool CellSystem::process()
{
//    qDebug()<<"Start process";
    for(unsigned int i=0; i< rulesMutationManual.size(); i++) {
        if(generation == rulesMutationManual[i]->getGeneration()){
            processRuleOfManualMutations(i);
        }
    }

    int tmpCurrentDeathHistory=0;
    for(unsigned int i=0; i<generationLastPositionCell;i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
//            if(generation>520)
//                qDebug()<<"b4 reproduce"<<i;
            if(!cellVec[i]->reproduce())
                tmpCurrentDeathHistory++;
//            if(generation>520)
//                qDebug()<<"after reproduce"<<i;
        }
        if(!(cellVec[i]->isAlive())){
            if (numberOfCells>0)
                numberOfCells--;
            delete cellVec[i];
            cellVec[i]=NULL;
            freePositions.push(i);
        }
    }
//    qDebug()<<"minorDebug";
    generation++;
    generationLastPositionCell = lastPositionCell;


    //Check Tumor Environment - if enabled
    if (stdGrowthRate>0){
        float cur_rate;
        if (populationHistory.empty())
            cur_rate = 1;
        else
            cur_rate = (float)numberOfCells/populationHistory.back();
        //rate is more than expected

        //qDebug()<<"num_cell: "<<numberOfCells<<"   |last pop  "<<populationHistory.back()<<"   Cur_RATE: "<<cur_rate;
        if(cur_rate>stdGrowthRate){
            int numberOfCellsThatMustDie = (int)numberOfCells- (int)populationHistory.back()*stdGrowthRate;
            //qDebug()<<"num_cell: "<<numberOfCells<<"  cells that must die: "<<numberOfCellsThatMustDie;
            int killed = killAtRandom(numberOfCellsThatMustDie);
            //qDebug()<<"KILLED: "<<killed;
            //qDebug()<<"num_cell after killing: "<<numberOfCells;
        }

    }

    //Logs recording
    std::vector<unsigned int> tmpMuts, tmpTape;
    int affectedCells=0;

    historyMutationGeneHistogram.push_back(MutationHistogram(myMutationTable->getMaxIdMutationGene()));
    historyMutationGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));
    historyFirstTapeGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));
    historySecondTapeGeneCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationGene()));

    historyMutationEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));
    historyFirstTapeEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));
    historySecondTapeEventCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutationEvent()));
    historyAncestorsCounters.push_back(AncestorHistogram(stdNumberOfAncestors));
    int tmpcounter=0;
    int affectedFirstStrand = 0;
    int affectedSecondStrand = 0;
    for(unsigned int i=0;i<generationLastPositionCell;i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            tmpcounter++;
            tmpMuts = cellVec[i]->getMutationsInOrder();
            if(tmpMuts.size()>0){
                //affectedCells++;
                for(unsigned int mutI=0;mutI<tmpMuts.size();mutI++)
                    historyMutationGeneCounters.back().incrementMutCounter(tmpMuts[mutI]);
                historyMutationGeneHistogram.back().incrementBin(tmpMuts.size());
            }

            tmpTape = cellVec[i]->getFirstTapeMutationGenes();
            if(tmpTape.size()>0){
                //affectedFirstStrand++;
                for(unsigned int mutI=0;mutI<tmpTape.size();mutI++)
                    historyFirstTapeGeneCounters.back().incrementMutCounter(tmpTape[mutI]);
            }

            tmpTape = cellVec[i]->getSecondTapeMutationGenes();
            if(tmpTape.size()>0){
                //affectedSecondStrand++;
                for(unsigned int mutI=0;mutI<tmpTape.size();mutI++)
                    historySecondTapeGeneCounters.back().incrementMutCounter(tmpTape[mutI]);
            }
            tmpMuts = cellVec[i]->getMutationEventsInOrder();
            if(tmpMuts.size()>0){
                affectedCells++;
                for(unsigned int mutI=0;mutI<tmpMuts.size();mutI++)
                    historyMutationEventCounters.back().incrementMutCounter(tmpMuts[mutI]);
            }

            tmpTape = cellVec[i]->getFirstTapeMutationEvents();
            if(tmpTape.size()>0){
                affectedFirstStrand++;
                for(unsigned int mutI=0;mutI<tmpTape.size();mutI++)
                    historyFirstTapeEventCounters.back().incrementMutCounter(tmpTape[mutI]);
            }
            tmpTape = cellVec[i]->getSecondTapeMutationEvents();
            if(tmpTape.size()>0){
                affectedSecondStrand++;
                for(unsigned int mutI=0;mutI<tmpTape.size();mutI++)
                    historySecondTapeEventCounters.back().incrementMutCounter(tmpTape[mutI]);
            }
            historyAncestorsCounters.back().incrementAncestorCounter(cellVec[i]->getAncestor());
            /**/
        }
    }


    qDebug()<<"Generation: "<< generation <<"Population: "<<numberOfCells<<" Affected: "<<affectedCells<<"   Counter: "<<tmpcounter<<"   Seed: "<<seed;
    //qDebug()<<"      BothAllels: "<<affectedCells<<"   1st: "<<affectedFirstStrand<<"   2nd: "<<affectedSecondStrand;
    populationHistory.push_back(numberOfCells);
    historyNumberOfAffectedCells.push_back(affectedCells);


    int tmpMonoallelicMutatedCells = affectedFirstStrand+affectedSecondStrand -2*affectedCells;
    //Stop Conditions
    if (stopGenerations>-1 && generation >= (unsigned int)stopGenerations)
        return false;
    if (stopNumberOfCells>-1 && numberOfCells>=(unsigned int)stopNumberOfCells)
        return false;
    if (stopMutatedCells>-1 && tmpMonoallelicMutatedCells>=stopMutatedCells)
        return false;
    if (numberOfCells==0)
        return false;
//    qDebug()<<"End process";
    return true;
}

void CellSystem::startSimulation()
{
    reset();

    lastPositionCell = startingNumberOfCells;
    generationLastPositionCell = startingNumberOfCells;
    cellVec.resize(30*startingNumberOfCells,NULL);
    for(unsigned int i=0;i<startingNumberOfCells;i++){
        cellVec[i] = new Cell((long)i,stdProliferation,stdDeath,stdMutationsPerDivision,stdTelomeres,myMutationTable->getMaxIdMutationGene());
        cellVec[i]->setAncestor(i);
        cellVec[i]->registerCellSystem(this);
    }
    counterId = startingNumberOfCells;
    numberOfCells = startingNumberOfCells;

    //logs
    qDebug()<<"Generation: "<< generation <<"Population: "<<numberOfCells;
    populationHistory.push_back(numberOfCells);
    //historyAncestorsCounters.push_back(AncestorHistogram(numberOfCells,1));
}

//return how much where killed
int CellSystem::killAtRandom(int numberOfCellsToDie)
{
    int numberOfDead=0;
    int tmpNumber;
    int i;
    for (i=0; i<numberOfCellsToDie;){
        tmpNumber = this->getRandomIntegerInRange(generationLastPositionCell);
        if(cellVec[tmpNumber]==NULL)
            continue;
        if(cellVec[tmpNumber]->isAlive()){
            i++;
            if ((cellVec[tmpNumber]->getFirstTapeMutationGenes().empty()) && (cellVec[tmpNumber]->getSecondTapeMutationGenes().empty()))
                continue;
            else{
                delete cellVec[tmpNumber];
                cellVec[tmpNumber]=NULL;
                freePositions.push(tmpNumber);
                numberOfDead++;
                numberOfCells--;
            }
        } else{
            continue;
        }
    }
    return numberOfDead;
}

void CellSystem::loadMutationTableFile(QString filename)
{
    myMutationTable->loadMutations(filename);
    genome_size = myMutationTable->getGenomeSize();
}

void CellSystem::loadSinergyTableFile(QString filename)
{
    myMutationTable->loadSinergyPairs(filename);
}

void CellSystem::loadDefaultValues(float tprobProl, float tprobDeath, int tTelomeres, int mutPerDiv)
{
    stdProliferation = tprobProl;
    stdDeath = tprobDeath;
    stdTelomeres = tTelomeres;
    stdMutationsPerDivision = mutPerDiv;
}

void CellSystem::loadSimulationParameters(int tseed, int tnumberOfCells)
{
    seed = tseed;
    startingNumberOfCells = tnumberOfCells;
}

void CellSystem::loadStopConditions(int tgenerations, int tCells, int tMutatedCells)
{
    stopGenerations = tgenerations;
    stopNumberOfCells = tCells;
    stopMutatedCells = tMutatedCells;
}

void CellSystem::addMutationRule(MutationManual *tMut)
{
    rulesMutationManual.push_back(tMut);
    for(unsigned int i=0;i<rulesMutationManual.size(); i++)
        qDebug()<<rulesMutationManual[i]->getMutationGeneName()<<" - "<<rulesMutationManual[i]->getMutationEventName()<<" - "<<rulesMutationManual[i]->getGeneration()<<" - "<<rulesMutationManual[i]->getPercentage()<<" - "<<rulesMutationManual[i]->getMutationAllelsChanged();
}

void CellSystem::updateMutationRule(int id, unsigned int generation, int percentage,QString mutationFullname, int mutationId, manualMutationAllelsChanged tallelschanged)
{
    if((unsigned int)id<rulesMutationManual.size()){
        rulesMutationManual[id]->setGeneration(generation);
        rulesMutationManual[id]->setPercentage(percentage);
        rulesMutationManual[id]->setMutationFullName(mutationFullname);
        rulesMutationManual[id]->setMutationEventId(mutationId);
        rulesMutationManual[id]->setMutationAllelsChanged(tallelschanged);
        for(unsigned int i=0;i<rulesMutationManual.size(); i++)
            qDebug()<<rulesMutationManual[i]->getMutationGeneName()<<" - "<<rulesMutationManual[i]->getMutationEventName()<<" - "<<rulesMutationManual[i]->getGeneration()<<" - "<<rulesMutationManual[i]->getPercentage()<<" - "<<rulesMutationManual[i]->getMutationAllelsChanged();
    }
}

void CellSystem::removeMutationRule(int index)
{
    if((unsigned int)index<rulesMutationManual.size()){
        delete(rulesMutationManual[index]);
        rulesMutationManual.erase(rulesMutationManual.begin()+index);
        for(unsigned int i=0;i<rulesMutationManual.size(); i++)
            qDebug()<<rulesMutationManual[i]->getMutationGeneName()<<rulesMutationManual[i]->getMutationEventName();
    }
}

void CellSystem::removeAllRules()
{
    for(unsigned int i=0; i<rulesMutationManual.size(); i++)
        delete(rulesMutationManual[i]);
    rulesMutationManual.clear();
}

MutationManual *CellSystem::getMutationManualFromIndex(int id) const
{
    if((unsigned int)id<rulesMutationManual.size())
        return rulesMutationManual[id];
    return NULL;
}

void CellSystem::processRuleOfManualMutations(unsigned int id)
{
    MutationManual *currentMutationEvent = rulesMutationManual[id];
    int tmpMutationEventId;
    float proliferationNew,deathNew,telomeresNew,tmpMutationsPerRepNew;
    int tapeChosen;
    tmpMutationEventId = currentMutationEvent->getMutationEventId();
    int tmpMutationGeneId = myMutationTable->getGeneOfEvent(tmpMutationEventId);
    for (unsigned int i=0; i<generationLastPositionCell; i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            if(getRandomProb()*100<=currentMutationEvent->getPercentage()){ //should receive the mutation
                std::vector<unsigned int> tmpfirstTapeEvents = cellVec[i]->getFirstTapeMutationEvents();
                std::vector<unsigned int> tmpsecondTapeEvents = cellVec[i]->getSecondTapeMutationEvents();
                std::vector<unsigned int> newMutationGenes = cellVec[i]->getMutationsInOrder();
                std::vector<unsigned int> newMutationEvents = cellVec[i]->getMutationEventsInOrder();
                std::vector<unsigned int> tmpfirstTapeGenes = cellVec[i]->getFirstTapeMutationGenes();
                std::vector<unsigned int> tmpsecondTapeGenes = cellVec[i]->getSecondTapeMutationGenes();
                proliferationNew=cellVec[i]->getProbReproduction();
                deathNew=cellVec[i]->getProbDeath();
                tmpMutationsPerRepNew=cellVec[i]->getMutationsPerRep();
                telomeresNew=(float)cellVec[i]->getTelomeres();
                float tmpMicroEnvironment = 1;

                if(currentMutationEvent->getMutationAllelsChanged() == MONOALLELIC){
                    addMutationEventToNewCell(tmpMutationEventId,tmpMutationGeneId,tmpfirstTapeGenes,tmpsecondTapeGenes,
                                              tmpfirstTapeEvents,tmpsecondTapeEvents,newMutationGenes,newMutationEvents,-1,proliferationNew,deathNew,
                                              tmpMutationsPerRepNew,telomeresNew,tmpMicroEnvironment);


                } else if(currentMutationEvent->getMutationAllelsChanged() == BIALLELIC){
                    addMutationEventToNewCell(tmpMutationEventId,tmpMutationGeneId,tmpfirstTapeGenes,tmpsecondTapeGenes,
                                              tmpfirstTapeEvents,tmpsecondTapeEvents,newMutationGenes,newMutationEvents,3,proliferationNew,deathNew,
                                              tmpMutationsPerRepNew,telomeresNew,tmpMicroEnvironment); //tapeChosen=3 for both allelles
                }
                cellVec[i]->setFirstTapeMutationGenes(tmpfirstTapeGenes);
                cellVec[i]->setSecondTapeMutationGenes(tmpsecondTapeGenes);
                cellVec[i]->setFirstTapeMutationEvents(tmpfirstTapeEvents);
                cellVec[i]->setSecondTapeMutationEvents(tmpsecondTapeEvents);
                cellVec[i]->setMutationsInOrder(newMutationGenes);
                cellVec[i]->setMutationEventsInOrder(newMutationEvents);
                cellVec[i]->setProbReproduction(proliferationNew);
                cellVec[i]->setProbDeath(deathNew);
                cellVec[i]->setTelomeres(telomeresNew);
                cellVec[i]->setMutationsPerRep(tmpMutationsPerRepNew);
                cellVec[i]->setMutationEventsInOrder(newMutationEvents);
                if (stdGrowthRate>0){
                    stdGrowthRate *= tmpMicroEnvironment;
                }
            }
        }
    }
}

QStringList CellSystem::getMutationListNames() const
{
    return myMutationTable->getMutationListNames();
}

void CellSystem::exportInfo(QString tname)
{
    exportInputParams(tname+"_parameters.txt");
    exportMutationGeneCounters(tname+"_genesMutationResults.csv");
    exportMutationEventCounters(tname+"_eventsMutationResults.csv");
    //exportMutationHistograms(tname+"_numberOfCellsWithNMutations.csv");
    exportAncestralResults(tname+"_ancestralResults.csv");
}

void CellSystem::exportLastLine(QString tname)
{
    exportMutationGeneCountersLastLine(tname+"_genesMutationResults.csv");
    exportMutationEventCountersLastLine(tname+"_eventsMutationResults.csv");
    //exportMutationHistogramsLastLine(tname+"_numberOfCellsWithNMutations.csv");
    exportAncestralResultsLastLine(tname+"_ancestralResults.csv");
    exportMutsEachCellLastLine(tname+"_sequenceEachCell.csv");
}

void CellSystem::exportInputParams(QString tname)
{
    QFile exportParameters(tname);
    if(!exportParameters.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamOutParam(&exportParameters);
    QString textInformation;
    textInformation.append("Parameters\n");
    textInformation.append("Seed: "+QString::number(seed)+"\n");
    textInformation.append("Number of cells (started): "+QString::number(startingNumberOfCells)+"\n");
    textInformation.append("StdProbOfProliferate: "+QString::number(stdProliferation)+"\n");
    textInformation.append("StdProbOfDeath: "+QString::number(stdDeath)+"\n");
    textInformation.append("Telomeres: "+QString::number(stdTelomeres)+"\n");
    textInformation.append("Mutations per div: "+QString::number(stdMutationsPerDivision)+"\n");
    streamOutParam<<textInformation;
    exportParameters.close();
}

void CellSystem::exportMutationGeneCounters(QString tname)
{
    QFile exportintermediateResults(tname);
    if(!exportintermediateResults.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamIntermediateResults(&exportintermediateResults);
    streamIntermediateResults<<"Generation;";
    QStringList mutationGeneNames = myMutationTable->getMutationGenesListNames();
    int tmpBiAllellic, tmpMonoAllellic;
    for(int i=0;i<mutationGeneNames.size();i++){
        //streamIntermediateResults<<mutationNames[i]<<" FT; ";
        //streamIntermediateResults<<mutationNames[i]<<" ST; ";
        streamIntermediateResults<<mutationGeneNames[i]<<" monoAl; ";
        streamIntermediateResults<<mutationGeneNames[i]<<" biAl; ";
    }
    //streamIntermediateResults<<"Affected Cells ; ";
    streamIntermediateResults<<"Total Population ; \n";
    if(historyMutationGeneCounters.size()>0){
        for(unsigned int i=0;i<historyMutationGeneCounters.size();i++){
            streamIntermediateResults<<QString::number(i)+" ; ";
            for(unsigned int j=0;j<myMutationTable->getMaxIdMutationGene();j++) {

                //using union and disjoint logic
                tmpBiAllellic = historyMutationGeneCounters[i].getNumberGivenMut(j);
                tmpMonoAllellic = historyFirstTapeGeneCounters[i].getNumberGivenMut(j) + historySecondTapeGeneCounters[i].getNumberGivenMut(j) - 2*tmpBiAllellic;
                streamIntermediateResults<<QString::number(tmpMonoAllellic)+" ; ";
                streamIntermediateResults<<QString::number(tmpBiAllellic)+" ; ";
            }
            //streamIntermediateResults<<QString::number(historyNumberOfAffectedCells[i])+" ; ";
            streamIntermediateResults<<QString::number(populationHistory[i])+" ; \n";
        }
    }
    exportintermediateResults.close();
}

void CellSystem::exportMutationEventCounters(QString tname)
{
    QFile exportintermediateResults(tname);
    if(!exportintermediateResults.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamIntermediateResults(&exportintermediateResults);
    streamIntermediateResults<<"Generation;";
    QStringList mutationEventNames = myMutationTable->getMutationListNames();
    for(int i=0;i<mutationEventNames.size();i++){
        //streamIntermediateResults<<mutationNames[i]<<" FT; ";
        //streamIntermediateResults<<mutationNames[i]<<" ST; ";
        //streamIntermediateResults<<mutationEventNames[i]<<" monoAl; ";
        //streamIntermediateResults<<mutationEventNames[i]<<" biAl; ";
        streamIntermediateResults<<mutationEventNames[i]<<+" ; ";
    }
    //streamIntermediateResults<<"Affected Cells ; ";
    streamIntermediateResults<<"Total Population ; \n";
    if(historyFirstTapeEventCounters.size()>0){
        for(unsigned int i=0;i<historyFirstTapeEventCounters.size();i++){
            streamIntermediateResults<<QString::number(i)+" ; ";
            for(unsigned int j=0;j<myMutationTable->getMaxIdMutationEvent();j++) {

                //using union and disjoint logic
                //tmpBiAllellic = historyMutationEventCounters[i].getNumberGivenMut(j);
                //tmpMonoAllellic = historyFirstTapeEventCounters[i].getNumberGivenMut(j) + historySecondTapeEventCounters[i].getNumberGivenMut(j) - 2*tmpBiAllellic;
                //streamIntermediateResults<<QString::number(tmpMonoAllellic)+" ; ";
                //summing individual events subtracting intersection
                streamIntermediateResults<<QString::number(historyFirstTapeEventCounters[i].getNumberGivenMut(j) + historySecondTapeEventCounters[i].getNumberGivenMut(j)-historyMutationEventCounters[i].getNumberGivenMut(j))+" ; ";
            }
            //streamIntermediateResults<<QString::number(historyNumberOfAffectedCells[i])+" ; ";
            streamIntermediateResults<<QString::number(populationHistory[i])+" ; \n";
        }
    }
    exportintermediateResults.close();
}



void CellSystem::exportMutationHistograms(QString tname)
{
    QFile exportmutHistResults(tname);
    if(!exportmutHistResults.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streammutHistResults(&exportmutHistResults);
    streammutHistResults<<"Generation;";
    streammutHistResults<<"Population;";
    streammutHistResults<<"Affected cells;";
    for (unsigned int i=1;i<myMutationTable->getMaxIdMutationGene();i++)
        streammutHistResults<<QString::number(i)+" mutation(s);";
    streammutHistResults<<"\n";
    for(unsigned int i=0; i<historyMutationGeneHistogram.size();i++){
        streammutHistResults<<QString::number(i)+";";
        streammutHistResults<<QString::number(populationHistory[i])+";";
        streammutHistResults<<QString::number(historyNumberOfAffectedCells[i])+";";
        for (unsigned int j=1; j<myMutationTable->getMaxIdMutationGene(); j++){
            streammutHistResults<<QString::number(historyMutationGeneHistogram[i].getNumberOfSamples(j))+";";
        }
        streammutHistResults<<"\n";
    }
    exportmutHistResults.close();
}

void CellSystem::exportAncestralResults(QString tname)
{
    QFile exportAncResults(tname);
    if(!exportAncResults.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamAncestralResults(&exportAncResults);
    streamAncestralResults<<"Generation;";
    streamAncestralResults<<"Population;";
    //streamAncestralResults<<"Affected cells;";
    for (unsigned int i=0;i<startingNumberOfCells;i++)
        streamAncestralResults<<"A"+QString::number(i)+" ;";
    streamAncestralResults<<"\n";
    for(unsigned int i=0;i<historyAncestorsCounters.size();i++){
        streamAncestralResults<<QString::number(i)+";";
        streamAncestralResults<<QString::number(populationHistory[i])+";";
        //streamAncestralResults<<QString::number(historyNumberOfAffectedCells[i])+";";
        for (unsigned int j=0; j<startingNumberOfCells; j++){
            streamAncestralResults<<QString::number(historyAncestorsCounters[i].getNumberofGivenAncestor(j))+";";
        }
        streamAncestralResults<<"\n";
    }
    exportAncResults.close();
}

void CellSystem::exportMutationGeneCountersLastLine(QString tname)
{
    bool first = !(fileExists(tname));
    QFile exportintermediateResults(tname);
    if(!exportintermediateResults.open(QIODevice::WriteOnly | QIODevice::Append) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamIntermediateResults(&exportintermediateResults);
    if (first){
        streamIntermediateResults<<"Seed ; ";
        streamIntermediateResults<<"Generation ; ";
        QStringList mutationGeneNames = myMutationTable->getMutationGenesListNames();
        for(int i=0;i<mutationGeneNames.size();i++){
            //streamIntermediateResults<<mutationNames[i]<<" 1s; ";
            //streamIntermediateResults<<mutationNames[i]<<" 2s; ";
            streamIntermediateResults<<mutationGeneNames[i]<<" monoAl; ";
            streamIntermediateResults<<mutationGeneNames[i]<<" biAl; ";
        }
        //streamIntermediateResults<<"Affected Cells ; ";
        streamIntermediateResults<<"Total Population ; \n";
    }
    streamIntermediateResults<<QString::number(seed)+" ; ";
    streamIntermediateResults<<QString::number(historyMutationGeneCounters.size()-1)+" ; ";
    int tmpMonoAllellic, tmpBiAllellic;
    for(unsigned int j=0;j<myMutationTable->getMaxIdMutationGene();j++){

        //using union and disjoint logic
        tmpBiAllellic = historyMutationGeneCounters.back().getNumberGivenMut(j);
        tmpMonoAllellic = historyFirstTapeGeneCounters.back().getNumberGivenMut(j) + historySecondTapeGeneCounters.back().getNumberGivenMut(j) - 2*tmpBiAllellic;
        streamIntermediateResults<<QString::number(tmpMonoAllellic)+" ; ";
        streamIntermediateResults<<QString::number(tmpBiAllellic)+" ; ";
    }
    //streamIntermediateResults<<QString::number(historyNumberOfAffectedCells.back())+" ; ";
    streamIntermediateResults<<QString::number(populationHistory.back())+" ; \n";
    exportintermediateResults.close();
}


void CellSystem::exportMutationEventCountersLastLine(QString tname)
{
    bool first = !(fileExists(tname));
    QFile exportintermediateResults(tname);
    if(!exportintermediateResults.open(QIODevice::WriteOnly | QIODevice::Append) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamIntermediateResults(&exportintermediateResults);
    if (first){
        streamIntermediateResults<<"Seed ; ";
        streamIntermediateResults<<"Generation ; ";
        QStringList mutationEventNames = myMutationTable->getMutationListNames();
        for(int i=0;i<mutationEventNames.size();i++){
            //streamIntermediateResults<<mutationNames[i]<<" 1s; ";
            //streamIntermediateResults<<mutationNames[i]<<" 2s; ";
            //streamIntermediateResults<<mutationEventNames[i]<<" monoAl; ";
            //streamIntermediateResults<<mutationEventNames[i]<<" biAl; ";

            streamIntermediateResults<<mutationEventNames[i]<<" ; ";
        }
        //streamIntermediateResults<<"Affected Cells ; ";
        streamIntermediateResults<<"Total Population ; \n";
    }
    streamIntermediateResults<<QString::number(seed)+" ; ";
    streamIntermediateResults<<QString::number(historyFirstTapeEventCounters.size()-1)+" ; ";
    for(unsigned int j=0;j<myMutationTable->getMaxIdMutationEvent();j++){

        //using union and disjoint logic
        //tmpBiAllellic = historyMutationEventCounters.back().getNumberGivenMut(j);
        //tmpMonoAllellic = historyFirstTapeEventCounters.back().getNumberGivenMut(j) + historySecondTapeEventCounters.back().getNumberGivenMut(j) - 2*tmpBiAllellic;
        //streamIntermediateResults<<QString::number(tmpMonoAllellic)+" ; ";
        //summing individual events subtracting the intersection
        streamIntermediateResults<<QString::number(historyFirstTapeEventCounters.back().getNumberGivenMut(j) + historySecondTapeEventCounters.back().getNumberGivenMut(j) - historyMutationEventCounters.back().getNumberGivenMut(j))+" ; ";
    }
    //streamIntermediateResults<<QString::number(historyNumberOfAffectedCells.back())+" ; ";
    streamIntermediateResults<<QString::number(populationHistory.back())+" ; \n";
    exportintermediateResults.close();
}




void CellSystem::exportMutationHistogramsLastLine(QString tname)
{
    bool first = !(fileExists(tname));
    QFile exportmutHistResults(tname);
    if(!exportmutHistResults.open(QIODevice::WriteOnly | QIODevice::Append) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streammutHistResults(&exportmutHistResults);
    if (first){
        streammutHistResults<<"Seed;";
        streammutHistResults<<"Generation;";
        streammutHistResults<<"Population;";
        //streammutHistResults<<"Affected cells;";
        for (unsigned int i=1;i<myMutationTable->getMaxIdMutationGene();i++)
            streammutHistResults<<QString::number(i)+" mutation(s);";
        streammutHistResults<<"\n";
    }
    streammutHistResults<<QString::number(seed)+";";
    streammutHistResults<<QString::number(historyMutationGeneHistogram.size()-1)+";";
    streammutHistResults<<QString::number(populationHistory.back())+";";
    //streammutHistResults<<QString::number(historyNumberOfAffectedCells.back())+";";
    for (unsigned int j=1; j<myMutationTable->getMaxIdMutationGene(); j++){
        streammutHistResults<<QString::number(historyMutationGeneHistogram.back().getNumberOfSamples(j))+";";
    }
    streammutHistResults<<"\n";
    exportmutHistResults.close();
}

void CellSystem::exportAncestralResultsLastLine(QString tname)
{
    bool first = !(fileExists(tname));
    QFile exportAncResults(tname);
    if(!exportAncResults.open(QIODevice::WriteOnly | QIODevice::Append) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamAncestralResults(&exportAncResults);
    if (first){
        streamAncestralResults<<"Seed;";
        streamAncestralResults<<"Generation;";
        streamAncestralResults<<"NumberOfDivisions;";
        streamAncestralResults<<"Population;";
        //streamAncestralResults<<"Affected cells;";
        for (unsigned int i=0;i<startingNumberOfCells;i++)
            streamAncestralResults<<"A"+QString::number(i)+" ;";
        streamAncestralResults<<"\n";
    }
    streamAncestralResults<<QString::number(seed)+";";
    streamAncestralResults<<QString::number(historyAncestorsCounters.size()-1)+";";
    streamAncestralResults<<QString::number(numberOfDivisions)+";";
    streamAncestralResults<<QString::number(populationHistory.back())+";";
    //streamAncestralResults<<QString::number(historyNumberOfAffectedCells.back())+";";
    for (unsigned int j=0; j<startingNumberOfCells; j++){
        streamAncestralResults<<QString::number(historyAncestorsCounters.back().getNumberofGivenAncestor(j))+";";
    }
    streamAncestralResults<<"\n";

    exportAncResults.close();
}

void CellSystem::exportMutsEachCellLastLine(QString tname)
{
    bool first = !(fileExists(tname));
    QFile exportMutsEachCellResults(tname);
    if(!exportMutsEachCellResults.open(QIODevice::WriteOnly | QIODevice::Append) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamMutsEachCellResults(&exportMutsEachCellResults);
    if (first){
        streamMutsEachCellResults<<"Seed;";
        streamMutsEachCellResults<<"Generation;";
        streamMutsEachCellResults<<"Population;";
        //streamMutsEachCellResults<<"Affected cells;";
        for (unsigned int i=0;i<numberOfCells;i++){
            //streamMutsEachCellResults<<"Cell"+QString::number(i)+" ;";
            streamMutsEachCellResults<<"Cell"+QString::number(i)+" firstStrand;";
            streamMutsEachCellResults<<"Cell"+QString::number(i)+" secondStrand;";
        }
        streamMutsEachCellResults<<"\n";
    }
    streamMutsEachCellResults<<QString::number(seed)+";";
    streamMutsEachCellResults<<QString::number(generation)+";";
    streamMutsEachCellResults<<QString::number(populationHistory.back())+";";
    //streamMutsEachCellResults<<QString::number(historyNumberOfAffectedCells.back())+";";
    QString tmpString = "";
    QString tmpChar;
    QString tmpStringFirstTape = "";
    QString tmpStringSecondTape = "";
    for(unsigned int i=0; i<generationLastPositionCell;i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            tmpString="m";
            tmpStringFirstTape="1s";
            tmpStringSecondTape="2s";
            for (unsigned int iMut=0; iMut<myMutationTable->getMaxIdMutationGene(); iMut++){
                if(cellVec[i]->hasThisMut(iMut))
                    tmpChar = "1";
                else
                    tmpChar = "0";
                tmpString = tmpString + tmpChar;

                if(cellVec[i]->hasThisMutGeneOnTape(iMut,0)) // firstTape
                    tmpChar = "1";
                else
                    tmpChar = "0";
                tmpStringFirstTape =tmpStringFirstTape+tmpChar;

                if(cellVec[i]->hasThisMutGeneOnTape(iMut,1)) // secondTape
                    tmpChar = "1";
                else
                    tmpChar = "0";
                tmpStringSecondTape =tmpStringSecondTape+tmpChar;
            }
            //streamMutsEachCellResults<<tmpString+";";
            streamMutsEachCellResults<<tmpStringFirstTape+";";
            streamMutsEachCellResults<<tmpStringSecondTape+";";
        }
    }
    streamMutsEachCellResults<<"\n";

    exportMutsEachCellResults.close();
}

void CellSystem::show()
{
    //qDebug()<<"Size: "<<historyMutationEventCounters.size();
    //for(unsigned int i=0;i<myMutationTable->getMaxIdMutationEvent();i++)
    //    qDebug()<<qPrintable(QString::fromStdString(myMutationTable->getMutationEvent(i)->mName))<<": "<<historyMutationEventCounters.back().getNumberGivenMut(i);
    //qDebug()<<"1 mutations: "<<historyMutationGeneHistogram.back().getNumberOfSamples(1);
}
