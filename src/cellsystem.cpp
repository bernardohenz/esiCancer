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
    historyMutationHistogram.clear();
    historyMutationHistogram.push_back(MutationHistogram(myMutationTable->getMaxIdMutation()));
    historyMutationCounters.clear();
    historyMutationCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutation()));

    stdNumberOfAncestors = startingNumberOfCells;
    historyAncestorsCounters.clear();
    historyAncestorsCounters.push_back(AncestorHistogram(stdNumberOfAncestors));
}

Cell *CellSystem::addNewCell(Cell *father){
    unsigned int tmpMutatedPosition;
    int tmpMutationId;
    unsigned int tapeChosen;
    //copy from father
    std::vector<unsigned int> newMutations = father->getMutationsInOrder();
    std::vector<unsigned int> newFirstTapeMutations = father->getFirstTapeMutations();
    std::vector<unsigned int> newSecondTapeMutations = father->getSecondTapeMutations();
    std::vector<long> newParentsOrder = father->getParentsInOrder();

    float proliferationChild = father->getProbReproduction();
    float deathChild = father->getProbDeath();
    float telomeresChild = father->getTelomeres();
    float tmpMutationsPerRepChild = father->getMutationsPerRep();
    float tmpMutationsPerRepFather = father->getMutationsPerRep();
    bool addedAnySuppressor=false;
    bool workSynergy = false;
//    if (generation>520){
//        qDebug()<<"Adding New Cell:  mutPerRep: "<<tmpMutationsPerRepFather<<"    mutations: "<<newMutations.size()<<"    Muts:";
//        for (int i=0;i<newMutations.size();i++)
//            qDebug()<<"   Mut "<<newMutations[i];
//    }
    for(unsigned int i=0; i<tmpMutationsPerRepFather; i++){
        //get a random basis position
        tmpMutatedPosition = getRandomInteger();

        //get the mutation id corresponding to the mutated base
        tmpMutationId = myMutationTable->getNumberofMutation(tmpMutatedPosition);
        //check if it is a valid mutation
        if(tmpMutationId>-1){

            //tmpMutationType = myMutationTable->getTypeOfMutationGivenAMutation(tmpMutationId);
            bool activateMutation = false;
            if (getRandomProb()>0.5){
                tapeChosen = 1;
                activateMutation = myMutationTable->checkIfAddMutationByReference(tmpMutationId,newFirstTapeMutations,newSecondTapeMutations,tapeChosen);
                if(!Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
                    newFirstTapeMutations.push_back(tmpMutationId);
                //newFirstTapeMutations[tmpMutationId] = true;
            }else {
                tapeChosen = 2;
                activateMutation = myMutationTable->checkIfAddMutationByReference(tmpMutationId,newFirstTapeMutations,newSecondTapeMutations,tapeChosen);
                if(!Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
                    newSecondTapeMutations.push_back(tmpMutationId);
                //newSecondTapeMutations[tmpMutationId] = true;
            }

            if (activateMutation){
                //if(myMutationTable->getTypeOfMutationGivenAMutation(tmpMutationId) == MutationType::TUMORSUPPRESSOR)
                //    qDebug()<<"MutationID: "<<tmpMutationId<<"  ACTIVATED";
                Mutation *activatedMutation = myMutationTable->getMutation(tmpMutationId);
                if(father->getHasSuppressorMut() || addedAnySuppressor){
                    proliferationChild *= activatedMutation->prolRateMultAfterSuppressor;
                    proliferationChild += activatedMutation->prolRateAddAfterSuppressor;
                    deathChild *= activatedMutation->deathRateMultAfterSuppressor;
                    deathChild += activatedMutation->deathRateAddAfterSuppressor;
                    telomeresChild *= activatedMutation->telomeresRateMultAfterSuppressor;
                    telomeresChild += activatedMutation->telomeresRateAddAfterSuppressor;
                } else{
                    proliferationChild *= activatedMutation->prolRateMultBeforeSuppressor;
                    proliferationChild += activatedMutation->prolRateAddBeforeSuppressor;
                    deathChild *= activatedMutation->deathRateMultBeforeSuppressor;
                    deathChild += activatedMutation->deathRateAddBeforeSuppressor;
                    telomeresChild *= activatedMutation->telomeresRateMultBeforeSuppressor;
                    telomeresChild += activatedMutation->telomeresRateAddBeforeSuppressor;
                }
                tmpMutationsPerRepChild *= activatedMutation->moreMutRateMultMutations;
                tmpMutationsPerRepChild += activatedMutation->moreMutRateAddMutations;
                newMutations.push_back(tmpMutationId);
                workSynergy = true;
                if (activatedMutation->mType==MutationType::TUMORSUPPRESSOR)
                    addedAnySuppressor=true;
            }
        }
    }
    float sinergyTmp=1;
    if (workSynergy)
        for (unsigned int i=0; i<newMutations.size();i++)
            sinergyTmp *= myMutationTable->getSynergyOfMutation(newMutations[i]);
    if(sinergyTmp<=0)
        sinergyTmp=1;
    newParentsOrder.push_back(father->getId());
    tmpMutationsPerRepChild = std::min((unsigned int)tmpMutationsPerRepChild,maxMutationsPerDivision);
    proliferationChild = std::min(proliferationChild,maxProliferation);
    deathChild = std::min(deathChild,maxDeath);
    Cell *newCell = new Cell(counterCellId,proliferationChild, deathChild, tmpMutationsPerRepChild,
                                 telomeresChild, newMutations, newFirstTapeMutations, newSecondTapeMutations, newParentsOrder);
    newCell->setSinergy(sinergyTmp);
    newCell->setAncestor(father->getAncestor());
    if(father->getHasSuppressorMut()){
        newCell->setHasSuppressorMut(true);
    }else{
        if(addedAnySuppressor){
            updateOncsNowTumorSuppressor(newCell);
            newCell->setHasSuppressorMut(true);
        }
    }
    newCell->registerCellSystem(this);
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
    return newCell;
}

void CellSystem::updateOncsNowTumorSuppressor(Cell *tcell)
{
    float prolCell=stdProliferation;
    float deathCell=stdDeath;
    float telomeresCell=stdTelomeres;
    std::vector<unsigned int> mutationsOfCell = tcell->getMutationsInOrder();
    Mutation *tmpMutation;
    for (unsigned int i=0; i<mutationsOfCell.size(); i++){
        tmpMutation = myMutationTable->getMutation(mutationsOfCell[i]);
        prolCell *= tmpMutation->prolRateMultAfterSuppressor;
        prolCell += tmpMutation->prolRateAddAfterSuppressor;
        deathCell *= tmpMutation->deathRateMultAfterSuppressor;
        deathCell += tmpMutation->deathRateAddAfterSuppressor;
        telomeresCell *= tmpMutation->telomeresRateMultAfterSuppressor;
        telomeresCell += tmpMutation->telomeresRateAddAfterSuppressor;
    }
    tcell->setProbReproduction(prolCell);
    tcell->setProbDeath(deathCell);
    tcell->setTelomeres(telomeresCell - tcell->getNumberOfReproductions());
}

void CellSystem::informNaturalKilled(long id)
{

}

float CellSystem::getRandomProb()
{
    return (float)probabilityDistribution(probabilityGenerator);
    //float a = (float)probabilityDistribution(probabilityGenerator);
    //return a;
}

long CellSystem::getRandomInteger()
{
    return (randomGenerator)()%genome_size;
    //long randomNumber = (randomGenerator)()%genome_size;
    //qDebug()<<"GenomeSize: "<<QString::number(genome_size);
    //return randomNumber;
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
            numberOfCells--;
            delete cellVec[i];
            cellVec[i]=NULL;
            freePositions.push(i);
        }
    }
//    qDebug()<<"minorDebug";
    generation++;
    generationLastPositionCell = lastPositionCell;

    //Logs recording
    std::vector<unsigned int> tmpMuts;
    int affectedCells=0;

    historyMutationHistogram.push_back(MutationHistogram(myMutationTable->getMaxIdMutation()));
    historyMutationCounters.push_back(MutationCounter(myMutationTable->getMaxIdMutation()));
    historyAncestorsCounters.push_back(AncestorHistogram(stdNumberOfAncestors));
    for(unsigned int i=0;i<generationLastPositionCell;i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            tmpMuts = cellVec[i]->getMutationsInOrder();
            if(tmpMuts.size()>0){
                affectedCells++;
                for(unsigned int mutI=0;mutI<tmpMuts.size();mutI++)
                    historyMutationCounters.back().incrementMutCounter(tmpMuts[mutI]);
                //qDebug()<<"SIZE: "<<tmpMuts.size();
                historyMutationHistogram.back().incrementBin(tmpMuts.size());
            }
            historyAncestorsCounters.back().incrementAncestorCounter(cellVec[i]->getAncestor());
        }
    }
    qDebug()<<"Generation: "<< generation <<"Population: "<<numberOfCells<<" Affected: "<<affectedCells;
    populationHistory.push_back(numberOfCells);
    historyNumberOfAffectedCells.push_back(affectedCells);

    //Stop Conditions
    if (stopGenerations>-1 && generation >= (unsigned int)stopGenerations)
        return false;
    if (stopNumberOfCells>-1 && numberOfCells>=(unsigned int)stopNumberOfCells)
        return false;
    if (stopMutatedCells>-1 && affectedCells>=stopMutatedCells)
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
    unsigned int maxIdMut = myMutationTable->getMaxIdMutation();
    for(unsigned int i=0;i<startingNumberOfCells;i++){
        cellVec[i] = new Cell((long)i,stdProliferation,stdDeath,stdMutationsPerDivision,stdTelomeres,maxIdMut);
        cellVec[i]->setAncestor(i);
        cellVec[i]->registerCellSystem(this);
    }
    counterId = startingNumberOfCells;
    numberOfCells = startingNumberOfCells;

    //logs
    populationHistory.push_back(numberOfCells);
    //historyAncestorsCounters.push_back(AncestorHistogram(numberOfCells,1));
}

void CellSystem::loadMutationTableFile(QString filename)
{
    myMutationTable->loadMutations(filename);
    genome_size = myMutationTable->getGenomeSize();
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
        qDebug()<<rulesMutationManual[i]->getMutationName()<<" - "<<rulesMutationManual[i]->getGeneration()<<" - "<<rulesMutationManual[i]->getPercentage()<<" - "<<rulesMutationManual[i]->getMutationAllelsChanged();
}

void CellSystem::updateMutationRule(int id, unsigned int generation, int percentage, QString mutationname, int mutationId, manualMutationAllelsChanged tallelschanged)
{
    if((unsigned int)id<rulesMutationManual.size()){
        rulesMutationManual[id]->setGeneration(generation);
        rulesMutationManual[id]->setPercentage(percentage);
        rulesMutationManual[id]->setMutationName(mutationname);
        rulesMutationManual[id]->setMutationId(mutationId);
        rulesMutationManual[id]->setMutationAllelsChanged(tallelschanged);
        for(unsigned int i=0;i<rulesMutationManual.size(); i++)
            qDebug()<<rulesMutationManual[i]->getMutationName()<<" - "<<rulesMutationManual[i]->getGeneration()<<" - "<<rulesMutationManual[i]->getPercentage()<<" - "<<rulesMutationManual[i]->getMutationAllelsChanged();
    }
}

void CellSystem::removeMutationRule(int index)
{
    if((unsigned int)index<rulesMutationManual.size()){
        delete(rulesMutationManual[index]);
        rulesMutationManual.erase(rulesMutationManual.begin()+index);
        for(unsigned int i=0;i<rulesMutationManual.size(); i++)
            qDebug()<<rulesMutationManual[i]->getMutationName();
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
    MutationManual *currentMutation = rulesMutationManual[id];
    int tmpMutationId;
    float proliferationNew,deathNew,telomeresNew,tmpMutationsPerRepNew;
    float sinergyTmp;
    int tapeChosen;
    for (unsigned int i=0; i<generationLastPositionCell; i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            if(getRandomProb()*100<=currentMutation->getPercentage()){ //should receive the mutation
                tmpMutationId = currentMutation->getMutationId();
                std::vector<unsigned int> tmpfirstTape = cellVec[i]->getFirstTapeMutations();
                std::vector<unsigned int> tmpsecondTape = cellVec[i]->getSecondTapeMutations();
                bool activateMutation = false;
                if(currentMutation->getMutationAllelsChanged() == MONOALLELIC){
                    if (getRandomProb()>0.5){
                        tapeChosen = 1;
                        activateMutation = myMutationTable->checkIfAddMutationByReference(tmpMutationId,tmpfirstTape,tmpsecondTape,tapeChosen);
                        //if(!Contains(tmpfirstTape,(unsigned int)tmpMutationId))
                        //    tmpfirstTape.push_back(tmpMutationId);
                        //tmpfirstTape[tmpMutationId] = true;
                        cellVec[i]->setFirstTapeMutations(tmpfirstTape);
                    }else {
                        tapeChosen = 2;
                        activateMutation = myMutationTable->checkIfAddMutationByReference(tmpMutationId,tmpfirstTape,tmpsecondTape,tapeChosen);
                        //if(!Contains(tmpsecondTape,(unsigned int)tmpMutationId))
                        //    tmpsecondTape.push_back(tmpMutationId);
                        //tmpsecondTape[tmpMutationId] = true;
                        cellVec[i]->setSecondTapeMutations(tmpsecondTape);
                    }
                } else if(currentMutation->getMutationAllelsChanged() == BIALLELIC){
                    if(!(cellVec[i]->hasThisMut(tmpMutationId))){
                        activateMutation=true;
                        if (!Contains(tmpfirstTape,(unsigned int)tmpMutationId)) {
                            tmpfirstTape.push_back(tmpMutationId);
                            cellVec[i]->setFirstTapeMutations(tmpfirstTape);
                        }
                        if (!Contains(tmpsecondTape,(unsigned int)tmpMutationId)){
                            tmpsecondTape.push_back(tmpMutationId);
                            cellVec[i]->setSecondTapeMutations(tmpsecondTape);
                        }
                    }
                }

                if(activateMutation){
                    Mutation *activatedMutation = myMutationTable->getMutation(tmpMutationId);
                    std::vector<unsigned int> newMutations = cellVec[i]->getMutationsInOrder();
                    proliferationNew=cellVec[i]->getProbReproduction();
                    deathNew=cellVec[i]->getProbDeath();
                    telomeresNew=(float)cellVec[i]->getTelomeres();
                    tmpMutationsPerRepNew=cellVec[i]->getMutationsPerRep();
                    if(cellVec[i]->getHasSuppressorMut()){
                        proliferationNew *= activatedMutation->prolRateMultAfterSuppressor;
                        proliferationNew += activatedMutation->prolRateAddAfterSuppressor;
                        deathNew *= activatedMutation->deathRateMultAfterSuppressor;
                        deathNew += activatedMutation->deathRateAddAfterSuppressor;
                        telomeresNew *= activatedMutation->telomeresRateMultAfterSuppressor;
                        telomeresNew += activatedMutation->telomeresRateAddAfterSuppressor;
                    } else{
                        proliferationNew *= activatedMutation->prolRateMultBeforeSuppressor;
                        proliferationNew += activatedMutation->prolRateAddBeforeSuppressor;
                        deathNew *= activatedMutation->deathRateMultBeforeSuppressor;
                        deathNew += activatedMutation->deathRateAddBeforeSuppressor;
                        telomeresNew *= activatedMutation->telomeresRateMultBeforeSuppressor;
                        telomeresNew += activatedMutation->telomeresRateAddBeforeSuppressor;
                    }
                    tmpMutationsPerRepNew *= activatedMutation->moreMutRateMultMutations;
                    tmpMutationsPerRepNew += activatedMutation->moreMutRateAddMutations;
                    cellVec[i]->addMutation(tmpMutationId);
                    newMutations.push_back(tmpMutationId);
                    sinergyTmp=1;
                    if (newMutations.size()>1)
                        for (unsigned int i=0; i<newMutations.size();i++)
                            sinergyTmp *= myMutationTable->getSynergyOfMutation(newMutations[i]);
                    cellVec[i]->setProbReproduction(proliferationNew);
                    cellVec[i]->setProbDeath(deathNew);
                    cellVec[i]->setTelomeres(telomeresNew);
                    cellVec[i]->setMutationsPerRep(tmpMutationsPerRepNew);
                    cellVec[i]->setSinergy(sinergyTmp);
                    if ((!(cellVec[i]->getHasSuppressorMut()))&&(activatedMutation->mType==MutationType::TUMORSUPPRESSOR)){
                        updateOncsNowTumorSuppressor(cellVec[i]);
                        cellVec[i]->setHasSuppressorMut(true);
                    }
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
    exportMutationCounters(tname+"_GeneMutationResults.csv");
    exportMutationHistograms(tname+"_numberOfCellsWithNMutations.csv");
    exportAncestralResults(tname+"_ancestralResults.csv");
}

void CellSystem::exportLastLine(QString tname)
{
    exportMutationCountersLastLine(tname+"_GeneMutationResults.csv");
    exportMutationHistogramsLastLine(tname+"_numberOfCellsWithNMutations.csv");
    exportAncestralResultsLastLine(tname+"_ancestralResults.csv");
    exportMutsEachCellLastLine(tname+"_mutsEachCell.csv");
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

void CellSystem::exportMutationCounters(QString tname)
{
    QFile exportintermediateResults(tname);
    if(!exportintermediateResults.open(QIODevice::WriteOnly) ){
        qDebug()<<"Error opening: "<<qPrintable(tname);
        return;
    }
    QTextStream streamIntermediateResults(&exportintermediateResults);
    streamIntermediateResults<<"Generation;";
    QStringList mutationNames = myMutationTable->getMutationListNames();
    for(unsigned int i=0;i<myMutationTable->getMaxIdMutation();i++){
        streamIntermediateResults<<mutationNames[i]<<" ; ";
    }
    streamIntermediateResults<<"Affected Cells ; ";
    streamIntermediateResults<<"Total Population ; \n";
    if(historyMutationCounters.size()>0){
        for(unsigned int i=0;i<historyMutationCounters.size();i++){
            streamIntermediateResults<<QString::number(i)+" ; ";
            for(unsigned int j=0;j<myMutationTable->getMaxIdMutation();j++)
                streamIntermediateResults<<QString::number(historyMutationCounters[i].getNumberGivenMut(j))+" ; ";
            streamIntermediateResults<<QString::number(historyNumberOfAffectedCells[i])+" ; ";
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
    for (unsigned int i=1;i<myMutationTable->getMaxIdMutation();i++)
        streammutHistResults<<QString::number(i)+" mutation(s);";
    streammutHistResults<<"\n";
    for(unsigned int i=0; i<historyMutationHistogram.size();i++){
        streammutHistResults<<QString::number(i)+";";
        streammutHistResults<<QString::number(populationHistory[i])+";";
        streammutHistResults<<QString::number(historyNumberOfAffectedCells[i])+";";
        for (unsigned int j=1; j<myMutationTable->getMaxIdMutation(); j++){
            streammutHistResults<<QString::number(historyMutationHistogram[i].getNumberOfSamples(j))+";";
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
    streamAncestralResults<<"Affected cells;";
    for (unsigned int i=0;i<startingNumberOfCells;i++)
        streamAncestralResults<<"A"+QString::number(i)+" ;";
    streamAncestralResults<<"\n";
    for(unsigned int i=0;i<historyAncestorsCounters.size();i++){
        streamAncestralResults<<QString::number(i)+";";
        streamAncestralResults<<QString::number(populationHistory[i])+";";
        streamAncestralResults<<QString::number(historyNumberOfAffectedCells[i])+";";
        for (unsigned int j=0; j<startingNumberOfCells; j++){
            streamAncestralResults<<QString::number(historyAncestorsCounters[i].getNumberofGivenAncestor(j))+";";
        }
        streamAncestralResults<<"\n";
    }
    exportAncResults.close();
}

void CellSystem::exportMutationCountersLastLine(QString tname)
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
        QStringList mutationNames = myMutationTable->getMutationListNames();
        for(unsigned int i=0;i<myMutationTable->getMaxIdMutation();i++){
            streamIntermediateResults<<mutationNames[i]<<" ; ";
        }
        streamIntermediateResults<<"Affected Cells ; ";
        streamIntermediateResults<<"Total Population ; \n";
    }
    streamIntermediateResults<<QString::number(seed)+" ; ";
    streamIntermediateResults<<QString::number(historyMutationCounters.size()-1)+" ; ";
    for(unsigned int j=0;j<myMutationTable->getMaxIdMutation();j++)
        streamIntermediateResults<<QString::number(historyMutationCounters.back().getNumberGivenMut(j))+" ; ";
    streamIntermediateResults<<QString::number(historyNumberOfAffectedCells.back())+" ; ";
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
        streammutHistResults<<"Affected cells;";
        for (unsigned int i=1;i<myMutationTable->getMaxIdMutation();i++)
            streammutHistResults<<QString::number(i)+" mutation(s);";
        streammutHistResults<<"\n";
    }
    streammutHistResults<<QString::number(seed)+";";
    streammutHistResults<<QString::number(historyMutationHistogram.size()-1)+";";
    streammutHistResults<<QString::number(populationHistory.back())+";";
    streammutHistResults<<QString::number(historyNumberOfAffectedCells.back())+";";
    for (unsigned int j=1; j<myMutationTable->getMaxIdMutation(); j++){
        streammutHistResults<<QString::number(historyMutationHistogram.back().getNumberOfSamples(j))+";";
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
        streamAncestralResults<<"Affected cells;";
        for (unsigned int i=0;i<startingNumberOfCells;i++)
            streamAncestralResults<<"A"+QString::number(i)+" ;";
        streamAncestralResults<<"\n";
    }
    streamAncestralResults<<QString::number(seed)+";";
    streamAncestralResults<<QString::number(historyAncestorsCounters.size()-1)+";";
    streamAncestralResults<<QString::number(numberOfDivisions)+";";
    streamAncestralResults<<QString::number(populationHistory.back())+";";
    streamAncestralResults<<QString::number(historyNumberOfAffectedCells.back())+";";
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
        streamMutsEachCellResults<<"Affected cells;";
        for (unsigned int i=0;i<numberOfCells;i++)
            streamMutsEachCellResults<<"Cell"+QString::number(i)+" ;";
        streamMutsEachCellResults<<"\n";
    }
    streamMutsEachCellResults<<QString::number(seed)+";";
    streamMutsEachCellResults<<QString::number(generation)+";";
    streamMutsEachCellResults<<QString::number(populationHistory.back())+";";
    streamMutsEachCellResults<<QString::number(historyNumberOfAffectedCells.back())+";";
    QString tmpString = "";
    QString tmpChar;
    for(unsigned int i=0; i<generationLastPositionCell;i++){
        if(cellVec[i]==NULL)
            continue;
        if(cellVec[i]->isAlive()){
            tmpString="m";
            for (unsigned int iMut=0; iMut<myMutationTable->getMaxIdMutation(); iMut++){
                if(cellVec[i]->hasThisMut(iMut))
                    tmpChar = "1";
                else
                    tmpChar = "0";
                tmpString = tmpString + tmpChar;
            }
            streamMutsEachCellResults<<tmpString+";";
        }
    }
    streamMutsEachCellResults<<"\n";

    exportMutsEachCellResults.close();
}

void CellSystem::show()
{
    qDebug()<<"Size: "<<historyMutationCounters.size();
    for(unsigned int i=0;i<myMutationTable->getMaxIdMutation();i++)
        qDebug()<<qPrintable(QString::fromStdString(myMutationTable->getMutation(i)->mName))<<": "<<historyMutationCounters.back().getNumberGivenMut(i);
    qDebug()<<"1 mutations: "<<historyMutationHistogram.back().getNumberOfSamples(1);
}
