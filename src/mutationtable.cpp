#include "mutationtable.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>

int MutationGene::numberOfMutationGenes=0;

MutationTable::MutationTable()
{
    genome_size = 4568;
    maxIdMutationEvent=0;
    maxIdMutationGene=0;
    total_oncogenic_size=0;
}

void MutationTable::loadMutations(QString fileName)
{
    if (mMutationEvents.size()>0){
        mMutationEvents.clear();
        maxIdMutationEvent = 0;
        total_oncogenic_size=0;
    }
    if(mMutationGenes.size()>0){
        MutationGene::numberOfMutationGenes=0;
        mMutationGenes.clear();
        maxIdMutationGene = 0;
    }
    QFile fileTable(fileName);
    if(!(fileTable.open(QIODevice::ReadOnly))){
        qDebug()<<"Error opening the files";
        return;
    }
    QTextStream fileTextStream(&fileTable);
    QString tmpTxt, tmpTxtMut;
    QStringList tmpstringList;
    bool flag_first_line_read = false;
    unsigned int id_gene_name = 0, id_event_name = 0;
    unsigned int id_proliferate_rate=0,id_proliferate_func=0;
    unsigned int id_death_rate=0, id_death_func=0;
    unsigned int id_mut=0, id_mut_func=0;
    unsigned int id_oncogenic_size,id_genome_size=0;
    unsigned int id_telomeres_rate=0,id_telomeres_func=0;
    //unsigned int id_prolOnTreat,id_deathOnTreat;
    unsigned int id_dominance=0,id_micro_environment=0,id_synergy;

    while (true) {
        if (fileTextStream.atEnd()){
            break;
        }
        tmpTxt = fileTextStream.readLine().toLower();
        tmpstringList = tmpTxt.replace(',',';').split(';');
        if (tmpstringList[0] == "basic"){
            for(int i=1; i<tmpstringList.size(); i++){
                if ((tmpstringList[i]=="gene name")||(tmpstringList[i]=="genes"))
                    id_gene_name=i;
                else if (tmpstringList[i]=="events")
                    id_event_name=i;
                else if ((tmpstringList[i]=="prorate")||(tmpstringList[i]=="divrate"))
                    id_proliferate_rate=i;
                else if ((tmpstringList[i]=="profunc")||(tmpstringList[i]=="divfunc"))
                    id_proliferate_func=i;
                else if(tmpstringList[i]=="synergy")
                    id_synergy=i;
                else if(tmpstringList[i]=="deathrate")
                    id_death_rate=i;
                else if(tmpstringList[i]=="defunc")
                    id_death_func=i;
                else if(tmpstringList[i]=="mutfunc")
                    id_mut_func=i;
                else if(tmpstringList[i]=="mut")
                    id_mut=i;
                else if((tmpstringList[i]=="telfunc")||(tmpstringList[i]=="maxdivfunc"))
                    id_telomeres_func=i;
                else if((tmpstringList[i]=="telrate")||(tmpstringList[i]=="maxdivrate"))
                    id_telomeres_rate=i;
                else if((tmpstringList[i]=="oncbase")||(tmpstringList[i]=="oncevents")||(tmpstringList[i]=="probevent"))
                    id_oncogenic_size=i;
                else if(tmpstringList[i]=="dominance")
                    id_dominance=i;
                else if(tmpstringList[i]=="microenvironment")
                    id_micro_environment=i;
                else if(tmpstringList[i]=="genome size")
                    id_genome_size=i;
            }
            continue;
        }
        if(!flag_first_line_read){
            if (id_genome_size>0){
                genome_size = tmpstringList[id_genome_size].toInt();
            }
            flag_first_line_read = true;
        }
        MutationEvent *tmpMutationEvent = new MutationEvent();
        if(id_proliferate_func>0){
            if (tmpstringList[id_proliferate_func] == "mult")
                tmpMutationEvent->prolModifierType = ModifierType::MULT;
            else if (tmpstringList[id_proliferate_func] == "add")
                tmpMutationEvent->prolModifierType = ModifierType::ADD;
            tmpMutationEvent->prolRate = tmpstringList[id_proliferate_rate].toFloat();
        }

        if(id_death_func>0){
            if (tmpstringList[id_death_func] == "mult")
                tmpMutationEvent->deathModifierType = ModifierType::MULT;
            else if (tmpstringList[id_death_func] == "add")
                tmpMutationEvent->deathModifierType = ModifierType::ADD;
            tmpMutationEvent->deathRate = tmpstringList[id_death_rate].toFloat();
        }
        if(id_mut>0){
            if(tmpstringList[id_mut_func]=="mult")
                tmpMutationEvent->moreMutModifierType = ModifierType::MULT;
            else if (tmpstringList[id_mut_func]=="add")
                tmpMutationEvent->moreMutModifierType = ModifierType::ADD;
            tmpMutationEvent->moreMutRateMutations = tmpstringList[id_mut].toFloat();
        }

        if(id_telomeres_func>0){
            if(tmpstringList[id_telomeres_func]=="mult")
                tmpMutationEvent->telomeresModifierType = ModifierType::MULT;
            else if (tmpstringList[id_telomeres_func]=="add")
                tmpMutationEvent->telomeresModifierType = ModifierType::ADD;
            tmpMutationEvent->telomeresRate = tmpstringList[id_telomeres_rate].toFloat();
        }
        if(id_gene_name>0){
            tmpTxtMut = tmpstringList[id_gene_name];
            bool tmpExistGene = false;
            unsigned int tmpIndex = 0;
            for (unsigned int i=0; i<mMutationGenes.size();i++){
                if (mMutationGenes[i]->mName == tmpTxtMut.toStdString()){
                    tmpExistGene=true;
                    tmpIndex=i;
                    break;
                }
            }
            if(tmpExistGene){
                tmpMutationEvent->mMutationGene = mMutationGenes[tmpIndex];
            }
            else{
                MutationGene *tmpMutationGene = new MutationGene(tmpTxtMut.toStdString());
                tmpMutationEvent->mMutationGene = tmpMutationGene;
                mMutationGenes.push_back(tmpMutationGene);
                maxIdMutationGene++;
            }
        }
        if(id_event_name>0){
            tmpMutationEvent->mName = tmpstringList[id_event_name].toStdString();
        }

        if(id_dominance>0){
            tmpMutationEvent->dominance = tmpstringList[id_dominance].toFloat();
        }

        if(id_micro_environment>0){
            tmpMutationEvent->microEnvironment = tmpstringList[id_micro_environment].toFloat();
        }

        if(id_oncogenic_size>0){
            tmpMutationEvent->oncogenic_size = tmpstringList[id_oncogenic_size].toFloat();
        }
        total_oncogenic_size += tmpMutationEvent->oncogenic_size;
        maxIdMutationEvent++;
        mMutationEvents.push_back(tmpMutationEvent);
    }
    fileTable.close();
    //initialize oncogenic tape
    oncogenic_tape = std::vector<int>(total_oncogenic_size,0);
    oncogenic_tape_genes = std::vector<int>(total_oncogenic_size,0);
    int tmpcounter=0;
    for (unsigned int i=0;i<mMutationEvents.size();i++){
        for(unsigned int j=0;j<mMutationEvents[i]->oncogenic_size;j++){
            oncogenic_tape[tmpcounter+j] = i;
            oncogenic_tape_genes[tmpcounter+j] = mMutationEvents[i]->mMutationGene->mId;
        }
        tmpcounter+=mMutationEvents[i]->oncogenic_size;
        //qDebug()<<mMutations[i]->oncogenic_size<<": ";
    }

    //Resetting sinergy vectors
    for (unsigned int i=0;i<mMutationEvents.size();i++){
        mMutationEvents[i]->resetSinergyVector(maxIdMutationEvent);
    }

    qDebug()<<"Number of mutations: "<<mMutationEvents.size();
    qDebug()<<"genome: "<<genome_size;
    qDebug()<<"oncogenic_size: "<<total_oncogenic_size;
    //for (unsigned int i=0;i<total_oncogenic_size;i++)
    //    qDebug()<<oncogenic_tape[i];

    for (unsigned int i=0;i<mMutationEvents.size();i++)
        qDebug()<<"Gene: "<<QString::fromStdString(mMutationEvents[i]->mMutationGene->mName)<<" - "<<QString::fromStdString(mMutationEvents[i]->mName);
}

void MutationTable::loadSinergyPairs(QString fileName)
{
    QFile fileTable(fileName);
    if(!(fileTable.open(QIODevice::ReadOnly))){
        qDebug()<<"Error opening the file";
        return;
    }
    qDebug()<<"Loading sinergyTable";
    for (unsigned int i=0;i<mMutationEvents.size();i++){
        mMutationEvents[i]->resetSinergyVector(maxIdMutationEvent);
    }
    QTextStream fileTextStream(&fileTable);
    QString tmpTxt;
    QStringList tmpstringList;
    unsigned int id_csv_before_gene,id_csv_before_event, id_csv_after_gene, id_csv_after_event, id_csv_sinergy_prol_value,id_csv_sinergy_death_value,id_csv_sinergy_tel_value,id_csv_chain_reaction;
    int beforeMutEventId, afterMutEventId;
    bool flag_first_line_read = false;
    float tmp_prol_sinergy_value,tmp_death_sinergy_value,tmp_tel_sinergy_value;
    while (true) {
        if (fileTextStream.atEnd()){
            break;
        }
        tmpTxt = fileTextStream.readLine().toLower();
        tmpstringList = tmpTxt.replace(',',';').split(';');
        if(!flag_first_line_read){
            for(int i=0; i<tmpstringList.size(); i++){
                if (tmpstringList[i]=="gene_before")
                    id_csv_before_gene=i;
                else if (tmpstringList[i]=="event_before")
                    id_csv_before_event=i;
                else if (tmpstringList[i]=="gene_after")
                    id_csv_after_gene=i;
                else if (tmpstringList[i]=="event_after")
                    id_csv_after_event=i;
                else if ((tmpstringList[i]=="promod")||(tmpstringList[i]=="divmod"))
                    id_csv_sinergy_prol_value=i;
                else if (tmpstringList[i]=="demod")
                    id_csv_sinergy_death_value=i;
                else if ((tmpstringList[i]=="telmod")||(tmpstringList[i]=="maxdivmod"))
                    id_csv_sinergy_tel_value=i;
                else if ((tmpstringList[i]=="chain_reaction")||(tmpstringList[i]=="link"))
                    id_csv_chain_reaction = i;
            }
            flag_first_line_read = true;
            continue;
        }
        beforeMutEventId = getEventIdFromFullName(tmpstringList[id_csv_before_gene].toStdString(),tmpstringList[id_csv_before_event].toStdString()  );
        afterMutEventId = getEventIdFromFullName(tmpstringList[id_csv_after_gene].toStdString(),tmpstringList[id_csv_after_event].toStdString()  );
        tmp_prol_sinergy_value = tmpstringList[id_csv_sinergy_prol_value].toFloat();
        tmp_death_sinergy_value = tmpstringList[id_csv_sinergy_death_value].toFloat();
        tmp_tel_sinergy_value = tmpstringList[id_csv_sinergy_tel_value].toFloat();
        /*qDebug()<<"trying: "<<tmpstringList[id_csv_before_gene]<<"-"<<tmpstringList[id_csv_before_event]<<
                  " to "<<tmpstringList[id_csv_after_gene]<<"-"<<tmpstringList[id_csv_after_event];
        qDebug()<<"IDS: "<<beforeMutEventId<<"  "<<afterMutEventId;*/
        if ((beforeMutEventId>=0)&&(beforeMutEventId<maxIdMutationEvent)&&(afterMutEventId>=0)&&(afterMutEventId<maxIdMutationEvent)){
            qDebug()<<"Adding sinergy between "<<beforeMutEventId<<" and "<<afterMutEventId<<": "<<tmp_prol_sinergy_value<<" ; "<<tmp_death_sinergy_value<<" ; "<<tmp_tel_sinergy_value;
            mMutationEvents[afterMutEventId]->setSinergyValue(beforeMutEventId,tmp_prol_sinergy_value,tmp_death_sinergy_value,tmp_tel_sinergy_value);
            if ( tmpstringList[id_csv_chain_reaction].toFloat()>0.5 ){
                qDebug()<<"Adding chainReaction between "<<beforeMutEventId<<" and "<<afterMutEventId;
                mMutationEvents[beforeMutEventId]->addChainReactionMutationEvent(afterMutEventId);
            }
        }
    }
    fileTable.close();
    qDebug()<<"Finished loading sinergyTable";
}

int MutationTable::getMutGeneIdFromName(std::string tname)
{
    int id = -1;
    for (unsigned int i=0;i<mMutationGenes.size();i++){
        if ( mMutationGenes[i]->mName == tname ){
            id = i;
            break;
        }
    }
    return id;
}

int MutationTable::getEventIdFromFullName(std::string tGeneName, std::string tEventName)
{
    int id = -1;
    for (unsigned int i=0;i<mMutationEvents.size();i++){
        if ( mMutationEvents[i]->mName == tEventName ){
            if (mMutationEvents[i]->mMutationGene != NULL){
                if((mMutationEvents[i]->mMutationGene->mName == tGeneName)){
                    id = i;
                    break;
                }
            }
        }
    }
    return id;
}


int MutationTable::getNumberofMutationEvent(unsigned int tbasis)
{
    if(tbasis<0)
        qDebug()<<"ERROR:--------- TBASIS NEGATIVE =====";
    if(tbasis<total_oncogenic_size){
        return oncogenic_tape[tbasis];
    }
    return -1;
}

int MutationTable::getNumberOfMutationGene(unsigned int tbasis)
{
    if(tbasis<0)
        qDebug()<<"ERROR:--------- TBASIS NEGATIVE =====";
    if(tbasis<total_oncogenic_size){
        return oncogenic_tape_genes[tbasis];
    }
    return -1;
}

MutationGene *MutationTable::getMutationGene(int tid) const
{
    if((unsigned int)tid<mMutationGenes.size())
        return mMutationGenes[tid];
    return NULL;
}

MutationEvent *MutationTable::getMutationEvent(int tid) const
{
    if((unsigned int)tid<mMutationEvents.size())
        return mMutationEvents[tid];
    return NULL;
}

QStringList MutationTable::getMutationListNames()
{
    QStringList tmpList;
    for (unsigned int i=0;i<mMutationEvents.size();i++)
        tmpList.push_back(QString::fromStdString(mMutationEvents[i]->mMutationGene->mName)+" - "+QString::fromStdString(mMutationEvents[i]->mName));
    return tmpList;
}

QStringList MutationTable::getMutationGenesListNames()
{
    QStringList tmpList;
    for (unsigned int i=0;i<mMutationGenes.size();i++)
        tmpList.push_back(QString::fromStdString(mMutationGenes[i]->mName));
    return tmpList;
}

int MutationTable::getGeneOfEvent(int teventId)
{
    if((unsigned int)teventId<mMutationEvents.size())
        return mMutationEvents[teventId]->mMutationGene->mId;
    return -1;
}
/*
bool MutationTable::checkIfAddMutationArrayOfBool(int tmpMutationId, std::vector<bool> newFirstTapeMutations, std::vector<bool> newSecondTapeMutations, unsigned int tapeChosen)
{
    if((unsigned int)tmpMutationId>=mMutations.size())
        return false;
    MutationType tmpType = mMutations[tmpMutationId]->mType;
    if(tapeChosen==1){
        if(newFirstTapeMutations[tmpMutationId]==true)
            return false;
        else{
            if(tmpType==MutationType::ONCOGENE){
                if(newSecondTapeMutations[tmpMutationId]==false)
                    return true;
            } else if (tmpType==MutationType::TUMORSUPPRESSOR){
                if(newSecondTapeMutations[tmpMutationId]==true)
                    return true;
            }
        }
    }else if(tapeChosen==2){
        if(newSecondTapeMutations[tmpMutationId]==true)
            return false;
        else{
            if(tmpType==MutationType::ONCOGENE){
                if(newFirstTapeMutations[tmpMutationId]==false)
                    return true;
            } else if (tmpType==MutationType::TUMORSUPPRESSOR){
                if(newFirstTapeMutations[tmpMutationId]==true)
                    return true;
            }
        }
    }
    return false;
}

bool MutationTable::checkIfAddMutation(int tmpMutationId, std::vector<unsigned int> newFirstTapeMutations, std::vector<unsigned int> newSecondTapeMutations, unsigned int tapeChosen)
{
    if((unsigned int)tmpMutationId>=mMutations.size())
        return false;
    MutationType tmpType = mMutations[tmpMutationId]->mType;
    if(tapeChosen==1){
        if (Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
            return false;
        else{
            if(tmpType==MutationType::ONCOGENE){
                if(!Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
                    return true;
            } else if (tmpType==MutationType::TUMORSUPPRESSOR){
                if (Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
                    return true;
            }
        }
    }else if(tapeChosen==2){
        if (Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
            return false;
        else{
            if(tmpType==MutationType::ONCOGENE){
                if(!Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
                    return true;
            } else if (tmpType==MutationType::TUMORSUPPRESSOR){
                if (Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
                    return true;
            }
        }
    }
    return false;
}
*/

int MutationTable::checkIfAddMutationByReference(int tmpMutationId, std::vector<unsigned int> &newFirstTapeMutations, std::vector<unsigned int> &newSecondTapeMutations, unsigned int tapeChosen)
{
    if((unsigned int)tmpMutationId>=mMutationGenes.size())
        return 0;
    if(tapeChosen==1){
        if (Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
            return 0;
        else{
            newFirstTapeMutations.push_back(tmpMutationId);
            if(Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
                return 2;
            else
                return 1;
        }
    }else if(tapeChosen==2){
        if (Contains(newSecondTapeMutations,(unsigned int)tmpMutationId))
            return 0;
        else{
            newSecondTapeMutations.push_back(tmpMutationId);
            if(Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
                return 2;
            else
                return 1;
        }
    }
    return 0;
}





float MutationEvent::computeProlSinergyModifier(std::vector<unsigned int> listOfMutations)
{
    float modifier=1;
    for (unsigned int i=0; i<listOfMutations.size(); i++){
        modifier*= sinergyProlValues[listOfMutations[i]];
    }
    return modifier;
}

float MutationEvent::computeProlSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape)
{
    float modifier=1;
    for (unsigned int i=0; i<tFirstTape.size(); i++){
        modifier*= sinergyProlValues[tFirstTape[i]];
    }
    for (unsigned int i=0; i<tSecondTape.size(); i++){
        if (!Contains(tFirstTape,tSecondTape[i]))
            modifier*= sinergyProlValues[tSecondTape[i]];
    }
    return modifier;
}

float MutationEvent::computeDeathSinergyModifier(std::vector<unsigned int> listOfMutations)
{
    float modifier=1;
    for (unsigned int i=0; i<listOfMutations.size(); i++){
        modifier*= sinergyDeathValues[listOfMutations[i]];
    }
    return modifier;
}

float MutationEvent::computeDeathSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape)
{
    float modifier=1;
    for (unsigned int i=0; i<tFirstTape.size(); i++){
        modifier*= sinergyDeathValues[tFirstTape[i]];
    }
    for (unsigned int i=0; i<tSecondTape.size(); i++){
        if (!Contains(tFirstTape,tSecondTape[i]))
            modifier*= sinergyDeathValues[tSecondTape[i]];
    }
    return modifier;
}

float MutationEvent::computeTelSinergyModifier(std::vector<unsigned int> listOfMutations)
{
    float modifier=1;
    for (unsigned int i=0; i<listOfMutations.size(); i++){
        modifier*= sinergyTelValues[listOfMutations[i]];
    }
    return modifier;
}

float MutationEvent::computeTelSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape)
{
    float modifier=1;
    for (unsigned int i=0; i<tFirstTape.size(); i++){
        modifier*= sinergyTelValues[tFirstTape[i]];
    }
    for (unsigned int i=0; i<tSecondTape.size(); i++){
        if (!Contains(tFirstTape,tSecondTape[i]))
            modifier*= sinergyTelValues[tSecondTape[i]];
    }
    return modifier;
}


float MutationEvent::computeNewProlRate(float celProlRate, float sinergyModifier, int activationLevel)
{
    float newProlRate;
    if (activationLevel==1)
        newProlRate = (prolRate-1)*dominance+1;
    else if (activationLevel==2)
        newProlRate = (prolRate-1)*(1-dominance)+1;
    else if (activationLevel==3)
        newProlRate = prolRate;
    else
        return celProlRate;

    if (prolModifierType == ModifierType::ADD){
        newProlRate = ((newProlRate-1)*sinergyModifier + 1) + (celProlRate-1);
    } else if (prolModifierType == ModifierType::MULT){
        newProlRate = ((newProlRate-1)*sinergyModifier + 1) * celProlRate;
    }else
        newProlRate = celProlRate;

    /*if (sinergyModifier>1){
        qDebug()<<"sinergyModifier: "<<sinergyModifier;
        qDebug()<<"celProlRate: "<<celProlRate;
        qDebug()<<"resulting prol rate: "<<newProlRate;
    }*/


    return newProlRate;
}

float MutationEvent::computeNewDeathRate(float celDeathRate, float sinergyModifier, int activationLevel)
{
    float newDeathRate;

    if (activationLevel==1)
        newDeathRate = (deathRate-1)*dominance+1;
    else if (activationLevel==2)
        newDeathRate = (deathRate-1)*(1-dominance)+1;
    else if (activationLevel==3)
        newDeathRate = deathRate;
    else
        return celDeathRate;


    if (deathModifierType == ModifierType::ADD){
        newDeathRate = ((newDeathRate-1)*sinergyModifier+1) + (celDeathRate-1);
    } else if (deathModifierType == ModifierType::MULT){
        newDeathRate = ((newDeathRate-1)*sinergyModifier+1) * celDeathRate;
    }else
        newDeathRate = celDeathRate;


    return newDeathRate;
}

float MutationEvent::computeNewMoreMutRate(float celMoreMutRate, float sinergyModifier)
{
    float newMoreMutRate = moreMutRateMutations;

    if (moreMutModifierType == ModifierType::ADD){
        newMoreMutRate = ((newMoreMutRate-1)*sinergyModifier+1) + (celMoreMutRate-1);
    } else if (moreMutModifierType == ModifierType::MULT){
        newMoreMutRate = ((newMoreMutRate-1)*sinergyModifier+1) * celMoreMutRate;
    }else
        newMoreMutRate = celMoreMutRate;

    return newMoreMutRate;
}

float MutationEvent::computeNewTelomeresRate(float celTelomeresRate, float sinergyModifier, int activationLevel)
{
    float newTelomeresRate;

    if (activationLevel==1)
        newTelomeresRate = (telomeresRate-1)*dominance +1;
    else if (activationLevel==2)
        newTelomeresRate = (telomeresRate-1)*(1-dominance) +1;
    else if (activationLevel==3)
        newTelomeresRate = telomeresRate;
    else
        return celTelomeresRate;

    if (telomeresModifierType == ModifierType::ADD){
        newTelomeresRate = ((newTelomeresRate-1)*sinergyModifier+1) + (celTelomeresRate-1);
    } else if (telomeresModifierType == ModifierType::MULT){
        newTelomeresRate = ((newTelomeresRate-1)*sinergyModifier+1) * celTelomeresRate;
    }else
        newTelomeresRate = celTelomeresRate;

    return newTelomeresRate;
}

float MutationEvent::computeNewMicroEnvironmentModifier(float curMicroEnvironment,int activationLevel)
{
    float microEnvironModifier=1;
    if (activationLevel==1)
        microEnvironModifier = (microEnvironment-1)*dominance +1;
    else if (activationLevel==2)
        microEnvironModifier = (microEnvironment-1)*(1-dominance) +1;
    else if (activationLevel==3)
        microEnvironModifier = microEnvironment;

    return microEnvironModifier*curMicroEnvironment;
}
