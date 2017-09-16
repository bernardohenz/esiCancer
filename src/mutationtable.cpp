#include "mutationtable.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>


MutationTable::MutationTable()
{
    genome_size = 4568;
    maxIdMutation=0;
    total_oncogenic_size=0;
}

void MutationTable::loadMutations(QString fileName)
{
    if (mMutations.size()>0){
        mMutations.clear();
        maxIdMutation = 0;
        total_oncogenic_size=0;
    }
    QFile fileTable(fileName);
    if(!(fileTable.open(QIODevice::ReadOnly))){
        qDebug()<<"Error opening the files";
        return;
    }
    QTextStream fileTextStream(&fileTable);
    QString tmpTxt, tmpTxtId;
    QStringList tmpstringList;
    bool flag_first_line_read = false;
    unsigned int id_name = 0;
    unsigned int id_proliferate_without_sup=0, id_proliferate_with_sup=0,id_proliferate_func=0;
    unsigned int id_die_without_sup=0, id_die_with_sup=0, id_death_func=0;
    unsigned int id_mut=0, id_mut_func=0;
    unsigned int id_oncogenic_size,id_genome_size=0;
    unsigned int id_mutation_type=0,id_telomeres_without_sup=0,id_telomeres_with_sup=0,id_telomeres_func=0;
    //unsigned int id_prolOnTreat,id_deathOnTreat;
    unsigned int id_synergy;

    while (true) {
        if (fileTextStream.atEnd()){
            break;
        }
        tmpTxt = fileTextStream.readLine().toLower();
        tmpstringList = tmpTxt.replace(',',';').split(';');
        if (tmpstringList[0] == "basic"){
            for(int i=1; i<tmpstringList.size(); i++){
                if ((tmpstringList[i]=="gene name")||(tmpstringList[i]=="gene"))
                    id_name=i;
                else if (tmpstringList[i]=="prowotsg")
                    id_proliferate_without_sup=i;
                else if(tmpstringList[i]=="prowitsg")
                    id_proliferate_with_sup=i;
                else if(tmpstringList[i]=="profunc")
                    id_proliferate_func=i;
                else if(tmpstringList[i]=="synergy")
                    id_synergy=i;
                else if(tmpstringList[i]=="dewotsg")
                    id_die_without_sup=i;
                else if(tmpstringList[i]=="dewitsg")
                    id_die_with_sup=i;
                else if(tmpstringList[i]=="defunc")
                    id_death_func=i;
                else if(tmpstringList[i]=="mutfunc")
                    id_mut_func=i;
                else if(tmpstringList[i]=="mut")
                    id_mut=i;
                else if(tmpstringList[i]=="telfunc")
                    id_telomeres_func=i;
                else if(tmpstringList[i]=="telwotsg")
                    id_telomeres_without_sup=i;
                else if(tmpstringList[i]=="telwitsg")
                    id_telomeres_with_sup=i;
                else if(tmpstringList[i]=="type")
                    id_mutation_type=i;
                else if(tmpstringList[i]=="oncbase")
                    id_oncogenic_size=i;
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
        Mutation *tmpMutation = new Mutation();
        if(id_proliferate_func>0){
            if (tmpstringList[id_proliferate_func] == "mult") {
                tmpMutation->prolRateMultBeforeSuppressor = tmpstringList[id_proliferate_without_sup].toFloat();
                tmpMutation->prolRateMultAfterSuppressor = tmpstringList[id_proliferate_with_sup].toFloat();
            } else if (tmpstringList[id_proliferate_func] == "add"){
                tmpMutation->prolRateAddBeforeSuppressor = tmpstringList[id_proliferate_without_sup].toFloat();
                tmpMutation->prolRateAddAfterSuppressor = tmpstringList[id_proliferate_with_sup].toFloat();
            }
        }
        if(id_synergy>0)
            tmpMutation->synergy = tmpstringList[id_synergy].toFloat();

        if(id_death_func>0){
            if (tmpstringList[id_death_func] == "mult") {
                tmpMutation->deathRateMultBeforeSuppressor = tmpstringList[id_die_without_sup].toFloat();
                tmpMutation->deathRateMultAfterSuppressor = tmpstringList[id_die_with_sup].toFloat();
            } else if (tmpstringList[id_death_func] == "add"){
                tmpMutation->deathRateAddBeforeSuppressor = tmpstringList[id_die_without_sup].toFloat();
                tmpMutation->deathRateAddAfterSuppressor = tmpstringList[id_die_with_sup].toFloat();
            }
        }
        if(id_mut>0){
            if(tmpstringList[id_mut_func]=="mult")
                tmpMutation->moreMutRateMultMutations = tmpstringList[id_mut].toFloat();
            else if (tmpstringList[id_mut_func]=="add")
                tmpMutation->moreMutRateAddMutations = tmpstringList[id_mut].toFloat();
        }

        if(id_telomeres_func>0){
            if(tmpstringList[id_telomeres_func]=="mult"){
                tmpMutation->telomeresRateMultBeforeSuppressor = tmpstringList[id_telomeres_without_sup].toFloat();
                tmpMutation->telomeresRateMultAfterSuppressor = tmpstringList[id_telomeres_with_sup].toFloat();
            } else if(tmpstringList[id_telomeres_func]=="add"){
                tmpMutation->telomeresRateAddBeforeSuppressor = tmpstringList[id_telomeres_without_sup].toFloat();
                tmpMutation->telomeresRateAddAfterSuppressor = tmpstringList[id_telomeres_with_sup].toFloat();
            }
        }
        if(id_name>0){
            tmpMutation->mName = tmpstringList[id_name].toStdString();
        }
        if(id_mutation_type>0){
            if ((tmpstringList[id_mutation_type]=="onc")||(tmpstringList[id_mutation_type]=="oncogene")||(tmpstringList[id_mutation_type]=="normd"))
                tmpMutation->mType = MutationType::ONCOGENE;
            else if ((tmpstringList[id_mutation_type]=="tsg")||(tmpstringList[id_mutation_type]=="tumor_suppressor_gene")||(tmpstringList[id_mutation_type]=="normr"))
                tmpMutation->mType = MutationType::TUMORSUPPRESSOR;
        }
        if(id_oncogenic_size>0){

            tmpMutation->oncogenic_size = tmpstringList[id_oncogenic_size].toFloat();
        }
        total_oncogenic_size += tmpMutation->oncogenic_size;
        maxIdMutation++;
        mMutations.push_back(tmpMutation);
    }
    fileTable.close();
    //initialize oncogenic tape
    oncogenic_tape = std::vector<int>(total_oncogenic_size,0);
    int tmpcounter=0;
    for (unsigned int i=0;i<mMutations.size();i++){
        for(unsigned int j=0;j<mMutations[i]->oncogenic_size;j++)
            oncogenic_tape[tmpcounter+j] = i;
        tmpcounter+=mMutations[i]->oncogenic_size;
        //qDebug()<<mMutations[i]->oncogenic_size<<": ";
    }
    qDebug()<<"genome: "<<genome_size;
    qDebug()<<"oncogenic_size: "<<total_oncogenic_size;
    for (unsigned int i=0;i<total_oncogenic_size;i++)
        qDebug()<<oncogenic_tape[i];

}

int MutationTable::getNumberofMutation(unsigned int tbasis)
{
    if(tbasis<total_oncogenic_size){
        if(oncogenic_tape[tbasis]>=2)
            qDebug()<<"fqfdqwdwq";
        return oncogenic_tape[tbasis];
    }
    return -1;
}

MutationType MutationTable::getTypeOfMutationGivenAMutation(int tmut)
{
    if ((unsigned int) tmut<mMutations.size()){
        return mMutations[tmut]->mType;
    }
    return MutationType::NOMUTATION;
}

Mutation *MutationTable::getMutation(int tid) const
{
    if((unsigned int)tid<mMutations.size())
        return mMutations[tid];
    return NULL;
}

QStringList MutationTable::getMutationListNames()
{
    QStringList tmpList;
    for (unsigned int i=0;i<mMutations.size();i++)
        tmpList.push_back(QString::fromStdString(mMutations[i]->mName));
    return tmpList;
}

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

bool MutationTable::checkIfAddMutationByReference(int tmpMutationId, std::vector<unsigned int> &newFirstTapeMutations, std::vector<unsigned int> &newSecondTapeMutations, unsigned int tapeChosen)
{
    if((unsigned int)tmpMutationId>=mMutations.size())
        return false;
    MutationType tmpType = mMutations[tmpMutationId]->mType;
    if(tapeChosen==1){
        if (Contains(newFirstTapeMutations,(unsigned int)tmpMutationId))
            return false;
        else{
            newFirstTapeMutations.push_back(tmpMutationId);
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
            newSecondTapeMutations.push_back(tmpMutationId);
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



float MutationTable::getSynergyOfMutation(int tmpMutationId)
{
    if((unsigned int)tmpMutationId<mMutations.size())
        return mMutations[tmpMutationId]->synergy;
    return 1;
}

