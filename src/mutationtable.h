#ifndef MUTATIONTABLE_H
#define MUTATIONTABLE_H

#include <vector>
#include <string>
#include <QString>
#include <algorithm>
#include <QDebug>


enum ModifierType{
    ADD,
    MULT,
    NONE
};

template <typename T>
const bool Contains( std::vector<T>& Vec, const T& Element )
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
        return true;

    return false;
}

class MutationGene{
public:
   static int numberOfMutationGenes;
   std::string mName;
   int mId;
   MutationGene(){
        mId = numberOfMutationGenes;
        numberOfMutationGenes++;
   }
   MutationGene(std::string tname){
       mName = tname;
       mId = numberOfMutationGenes;
       numberOfMutationGenes++;
   }
};

class MutationEvent{
public:

   ModifierType prolModifierType;
   float prolRate;

   ModifierType deathModifierType;
   float deathRate;

   ModifierType moreMutModifierType;
   float moreMutRateMutations;

   ModifierType telomeresModifierType;
   float telomeresRate;

   float dominance; //if 0, it's oncogene;  if 1, it's TSG
   float microEnvironment;

   std::string mName;
   unsigned int oncogenic_size;
   std::vector<float> sinergyProlValues; //key are before mutation events
   std::vector<float> sinergyDeathValues;
   std::vector<float> sinergyTelValues;
   std::vector<unsigned int> listOfChainReactionlMutationEvents;
   MutationGene *mMutationGene;

   //Initializing with default values
   MutationEvent(){
       prolModifierType = ModifierType::MULT;
       prolRate=1;

       deathModifierType = ModifierType::MULT;
       deathRate=1;

       moreMutModifierType = ModifierType::MULT;
       moreMutRateMutations=1;

       telomeresModifierType = ModifierType::MULT;
       telomeresRate=1;

       dominance = 0;
       microEnvironment = 1;

       oncogenic_size=0;
       mName="";
       sinergyProlValues = std::vector<float>();
       sinergyDeathValues = std::vector<float>();
       sinergyTelValues = std::vector<float>();
       listOfChainReactionlMutationEvents = std::vector<unsigned int>();
       mMutationGene = NULL;
   }

   void resetSinergyVector(unsigned int numberOfTotalEvents){
       sinergyProlValues.clear();
       sinergyDeathValues.clear();
       sinergyTelValues.clear();
       sinergyProlValues = std::vector<float>(numberOfTotalEvents,1);
       sinergyDeathValues = std::vector<float>(numberOfTotalEvents,1);
       sinergyTelValues = std::vector<float>(numberOfTotalEvents,1);
   }

   void resetChainReactionVector(){
       listOfChainReactionlMutationEvents.clear();
   }


   void setSinergyValue(unsigned int tBeforeMutEvent, float tProlSinergy, float tDeathSinergy, float tTelSinergy){
       if (tBeforeMutEvent<sinergyProlValues.size()){
           sinergyProlValues[tBeforeMutEvent] = tProlSinergy;
           sinergyDeathValues[tBeforeMutEvent] = tDeathSinergy;
           sinergyTelValues[tBeforeMutEvent] = tTelSinergy;
       }
   }
   void addChainReactionMutationEvent(unsigned int tAfterMutEvent){
       listOfChainReactionlMutationEvents.push_back(tAfterMutEvent);
   }


   float computeProlSinergyModifier(std::vector<unsigned int> listOfMutations);
   float computeProlSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape);
   float computeDeathSinergyModifier(std::vector<unsigned int> listOfMutations);
   float computeDeathSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape);
   float computeTelSinergyModifier(std::vector<unsigned int> listOfMutations);
   float computeTelSinergyModifier(std::vector<unsigned int> tFirstTape,std::vector<unsigned int> tSecondTape);
   float computeNewProlRate(float celProlRate, float sinergyModifier, int activationLevel);
   float computeNewDeathRate(float celDeathRate, float sinergyModifier,int activationLevel);
   float computeNewMoreMutRate(float celMoreMutRate, float sinergyModifier);
   float computeNewTelomeresRate(float celTelomeresRate, float sinergyModifier, int activationLevel);
   float computeNewMicroEnvironmentModifier(float curMicroEnvironment, int activationLevel);

};

class MutationTable
{
private:
    std::vector<MutationEvent *> mMutationEvents;
    std::vector<MutationGene *> mMutationGenes;

    //Classification
    unsigned int total_oncogenic_size;
    unsigned int maxIdMutationEvent;
    unsigned int maxIdMutationGene;
    std::vector<int> oncogenic_tape;
    std::vector<int> oncogenic_tape_genes;
    int genome_size;
public:
    MutationTable();
    void loadMutations(QString fileName);
    void loadSinergyPairs(QString fileName);
    int getMutGeneIdFromName(std::string tname);
    int getEventIdFromFullName(std::string tGeneName, std::string tEventName);

    unsigned int getMaxIdMutationEvent(){ return maxIdMutationEvent; }
    unsigned int getMaxIdMutationGene(){ return maxIdMutationGene; }
    int getNumberofMutationEvent(unsigned int tbasis);
    int getNumberOfMutationGene(unsigned int tbasis);
    MutationGene *getMutationGene(int tid) const;
    MutationEvent *getMutationEvent(int tid) const;
    int getGenomeSize() const {return genome_size;}
    QStringList getMutationListNames();
    QStringList getMutationGenesListNames();
    int getGeneOfEvent(int teventId);

    //bool checkIfAddMutation(int tmpMutationId,std::vector<unsigned int> newFirstTapeMutations,std::vector<unsigned int> newSecondTapeMutations,unsigned int tapeChosen);
    int checkIfAddMutationByReference(int tmpMutationId,std::vector<unsigned int> &newFirstTapeMutations,std::vector<unsigned int> &newSecondTapeMutations,unsigned int tapeChosen); //0=noActivation 1=partialActivation 2=fullActivation
    //bool checkIfAddMutationArrayOfBool(int tmpMutationId,std::vector<bool> newFirstTapeMutations,std::vector<bool> newSecondTapeMutations,unsigned int tapeChosen);

};

#endif // MUTATIONTABLE_H
