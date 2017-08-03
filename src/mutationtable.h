#ifndef MUTATIONTABLE_H
#define MUTATIONTABLE_H

#include <vector>
#include <string>
#include <QString>
#include <algorithm>



enum MutationType{
    ONCOGENE,
    TUMORSUPPRESSOR,
    NOMUTATION
};

template <typename T>
const bool Contains( std::vector<T>& Vec, const T& Element )
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
        return true;

    return false;
}

class Mutation{
public:

   float prolRateMultBeforeSuppressor;
   float prolRateMultAfterSuppressor;
   float prolRateAddBeforeSuppressor;
   float prolRateAddAfterSuppressor;

   float deathRateMultBeforeSuppressor;
   float deathRateMultAfterSuppressor;
   float deathRateAddBeforeSuppressor;
   float deathRateAddAfterSuppressor;

   float moreMutRateMultMutations;
   float moreMutRateAddMutations;

   float telomeresRateMultBeforeSuppressor;
   float telomeresRateMultAfterSuppressor;
   float telomeresRateAddBeforeSuppressor;
   float telomeresRateAddAfterSuppressor;

   float synergy;

   MutationType mType;
   std::string mName;
   unsigned int oncogenic_size;

   //Initializing with default values
   Mutation(){
       prolRateMultBeforeSuppressor=1;
       prolRateMultAfterSuppressor=1;
       prolRateAddBeforeSuppressor=0;
       prolRateAddAfterSuppressor=0;
       deathRateMultBeforeSuppressor=1;
       deathRateMultAfterSuppressor=1;
       deathRateAddBeforeSuppressor=0;
       deathRateAddAfterSuppressor=0;
       moreMutRateMultMutations=1;
       moreMutRateAddMutations=0;
       telomeresRateMultBeforeSuppressor=1;
       telomeresRateMultAfterSuppressor=1;
       telomeresRateAddBeforeSuppressor=0;
       telomeresRateAddAfterSuppressor=0;
       synergy=1;
       oncogenic_size=0;
       mName="";
       mType=MutationType::NOMUTATION;
   }
};

class MutationTable
{
private:
    std::vector<Mutation *> mMutations;

    //Classification
    unsigned int total_oncogenic_size;
    unsigned int maxIdMutation;
    std::vector<int> oncogenic_tape;
    int genome_size;
public:
    MutationTable();
    void loadMutations(QString fileName);

    unsigned int getMaxIdMutation(){ return maxIdMutation; }
    int getNumberofMutation(unsigned int tbasis);
    MutationType getTypeOfMutationGivenAMutation(int tmut);
    Mutation *getMutation(int tid) const;
    int getGenomeSize() const {return genome_size;}
    QStringList getMutationListNames();

    bool checkIfAddMutation(int tmpMutationId,std::vector<unsigned int> newFirstTapeMutations,std::vector<unsigned int> newSecondTapeMutations,unsigned int tapeChosen);
    bool checkIfAddMutationByReference(int tmpMutationId,std::vector<unsigned int> &newFirstTapeMutations,std::vector<unsigned int> &newSecondTapeMutations,unsigned int tapeChosen);
    bool checkIfAddMutationArrayOfBool(int tmpMutationId,std::vector<bool> newFirstTapeMutations,std::vector<bool> newSecondTapeMutations,unsigned int tapeChosen);
    float getSynergyOfMutation(int tmpMutationId);
};

#endif // MUTATIONTABLE_H
