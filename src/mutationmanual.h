#ifndef MUTATIONMANUAL_H
#define MUTATIONMANUAL_H

#include <QString>

enum manualMutationAllelsChanged{
    MONOALLELIC,
    BIALLELIC
};

class MutationManual
{
private:
    unsigned int generation;
    float percentage;
    unsigned int mutationEventId;
    QString mutationGeneName;
    QString mutationEventName;
    QString mutationFullName;
    manualMutationAllelsChanged mutationAllelsChanged;
public:
    MutationManual(unsigned int tgeneration, float tpercentage, unsigned int tmutationEventId,QString tmutationFullName,manualMutationAllelsChanged tmutationAllelsChanged = MONOALLELIC);

    unsigned int getGeneration() const;
    void setGeneration(unsigned int value);
    float getPercentage() const;
    void setPercentage(float value);
    unsigned int getMutationEventId() const;
    void setMutationEventId(int tvalue);
    void setMutationFullName(QString tname);
    QString getMutationGeneName() const;
    QString getMutationEventName() const;
    manualMutationAllelsChanged getMutationAllelsChanged() const;
    void setMutationAllelsChanged(const manualMutationAllelsChanged &value);
};

#endif // MUTATIONMANUAL_H
