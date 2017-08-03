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
    unsigned int mutationId;
    QString mutationName;
    manualMutationAllelsChanged mutationAllelsChanged;
public:
    MutationManual(unsigned int tgeneration, float tpercentage, unsigned int tmutationId, QString tmutationName,manualMutationAllelsChanged tmutationAllelsChanged = MONOALLELIC);

    unsigned int getGeneration() const;
    void setGeneration(unsigned int value);
    float getPercentage() const;
    void setPercentage(float value);
    unsigned int getMutationId() const;
    void setMutationId(unsigned int value);
    QString getMutationName() const;
    void setMutationName(const QString &value);
    manualMutationAllelsChanged getMutationAllelsChanged() const;
    void setMutationAllelsChanged(const manualMutationAllelsChanged &value);
};

#endif // MUTATIONMANUAL_H
