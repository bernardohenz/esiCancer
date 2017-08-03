#include "mutationmanual.h"

unsigned int MutationManual::getGeneration() const
{
    return generation;
}

void MutationManual::setGeneration(unsigned int value)
{
    generation = value;
}

float MutationManual::getPercentage() const
{
    return percentage;
}

void MutationManual::setPercentage(float value)
{
    percentage = value;
}

unsigned int MutationManual::getMutationId() const
{
    return mutationId;
}

void MutationManual::setMutationId(unsigned int value)
{
    mutationId = value;
}

QString MutationManual::getMutationName() const
{
    return mutationName;
}

void MutationManual::setMutationName(const QString &value)
{
    mutationName = value;
}

manualMutationAllelsChanged MutationManual::getMutationAllelsChanged() const
{
    return mutationAllelsChanged;
}

void MutationManual::setMutationAllelsChanged(const manualMutationAllelsChanged &value)
{
    mutationAllelsChanged = value;
}

MutationManual::MutationManual(unsigned int tgeneration, float tpercentage, unsigned int tmutationId, QString tmutationName,manualMutationAllelsChanged tmutationAllelsChanged)
{
    generation = tgeneration;
    percentage = tpercentage;
    mutationId = tmutationId;
    mutationName = tmutationName;
    mutationAllelsChanged = tmutationAllelsChanged;
}
