#include "mutationmanual.h"
#include <QStringList>
#include <QDebug>

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

unsigned int MutationManual::getMutationEventId() const
{
    return mutationEventId;
}

void MutationManual::setMutationEventId(int tvalue)
{
    mutationEventId = tvalue;
}

void MutationManual::setMutationFullName(QString tmutationFullName)
{
    mutationGeneName = tmutationFullName.split('-')[0];
    mutationEventName = tmutationFullName.split('-')[1];
}



QString MutationManual::getMutationEventName() const
{
    return mutationEventName;
}


QString MutationManual::getMutationGeneName() const
{
    return mutationGeneName;
}


manualMutationAllelsChanged MutationManual::getMutationAllelsChanged() const
{
    return mutationAllelsChanged;
}

void MutationManual::setMutationAllelsChanged(const manualMutationAllelsChanged &value)
{
    mutationAllelsChanged = value;
}

MutationManual::MutationManual(unsigned int tgeneration, float tpercentage, unsigned int tmutationEventId, QString tmutationFullName, manualMutationAllelsChanged tmutationAllelsChanged)
{
    generation = tgeneration;
    percentage = tpercentage;
    mutationEventId = tmutationEventId;
    mutationGeneName = tmutationFullName.split('-')[0];
    mutationEventName = tmutationFullName.split('-')[1];
    mutationFullName = tmutationFullName;
    mutationAllelsChanged = tmutationAllelsChanged;
    qDebug()<<"Created: "<<mutationGeneName<<"     event: "<<mutationEventName;
}
