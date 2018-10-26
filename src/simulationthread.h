#ifndef SIMULATIONTHREAD_H
#define SIMULATIONTHREAD_H


#include <QtCore>
#include <QTime>

class CellSystem;

class SimulationThread:public QThread
{
private:
    QString exportingFileName;
    QString exportingAutorunFileName;
    unsigned int outputFilesMode;
    QTime mTime;
public:
    SimulationThread();
    CellSystem *myCellSystem;
    void registerCellSystem(CellSystem *tsystem);
    void run();
    void setExportingFilename(QString tname);
    void setExportingAutorunFilename(QString tname);
    void setOutputFilesMode(unsigned int mode);
};

#endif // SIMULATIONTHREAD_H
