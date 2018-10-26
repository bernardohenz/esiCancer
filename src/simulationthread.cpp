#include "simulationthread.h"
#include "cellsystem.h"
#include <QDebug>

SimulationThread::SimulationThread()
{

}

void SimulationThread::registerCellSystem(CellSystem *tsystem)
{
    myCellSystem = tsystem;
}

void SimulationThread::run()
{
    mTime.restart();
    for(int i=0;i<myCellSystem->getStopGenerations();i++)
        if(!myCellSystem->process())
            break;
    qDebug()<<"Finished: "<<mTime.elapsed()<<"ms";
    if (outputFilesMode==1){
        myCellSystem->exportInfo(exportingFileName);
    }
    else if (outputFilesMode==2)
        myCellSystem->exportLastLine(exportingAutorunFileName);
    else if (outputFilesMode==3){
        myCellSystem->exportInfo(exportingFileName);
        myCellSystem->exportLastLine(exportingAutorunFileName);
    }
    //myCellSystem->show();
}

void SimulationThread::setExportingFilename(QString tname)
{
    exportingFileName = tname;
}

void SimulationThread::setExportingAutorunFilename(QString tname)
{
    exportingAutorunFileName = tname;
}

void SimulationThread::setOutputFilesMode(unsigned int mode)
{
    outputFilesMode = mode;
}
