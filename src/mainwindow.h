#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "simulationthread.h"

class CellSystem;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    CellSystem *mySystem;
    SimulationThread *mSimulationThread;

    //Automatics
    unsigned int parameterToIterate;
    bool runningAutomaticsBool;
    QString automaticsFilename;

public slots:
    void loadConfigFile();

    void loadMutationTable();
    void loadSinergyTable();
    void startSimulation();
    void forceStop();

    void changeMicroAmbient(bool state);
    //Manual Mutations
    void addMutationManualToTable();
    void editMutationManualFromTable();
    void deleteMutationTableCurrentRule();
    void clearAllMutationsRules();

    //Automatic Runs
    void startAutomaticRuns();
    void nextAutomatic();

};

#endif // MAINWINDOW_H
