#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "dialogmanualmutation.h"
#include "mutationmanual.h"
#include <QFileDialog>
#include "cellsystem.h"
#include <QDebug>
#include <QThread>
#include <QFile>
#include <QString>
#include <QXmlStreamReader>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->LoadTableButton,SIGNAL(pressed()),this,SLOT(loadMutationTable()));
    connect(ui->startSimulationButton,SIGNAL(pressed()),this,SLOT(startSimulation()));
    connect(ui->forceStopButton,SIGNAL(pressed()),this,SLOT(forceStop()));
    connect(ui->runAutomaticsButton,SIGNAL(pressed()),this,SLOT(startAutomaticRuns()));

    connect(ui->addNewRuleButton,SIGNAL(pressed()),this,SLOT(addMutationManualToTable()));
    connect(ui->editRuleButton,SIGNAL(pressed()),this,SLOT(editMutationManualFromTable()));
    connect(ui->deleteRuleButton,SIGNAL(pressed()),this,SLOT(deleteMutationTableCurrentRule()));
    connect(ui->clearAllMutationManual,SIGNAL(pressed()),this,SLOT(clearAllMutationsRules()));


    mySystem = new CellSystem();
    mSimulationThread = new SimulationThread();
    mSimulationThread->registerCellSystem(mySystem);
    connect(mSimulationThread,SIGNAL(finished()),this,SLOT(nextAutomatic()));
    runningAutomaticsBool = false;
    startSimulationTimer = new QTimer(this);
    connect(startSimulationTimer, SIGNAL(timeout()), this, SLOT(setTextStartSimulationLabel()));
    loadConfigFile();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::loadConfigFile()
{
    QFile xmlFile("configs.xml");
    if(!(xmlFile.open(QIODevice::ReadOnly))){
        qDebug()<<"Error opening the Config file";
        return;
    }
    qDebug()<<"Config file loaded";
    QXmlStreamReader *xmlReader = new QXmlStreamReader(&xmlFile);

    int i=0;
    QStringRef currentState;
    int tmpInt;
    float tmpFloat;
    while(!xmlReader->atEnd() && !xmlReader->hasError()) {
        QXmlStreamReader::TokenType token = xmlReader->readNext();
        if(token == QXmlStreamReader::StartDocument) {
            continue;
        }
        if(xmlReader->name()=="esiCancer")
            continue;
        if (currentState.isEmpty()){
            currentState = xmlReader->name();
        } else{
            if(currentState == "default_values"){
                if (xmlReader->name() == "default_values"){
                    currentState.clear();
                    continue;
                } else if (xmlReader->name() == "number_of_cells"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"Number of cells: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->numberOfCellsSpin->setValue(tmpInt);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "prob_prol") {
                    xmlReader->readNext();
                    tmpFloat = xmlReader->text().toFloat();
                    //qDebug()<<"prob_prol: "<<tmpFloat;
                    if ((tmpFloat>0)&&(tmpFloat<1))
                        ui->probProlSpin->setValue(tmpFloat);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "prob_death") {
                    xmlReader->readNext();
                    tmpFloat = xmlReader->text().toFloat();
                    //qDebug()<<"prob_death: "<<tmpFloat;
                    if ((tmpFloat>0)&&(tmpFloat<1))
                        ui->probDeathSpin->setValue(tmpFloat);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "max_divisions"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"Max divs: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->maxDivisionsSpin->setValue(tmpInt);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "muts_per_division"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"Muts per div: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->mutationsPerDivisionSpin->setValue(tmpInt);
                    xmlReader->readNext();
                }
            } else if (currentState == "stop_conditions"){
                if (xmlReader->name() == "stop_conditions"){
                    currentState.clear();
                    continue;
                } else if (xmlReader->name() == "num_of_generations"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"Stop num gen: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->numberOfGenerationsStopSpin->setValue(tmpInt);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "stop_num_of_cells"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"Stop num cells: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->numberOfCellsStopSpin->setValue(tmpInt);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "stop_num_of_cells"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"stop mut cells: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        ui->mutatedCellsStopSpin->setValue(tmpInt);
                    xmlReader->readNext();
                }
            } else if (currentState == "max_values"){
                if (xmlReader->name() == "max_values"){
                    currentState.clear();
                    continue;
                } else if (xmlReader->name() == "muts_per_division"){
                    xmlReader->readNext();
                    tmpInt = xmlReader->text().toInt();
                    //qDebug()<<"MaxMuts per division: "<<tmpInt;
                    if ((tmpInt>0)&&(tmpInt<9999999))
                        mySystem->setMaxMutationsPerDivision(tmpInt);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "prob_prol") {
                    xmlReader->readNext();
                    tmpFloat = xmlReader->text().toFloat();
                    //qDebug()<<"max prob_prol: "<<tmpFloat;
                    if ((tmpFloat>0)&&(tmpFloat<1))
                        mySystem->setMaxProliferation(tmpFloat);
                    xmlReader->readNext();
                } else if (xmlReader->name() == "prob_death") {
                    xmlReader->readNext();
                    tmpFloat = xmlReader->text().toFloat();
                    //qDebug()<<"max prob_death: "<<tmpFloat;
                    if ((tmpFloat>0)&&(tmpFloat<1))
                        mySystem->setMaxDeath(tmpFloat);
                    xmlReader->readNext();
                }
            }
        }

        //qDebug()<<"Current state: "<<currentState;
        //qDebug()<<"iteration: "<<i<<"  xml_name: "<<xmlReader->name();
        i++;
    }

    if(xmlReader->hasError()) {
        qDebug()<<"Error parsing the XML file";
        return;
    }
}

void MainWindow::loadMutationTable()
{
    QString filename = QFileDialog::getOpenFileName(this,"Open mutation table ...");
    if(!filename.isEmpty()){
        mySystem->loadMutationTableFile(filename);
        /*QStringList tmpList = mySystem->getMutationListNames();
        tmpList.append("ALL");
        ui->comboMutationList1->clear();
        ui->comboMutationList1->addItems(tmpList);
        ui->comboMutationList1->setCurrentIndex(tmpList.size()-1);
        ui->comboMutationList2->clear();
        ui->comboMutationList2->addItems(tmpList);
        ui->comboMutationList2->setCurrentIndex(tmpList.size()-1);
        ui->comboMutationList3->clear();
        ui->comboMutationList3->addItems(tmpList);
        ui->comboMutationList3->setCurrentIndex(tmpList.size()-1);
        ui->frame_MutationAddedManual->setEnabled(true);
        ui->tableMutationsManual->setEnabled(true);*/
        qDebug()<<"Loaded Mutation Table";
        ui->frame_ManualMutations->setEnabled(true);
    }
}

void MainWindow::startSimulation()
{
    runningAutomaticsBool = false;
    mySystem->loadDefaultValues(ui->probProlSpin->value(),ui->probDeathSpin->value(),ui->maxDivisionsSpin->value(),ui->mutationsPerDivisionSpin->value());
    mySystem->loadSimulationParameters(ui->seedSpin->value(),ui->numberOfCellsSpin->value());
    mySystem->loadStopConditions(ui->numberOfGenerationsStopSpin->value(),ui->numberOfCellsStopSpin->value(),ui->mutatedCellsStopSpin->value());
    mySystem->startSimulation();
    mSimulationThread->setExportingFilename("output");
    mSimulationThread->setOutputFilesMode(1);
    mSimulationThread->start();
    ui->start_simulation_label_finished->setText("Exported files to output_.csv");

    startSimulationTimer->start(1000);
    //for(int i=0;i<500;i++)
    //    mySystem->process();
}

void MainWindow::forceStop()
{
    mSimulationThread->terminate();
}

void MainWindow::setTextStartSimulationLabel()
{
    ui->start_simulation_label_finished->setText("");
}

void MainWindow::addMutationManualToTable()
{
    DialogManualMutation tmp(this,mySystem->getMutationListNames());
    if(tmp.exec()==QDialog::Accepted) {
        qDebug()<<"ALOHA";
        ui->tableMutationsManual->setRowCount(ui->tableMutationsManual->rowCount()+1);
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->rowCount()-1,0,new QTableWidgetItem(QString::number(tmp.getGeneration())));
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->rowCount()-1,1,new QTableWidgetItem(QString::number(tmp.getPercentage())));
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->rowCount()-1,2,new QTableWidgetItem(tmp.getMutationName()));
        if (tmp.getAllelsChanged()==MONOALLELIC){
            ui->tableMutationsManual->setItem(ui->tableMutationsManual->rowCount()-1,3,new QTableWidgetItem("Monoallelic"));
        } else{
            ui->tableMutationsManual->setItem(ui->tableMutationsManual->rowCount()-1,3,new QTableWidgetItem("Biallelic"));
        }

        ui->tableMutationsManual->resizeColumnsToContents();
        mySystem->addMutationRule(new MutationManual(tmp.getGeneration(),tmp.getPercentage(),tmp.getIndex(),tmp.getMutationName(),tmp.getAllelsChanged()));
    }
}

void MainWindow::editMutationManualFromTable()
{
    if(ui->tableMutationsManual->currentRow()==-1)
        return;
    DialogManualMutation tmp(this,mySystem->getMutationListNames());
    MutationManual *tmpmutation = mySystem->getMutationManualFromIndex(ui->tableMutationsManual->currentRow());
    tmp.setValues(tmpmutation->getGeneration(),tmpmutation->getPercentage(),tmpmutation->getMutationId(),tmpmutation->getMutationAllelsChanged());
    if(tmp.exec()==QDialog::Accepted) {
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->currentRow(),0,new QTableWidgetItem(QString::number(tmp.getGeneration())));
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->currentRow(),1,new QTableWidgetItem(QString::number(tmp.getPercentage())));
        ui->tableMutationsManual->setItem(ui->tableMutationsManual->currentRow(),2,new QTableWidgetItem(tmp.getMutationName()));
        if (tmp.getAllelsChanged()==MONOALLELIC){
            ui->tableMutationsManual->setItem(ui->tableMutationsManual->currentRow(),3,new QTableWidgetItem("Monoallelic"));
        } else{
            ui->tableMutationsManual->setItem(ui->tableMutationsManual->currentRow(),3,new QTableWidgetItem("Biallelic"));
        }
        //ui->tableMutationsManual->resizeColumnsToContents();
        mySystem->updateMutationRule(ui->tableMutationsManual->currentRow(),tmp.getGeneration(),tmp.getPercentage(),tmp.getMutationName(),tmp.getIndex(),tmp.getAllelsChanged());
    }
}

void MainWindow::deleteMutationTableCurrentRule()
{
    if(ui->tableMutationsManual->currentRow()==-1)
        return;
    mySystem->removeMutationRule(ui->tableMutationsManual->currentRow());
    ui->tableMutationsManual->removeRow(ui->tableMutationsManual->currentRow());
}

void MainWindow::clearAllMutationsRules()
{
    ui->tableMutationsManual->clearContents();
    ui->tableMutationsManual->setRowCount(0);
    mySystem->removeAllRules();
}

void MainWindow::startAutomaticRuns()
{
    automaticsFilename = QFileDialog::getSaveFileName(this,"Exporting path");
    if(automaticsFilename=="")
        return;
    parameterToIterate = ui->minValueAutomaticSpin->value();
    runningAutomaticsBool = true;
    nextAutomatic();
}

void MainWindow::nextAutomatic()
{
    //working only for seed
    if(!runningAutomaticsBool)
        return;
    if(parameterToIterate>(unsigned int)ui->maxValueAutomaticSpin->value()){
        runningAutomaticsBool = false;
        ui->automaticsProgressBar->setValue(100);
        return;
    }
    mySystem->loadDefaultValues(ui->probProlSpin->value(),ui->probDeathSpin->value(),ui->maxDivisionsSpin->value(),ui->mutationsPerDivisionSpin->value());
    mySystem->loadSimulationParameters(parameterToIterate,ui->numberOfCellsSpin->value());
    mySystem->loadStopConditions(ui->numberOfGenerationsStopSpin->value(),ui->numberOfCellsStopSpin->value(),ui->mutatedCellsStopSpin->value());
    mySystem->startSimulation();
    mSimulationThread->setExportingFilename(automaticsFilename+"_seed"+QString::number(parameterToIterate));
    mSimulationThread->setExportingAutorunFilename(automaticsFilename+"_automatics");
    mSimulationThread->setOutputFilesMode(ui->comboBoxOutputFiles->currentIndex()+1);
    mSimulationThread->start();


    int tmpProgress = 100*(parameterToIterate - ui->minValueAutomaticSpin->value())/(1+ui->maxValueAutomaticSpin->value() - ui->minValueAutomaticSpin->value());
    ui->automaticsProgressBar->setValue(tmpProgress);

    parameterToIterate+=ui->incrementAutomaticSpin_2->value();
}
