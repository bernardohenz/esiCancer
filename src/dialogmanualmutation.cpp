#include "dialogmanualmutation.h"
#include "ui_dialogmanualmutation.h"
#include <QStringList>

DialogManualMutation::DialogManualMutation(QWidget *parent, QStringList tList) :
    QDialog(parent),
    ui(new Ui::DialogManualMutation)
{
    ui->setupUi(this);
    connect(ui->buttonAccept,SIGNAL(pressed()),this,SLOT(accept()));
    connect(ui->buttonReject,SIGNAL(pressed()),this,SLOT(reject()));
    QStringList mutationNames = tList;
    for (int i=0; i<mutationNames.size(); i++){
        ui->comboBoxMutationFullName->addItem(mutationNames.at(i));
    }

    ui->buttonAccept->setFocus();
}

void DialogManualMutation::setValues(int generation, int percentage, int selectedIndexEvent ,manualMutationAllelsChanged tAllelsChanged)
{
    ui->spinMutationGeneration->setValue(generation);
    ui->spinMutationPercentage->setValue(percentage);
    ui->comboBoxMutationFullName->setCurrentIndex(selectedIndexEvent);
    if(tAllelsChanged==MONOALLELIC)
        ui->comboBoxAllelsChanged->setCurrentIndex(0);
    else
        ui->comboBoxAllelsChanged->setCurrentIndex(1);
}

int DialogManualMutation::getGeneration()
{
    return ui->spinMutationGeneration->value();
}

int DialogManualMutation::getPercentage()
{
    return ui->spinMutationPercentage->value();
}

QString DialogManualMutation::getMutationFullName()
{
    return ui->comboBoxMutationFullName->currentText();
}


int DialogManualMutation::getIndex()
{
    return ui->comboBoxMutationFullName->currentIndex();
}

manualMutationAllelsChanged DialogManualMutation::getAllelsChanged()
{
    if (ui->comboBoxAllelsChanged->currentIndex()==0)
        return MONOALLELIC;
    else
        return BIALLELIC;
}

DialogManualMutation::~DialogManualMutation()
{
    delete ui;
}
