#ifndef DIALOGMANUALMUTATION_H
#define DIALOGMANUALMUTATION_H

#include <QDialog>
#include "mutationmanual.h"

namespace Ui {
class DialogManualMutation;
}

class DialogManualMutation : public QDialog
{
    Q_OBJECT
    

public:
    explicit DialogManualMutation(QWidget *parent, QStringList tList);

    void setValues(int generation, int percentage, int selectedIndex, manualMutationAllelsChanged tAllelsChanged);
    int getGeneration();
    int getPercentage();
    QString getMutationName();
    int getIndex();
    manualMutationAllelsChanged getAllelsChanged();

    ~DialogManualMutation();
    
private:
    Ui::DialogManualMutation *ui;

};

#endif // DIALOGMANUALMUTATION_H
