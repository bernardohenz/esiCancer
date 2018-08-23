#ifndef MUTATIONCOUNTER_H
#define MUTATIONCOUNTER_H
#include <vector>
#include <QDebug>

class MutationCounter{
private:
    std::vector<unsigned int> vMutations;
public:
    MutationCounter(){
        vMutations = std::vector<unsigned int>();
    }
    MutationCounter(unsigned int tnumberOfMuts){
        vMutations = std::vector<unsigned int>(tnumberOfMuts,0);
    }
    void reset(){
        for (unsigned int i=0;i<vMutations.size(); i++)
            vMutations[i] = 0;
    }
    void incrementMutCounter(unsigned int tmut) {
        if (tmut>=vMutations.size()){
            return;
        }
        vMutations[tmut]++;
    }

    unsigned int getNumberGivenMut(unsigned int tmut){
        return vMutations[tmut];
    }
};

#endif // MUTATIONCOUNTER_H
