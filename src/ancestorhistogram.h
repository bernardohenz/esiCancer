#ifndef ANCESTORHISTOGRAM_H
#define ANCESTORHISTOGRAM_H

#include <vector>

class AncestorHistogram{
private:
    std::vector<unsigned int> vAncestors;
public:
    AncestorHistogram(unsigned int tnumberOfAncestors){
        vAncestors = std::vector<unsigned int>(tnumberOfAncestors,0);
    }

    AncestorHistogram(unsigned int tnumberOfAncestors,unsigned int counter){
        vAncestors = std::vector<unsigned int>(tnumberOfAncestors,counter);
    }

    void incrementAncestorCounter(unsigned int tanc) {
        if(tanc<vAncestors.size())
            vAncestors[tanc]++;
    }

    unsigned int getNumberofGivenAncestor(unsigned int tanc){
        return vAncestors[tanc];
    }


};

#endif // ANCESTORHISTOGRAM_H
