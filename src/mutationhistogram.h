#ifndef MUTATIONHISTOGRAM_H
#define MUTATIONHISTOGRAM_H
#include <vector>


class MutationHistogram
{
private:
    std::vector<unsigned int> vBins;

public:
    MutationHistogram(){
        vBins = std::vector<unsigned int>();
    }
    MutationHistogram(unsigned int tnumberOfBins){
        vBins = std::vector<unsigned int>(tnumberOfBins,0);
    }

    void initialize(unsigned int tnumberOfBins){
        vBins.clear();
        vBins = std::vector<unsigned int>(tnumberOfBins,0);
    }

    void resetBins(){
        for (unsigned int i=0;i<vBins.size(); i++)
            vBins[i] = 0;
    }

    void incrementBin(unsigned int tbin) {
        vBins[tbin]++;
    }

    unsigned int getNumberOfSamples(unsigned int tbin){
        return vBins[tbin];
    }

};

#endif // MUTATIONHISTOGRAM_H
