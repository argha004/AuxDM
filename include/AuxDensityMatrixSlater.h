//
// Created by Arghadwip Paul.
//

#ifndef AUXDM_AUXDENSITYMATRIXSLATER_H
#define AUXDM_AUXDENSITYMATRIXSLATER_H

#include "AuxDensityMatrix.h"
#include "SlaterBasisSet.h"
#include "SlaterBasisData.h"
#include "AtomInfo.h"
#include <vector>
#include <utility>
#include <map>

class AuxDensityMatrixSlater : public AuxDensityMatrix {

private:
    std::vector<double> quadpts;
    std::vector<double> quadWt;
    std::string coordFile;
    int nQuad;
    int nSpin;

    std::vector<Atom> atoms;
    SlaterBasisSet sbs;
    SlaterBasisData sbd;
    int nBasis;

    std::vector<double> DM;
    std::unordered_map<DensityDescriptorDataAttributes, std::vector<double>> attributeData;


public:
    // Constructor
    AuxDensityMatrixSlater(const std::vector<double>& quadpts,
                           const std::vector<double>& quadWt,
                           const std::string coordFile,
                           const int nQuad,
                           const int nSpin
                           );

    // Destructor
    virtual ~AuxDensityMatrixSlater();

    // Implementations of pure virtual functions from the base class
    virtual void setDMzero() override;
    virtual void initializeLocalDescriptors() override;
    virtual void setLocalDescriptors() override;
    virtual void getLocalDescriptors() override;

    virtual void evalLocalDescriptors() override;

    virtual void projectDensityMatrix() override;

    void setLocalDescriptors(DensityDescriptorDataAttributes attr,
                             std::pair<int, int> indexRange,
                             std::vector<double> values);

    std::vector<double> getLocalDescriptors(DensityDescriptorDataAttributes attr,
                                            std::pair<int, int> indexRange) const;


    void projectDensityMatrix(const std::vector<double>& Qpts,
                              const std::vector<double>& QWt,
                              const int nQ,
                              const std::vector<double>& psiFunc,
                              const std::vector<double>& fValues,
                              const std::pair<int, int> nPsi);

    virtual void projectDensity() override;
};

#endif //AUXDM_AUXDENSITYMATRIXSLATER_H
