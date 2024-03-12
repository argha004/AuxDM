//
// Created by Arghadwip Paul.
//

#ifndef AUXDM_AUXDENSITYMATRIX_H
#define AUXDM_AUXDENSITYMATRIX_H

enum class DensityDescriptorDataAttributes
{
    valuesTotal,
    valuesSpinUp,
    valuesSpinDown,
    gradValueSpinUp,
    gradValueSpinDown,
    sigma, // contracted gradients of the density using libxc format
    hessianSpinUp,
    hessianSpinDown,
    laplacianSpinUp,
    laplacianSpinDown,
    rhoStencilData,
    FDSpacing
};

class AuxDensityMatrix {
public:
    // Constructor
    AuxDensityMatrix();

    // Virtual destructor to ensure proper cleanup of derived classes
    virtual ~AuxDensityMatrix();

    // Pure virtual functions
    virtual void setDMzero() = 0;
    virtual void initializeLocalDescriptors() = 0;
    virtual void setLocalDescriptors() = 0;
    virtual void getLocalDescriptors() = 0;
    virtual void evalLocalDescriptors() = 0;
    virtual void projectDensityMatrix() = 0;
    virtual void projectDensity() = 0;
};

#endif //AUXDM_AUXDENSITYMATRIX_H
