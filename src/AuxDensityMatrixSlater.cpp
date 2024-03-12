//
// Created by Arghadwip Paul.
//

#include "AuxDensityMatrixSlater.h"

AuxDensityMatrixSlater::AuxDensityMatrixSlater(
        const std::vector<double>& quadpts,
        const std::vector<double>& quadWt,
        const std::string coordFile,
        const int nQuad,
        const int nSpin
        ):
        AuxDensityMatrix(),  // call the base class constructor
        quadpts(quadpts),
        quadWt(quadWt),
        coordFile(coordFile),
        nQuad(nQuad),
        nSpin(nSpin)
{
    // ------------------ Read AtomicCoords_Slater --------------
    std::vector<Atom> atoms = Atom::readCoordFile(coordFile);

    // ------------------ SlaterBasisSets -------------
    sbs.constructBasisSet(atoms);
    nBasis = sbs.getTotalBasisSize(atoms);
    std::cout << "nBasis : " << nBasis << std::endl;

    DM.assign(nSpin * nBasis * nBasis, 0.0);
    this -> initializeLocalDescriptors();

    int maxDerOrder = 2;
    sbd.evalBasisData(atoms, quadpts, sbs, nQuad, nBasis, maxDerOrder);
    sbd.printBasisInfo(nQuad, nBasis);
}

AuxDensityMatrixSlater::~AuxDensityMatrixSlater() {
    // Destructor implementation
}

void AuxDensityMatrixSlater::setDMzero() {
    // setDMzero implementation
    for (auto& element : DM) {
        element = 0.0;
    }
}

void
AuxDensityMatrixSlater::initializeLocalDescriptors(){
    attributeData[DensityDescriptorDataAttributes::valuesTotal] = std::vector<double>(nQuad, 0.0);
    attributeData[DensityDescriptorDataAttributes::valuesSpinUp] = std::vector<double>(nQuad, 0.0);
    attributeData[DensityDescriptorDataAttributes::valuesSpinDown] = std::vector<double>(nQuad, 0.0);
    attributeData[DensityDescriptorDataAttributes::gradValueSpinUp] = std::vector<double>(nQuad * 3, 0.0);
    attributeData[DensityDescriptorDataAttributes::gradValueSpinDown] = std::vector<double>(nQuad * 3, 0.0);
    attributeData[DensityDescriptorDataAttributes::sigma] = std::vector<double>(nQuad * 3, 0.0);
    attributeData[DensityDescriptorDataAttributes::hessianSpinUp] = std::vector<double>(nQuad * 9, 0.0);
    attributeData[DensityDescriptorDataAttributes::hessianSpinDown] = std::vector<double>(nQuad * 9, 0.0);
    attributeData[DensityDescriptorDataAttributes::laplacianSpinUp] = std::vector<double>(nQuad, 0.0);
    attributeData[DensityDescriptorDataAttributes::laplacianSpinDown] = std::vector<double>(nQuad, 0.0);
}

void
AuxDensityMatrixSlater::setLocalDescriptors(){

}

void
AuxDensityMatrixSlater::setLocalDescriptors(
        DensityDescriptorDataAttributes attr,
        std::pair<int, int> indexRange,
        std::vector<double> values) {

    if (indexRange.first >= indexRange.second) {
        throw std::invalid_argument("Invalid index range: The start index must be less than the end index.");
    }

    auto it = attributeData.find(attr);
    if (it == attributeData.end()) {
        throw std::invalid_argument("Specified attribute does not exist.");
    }

    size_t expectedRangeSize = indexRange.second - indexRange.first;
    if (values.size() != expectedRangeSize) {
        throw std::invalid_argument("Values vector size does not match the specified index range.");
    }

    if (indexRange.second > it->second.size()) {
        throw std::out_of_range("Index range is out of bounds for the specified attribute.");
    }

    std::copy(values.begin(), values.end(), it->second.begin() + indexRange.first);

}

std::vector<double> AuxDensityMatrixSlater::getLocalDescriptors(
        DensityDescriptorDataAttributes attr,
        std::pair<int, int> indexRange) const {

    if (indexRange.first >= indexRange.second) {
        throw std::invalid_argument("Invalid index range: The start index must be less than the end index.");
    }

    auto it = attributeData.find(attr);
    if (it == attributeData.end()) {
        throw std::invalid_argument("Specified attribute does not exist.");
    }

    if (indexRange.second > it->second.size()) {
        throw std::out_of_range("Index range is out of bounds for the specified attribute.");
    }

    // Extract and return the values vector for the specified range
    std::vector<double> values(it->second.begin() + indexRange.first, it->second.begin() + indexRange.second);

    return values;
}



void
AuxDensityMatrixSlater::evalLocalDescriptors() {
    // evalLocalDescriptors implementation
    std::cout << "eval local descriptors" << std::endl;

    int DMSpinOffset = nBasis * nBasis;
    std::pair<int, int> indexRange;
    for(int iQuad = 0; iQuad < nQuad; iQuad++) {

        std::vector<double> rhoUp(1, 0.0);
        std::vector<double> rhoDown(1, 0.0);
        std::vector<double> rhoTotal(1, 0.0);
        std::vector<double> gradrhoUp(3, 0.0);
        std::vector<double> gradrhoDown(3, 0.0);
        std::vector<double> sigma(3, 0.0);
        std::vector<double> HessianrhoUp(9, 0.0);
        std::vector<double> HessianrhoDown(9, 0.0);
        std::vector<double> LaplacianrhoUp(1, 0.0);
        std::vector<double> LaplacianrhoDown(1, 0.0);

        for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
                for(int iSpin = 0; iSpin < nSpin; iSpin++) {
                    if (iSpin == 0) {
                        rhoUp[0] += DM[i * nBasis + j] * sbd.getBasisValues(iQuad * nBasis + i) *
                                 sbd.getBasisValues(iQuad * nBasis + j);

                        for(int derIndex = 0; derIndex < 3; derIndex++) {
                            gradrhoUp[derIndex] += DM[i * nBasis + j] *
                                                    (sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex) *
                                                     sbd.getBasisValues(iQuad * nBasis + j) +
                                                     sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex) *
                                                     sbd.getBasisValues(iQuad * nBasis + i)
                                                     );
                        }

                        for(int derIndex1 = 0; derIndex1 < 3; derIndex1++) {
                            for(int derIndex2 = 0; derIndex2 < 3; derIndex2++){
                                HessianrhoUp[derIndex1 * 3 + derIndex2] += DM[i * nBasis + j] *
                                        (
                                        sbd.getBasisHessianValues(iQuad * nBasis + 9 * i + 3 * derIndex1 + derIndex2) *
                                        sbd.getBasisValues(iQuad * nBasis + j) +
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex1) *
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex2) +
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex2) *
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex1) +
                                        sbd.getBasisValues(iQuad * nBasis + i) *
                                        sbd.getBasisHessianValues(iQuad * nBasis + 9 * j + 3 * derIndex1 + derIndex2)
                                        );
                            }
                        }
                    }

                    if (iSpin == 1) {
                        rhoDown[0] += DM[DMSpinOffset + i * nBasis + j] * sbd.getBasisValues(iQuad * nBasis + i) *
                                   sbd.getBasisValues(iQuad * nBasis + j);

                        for(int derIndex = 0; derIndex < 3; derIndex++) {
                            gradrhoDown[derIndex] += DM[DMSpinOffset + i * nBasis + j] *
                                                   (sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex) *
                                                    sbd.getBasisValues(iQuad * nBasis + j) +
                                                    sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex) *
                                                    sbd.getBasisValues(iQuad * nBasis + i)
                                                    );
                        }

                        for(int derIndex1 = 0; derIndex1 < 3; derIndex1++) {
                            for(int derIndex2 = 0; derIndex2 < 3; derIndex2++){
                                HessianrhoDown[derIndex1 * 3 + derIndex2] += DM[DMSpinOffset + i * nBasis + j] *
                                        (
                                        sbd.getBasisHessianValues(iQuad * nBasis + 9 * i + 3 * derIndex1 + derIndex2) *
                                        sbd.getBasisValues(iQuad * nBasis + j) +
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex1) *
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex2) +
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * i + derIndex2) *
                                        sbd.getBasisGradValues(iQuad * nBasis + 3 * j + derIndex1) +
                                        sbd.getBasisValues(iQuad * nBasis + i) *
                                        sbd.getBasisHessianValues(iQuad * nBasis + 9 * j + 3 * derIndex1 + derIndex2)
                                        );
                            }
                        }
                    }
                }
            }
        }

        rhoTotal[0] = rhoUp[0] + rhoDown[0];
        sigma[0] = gradrhoUp[0] * gradrhoUp[0] + gradrhoUp[1] * gradrhoUp[1] + gradrhoUp[2] * gradrhoUp[2];
        sigma[1] = gradrhoUp[0] * gradrhoDown[0] + gradrhoUp[1] * gradrhoDown[1] + gradrhoUp[2] * gradrhoDown[2];
        sigma[2] = gradrhoDown[0] * gradrhoDown[0] + gradrhoDown[1] * gradrhoDown[1] + gradrhoDown[2] * gradrhoDown[2];
        LaplacianrhoUp[0] = HessianrhoUp[0] + HessianrhoUp[4] + HessianrhoUp[8];
        LaplacianrhoDown[0] = HessianrhoDown[0] + HessianrhoDown[4] + HessianrhoDown[8];

        indexRange = std::make_pair(iQuad, iQuad + 1);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::valuesSpinUp, indexRange, rhoUp);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::valuesSpinDown, indexRange, rhoDown);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::valuesTotal, indexRange, rhoTotal);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::laplacianSpinUp, indexRange, LaplacianrhoUp);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::laplacianSpinDown, indexRange, LaplacianrhoDown);

        indexRange = std::make_pair(iQuad * 3, iQuad * 3 + 3);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::gradValueSpinUp, indexRange, gradrhoUp);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::gradValueSpinDown, indexRange, gradrhoDown);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::sigma, indexRange, sigma);

        indexRange = std::make_pair(iQuad * 9, iQuad * 9 + 9);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::hessianSpinUp, indexRange, HessianrhoUp);
        this -> setLocalDescriptors(DensityDescriptorDataAttributes::hessianSpinDown, indexRange, HessianrhoDown);
    }
}

void
AuxDensityMatrixSlater::projectDensityMatrix() {

}

void
AuxDensityMatrixSlater::projectDensityMatrix(
        const std::vector<double>& Qpts,
        const std::vector<double>& QWt,
        const int nQ,
        const std::vector<double>& psiFunc,
        const std::vector<double>& fValues,
        const std::pair<int, int> nPsi) {


}


void AuxDensityMatrixSlater::projectDensity() {
    // projectDensity implementation
    std::cout << "Error : No implementation yet" << std::endl;
}
