#include <iostream>
#include <string>
#include <vector>
#include <torch/torch.h>

#include "include/SlaterBasisData.h"
#include "include/readQuadWeightsFile.h"
#include "include/AtomInfo.h"
#include "include/SlaterPrimitive.h"
#include "SlaterBasisSet.h"

#include "AuxDensityMatrixSlater.h"

// const std::string DataDirPath = "../Data/";

int main() {
    std::cout << "Hello, World!" << std::endl;

    // ------------------ Read QuadWeights --------------
    //

    std::string quadweightsFile = DataDirPath + "QWt";
    auto [quadpts, quadWt] = readQuadWeightsFile(quadweightsFile);

    std::cout << quadpts.size() << std::endl;
    int nQuad = quadpts.size()/3;

    std::cout << nQuad << std::endl;

    // ------------------ Read AtomicCoords_Slater --------------

    std::string coordFile = DataDirPath + "AtomicCoords_Slater";

    /*
    std::vector<Atom> atoms = Atom::readCoordFile(coordFile);

    // ------------------ SlaterBasisSets -------------
    SlaterBasisSet sbs;
    sbs.constructBasisSet(atoms);
    int nBasis = sbs.getTotalBasisSize(atoms);
    std::cout << "nBasis : " << nBasis << std::endl;

    // ----------------- SlaterBasisData----------------

    SlaterBasisData sbd;
    int maxDerOrder = 2;
    sbd.evalBasisData(atoms, quadpts, sbs, nQuad, nBasis, maxDerOrder);
    sbd.printBasisInfo(nQuad, nBasis);*/

    AuxDensityMatrixSlater auxDMSlater(quadpts, quadWt, coordFile, nQuad, 2 );

    return 0;
}
