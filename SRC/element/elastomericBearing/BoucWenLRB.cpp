BoucWenIM/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Yuhang Lu (ll05734823@gmail.com)
// Created: 2024/04/11
// Description: This file contains the implementation of the
// BoucWenLRB class.

#include "BoucWenLRB.h"

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>
#include <elementAPI.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <G3Globals.h>
#include <Message.h>
// using namespace std;
#include <iostream>

constexpr long double PI = 3.14159l;


// initialize the class wide variables
Matrix BoucWenLRB::theMatrix(12, 12);
Vector BoucWenLRB::theVector(12);


static int tag = 0;  // Tag to identify if bearing has failed in buckling
void* OPS_BoucWenLRB()
{
    int ndm = OPS_GetNDM();
    if (ndm != 3) {
        opserr << "BoucWenLRB is only available for ndm=3 three-dimensional models\n";
        return 0;
    }

    int ndf = OPS_GetNDF();
    if (ndf != 6) {
        opserr << "WARNING invalid ndf: " << ndf;
        opserr << ", for space problem need 6 - BoucWenLRB\n";
        return 0;
    }

    Element* theEle = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs == 0) { // parallel processing
        theEle = new BoucWenLRB();
        return theEle;
    }

    if (numArgs < 18) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: BoucWenLRB eleTag Nd1 Nd2 kInit fy alphaL dL ts tr n -P matTag -T matTag -My matTag -Mz matTag <-alpha1 alpha1> <-beta1 beta1> <-gamma1 gamma1> <-eta eta> <-alpha2 alpha2> <-beta2 beta2> <-gamma2 gamma2> <-c1 c1> <-tagLH> <-E2 E2> <-rhoL rhoL> <-cL cL> <-kS kS> <-alphaS alphaS> <-TL0 TL0> <-tagIH> <-c2 c2> <-c3 c3> <-orient x1 x2 x3 y1 y2 y3> <-sDratio sDratio> <-mass mass> <-maxIter maxIter> <-tol tol>\n";
        return 0;
    }

    // get the id, end nodes, and element properties
    int iData[3];
    double dData[7];
    int numData;

    numData = 3;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING: invalid element data\n";
        return 0;
    }

    int eleTag = iData[0];

    numData = 7;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING: error reading element properties for element" << eleTag << endln;
        return 0;
    }

    // materials
    UniaxialMaterial* mats[4] = { 0,0,0,0 };
    int matTag;

    const char* type = OPS_GetString();
    if (strcmp(type, "-P") != 0) {
        opserr << "WARNING: want -P\n";
        return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) < 0) {
        opserr << "WARNING: invalid -P matTag\n";
        return 0;
    }
    mats[0] = OPS_getUniaxialMaterial(matTag);
    if (mats[0] == 0) {
        opserr << "WARNING: -P material not found\n";
        return 0;
    }

    type = OPS_GetString();
    if (strcmp(type, "-T") != 0) {
        opserr << "WARNING: want -T\n";
        return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) < 0) {
        opserr << "WARNING: invalid -T matTag\n";
        return 0;
    }
    mats[1] = OPS_getUniaxialMaterial(matTag);
    if (mats[1] == 0) {
        opserr << "WARNING: -T material not found\n";
        return 0;
    }

    type = OPS_GetString();
    if (strcmp(type, "-My") != 0) {
        opserr << "WARNING: want -My\n";
        return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) < 0) {
        opserr << "WARNING: invalid -My matTag\n";
        return 0;
    }
    mats[2] = OPS_getUniaxialMaterial(matTag);
    if (mats[2] == 0) {
        opserr << "WARNING: -My material not found\n";
        return 0;
    }

    type = OPS_GetString();
    if (strcmp(type, "-Mz") != 0) {
        opserr << "WARNING: want -Mz\n";
        return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) < 0) {
        opserr << "WARNING: invalid -Mz matTag\n";
        return 0;
    }
    mats[3] = OPS_getUniaxialMaterial(matTag);
    if (mats[3] == 0) {
        opserr << "WARNING: -Mz material not found\n";
        return 0;
    }

    // The default values of the parameters
    double alpha1 = 1.0;                    // shape parameter of Z1
    double beta1 = 0.5;                     // shape parameter of Z1
    double gamma1 = 0.5;                    // shape parameter of Z1
    double eta = 1.0;                       // yielding exponent (sharpness of hysteresis loop corners)
    double alpha2 = 0.01;                   // shape parameter of Z2
    double beta2 = 0.08;                    // shape parameter of Z2
    double gamma2 = -0.16;                  // shape parameter of Z2
    double c1 = 1.8E-2;                     // parameter controls degradation
    int tagLH = 0;                          // boolean flag to control whether lead core heating is enabled (1) or disabled (0)
    double E2 = 0.0069;                     // shape parameter of lead core heating
    double rhoL = 11200.0;                  // density of lead
    double cL = 130.0;                      // specific heat of lead
    double kS = 50.0;                       // thermal conductivity of steel
    double alphaS = 1.41E-05;               // thermal diffusivity of steel
    double TL0 = 20.0;                      // initial temperature of lead core
    int tagIH = 0;                          // boolean flag to control whether initial hardening is enabled (1) or disabled (0)
    double c2 = 0.5;                        // parameter controls initial hardening
    double c3 = 1;                          // parameter controls initial hardening
    Vector x(0), y(3);                      // the orientation vector
    y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
    double sDratio = 0.5;                   // shear distance ratio from node I
    double mass = 0.0;                      // mass of element
    int maxIter = 500;                      // maximum number of iterations
    double tol = 1E-12;                     // tolerance for convergence criterion

    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        if (strcmp(type, "-alpha1") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &alpha1) != 0) {
                opserr << "WARNING: invalid alpha1 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-beta1") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &beta1) != 0) {
                opserr << "WARNING: invalid beta1 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-gamma1") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &gamma1) != 0) {
                opserr << "WARNING: invalid gamma1 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-eta") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &eta) != 0) {
                opserr << "WARNING: invalid eta for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-alpha2") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &alpha2) != 0) {
                opserr << "WARNING: invalid alpha2 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-beta2") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &beta2) != 0) {
                opserr << "WARNING: invalid beta2 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-gamma2") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &gamma2) != 0) {
                opserr << "WARNING: invalid gamma2 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-c1") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &c1) != 0) {
                opserr << "WARNING: invalid c1 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-tagLH") == 0) {
            tagLH = 1;
        }
        else if (strcmp(type, "-E2") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &E2) != 0) {
                opserr << "WARNING: invalid E2 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-rhoL") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &rhoL) != 0) {
                opserr << "WARNING: invalid rhoL for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-cL") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &cL) != 0) {
                opserr << "WARNING: invalid cL for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-kS") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &kS) != 0) {
                opserr << "WARNING: invalid kS for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-alphaS") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &alphaS) != 0) {
                opserr << "WARNING: invalid alphaS for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-TL0") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &TL0) != 0) {
                opserr << "WARNING: invalid TL0 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-tagIH") == 0) {
            tagIH = 1;
        }
        else if (strcmp(type, "-c2") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &c2) != 0) {
                opserr << "WARNING: invalid c2 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-c3") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &c3) != 0) {
                opserr << "WARNING: invalid c3 for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-orient") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 6) {
                opserr << "WARNING: insufficient args after -orient\n";
                return 0;
            }
            x.resize(3);
            numData = 3;
            if (OPS_GetDoubleInput(&numData, &x(0)) != 0) {
                opserr << "WARNING: invalid x vector for element" << eleTag << endln;
                return 0;
            }
            y.resize(3);
            numData = 3;
            if (OPS_GetDoubleInput(&numData, &y(0)) != 0) {
                opserr << "WARNING: invalid y vector for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-sDratio") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &sDratio) != 0) {
                opserr << "WARNING: invalid sDratio for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &mass) != 0) {
                opserr << "WARNING: invalid mass for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-maxIter") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetIntInput(&numData, &maxIter) != 0) {
                opserr << "WARNING: invalid maxIter for element" << eleTag << endln;
                return 0;
            }
        }
        else if (strcmp(type, "-tol") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: insufficient args\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &tol) != 0) {
                opserr << "WARNING: invalid tol for element" << eleTag << endln;
                return 0;
            }
        }
    }

    theEle = new BoucWenLRB(iData[0], iData[1], iData[2],
        dData[0], dData[1], dData[2],
        dData[3], dData[4], dData[5], dData[6],
        mats,
        alpha1, beta1,
        gamma1, eta,
        alpha2, beta2,
        gamma2, c1,
        tagLH,
        E2, rhoL,
        cL, kS,
        alphaS, TL0,
        tagIH,
        c2, c3,
        x, y,
        sDratio, mass,
        maxIter, tol);

    if (theEle == 0) {
        opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
        return 0;
    }

    return theEle;
}


BoucWenLRB::BoucWenLRB(int eleTag, int Nd1, int Nd2,
    double _kInit, double _fy, double _alphaL,
    double _dL, double _ts, double _tr, double _n,
    UniaxialMaterial** materials,
    double _alpha1, double _beta1,
    double _gamma1, double _eta,
    double _alpha2, double _beta2,
    double _gamma2, double _c1,
    int _tagLH,
    double _E2, double _rhoL,
    double _cL, double _kS,
    double _alphaS, double _TL0,
    int _tagIH,
    double _c2, double _c3,
    const Vector _x, const Vector _y,
    double _sDratio, double _mass,
    int _maxIter, double _tol)
    : Element(eleTag, ELE_TAG_BoucWenLRB), connectedExternalNodes(2),
    kInit(_kInit), fy(_fy), alphaL(_alphaL),
    dL(_dL), ts(_ts), tr(_tr), n(_n),
    alpha1(_alpha1), beta1(_beta1),
    gamma1(_gamma1), eta(_eta),
    alpha2(_alpha2), beta2(_beta2),
    gamma2(_gamma2), c1(_c1),
    tagLH(_tagLH),
    E2(_E2), rhoL(_rhoL),
    cL(_cL), kS(_kS),
    alphaS(_alphaS), TL0(_TL0),
    tagIH(_tagIH),
    c2(_c2), c3(_c3),
    x(_x), y(_y),
    sDratio(_sDratio), mass(_mass),
    maxIter(_maxIter), tol(_tol),
    L(0.0),
    ub(6), z1(2), z2(2), dz1du(2, 2), dz2du(2, 2), qb(6), Kb(6, 6), ul(12),
    Tgl(12, 12), Tlb(6, 12),
    ubC(6), z1C(2), z2C(2),
    KbInit(6, 6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
        opserr << "BoucWenLRB::BoucWenLRB() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

    // set node pointers to NULL
    for (int i = 0; i < 2; i++)
        theNodes[i] = 0;

    // horizontal motion
    rL = dL / 2.0;                                                                  // radius of lead core
    ALead = PI * rL * rL;                                                           // cross-sectional area of lead core
    Tr = n * tr;                                                                    // total thickness of rubber layers
    Ts = (n - 1) * ts;                                                              // total thickness of steel shims
    h = Tr + Ts;                                                                    // height of rubber + shims

    uy = fy / kInit;                                                                // yield displacement of bearing
    ke = alphaL * kInit;                                                            // stiffness of elastic component (due to rubber)
    k0 = (1 - alphaL) * kInit;                                                      // initial stiffness of hysteretic component (due to lead)

    // check material input
    if (materials == 0) {
        opserr << "BoucWenLRB::BoucWenLRB() - "
            << "null material array passed.\n";
        exit(-1);
    }

    // get copies of the uniaxial materials
    for (int i = 0; i < 4; i++) {
        if (materials[i] == 0) {
            opserr << "BoucWenLRB::BoucWenLRB() - "
                "null uniaxial material pointer passed.\n";
            exit(-1);
        }
        theMaterials[i] = materials[i]->getCopy();
        if (theMaterials[i] == 0) {
            opserr << "BoucWenLRB::BoucWenLRB() - "
                << "failed to copy uniaxial material.\n";
            exit(-1);
        }
    }

    // initialize initial stiffness matrix
    KbInit.Zero();
    KbInit(0, 0) = theMaterials[0]->getInitialTangent();
    KbInit(1, 1) = kInit;
    KbInit(2, 2) = kInit;
    KbInit(3, 3) = theMaterials[1]->getInitialTangent();
    KbInit(4, 4) = theMaterials[2]->getInitialTangent();
    KbInit(5, 5) = theMaterials[3]->getInitialTangent();

    // initialize variables
    this->revertToStart();
}


BoucWenLRB::BoucWenLRB()
    : Element(0, ELE_TAG_BoucWenLRB), connectedExternalNodes(2),
    kInit(0), fy(0), alphaL(0),
    dL(0), ts(0), tr(0), n(0),
    alpha1(1.0), beta1(0.5),
    gamma1(0.5), eta(1.0),
    alpha2(0.01), beta2(0.08),
    gamma2(-0.16), c1(1.8E-11),
    tagLH(0),
    E2(0.0069), rhoL(11200),
    cL(130), kS(50),
    alphaS(1.41E-05), TL0(20),
    tagIH(0),
    c2(0.5), c3(0.5),
    x(0), y(0),
    sDratio(0.5), mass(0),
    maxIter(500), tol(1E-12),
    L(0.0),
    ub(6), z1(2), z2(2), dz1du(2, 2), dz2du(2, 3), qb(6), Kb(6, 6), ul(12),
    Tgl(12, 12), Tlb(6, 12),
    ubC(6), z1C(2), z2C(2),
    KbInit(6, 6), theLoad(12)
{
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
        opserr << "BoucWenLRB::BoucWenLRB() - element: "
            << this->getTag() << " - failed to create an ID of size 2.\n";
        exit(-1);
    }

    // set node pointers to NULL
    for (int i = 0; i < 2; i++)
        theNodes[i] = 0;

    // initialize variables
    rL = 0.0;
    ALead = 0.0;
    Tr = 0.0;
    Ts = 0.0;
    h = 0.0;
    ke = 0.0;
    k0 = 0.0;
    uy = 0.0;

    // set material pointers to NULL
    for (int i = 0; i < 4; i++)
        theMaterials[i] = 0;

    this->revertToStart();
}


BoucWenLRB::~BoucWenLRB()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    for (int i = 0; i < 4; i++)
        if (theMaterials[i] != 0)
            delete theMaterials[i];
}


int BoucWenLRB::getNumExternalNodes() const
{
    return 2;
}


const ID& BoucWenLRB::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** BoucWenLRB::getNodePtrs()
{
    return theNodes;
}


int BoucWenLRB::getNumDOF()
{
    return 12;
}


void BoucWenLRB::setDomain(Domain* theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain) {
        theNodes[0] = 0;
        theNodes[1] = 0;

        return;
    }

    // first set the node pointers
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1]) {
        if (!theNodes[0]) {
            opserr << "WARNING BoucWenLRB::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        }
        else {
            opserr << "WARNING BoucWenLRB::setDomain() - Nd2: "
                << connectedExternalNodes(1)
                << " does not exist in the model for";
        }
        opserr << " element: " << this->getTag() << ".\n";

        return;
    }

    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    // if differing dof at the ends - print a warning message
    if (dofNd1 != 6) {
        opserr << "BoucWenLRB::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6) {
        opserr << "BoucWenLRB::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // set up the transformation matrix for orientation
    this->setUp();
}


int BoucWenLRB::commitState()
{
    int errCode = 0;

    // commit trial history variables
    ubC = ub;
    z1C = z1;
    z2C = z2;

    // commit trial history variables for degradation
    DSplusC = DSplus;
    DSminusC = DSminus;
    DSC = DS;

    // commit trial history variables for lead core heating
    dTLC = dTL;

    // commit trial history variables for initial hardening
    DIC = DI;

    // commit material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->commitState();

    // commit the base class
    errCode += this->Element::commitState();

    return errCode;
}


int BoucWenLRB::revertToLastCommit()
{
    int errCode = 0;

    // revert trial history variables
    ub = ubC;
    z1 = z1C;
    z2 = z2C;

    // revert trial history variables for degradation
    DSplus = DSplusC;
    DSminus = DSminusC;
    DS = DSC;

    // revert trial history variables for lead core heating
    dTL = dTLC;

    // revert trial history variables for initial hardening
    DI = DIC;

    // revert material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->revertToLastCommit();

    return errCode;
}


int BoucWenLRB::revertToStart()
{
    int errCode = 0;

    // reset trial history variables
    ub.Zero();
    z1.Zero();
    z2.Zero();
    qb.Zero();
    DSplus = 0.0;
    DSminus = 0.0;
    DS = 0.0;
    dTL = 0.0;
    DI = 0.0;

    // reset committed history variables
    ubC.Zero();
    z1C.Zero();
    z2C.Zero();
    DSplusC = 0.0;
    DSminusC = 0.0;
    DSC = 0.0;
    dTLC = 0.0;
    DIC = 0.0;


    // reset tangent of hysteretic evolution parameters
    dz1du(0, 0) = dz1du(1, 1) = alpha1;
    dz1du(1, 0) = dz1du(0, 1) = 0.0;
    dz2du(0, 0) = dz2du(1, 1) = alpha2;
    dz2du(1, 0) = dz2du(0, 1) = 0.0;

    // reset stiffness matrix in basic system
    Kb = KbInit;

    // revert material models
    for (int i = 0; i < 4; i++)
        errCode += theMaterials[i]->revertToStart();

    return errCode;
}


int BoucWenLRB::update()
{
    // get global trial displacements and velocities
    const Vector& dsp1 = theNodes[0]->getTrialDisp();
    const Vector& dsp2 = theNodes[1]->getTrialDisp();
    const Vector& vel1 = theNodes[0]->getTrialVel();
    const Vector& vel2 = theNodes[1]->getTrialVel();

    static Vector ug(12), ugdot(12), uldot(12), ubdot(6);
    for (int i = 0; i < 6; i++) {
        ug(i) = dsp1(i);  ugdot(i) = vel1(i);
        ug(i + 6) = dsp2(i);  ugdot(i + 6) = vel2(i);
    }

    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);

    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);

    // 1) get axial force and stiffess in basic x-direction
    theMaterials[0]->setTrialStrain(ub(0), ubdot(0));
    Kb(0, 0) = theMaterials[0]->getTangent();
    qb(0) = theMaterials[0]->getStress();

    // 2) calculate shear forces and stiffnesses in basic y- and z-direction

    // get displacement increments (trial - committed)
    Vector delta_ub = ub - ubC;
    // get horizontal displacement
    double uh = sqrt(pow(ub(1), 2) + pow(ub(2), 2));
    // get horontal velocity
    double vh = sqrt(pow(ubdot(1), 2) + pow(ubdot(2), 2));

    // Implement the degrading procedure here
    if (uh > DSplusC) {
        DSplus = uh;
        DSminus = DSminusC + (DSplus - DSplusC);
        DS = DSC;
    }
    if (uh < DSminusC) {
        DSplus = DSplusC;
        DSminus = uh;
        DS = DSC - (DSminus - DSminusC);
    }
    // get KL
    KL = exp(-c1 * pow(DS / Tr, 3));
    ke = KL * alphaL * kInit;
    k0 = (1 - KL * alphaL) * kInit;
    double qYield = (1 - KL * alphaL) * fy;

    // get KT1
    KT1 = 1.0;
    if (tagLH == 1) {
        updateCurrentDeltaTemp(vh);
        double TL = TL0 + dTL;
        if (TL <= 250) {
            KT1 = exp(-E2 * dTL);
        }
        else if (TL <= 327) {
            KT1 = exp(-E2 * (250 - TL0)) * (327 - TL) / (327 - 250);
        }
        else {
            KT1 = 0;
        }
    }

    // get KI
    KI = 1.0;
    if (tagIH == 1) {
        DI = DIC + sqrt(pow(delta_ub(1), 2) + pow(delta_ub(2), 2));
        KI = 1 + c2 * (DI / Tr) * exp(-c3 * DI / Tr);
    }

    if (sqrt(pow(delta_ub(1), 2) + pow(delta_ub(2), 2)) > DBL_EPSILON) {
        // calculate hysteretic evolution parameter z1 using Newton-Raphson
        int iter = 0;
        double z1Nrm, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
        Vector f1(2), delta_z1(2);
        Matrix Df1(2, 2);
        do {
            z1Nrm = z1.Norm();
            if (z1Nrm == 0.0)  // check because of negative exponents
                z1Nrm = DBL_EPSILON;
            tmp11 = gamma1 + beta1 * sgn(delta_ub(1) * z1(0));
            tmp12 = gamma1 + beta1 * sgn(delta_ub(2) * z1(1));
            tmp13 = delta_ub(1) * z1(0) * tmp11;
            tmp14 = delta_ub(2) * z1(1) * tmp12;
            tmp15 = pow(z1Nrm, eta - 2.0) * (tmp13 + tmp14);
            tmp16 = pow(z1Nrm, eta - 4.0) * (tmp13 + tmp14);

            // function and derivative
            f1(0) = z1(0) - z1C(0) - (1.0 / uy) * (alpha1 * delta_ub(1) - z1(0) * tmp15);
            f1(1) = z1(1) - z1C(1) - (1.0 / uy) * (alpha1 * delta_ub(2) - z1(1) * tmp15);

            Df1(0, 0) = 1.0 + (1 / uy) * (tmp15 + (eta - 2) * pow(z1(0), 2) * tmp16 + pow(z1Nrm, eta - 2) * tmp13);
            Df1(0, 1) = (z1(0) / uy) * ((eta - 2) * z1(1) * tmp16 + pow(z1Nrm, eta - 2) * delta_ub(2) * tmp12);
            Df1(1, 1) = 1.0 + (1 / uy) * (tmp15 + (eta - 2) * pow(z1(1), 2) * tmp16 + pow(z1Nrm, eta - 2) * tmp14);
            Df1(1, 0) = (z1(1) / uy) * ((eta - 2) * z1(0) * tmp16 + pow(z1Nrm, eta - 2) * delta_ub(1) * tmp11);

            // issue warning if diagonal of derivative Df1 is zero
            if ((fabs(Df1(0, 0)) <= DBL_EPSILON) || (fabs(Df1(1, 1)) <= DBL_EPSILON)) {
                opserr << "WARNING: BoucWenLRB::update() - "
                    << "zero Jacobian in Newton-Raphson scheme for hysteretic "
                    << "evolution parameter z1.\n";
                return -1;
            }

            // advance one step
            // delta_z1 = f1 / Df1;
            if (Df1.Solve(f1, delta_z1) < 0) {
                opserr << "WARNING: BoucWenLRB::update() - "
                    << "failed in solving for hysteretic evolution parameter z1.\n";
                return -1;
            }
            z1 -= delta_z1;
            iter++;
        } while ((delta_z1.Norm() >= tol) && (iter < maxIter));

        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter) {
            opserr << "WARNING: BoucWenLRB::update() - "
                << "did not find the hysteretic evolution parameters z1 after "
                << iter << " iterations and norm: " << delta_z1.Norm() << endln;
            return -2;
        }

        // calculate hysteretic evolution parameter z2 using Newton-Raphson
        iter = 0;
        double eta2 = 1;
        double z2Nrm, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26;
        Vector f2(2), delta_z2(2);
        Matrix Df2(2, 2);
        do {
            z2Nrm = z2.Norm();
            if (z2Nrm == 0.0)  // check because of negative exponents
                z2Nrm = DBL_EPSILON;
            tmp21 = gamma2 + beta2 * sgn(delta_ub(1) * z2(0));
            tmp22 = gamma2 + beta1 * sgn(delta_ub(2) * z2(1));
            tmp23 = delta_ub(1) * z2(0) * tmp21;
            tmp24 = delta_ub(2) * z2(1) * tmp22;
            tmp25 = pow(z2Nrm, eta2 - 2.0) * (tmp23 + tmp24);
            tmp26 = pow(z2Nrm, eta2 - 4.0) * (tmp23 + tmp24);

            // function and derivative
            f2(0) = z2(0) - z2C(0) - (1.0 / uy) * (alpha2 * delta_ub(1) - z2(0) * tmp25);
            f2(1) = z2(1) - z2C(1) - (1.0 / uy) * (alpha2 * delta_ub(2) - z2(1) * tmp25);

            Df2(0, 0) = 1.0 + (1 / uy) * (tmp25 + (eta2 - 2) * pow(z2(0), 2) * tmp26 + pow(z2Nrm, eta2 - 2) * tmp23);
            Df2(0, 1) = (z2(0) / uy) * ((eta2 - 2) * z2(1) * tmp26 + pow(z2Nrm, eta2 - 2) * delta_ub(2) * tmp22);
            Df2(1, 1) = 1.0 + (1 / uy) * (tmp25 + (eta2 - 2) * pow(z2(1), 2) * tmp26 + pow(z2Nrm, eta2 - 2) * tmp24);
            Df2(1, 0) = (z2(1) / uy) * ((eta2 - 2) * z2(0) * tmp26 + pow(z2Nrm, eta2 - 2) * delta_ub(1) * tmp21);

            // issue warning if diagonal of derivative Df2 is zero
            if ((fabs(Df2(0, 0)) <= DBL_EPSILON) || (fabs(Df2(1, 1)) <= DBL_EPSILON)) {
                opserr << "WARNING: BoucWenLRB::update() - "
                    << "zero Jacobian in Newton-Raphson scheme for hysteretic "
                    << "evolution parameter z2.\n";
                return -1;
            }

            // advance one step
            // delta_z2 = f2 / Df2;
            if (Df2.Solve(f2, delta_z2) < 0) {
                opserr << "WARNING: BoucWenLRB::update() - "
                    << "failed in solving for hysteretic evolution parameter z2.\n";
                return -1;
            }
            z2 -= delta_z2;
            iter++;
        } while ((delta_z2.Norm() >= tol) && (iter < maxIter));

        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter) {
            opserr << "WARNING: BoucWenLRB::update() - "
                << "did not find the hysteretic evolution parameters z2 after "
                << iter << " iterations and norm: " << delta_z2.Norm() << endln;
            return -2;
        }

        // get derivatives of hysteretic evolution parameters
        double du1du2, du2du1;
        if (delta_ub(1) * delta_ub(2) == 0) {
            du1du2 = 0.0;
            du2du1 = 0.0;
        }
        else {
            du1du2 = delta_ub(1) / delta_ub(2);
            du2du1 = delta_ub(2) / delta_ub(1);
        }

        double tmp17 = pow(z1Nrm, eta - 2.0) * z1(0) * tmp11;
        double tmp18 = pow(z1Nrm, eta - 2.0) * z1(1) * tmp12;
        dz1du(0, 0) = (1 / uy) * (alpha1 - z1(0) * (tmp17 + du2du1 * tmp18));
        dz1du(0, 1) = (1 / uy) * (alpha1 * du1du2 - z1(0) * (du1du2 * tmp17 + tmp18));
        dz1du(1, 1) = (1 / uy) * (alpha1 - z1(1) * (du1du2 * tmp17 + tmp18));
        dz1du(1, 0) = (1 / uy) * (alpha1 * du2du1 - z1(1) * (tmp17 + du2du1 * tmp18));

        double tmp27 = pow(z2Nrm, eta2 - 2.0) * z2(0) * tmp21;
        double tmp28 = pow(z2Nrm, eta2 - 2.0) * z2(1) * tmp22;
        dz2du(0, 0) = (1 / uy) * (alpha2 - z2(0) * (tmp27 + du2du1 * tmp28));
        dz2du(0, 1) = (1 / uy) * (alpha2 * du1du2 - z2(0) * (du1du2 * tmp27 + tmp28));
        dz2du(1, 1) = (1 / uy) * (alpha2 - z2(1) * (du1du2 * tmp27 + tmp28));
        dz2du(1, 0) = (1 / uy) * (alpha2 * du2du1 - z2(1) * (tmp27 + du2du1 * tmp28));

        // set shear force
        qb(1) = ke * ub(1) + KT1 * KI * qYield * z1(0) + ke * uy * z2(0);
        qb(2) = ke * ub(2) + KT1 * KI * qYield * z1(1) + ke * uy * z2(1);
        // set tangent stiffness
        Kb(1, 1) = ke + KT1 * KI * qYield * dz1du(0, 0) + ke * uy * dz2du(0, 0);
        Kb(1, 2) = KT1 * KI * qYield * dz1du(0, 1) + ke * uy * dz2du(0, 1);
        Kb(2, 1) = KT1 * KI * qYield * dz1du(1, 0) + ke * uy * dz2du(1, 0);
        Kb(2, 2) = ke + KT1 * KI * qYield * dz1du(1, 1) + ke * uy * dz2du(1, 1);
    }

    // 3) get moment and stiffness about basic x-direction
    theMaterials[1]->setTrialStrain(ub(3), ubdot(3));
    qb(3) = theMaterials[1]->getStress();
    Kb(3, 3) = theMaterials[1]->getTangent();

    // 4) get moment and stiffness about basic y-direction
    theMaterials[2]->setTrialStrain(ub(4), ubdot(4));
    qb(4) = theMaterials[2]->getStress();
    Kb(4, 4) = theMaterials[2]->getTangent();

    // 5) get moment and stiffness about basic z-direction
    theMaterials[3]->setTrialStrain(ub(5), ubdot(5));
    qb(5) = theMaterials[3]->getStress();
    Kb(5, 5) = theMaterials[3]->getTangent();

    return 0;
}


const Matrix& BoucWenLRB::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();

    // transform from basic to local system
    static Matrix Kl(12, 12);
    Kl.addMatrixTripleProduct(0.0, Tlb, Kb, 1.0);

    // add geometric stiffness to local stiffness
    double kGeo1 = 0.5 * qb(0);
    Kl(5, 1) -= kGeo1;
    Kl(5, 7) += kGeo1;
    Kl(11, 1) -= kGeo1;
    Kl(11, 7) += kGeo1;
    Kl(4, 2) += kGeo1;
    Kl(4, 8) -= kGeo1;
    Kl(10, 2) += kGeo1;
    Kl(10, 8) -= kGeo1;
    double kGeo2 = kGeo1 * sDratio * L;
    Kl(5, 5) += kGeo2;
    Kl(11, 5) -= kGeo2;
    Kl(4, 4) += kGeo2;
    Kl(10, 4) -= kGeo2;
    double kGeo3 = kGeo1 * (1.0 - sDratio) * L;
    Kl(5, 11) -= kGeo3;
    Kl(11, 11) += kGeo3;
    Kl(4, 10) -= kGeo3;
    Kl(10, 10) += kGeo3;

    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, Kl, 1.0);

    return theMatrix;
}


const Matrix& BoucWenLRB::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();

    // transform from basic to local system
    static Matrix KlInit(12, 12);
    KlInit.addMatrixTripleProduct(0.0, Tlb, KbInit, 1.0);

    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, KlInit, 1.0);

    return theMatrix;
}


const Matrix& BoucWenLRB::getDamp()
{
    // zero the matrix
    theMatrix.Zero();

    return theMatrix;
}


const Matrix& BoucWenLRB::getMass()
{
    // zero the matrix
    theMatrix.Zero();

    // check for quick return
    if (mass == 0.0) {
        return theMatrix;
    }

    double m = 0.5 * mass;
    for (int i = 0; i < 3; i++) {
        theMatrix(i, i) = m;
        theMatrix(i + 6, i + 6) = m;
    }

    return theMatrix;
}


void BoucWenLRB::zeroLoad()
{
    theLoad.Zero();
}


int BoucWenLRB::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "BoucWenLRB::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << ".\n";

    return -1;
}


int BoucWenLRB::addInertiaLoadToUnbalance(const Vector& accel)
{
    // check for quick return
    if (mass == 0.0) {
        return 0;
    }

    // get R * accel from the nodes
    const Vector& Raccel1 = theNodes[0]->getRV(accel);
    const Vector& Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr << "BoucWenLRB::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible.\n";
        return -1;
    }

    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5 * mass;
    for (int i = 0; i < 3; i++) {
        theLoad(i) -= m * Raccel1(i);
        theLoad(i + 6) -= m * Raccel2(i);
    }

    return 0;
}


const Vector& BoucWenLRB::getResistingForce()
{
    // zero the residual
    theVector.Zero();

    // determine resisting forces in local system
    static Vector ql(12);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);

    // add P-Delta moments to local forces
    double kGeo1 = 0.5 * qb(0);
    double MpDelta1 = kGeo1 * (ul(7) - ul(1));
    ql(5) += MpDelta1;
    ql(11) += MpDelta1;
    double MpDelta2 = kGeo1 * sDratio * L * ul(5);
    ql(5) += MpDelta2;
    ql(11) -= MpDelta2;
    double MpDelta3 = kGeo1 * (1.0 - sDratio) * L * ul(11);
    ql(5) -= MpDelta3;
    ql(11) += MpDelta3;
    double MpDelta4 = kGeo1 * (ul(8) - ul(2));
    ql(4) -= MpDelta4;
    ql(10) -= MpDelta4;
    double MpDelta5 = kGeo1 * sDratio * L * ul(4);
    ql(4) += MpDelta5;
    ql(10) -= MpDelta5;
    double MpDelta6 = kGeo1 * (1.0 - sDratio) * L * ul(10);
    ql(4) -= MpDelta6;
    ql(10) += MpDelta6;

    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);

    return theVector;
}


const Vector& BoucWenLRB::getResistingForceIncInertia()
{
    // this already includes damping forces from materials
    theVector = this->getResistingForce();

    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);

    // add inertia forces from element mass
    if (mass != 0.0) {
        const Vector& accel1 = theNodes[0]->getTrialAccel();
        const Vector& accel2 = theNodes[1]->getTrialAccel();

        double m = 0.5 * mass;
        for (int i = 0; i < 3; i++) {
            theVector(i) += m * accel1(i);
            theVector(i + 6) += m * accel2(i);
        }
    }

    return theVector;
}


int BoucWenLRB::sendSelf(int commitTag, Channel& sChannel)
{
    // send element parameters
    static Vector data(32);
    data(0) = this->getTag();
    data(1) = kInit;
    data(2) = fy;
    data(3) = alphaL;
    data(4) = dL;
    data(5) = ts;
    data(6) = tr;
    data(7) = n;
    data(8) = alpha1;
    data(9) = beta1;
    data(10) = gamma1;
    data(11) = eta;
    data(12) = alpha2;
    data(13) = beta2;
    data(14) = gamma2;
    data(15) = c1;
    data(16) = tagLH;
    data(17) = E2;
    data(18) = rhoL;
    data(19) = cL;
    data(20) = kS;
    data(21) = alphaS;
    data(22) = TL0;
    data(23) = tagIH;
    data(24) = c2;
    data(25) = c3;
    data(26) = x.Size();
    data(27) = y.Size();
    data(28) = sDratio;
    data(29) = mass;
    data(30) = maxIter;
    data(31) = tol;

    sChannel.sendVector(0, commitTag, data);

    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);

    // send the material class tags
    ID matClassTags(4);
    for (int i = 0; i < 4; i++)
        matClassTags(i) = theMaterials[i]->getClassTag();
    sChannel.sendID(0, commitTag, matClassTags);

    // send the material models
    for (int i = 0; i < 4; i++)
        theMaterials[i]->sendSelf(commitTag, sChannel);

    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);

    return 0;
}


int BoucWenLRB::recvSelf(int commitTag, Channel& rChannel,
    FEM_ObjectBroker& theBroker)
{
    // receive element parameters
    static Vector data(32);
    rChannel.recvVector(0, commitTag, data);

    this->setTag((int)data(0));
    kInit = data(1);
    fy = data(2);
    alphaL = data(3);
    dL = data(4);
    ts = data(5);
    tr = data(6);
    n = data(7);
    alpha1 = data(8);
    beta1 = data(9);
    gamma1 = data(10);
    eta = data(11);
    alpha2 = data(12);
    beta2 = data(13);
    gamma2 = data(14);
    c1 = data(15);
    tagLH = (int)data(16);
    E2 = data(17);
    rhoL = data(18);
    cL = data(19);
    kS = data(20);
    alphaS = data(21);
    TL0 = data(22);
    tagIH = (int)data(23);
    c2 = data(24);
    c3 = data(25);
    sDratio = data(28);
    mass = data(29);
    maxIter = (int)data(30);
    tol = data(31);

    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);

    // receive the material class tags
    ID matClassTags(4);
    rChannel.recvID(0, commitTag, matClassTags);

    // receive the material models
    for (int i = 0; i < 4; i++) {
        theMaterials[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
        if (theMaterials[i] == 0) {
            opserr << "BoucWenLRB::recvSelf() - "
                << "failed to get blank uniaxial material.\n";
            return -2;
        }
        theMaterials[i]->recvSelf(commitTag, rChannel, theBroker);
    }

    // receive remaining data
    if ((int)data(26) == 3) {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(27) == 3) {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }

    // horizontal motion
    ALead = PI * rL * rL;                                                           // cross-sectional area of lead core
    Tr = n * tr;                                                                    // total thickness of rubber layers
    Ts = (n - 1) * ts;                                                              // total thickness of steel shims
    h = Tr + Ts;                                                                    // height of rubber + shims

    uy = fy / kInit;                                                                // yield displacement of bearing
    ke = alphaL * kInit;                                                            // stiffness of elastic component (due to rubber)
    k0 = (1 - alphaL) * kInit;                                                      // initial stiffness of hysteretic component (due to lead)

    // initialize initial stiffness matrix
    KbInit.Zero();
    KbInit(0, 0) = theMaterials[0]->getInitialTangent();
    KbInit(1, 1) = kInit;
    KbInit(2, 2) = kInit;
    KbInit(3, 3) = theMaterials[1]->getInitialTangent();
    KbInit(4, 4) = theMaterials[2]->getInitialTangent();
    KbInit(5, 5) = theMaterials[3]->getInitialTangent();

    // initialize variables
    this->revertToStart();

    return 0;
}


int BoucWenLRB::displaySelf(Renderer& theViewer,
    int displayMode, float fact, const char** argv, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void BoucWenLRB::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln;
        s << "  type: BoucWenLRB\n";
        s << "  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  kInit: " << kInit << "  fy: " << fy
            << "  alphaL: " << alphaL << endln;
        s << "  dL: " << dL << "  ts: " << ts
            << "  tr: " << tr << "  n: " << n << endln;
        s << "  Material ux: " << theMaterials[0]->getTag() << endln;
        s << "  Material rx: " << theMaterials[1]->getTag() << endln;
        s << "  Material ry: " << theMaterials[2]->getTag() << endln;
        s << "  Material rz: " << theMaterials[3]->getTag() << endln;
        s << "  alpha1: " << alpha1 << "  beta1: " << beta1
            << "  gamma1: " << gamma1 << "  eta: " << eta << endln;
        s << "  alpha2: " << alpha2 << "  beta2: " << beta2
            << "  gamma2: " << gamma2 << "  c1: " << c1 << endln;
        s << "  tagLH: " << tagLH << endln;
        s << "  E2: " << E2 << "  rhoL: " << rhoL << endln;
        s << "  cL: " << cL << "  kS: " << kS << endln;
        s << "  alphaS: " << alphaS << "  TL0: " << TL0 << endln;
        s << "  tagIH: " << tagIH << endln;
        s << "  c2: " << c2 << "  c3: " << c3 << endln;
        s << "  sDratio: " << sDratio << "  mass: " << mass << endln;
        s << "  maxIter: " << sDratio << "  tol: " << mass << endln;
        s << "  k0: " << k0 << "  ke: " << ke << endln;
        // mechanical properties: large deformation and lead core heating
        s << "  DSplus: " << DSplus << "  DSminus: " << DSminus
            << "  DS: " << DS << endln;
        s << "  dTL: " << dTL << "  DI: " << DI << endln;
        s << "  KL: " << KL << "  KT1: " << KT1
            << "  KI: " << KI << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BoucWenLRB\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"kInit\": " << kInit << ", ";
        s << "\"fy\": " << fy << ", ";
        s << "\"alphaL\": " << alphaL << ", ";
        s << "\"dL\": " << dL << ", ";
        s << "\"ts\": " << ts << ", ";
        s << "\"tr\": " << tr << ", ";
        s << "\"n\": " << n << ", ";
        s << "\"materials\": [\"";
        s << theMaterials[0]->getTag() << "\", \"";
        s << theMaterials[1]->getTag() << "\", \"";
        s << theMaterials[2]->getTag() << "\", \"";
        s << theMaterials[3]->getTag() << "\"], ";
        s << "\"alpha1\": " << alpha1 << ", ";
        s << "\"beta1\": " << beta1 << ", ";
        s << "\"gamma1\": " << gamma1 << ", ";
        s << "\"eta\": " << eta << ", ";
        s << "\"alpha2\": " << alpha2 << ", ";
        s << "\"beta2\": " << beta2 << ", ";
        s << "\"gamma2\": " << gamma2 << ", ";
        s << "\"c1\": " << c1 << ", ";
        s << "\"tagLH\": " << tagLH << ", ";
        s << "\"E2\": " << E2 << ", ";
        s << "\"rhoL\": " << rhoL << ", ";
        s << "\"cL\": " << cL << ", ";
        s << "\"kS\": " << kS << ", ";
        s << "\"alphaS\": " << alphaS << ", ";
        s << "\"TL0\": " << TL0 << ", ";
        s << "\"tagIH\": " << tagIH << ", ";
        s << "\"c2\": " << c2 << ", ";
        s << "\"c3\": " << c3 << ", ";
        s << "\"sDratio\": " << sDratio << ", ";
        s << "\"mass\": " << mass << ", ";
        s << "\"maxIter\": " << maxIter << ", ";
        s << "\"tol\": " << tol << ", ";
        s << "\"k0\": " << k0 << ", ";
        s << "\"ke\": " << ke << ", ";
        s << "\"DSplus\": " << DSplus << ", ";
        s << "\"DSminus\": " << DSminus << ", ";
        s << "\"DS\": " << DS << ", ";
        s << "\"dTL\": " << dTL << ", ";
        s << "\"DI\": " << DI << ", ";
        s << "\"KL\": " << KL << ", ";
        s << "\"KT1\": " << KT1 << ", ";
        s << "\"KI\": " << KI << "}";
    }
}


Response* BoucWenLRB::setResponse(const char** argv, int argc,
    OPS_Stream& output)
{
    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "BoucWenLRB");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    // global forces
    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0)
    {
        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0)
    {
        output.tag("ResponseType", "N_ 1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Tz_1");
        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 2, theVector);
    }
    // basic forces
    else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0)
    {
        output.tag("ResponseType", "qb1");
        output.tag("ResponseType", "qb2");
        output.tag("ResponseType", "qb3");
        output.tag("ResponseType", "qb4");
        output.tag("ResponseType", "qb5");
        output.tag("ResponseType", "qb6");

        theResponse = new ElementResponse(this, 3, Vector(6));
    }
    // local displacements
    else if (strcmp(argv[0], "localDisplacement") == 0 ||
        strcmp(argv[0], "localDisplacements") == 0)
    {
        output.tag("ResponseType", "ux_1");
        output.tag("ResponseType", "uy_1");
        output.tag("ResponseType", "uz_1");
        output.tag("ResponseType", "rx_1");
        output.tag("ResponseType", "ry_1");
        output.tag("ResponseType", "rz_1");
        output.tag("ResponseType", "ux_2");
        output.tag("ResponseType", "uy_2");
        output.tag("ResponseType", "uz_2");
        output.tag("ResponseType", "rx_2");
        output.tag("ResponseType", "ry_2");
        output.tag("ResponseType", "rz_2");

        theResponse = new ElementResponse(this, 4, theVector);
    }
    // basic displacements
    else if (strcmp(argv[0], "deformation") == 0 || strcmp(argv[0], "deformations") == 0 ||
        strcmp(argv[0], "basicDeformation") == 0 || strcmp(argv[0], "basicDeformations") == 0 ||
        strcmp(argv[0], "basicDisplacement") == 0 || strcmp(argv[0], "basicDisplacements") == 0)
    {
        output.tag("ResponseType", "ub1");
        output.tag("ResponseType", "ub2");
        output.tag("ResponseType", "ub3");
        output.tag("ResponseType", "ub4");
        output.tag("ResponseType", "ub5");
        output.tag("ResponseType", "ub6");

        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // hysteretic evolution parameter
    else if (strcmp(argv[0], "hystereticParameter") == 0 || strcmp(argv[0], "hystParameter") == 0 ||
        strcmp(argv[0], "hystereticParam") == 0 || strcmp(argv[0], "hystParam") == 0 ||
        strcmp(argv[0], "z") == 0)
    {
        output.tag("ResponseType", "z11");
        output.tag("ResponseType", "z12");
        output.tag("ResponseType", "z21");
        output.tag("ResponseType", "z22");

        theResponse = new ElementResponse(this, 6, Vector(4));
    }
    // dzdu
    else if (strcmp(argv[0], "dzdu") == 0)
    {
        output.tag("ResponseType", "dz11du1");
        output.tag("ResponseType", "dz11du2");
        output.tag("ResponseType", "dz12du1");
        output.tag("ResponseType", "dz12du2");
        output.tag("ResponseType", "dz21du1");
        output.tag("ResponseType", "dz21du2");
        output.tag("ResponseType", "dz22du1");
        output.tag("ResponseType", "dz22du2");

        theResponse = new ElementResponse(this, 7, Vector(8));
    }
    // basic stiffness
    else if (strcmp(argv[0], "kb") == 0 || strcmp(argv[0], "basicStiff") == 0 ||
        strcmp(argv[0], "basicStiffness") == 0)
    {
        output.tag("ResponseType", "Kb22");
        output.tag("ResponseType", "Kb23");
        output.tag("ResponseType", "Kb32");
        output.tag("ResponseType", "Kb33");

        theResponse = new ElementResponse(this, 8, Vector(4));
    }
    // parameters that vary with time
    else if (strcmp(argv[0], "params") == 0 || strcmp(argv[0], "Params") == 0 ||
        strcmp(argv[0], "parameters") == 0 || strcmp(argv[0], "Parameters") == 0)
    {
        output.tag("ResponseType", "DSplus");
        output.tag("ResponseType", "DSminus");
        output.tag("ResponseType", "DS");
        output.tag("ResponseType", "dTL");
        output.tag("ResponseType", "DI");
        output.tag("ResponseType", "KL");
        output.tag("ResponseType", "KT1");
        output.tag("ResponseType", "KI");

        theResponse = new ElementResponse(this, 9, Vector(8));
    }

    output.endTag(); // ElementOutput

    return theResponse;
}


int BoucWenLRB::getResponse(int responseID, Information& eleInfo)
{
    double kGeo1, MpDelta1, MpDelta2, MpDelta3, MpDelta4, MpDelta5, MpDelta6;
    Vector z(4), dzduVec(8), kbVec(4), params(8);

    switch (responseID) {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());

    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta moments
        kGeo1 = 0.5 * qb(0);
        MpDelta1 = kGeo1 * (ul(7) - ul(1));
        theVector(5) += MpDelta1;
        theVector(11) += MpDelta1;
        MpDelta2 = kGeo1 * sDratio * L * ul(5);
        theVector(5) += MpDelta2;
        theVector(11) -= MpDelta2;
        MpDelta3 = kGeo1 * (1.0 - sDratio) * L * ul(11);
        theVector(5) -= MpDelta3;
        theVector(11) += MpDelta3;
        MpDelta4 = kGeo1 * (ul(8) - ul(2));
        theVector(4) -= MpDelta4;
        theVector(10) -= MpDelta4;
        MpDelta5 = kGeo1 * sDratio * L * ul(4);
        theVector(4) += MpDelta5;
        theVector(10) -= MpDelta5;
        MpDelta6 = kGeo1 * (1.0 - sDratio) * L * ul(10);
        theVector(4) -= MpDelta6;
        theVector(10) += MpDelta6;
        return eleInfo.setVector(theVector);

    case 3:  // basic forces
        return eleInfo.setVector(qb);

    case 4:  // local displacements
        return eleInfo.setVector(ul);

    case 5:  // basic displacements
        return eleInfo.setVector(ub);

    case 6:  // hysteretic evolution parameter
        z(0) = z1(0); z(1) = z1(1);
        z(2) = z2(0); z(3) = z2(1);
        return eleInfo.setVector(z);

    case 7:  // dzdu
        dzduVec(0) = dz1du(0, 0); dzduVec(1) = dz1du(0, 1);
        dzduVec(2) = dz1du(1, 0); dzduVec(3) = dz1du(1, 1);
        dzduVec(4) = dz2du(0, 0); dzduVec(5) = dz2du(0, 1);
        dzduVec(6) = dz2du(1, 0); dzduVec(7) = dz2du(1, 1);
        return eleInfo.setVector(dzduVec);

    case 8:  // basic stiffness
        kbVec(0) = Kb(1, 1); kbVec(1) = Kb(1, 2);
        kbVec(2) = Kb(2, 1); kbVec(3) = Kb(2, 2);
        return eleInfo.setVector(kbVec);

    case 9:  // parameters that vary with time
        params(0) = DSplus / Tr; params(1) = DSminus / Tr;
        params(2) = DS / Tr; params(3) = dTL;
        params(4) = DI; params(5) = KL;
        params(6) = KT1; params(7) = KI;
        return eleInfo.setVector(params);

    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void BoucWenLRB::setUp()
{
    const Vector& end1Crd = theNodes[0]->getCrds();
    const Vector& end2Crd = theNodes[1]->getCrds();
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();

    if (L > DBL_EPSILON) {
        if (x.Size() == 0) {
            x.resize(3);
            x = xp;
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3) {
        opserr << "BoucWenLRB::setUp() - "
            << "element: " << this->getTag() << endln
            << " - incorrect dimension of orientation vectors.\n";
        exit(-1);
    }

    // establish orientation of element for the transformation matrix
    // z = x cross y
    static Vector z(3);
    z(0) = x(1) * y(2) - x(2) * y(1);
    z(1) = x(2) * y(0) - x(0) * y(2);
    z(2) = x(0) * y(1) - x(1) * y(0);

    // y = z cross x
    y(0) = z(1) * x(2) - z(2) * x(1);
    y(1) = z(2) * x(0) - z(0) * x(2);
    y(2) = z(0) * x(1) - z(1) * x(0);

    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();

    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0) {
        opserr << "BoucWenLRB::setUp() - "
            << "element: " << this->getTag() << endln
            << " - invalid orientation vectors.\n";
        exit(-1);
    }

    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0, 0) = Tgl(3, 3) = Tgl(6, 6) = Tgl(9, 9) = x(0) / xn;
    Tgl(0, 1) = Tgl(3, 4) = Tgl(6, 7) = Tgl(9, 10) = x(1) / xn;
    Tgl(0, 2) = Tgl(3, 5) = Tgl(6, 8) = Tgl(9, 11) = x(2) / xn;
    Tgl(1, 0) = Tgl(4, 3) = Tgl(7, 6) = Tgl(10, 9) = y(0) / yn;
    Tgl(1, 1) = Tgl(4, 4) = Tgl(7, 7) = Tgl(10, 10) = y(1) / yn;
    Tgl(1, 2) = Tgl(4, 5) = Tgl(7, 8) = Tgl(10, 11) = y(2) / yn;
    Tgl(2, 0) = Tgl(5, 3) = Tgl(8, 6) = Tgl(11, 9) = z(0) / zn;
    Tgl(2, 1) = Tgl(5, 4) = Tgl(8, 7) = Tgl(11, 10) = z(1) / zn;
    Tgl(2, 2) = Tgl(5, 5) = Tgl(8, 8) = Tgl(11, 11) = z(2) / zn;

    // create transformation matrix from local to basic system (linear)
    Tlb.Zero();
    Tlb(0, 0) = Tlb(1, 1) = Tlb(2, 2) = Tlb(3, 3) = Tlb(4, 4) = Tlb(5, 5) = -1.0;
    Tlb(0, 6) = Tlb(1, 7) = Tlb(2, 8) = Tlb(3, 9) = Tlb(4, 10) = Tlb(5, 11) = 1.0;
    Tlb(1, 5) = -sDratio * L;
    Tlb(1, 11) = -(1.0 - sDratio) * L;
    Tlb(2, 4) = -Tlb(1, 5);
    Tlb(2, 10) = -Tlb(1, 11);
}


double BoucWenLRB::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}


void BoucWenLRB::updateCurrentDeltaTemp(double vh)
{
    // lead core heating
    double tCurrent = (this->getDomain())->getCurrentTime();

    double dt = (this->getDomain())->getDT();
    // if (dT>1) dT = 0;
    double tau = (alphaS * tCurrent) / (pow(rL, 2));
    double F;
    if (tau < 0.6) {
        F = 2.0 * sqrt(tau / PI) - (tau / PI) * (2.0 - (tau / 4.0) - pow(tau / 4.0, 2) - (15.0 / 4.0) * pow(tau / 4.0, 3));
    }
    else {
        F = 8.0 / (3.0 * PI) - (1.0 / (2.0 * sqrt(PI * tau))) * (1.0 - (1.0 / (12.0 * tau)) + (1.0 / (6.0 * pow(4.0 * tau, 2))) - (1.0 / (12.0 * pow(4.0 * tau, 3))));
    }
    double deltaT1 = (dt / (rhoL * cL * h)) * ((KT1 * (1 - alphaL) * fy * vh) / ALead - (kS * dTLC / rL) * (1.0 / F + 1.274 * (Ts / rL) * pow(tau, -1.0 / 3.0)));

    // use improved euler method to obtain final temperature
    double dTL1 = dTLC + deltaT1;
    double tCurrent2 = tCurrent + dt;
    tau = (alphaS * tCurrent2) / (pow(rL, 2));
    if (tau < 0.6) {
        F = 2.0 * sqrt(tau / PI) - (tau / PI) * (2.0 - (tau / 4.0) - pow(tau / 4.0, 2) - (15.0 / 4.0) * pow(tau / 4.0, 3));
    }
    else {
        F = 8.0 / (3.0 * PI) - (1.0 / (2.0 * sqrt(PI * tau))) * (1.0 - (1.0 / (12.0 * tau)) + (1.0 / (6.0 * pow(4.0 * tau, 2))) - (1.0 / (12.0 * pow(4.0 * tau, 3))));
    }
    double deltaT2 = (dt / (rhoL * cL * h)) * ((KT1 * (1 - alphaL) * fy * vh) / ALead - (kS * dTL1 / rL) * (1.0 / F + 1.274 * (Ts / rL) * pow(tau, -1.0 / 3.0)));

    dTL = dTLC + 0.5 * (deltaT1 + deltaT2);
    if (dTL < 0) dTL = 0;
}