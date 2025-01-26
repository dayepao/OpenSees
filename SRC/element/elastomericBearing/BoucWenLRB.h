/* ****************************************************************** **
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

#ifndef BoucWenLRB_h
#define BoucWenLRB_h

// Written: Yuhang Lu (ll05734823@gmail.com)
// Created: 2024/04/11
// Description: This file contains the class definition for BoucWenLRB
// element. BoucWenLRB is created to describe the lead rubber bearings
// with initial hardening, large strain stiffening, degradation and 
// strength reduction due to lead core heating

#include <Element.h>
#include <Matrix.h>

class Channel;
class UniaxialMaterial;
class Response;

class BoucWenLRB : public Element
{
public:
    // constructor
    BoucWenLRB(int eleTag, int Nd1, int Nd2,
        double kInit, double fy, double alphaL,
        double dL, double ts, double tr, double n,
        UniaxialMaterial** theMaterials,
        double alpha1 = 1.0, double beta1 = 0.5,
        double gamma1 = 0.5, double eta = 1.0,
        double alpha2 = 0.01, double beta2 = 0.08,
        double gamma2 = -0.16, double c1 = 1.8E-2,
        int tagLH = 0,
        double E2 = 0.0069, double rhoL = 11200.0,
        double cL = 130.0, double kS = 50.0,
        double alphaS = 1.41E-05, double TL0 = 20.0,
        int tagIH = 0,
        double c2 = 0.5, double c3 = 1,
        const Vector x = 0, const Vector y = 0,
        double sDratio = 0.5, double mass = 0.0,
        int maxIter = 500, double tol = 1E-12);
    BoucWenLRB();

    // destructor
    ~BoucWenLRB();

    // method to get class type
    const char* getClassType() const { return "BoucWenLRB"; };

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getDamp();
    const Matrix& getMass();

    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode);
    void Print(OPS_Stream& s, int flag = 0);

    // public methods for element recorder
    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& eleInfo);

protected:

private:
    // private methods
    void setUp();
    double sgn(double x);
    void updateCurrentDeltaTemp(double vh);

    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;          // contains the tags of the end nodes
    Node* theNodes[2];                  // array of nodes
    UniaxialMaterial* theMaterials[4];  // array of uniaxial materials

    // parameters
    double kInit;                       // initial stiffness of bearing
    double k0;                          // initial stiffness of hysteretic component (due to lead)
    double ke;                          // stiffness of elastic component (due to rubber)
    double uy;                          // yield displacement of bearing
    double fy;                          // yield force of bearing
    double alphaL;                      // stiffness ratio of elastic component

    double dL;                          // lead core diameter of the bearing
    double ts;                          // single steel shim layer thickness
    double tr;                          // single rubber layer thickness
    double n;                           // number of rubber layers

    double rL;                          // radius of lead core
    double ALead;                       // cross-sectional area of lead core
    double Tr;                          // total thickness of rubber layers
    double Ts;                          // total thickness of steel shims
    double h;                           // height of rubber + shims

    Vector x;                           // local x direction
    Vector y;                           // local y direction

    double alpha1, beta1, gamma1;       // shape parameter of Z1
    double eta;                         // yielding exponent (sharpness of hysteresis loop corners)
    double alpha2, beta2, gamma2;       // shape parameter of Z2
    double c1;                          // parameter controls degradation

    int tagLH;                          // Boolean flag to control whether lead core heating is enabled (1) or disabled (0)
    double E2;                          // shape parameter of lead core heating
    double rhoL, cL, kS, alphaS;        // heating parameters for lead core
    double TL0;                         // initial temperature of lead core

    int tagIH;                          // Boolean flag to control whether initial hardening is enabled (1) or disabled (0)
    double c2, c3;                      // parameters of initial hardening

    double sDratio;                     // shear distance ratio from node I
    double mass;                        // mass of element
    int maxIter;                        // maximum number of iterations
    double tol;                         // tolerance for convergence criterion
    double L;                           // element length

    // state variables
    Vector ub;                          // displacements in basic system
    Vector z1;                          // hysteretic evolution parameter Z1
    Vector z2;                          // evolution parameter Z2
    Matrix dz1du;                       // tangent of hysteretic evolution parameters z1
    Matrix dz2du;                       // tangent of hysteretic evolution parameters z2
    Vector qb;                          // forces in basic system
    Matrix Kb;                          // stiffness matrix in basic system
    Vector ul;                          // displacements in local system
    Matrix Tgl;                         // transformation matrix from global to local system
    Matrix Tlb;                         // transformation matrix from local to basic system
    double DSplus, DSminus, DS;         // parameters of degradation
    double dTL;                         // change in temperature of lead core from the reference temperature
    double DI;                          // cumulative shear deformation
    double KL;                          // correction coefficient for large deformation
    double KT1;                         // correction coefficient for lead core heating
    double KI;                          // correction coefficient for initial hardening

    // committed history variables
    Vector ubC;                         // displacements in basic system
    Vector z1C;                         // hysteretic evolution parameters
    Vector z2C;                         // evolution parameters
    double DSplusC, DSminusC, DSC;      // parameters of degradation
    double dTLC;                        // Change in temperature of lead core from the reference temperature
    double DIC;                         // cumulative shear deformation

    // initial stiffness matrix in basic system
    Matrix KbInit;

    static Matrix theMatrix;            // a class wide Matrix
    static Vector theVector;            // a class wide Vector
    Vector theLoad;
};

#endif
