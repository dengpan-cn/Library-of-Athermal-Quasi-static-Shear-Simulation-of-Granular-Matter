#ifndef __FireFunc_ndim_dengPan_h__
#define __FireFunc_ndim_dengPan_h__

#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ConstBoxShapeFIRE {
    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet;
    double vdotf, vdotv, fdotf;
    double aveForce;
    double dt;
    
    bool isInit;
} ConstBoxShapeFIRE;
typedef struct ConstPressIsoFIRE {
    double Pset, ptol;
    double gVeloc, sfact, gForce;

    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet;
    double vdotf, vdotv, fdotf;
    double aveForce;
    double dt;
    
    bool isInit;
} ConstPressIsoFIRE;

#if defined(__triBox__)
typedef struct ConstPressTriFIRE {
    double Pset, ptol;
    uptriMat gVeloc, gForce;

    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet;
    double vdotf, vdotv, fdotf;
    double aveForce;
    double dt;
    
    bool isInit;
} ConstPressTriFIRE;

typedef struct ConstVolStressFIRE {
    int shearDim, gradDim;
    double stressSet, stressTol;
    
    double deltaStrain, maxDeltaStrain;
    double gVeloc, gForce;
    
    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet, dt;
    double vdotf, vdotv, fdotf;
    double aveForce;
    
    int rtype;//0: converge, -1: reaching maximum strain.
    bool isInit;
} ConstVolStressFIRE;
#endif

#define VelocVerletFIRE

#define maxDeltaVF 1E-4
#define maxStrainRate 1E-6
#define __Ptol__ 1E-4
#define __Stol__ 1E-3

#define __ZeroStress__ 1E-14
#define __ZeroForce__ 1E-14
#define __ZeroEnergy__ 1E-20
#define __ZeroOverlap__ 1E-8

#ifdef Hertzian
#define __MinDeltaE__ 1E-16
#endif

#ifdef VelocVerletFIRE
//=================const box shape========
ConstBoxShapeFIRE *addConstBoxShapeFIRE(Box *box, Particle *particle, Update *update);
int reInitConstBoxShapeFIRE(ConstBoxShapeFIRE *fire);
ConstBoxShapeFIRE *getConstBoxShapeFIRE(Update *update);
int minConstBoxShapeFIRE(Box *box, Particle *particle, Update *update);
int delConstBoxShapeFIRE(Update *update);

// eliminating residual overlapping in unjammed configuration
int eliminateResOverlap(Box *box, Particle *particle, Update *update);

//=================const iso pressure========
ConstPressIsoFIRE *addConstPressIsoFIRE(Box *box, Particle *particle, Update *update, Variable *var);
int reInitConstPressIsoFIRE(ConstPressIsoFIRE*fire, double pTarget);
ConstPressIsoFIRE *getConstPressIsoFIRE(Update *update);
int minConstPressIsoFIRE(Box *box, Particle *particle, Update *update);
int delConstPressIsoFIRE(Update *update);

#if defined(__triBox__)
//=================const tri pressure========
ConstPressTriFIRE *addConstPressTriFIRE(Box *box, Particle *particle, Update *update, Variable *var);
int reInitConstPressTriFIRE(ConstPressTriFIRE *fire, double pTarget);
ConstPressTriFIRE *getConstPressTriFIRE(Update *update);
int minConstPressTriFIRE(Box *box, Particle *particle, Update *update);
int delConstPressTriFIRE(Update *update);

#if (DIM == 3 || DIM == 2)
//================const vol and stress=======
ConstVolStressFIRE *addConstVolStressFIRE(Box *box, Particle *particle, Update *update, Variable *var);
int reInitConstVolStressFIRE(ConstVolStressFIRE *fire, double tStress, char *sType, double maxDeltaStrain);
ConstVolStressFIRE *getConstVolStressFIRE(Update *update);
int minConstVolStressFIRE(Box *box, Particle *particle, Update *update);
int delConstVolStressFIRE(Update *update);
#endif


#endif

#endif

#ifdef __cplusplus
}
#endif

#endif
