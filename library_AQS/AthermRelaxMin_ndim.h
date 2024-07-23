//
//  AthermalRelaxMin.h
//  TimeDrivenSimLib
//
//  Created by Deng Pan on 2023/9/7.
//

#ifndef AthermalRelaxMin_h
#define AthermalRelaxMin_h

#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"
#include "FireMinAthermal_ndim.h"
#include "GeoFamily_ndim.h"

typedef struct AthermRelaxMin {
    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet, dt;
    double vdotf, vdotv, fdotf;
    double aveForce;

    double Tlimit, scaleVeloc;
    
    bool isInit;
} AthermRelaxMin;

//=================const box shape========
AthermRelaxMin *addAthermRelaxMin(Box *box, Particle *particle, Update *update, Variable *var);
int reInitAthermRelaxMin(AthermRelaxMin *fire, double T_limit);
AthermRelaxMin *getAthermRelaxMin(Update *update);
int minAthermRelax(Box *box, Particle *particle, Update *update);
int delAthermRelaxMin(Update *update);

#endif /* AthermalRelaxMin_h */
