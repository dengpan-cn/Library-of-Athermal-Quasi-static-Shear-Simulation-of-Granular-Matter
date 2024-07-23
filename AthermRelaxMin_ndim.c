//
//  AthermalRelaxMin.c
//  TimeDrivenSimLib
//
//  Created by Deng Pan on 2023/9/7.
//

#include "AthermRelaxMin_ndim.h"

#define DELAYSTEP 5
#define DT_GROW 1.1
#define ALPHA_SHRINK 0.99
#define DT_SHRINK 0.5
#define ALPHA0 0.1

void mAR_Integrate(Box *box, Particle *particle, Update *update, AthermRelaxMin *fire) {
    // VV Method: do first half step
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScale(veloc, fire->scaleVeloc, veloc);
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    particle->isForceValid = false;
}
void mAR_DotProduct(Box *box, Particle *particle, Update *update, AthermRelaxMin *fire) {
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    // double dotPrd = 0.0;
    //__float128 sumForce = 0;
    double sumForce = 0;  //?
    // do second half step then calculate dotProduct
    //double totKin = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        //totKin += sNormP2(veloc);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    fire->aveForce = sumForce / particle->nAtom / update->forceUnits;
    
    double Tcur = vdotv / (DIM * particle->nAtom);
    double fact = pow((5E-3 / update->timeUnits) / fire->dt, 2);
    fact = (fact > 1.0 ? 1.0 : fact);
    fire->scaleVeloc = (Tcur < fire->Tlimit * fact ? 1.0 : sqrt(fire->Tlimit * fact / Tcur));
}

AthermRelaxMin *addAthermRelaxMin(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getAthermRelaxMin(update) != NULL) {
        return getAthermRelaxMin(update);
    }
    AthermRelaxMin *fire = (AthermRelaxMin *)calloc(sizeof(AthermRelaxMin), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__AthermRelaxMin__");
    
    if (var == NULL) {
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
    
    cmdArg *cmd = findVariable(var, "atherm");
    if (cmd == NULL) {
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
        
    if(cmd->cmdArgc == 0){
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
    if (cmd->cmdArgc != 1) Abort("--atherm 1E-8");
    
    fire->dtSet = 5E-3;
    fire->Tlimit = atof(cmd->cmdArgv[0]);
    if (fire->Tlimit < 1E-20)
        Abort("upper limits of temperature is too small!");
    
    fire->isInit = true;
    return fire;
}
int reInitAthermRelaxMin(AthermRelaxMin *fire, double T_limit) {
    memset(fire, '\0', sizeof(AthermRelaxMin));
    
    fire->dtSet = 5E-3;
    fire->Tlimit = T_limit;
    if (fire->Tlimit < 1E-20)
        Abort("upper limits of temperature is too small!");
    
    fire->isInit = true;
    return 1;
}
AthermRelaxMin *getAthermRelaxMin(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__AthermRelaxMin__");
    if (whichTool < 0)
        return NULL;
    return (AthermRelaxMin *)update->toolkit.toolkit[whichTool];
}
int minAthermRelax(Box *box, Particle *particle, Update *update) {
    AthermRelaxMin *fire = getAthermRelaxMin(update);
    if (!fire || !fire->isInit)
        Abort("Call addAthermRelaxMin(...) or reInitAthermRelaxMin(...)!");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    update->nebrList.nDelay = 0;
    update->nebrList.skinSet = 0.1;
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    
    double dtmax = 0.5 * update->timeUnits;
    int currStep = 0, last_negative = 0;
    double alpha = ALPHA0;
    fire->dt = fire->dtSet * update->timeUnits;
    double dt0 = fire->dt;
    double ePairLastNeg = 0;
    int rVal = 0;
    
    calcForce(box, particle, update);
    mAR_DotProduct(box, particle, update, fire);
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    fire->scaleVeloc = 1.0;
    ePairLastNeg = update->ePair;
    
    while (!(fire->aveForce <= __ZeroForce__)) {
        // update x
//        if (fire->scaleVeloc < 1.0)
//            printf("%d %g\n", currStep, fire->scaleVeloc);
        mAR_Integrate(box, particle, update, fire);
        
        // update v
        calcForce(box, particle, update);
        mAR_DotProduct(box, particle, update, fire);
        currStep++;
        
        // mixing v and update timestep
        if (fire->vdotf > 0) {
            if (currStep - last_negative > DELAYSTEP) {
                fire->dt *= DT_GROW;
                fire->dt = cpuMin(fire->dt, dtmax);
                alpha *= ALPHA_SHRINK;
            }
            
            double scale1 = 1.0 - alpha;
            double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                doubleVecPtr veloc = particle->veloc[iatom];
                doubleVecPtr force = particle->force[iatom];
                vScale(veloc, scale1, veloc);
                vScaleAdd(veloc, veloc, scale2, force);
            }
        } else {
            int deltaStep = currStep - last_negative;
            if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 0.9;
                dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * update->timeUnits);
            } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 1.1;
                dtmax = cpuMin(dtmax, 0.5 * update->timeUnits);
            }
            
            fire->dt *= DT_SHRINK;
            fire->dt = cpuMax(fire->dt, 0.1 * fire->dtSet * update->timeUnits);
            
#ifdef Hertzian
            if (ePairLastNeg - update->ePair < __MinDeltaE__) {
                rVal = 2;
                break;
            }
#endif
            
            dt0 = fire->dt;
            last_negative = currStep;
            ePairLastNeg = update->ePair;
            alpha = ALPHA0;
            memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
            fire->scaleVeloc = 1.0;
        }
        
        if (update->ePair / update->energyUnits < __ZeroEnergy__) {
            rVal = 1;
            break;
        }
        
        if (currStep >= 1E7) {
            rVal = -1;
            break;
        }
    }
    
    fire->isInit = false;
    return rVal;
}
int delAthermRelaxMin(Update *update) {
    AthermRelaxMin *fire = getAthermRelaxMin(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__AthermRelaxMin__");
    return 0;
}

