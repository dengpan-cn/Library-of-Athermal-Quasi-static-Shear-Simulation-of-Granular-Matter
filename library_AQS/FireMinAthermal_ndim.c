#include "FireMinAthermal_ndim.h"

#ifdef __cplusplus
#error "This file must NOT be compiled by C++ compliler!"
#endif

#define DELAYSTEP 5
#define DT_GROW 1.1
#define ALPHA_SHRINK 0.99
#define DT_SHRINK 0.5
#define ALPHA0 0.1

typedef struct pureShear {
    FILE *shearLog;
    
    // positive strain: extension in x direction and compression in z direction;
    double strain, maxStrain;
    // pure shear: Lx,t = x * Lx,0; Lz,t = 1/x * Lz,0;
    // strain = ln(Lx,t/Lz,t) - ln(Lx,0/Lz,0) = ln(x^2);
    double deltaStrain;
    
    bool isInit;
} pureShear;

#ifdef VelocVerletFIRE
//=================const box shape========

void cbsFire_Integrate(Box *box, Particle *particle, ConstBoxShapeFIRE *fire) {
    // VV Method: do first half step
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    particle->isForceValid = false;
}
void cbsFire_DotProduct(Box *box, Particle *particle, ConstBoxShapeFIRE *fire) {
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    // double dotPrd = 0.0;
    //__float128 sumForce = 0;
    double sumForce = 0;  //?
    // do second half step then calculate dotProduct
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    fire->aveForce = sumForce / particle->nAtom;
}

ConstBoxShapeFIRE *addConstBoxShapeFIRE(Box *box, Particle *particle, Update *update) {
    if (getConstBoxShapeFIRE(update) != NULL) {
        return getConstBoxShapeFIRE(update);
    }
    ConstBoxShapeFIRE *fire = (ConstBoxShapeFIRE *)calloc(sizeof(ConstBoxShapeFIRE), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__ConstBoxShapeFireRelax__");
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return fire;
}
int reInitConstBoxShapeFIRE(ConstBoxShapeFIRE *fire){
    memset(fire, '\0', sizeof(ConstBoxShapeFIRE));
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return 1;
}
ConstBoxShapeFIRE *getConstBoxShapeFIRE(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__ConstBoxShapeFireRelax__");
    if (whichTool < 0)
        return NULL;
    return (ConstBoxShapeFIRE *)update->toolkit.toolkit[whichTool];
}
int minConstBoxShapeFIRE(Box *box, Particle *particle, Update *update) {
    ConstBoxShapeFIRE *fire = getConstBoxShapeFIRE(update);
    if (!fire || !fire->isInit) Abort("Call addConstBoxShapeFireRelax(...) or reInitConstBoxShapeFireRelax(...)!");
    
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
    cbsFire_DotProduct(box, particle, fire);
    fire->aveForce /= update->forceUnits;
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    ePairLastNeg = update->ePair;
    while (!(fire->aveForce <= __ZeroForce__)) {
        // update x
        cbsFire_Integrate(box, particle, fire);
        
        // update v
        calcForce(box, particle, update);
        cbsFire_DotProduct(box, particle, fire);
        fire->aveForce /= update->forceUnits;
        currStep++;
        
        // mixing v and update timestep
        if (fire->vdotf > 0) {
            if (currStep - last_negative > DELAYSTEP) {
                fire->dt *= DT_GROW;
                fire->dt = fmin(fire->dt, dtmax);
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
                dtmax = fmax(dtmax, 5.0 * fire->dtSet * update->timeUnits);
            } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 1.1;
                dtmax = fmin(dtmax, 0.5 * update->timeUnits);
            }
            
            fire->dt *= DT_SHRINK;
            fire->dt = fmax(fire->dt, 0.1 * fire->dtSet * update->timeUnits);
            
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
        }
        
        if (update->ePair / update->energyUnits < __ZeroEnergy__) {
            rVal = 1;
            break;
        }
        
//        if (currStep % 1000==0) {
//            printf("%d %g %g\n",currStep,update->ePair,fire->aveForce);
//        }
        if (currStep >= 1E7) {
            rVal = -1;
            break;
        }
    }
    
    fire->isInit = false;
    return rVal;
    //0: force converge;
    //1: unjamming due to energy criteria；
    //2: deltaE between two successive over runing step is less than 1E-16;
    //-1: not converge and return due to runing step limits 1E7;
}
int delConstBoxShapeFIRE(Update *update) {
    ConstBoxShapeFIRE *fire = getConstBoxShapeFIRE(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__ConstBoxShapeFireRelax__");
    return 0;
}

//===============================
void eliminate_fire(Box *box, Particle *particle, Update *update) {
    ConstBoxShapeFIRE *fire = getConstBoxShapeFIRE(update);
    reInitConstBoxShapeFIRE(fire);
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    
    double dtmax = 0.5 * update->timeUnits;
    int currStep = 0, last_negative = 0;
    double alpha = ALPHA0;
    fire->dt = fire->dtSet * update->timeUnits;
    double dt0 = fire->dt;
    
    calcForce(box, particle, update);
    cbsFire_DotProduct(box, particle, fire);
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    while (true) {
        // update x
        cbsFire_Integrate(box, particle, fire);
        
        // update v
        calcForce(box, particle, update);
        cbsFire_DotProduct(box, particle, fire);
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
            
            dt0 = fire->dt;
            last_negative = currStep;
            alpha = ALPHA0;
            memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
        }
        
        if (update->maxOvlp <= __ZeroOverlap__)
            break;  // unjamming
        if (currStep > 1E6 && (fire->aveForce / update->forceUnits <= __ZeroForce__))
            break;
    }
}
int eliminateResOverlap(Box *box, Particle *particle, Update *update) {
    ConstBoxShapeFIRE *fire = getConstBoxShapeFIRE(update);
    if (!fire) {
        fire = addConstBoxShapeFIRE(box, particle, update);
    }else{
        reInitConstBoxShapeFIRE(fire);
    }
    minConstBoxShapeFIRE(box, particle, update);
    
    contactInfo *cinfo = getContactInfo(update);
    if (!cinfo) {
        cinfo = addContactInfo(box, particle, update);
    }
    computeContactInfo(box, particle, update);
    
    int nNR = particle->nAtom - cinfo->nRattler;
    if (nNR != 0) {
        if (cinfo->aveCoordNumExRattler - 2.0 * DIM + 2.0 * (DIM - 1.0) / (double)nNR > -1E-5)
            return -1;  // jammed, with finite size and boundary correction.
    }
    
    // expanding the particles
    double dVFmin = (1.0 - pow((1.0 - __ZeroOverlap__), DIM)) / (1.0 + 1.0 - pow((1.0 - __ZeroOverlap__), DIM)) * update->volFrac * 1.2;
    instant_inflate(box, particle, update, dVFmin);
    eliminate_fire(box, particle, update);
    instant_inflate(box, particle, update, -dVFmin);
    
    calcForce(box, particle, update);
    
    // check overlapping
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        for (int jatom = iatom + 1; jatom < particle->nAtom; jatom++) {
            doubleVector dRij;
            vSub(dRij, particle->pos[iatom], particle->pos[jatom]);
            PBC(dRij, box);
            double sRc = (particle->diameterScale[iatom] + particle->diameterScale[jatom]) * particle->meanDiameter * 0.5;
            double rij = sNorm(dRij);
            if (rij < sRc) {
                return -1;
            }
        }
    }
    
    return 0;
}

//=================const iso pressure========
void cpiFire_Integrate(Box *box, Particle *particle, Update *update, ConstPressIsoFIRE *fire) {
    // VV Method: do first half step
    //==============volume=============
    double relativeP = fabs(update->pVir / (fire->Pset * update->pressureUnits) - 1.0);
    double Qmass = particle->nAtom * (exp(-relativeP / fire->ptol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    fire->sfact = exp(fire->dt * fire->gVeloc);
    double phi = update->volFrac / pow(fire->sfact, DIM);
    if (update->volFrac - phi > maxDeltaVF) {
        fire->sfact =
        pow(update->volFrac / (update->volFrac - maxDeltaVF), 1.0 / DIM);
        fire->gVeloc = log(fire->sfact) / fire->dt;
    } else if (phi - update->volFrac > maxDeltaVF) {
        fire->sfact =
        pow(update->volFrac / (update->volFrac + maxDeltaVF), 1.0 / DIM);
        fire->gVeloc = log(fire->sfact) / fire->dt;
    }
    update->volFrac /= pow(fire->sfact, DIM);
    particle->meanDiameter /= fire->sfact;
    setUnits(update, particle->meanDiameter);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
}
void cpiFire_DotProduct(Box *box, Particle *particle, Update *update, ConstPressIsoFIRE *fire) {
    // do second half step then calculate dotProduct
    fire->gForce =
    (update->pVir - fire->Pset * update->pressureUnits) * box->volume;
    
    double relativeP =
    fabs(update->pVir / (fire->Pset * update->pressureUnits) - 1.0);
    double Qmass = particle->nAtom * (exp(-relativeP / fire->ptol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    // double dotPrd = 0.0;
    double sumForce = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf + fire->gVeloc * fire->gForce;
    fire->vdotv = vdotv + fire->gVeloc * fire->gVeloc;
    fire->fdotf = fdotf + fire->gForce * fire->gForce;
    fire->aveForce = sumForce / update->forceUnits / particle->nAtom;
}

ConstPressIsoFIRE *addConstPressIsoFIRE(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getConstPressIsoFIRE(update) != NULL) {
        return getConstPressIsoFIRE(update);
    }
    
    ConstPressIsoFIRE *fire = (ConstPressIsoFIRE *)calloc(sizeof(ConstPressIsoFIRE), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__ConstPressIsoFireRelax__");
    
    cmdArg *cmd = findVariable(var, "cpmin");
    if(!cmd) Abort("--cpmin pTarget or --cpmin");
    if(cmd->cmdArgc == 0){
        fire->isInit = false;
        return fire;
    }
    
    if (cmd->cmdArgc != 1) Abort("--cpmin pTarget or --cpmin");
    
    fire->Pset = atof(cmd->cmdArgv[0]);
    fire->ptol = __Ptol__;
    if (fire->Pset < 1E-8) {
        Abort("The Pset is too small!");
    }
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return fire;
}
int reInitConstPressIsoFIRE(ConstPressIsoFIRE*fire, double pTarget) {
    memset(fire, '\0', sizeof(ConstPressIsoFIRE));
    
    fire->Pset = pTarget;
    fire->ptol = __Ptol__;
    if (fire->Pset < 1E-8) {
        Abort("The Pset is too small!");
    }
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return 1;
}
ConstPressIsoFIRE *getConstPressIsoFIRE(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__ConstPressIsoFireRelax__");
    if (whichTool < 0)
        return NULL;
    return (ConstPressIsoFIRE *)update->toolkit.toolkit[whichTool];
}
int minConstPressIsoFIRE(Box *box, Particle *particle, Update *update) {
    ConstPressIsoFIRE *fire = getConstPressIsoFIRE(update);
    if (!fire) Abort("Call addConstPressIsoFireRelax(...) or reInitConstPressIsoFireRelax(...)");
    if(!fire->isInit) Abort("Call addConstPressIsoFireRelax(...) or reInitConstPressIsoFireRelax(...)");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = false;
    update->nebrList.nDelay = 0;
    update->nebrList.skinSet = 0.1;
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    
    double dtmax = 0.5 * update->timeUnits;
    int currStep = 0, last_negative = 0;
    double alpha = ALPHA0;
    double ePairLastNeg = 0.0;
    int rVal = 0;
    
    fire->dt = fire->dtSet * update->timeUnits;
    double dt0 = fire->dt;
    calcForce(box, particle, update);
    cpiFire_DotProduct(box, particle, update, fire);
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    fire->gVeloc = 0;
    ePairLastNeg = update->ePair;
    
    double relativeP = fabs(update->pVir / (fire->Pset * update->pressureUnits) - 1.0);
    while (!(relativeP <= fire->ptol && fire->aveForce <= __ZeroForce__)) {
        // update x and box
        cpiFire_Integrate(box, particle, update, fire);
        
        // calc force
        calcForce(box, particle, update);
        cpiFire_DotProduct(box, particle, update, fire);
        relativeP = fabs(update->pVir / (fire->Pset * update->pressureUnits) - 1.0);
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
            fire->gVeloc = fire->gVeloc * scale1 + fire->gForce * scale2;
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
            fire->gVeloc = 0;
        }
    }
    
    fire->isInit = false;
    return rVal;
}
int delConstPressIsoFIRE(Update *update) {
    ConstPressIsoFIRE *fire = getConstPressIsoFIRE(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__ConstPressIsoFireRelax__");
    return 0;
}

#if defined(__triBox__)
//=================const tri pressure========
void cptFire_Integrate(Box *box, Particle *particle, Update *update, ConstPressTriFIRE *fire) {
    uptriMat sfact;
    double phi = update->volFrac;
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            int sidx = spaceIdx2voigt(idim, jdim);
            if (idim == jdim) {
                double relativeP = fabs(update->pVirTens[sidx] / fire->Pset - 1.0);
                double Qmass = DIM * particle->nAtom * (exp(-relativeP / fire->ptol) + 1.0);
                
                fire->gVeloc[sidx] += 0.5 * fire->dt * fire->gForce[sidx] / Qmass;
                sfact[sidx] = exp(fire->dt * fire->gVeloc[sidx]);
                phi /= sfact[sidx];
            } else {
                double Qmass = DIM * particle->nAtom * (exp(fabs(update->pVirTens[sidx] / update->pressureUnits)) + 1.0);
                fire->gVeloc[sidx] += 0.5 * fire->dt * fire->gForce[sidx] / Qmass;
                sfact[sidx] = fire->dt * fire->gVeloc[sidx];
            }
        }
    }
    if (fabs(phi - update->volFrac) > maxDeltaVF) {
        double scale = 1.0;
        if (phi > update->volFrac) {
            scale = pow(phi / (update->volFrac + maxDeltaVF), 1.0 / DIM);
        } else {
            scale = pow(phi / (update->volFrac - maxDeltaVF), 1.0 / DIM);
        }
        for (int idim = 0; idim < DIM; idim++) {
            sfact[spaceIdx2voigt(idim, idim)] *= scale;
            fire->gVeloc[spaceIdx2voigt(idim, idim)] =
            log(sfact[spaceIdx2voigt(idim, idim)]) / fire->dt;
        }
    }
    
    //=============update box=========
    for (int iedge = 0; iedge < DIM; iedge++) {
        for (int idim = 0; idim < DIM; idim++) {
            box->boxEdge[iedge][idim] *= sfact[spaceIdx2voigt(idim, idim)];
            for (int jdim = idim + 1; jdim < DIM; jdim++) {
                box->boxEdge[iedge][idim] += sfact[spaceIdx2voigt(idim, jdim)] * box->boxEdge[iedge][jdim];
            }
        }
    }
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            box->boxH[spaceIdx2voigt(idim, jdim)] = box->boxEdge[jdim][idim];
        }
        update->volFrac = update->volFrac / sfact[spaceIdx2voigt(idim, idim)];
    }
    setBoxPara(box);
    //=============update particle=========
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        for (int idim = 0; idim < DIM; idim++) {
            pos[idim] *= sfact[spaceIdx2voigt(idim, idim)];
            for (int jdim = idim + 1; jdim < DIM; jdim++) {
                pos[idim] += sfact[spaceIdx2voigt(idim, jdim)] * pos[jdim];
            }
        }
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
}
void cptFire_DotProduct(Box *box, Particle *particle, Update *update, ConstPressTriFIRE *fire) {
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            int sidx = spaceIdx2voigt(idim, jdim);
            if (idim == jdim) {
                fire->gForce[sidx] = (update->pVirTens[sidx] - fire->Pset) * box->volume;
                
                double relativeP = fabs(update->pVirTens[sidx] / fire->Pset - 1.0);
                double Qmass = DIM * particle->nAtom * (exp(-relativeP / fire->ptol) + 1.0);
                
                fire->gVeloc[sidx] += 0.5 * fire->dt * fire->gForce[sidx] / Qmass;
            } else {
                fire->gForce[sidx] = update->pVirTens[sidx] * box->volume;
                
                double Qmass = DIM * particle->nAtom * (exp(fabs(update->pVirTens[sidx]/update->pressureUnits)) + 1.0);
                fire->gVeloc[sidx] += 0.5 * fire->dt * fire->gForce[sidx] / Qmass;
            }
        }
    }
    
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    // double dotPrd = 0.0;
    double sumForce = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    
    for (int sidx = 0; sidx < (DIM * (DIM + 1) / 2); sidx++) {
        fire->vdotf += fire->gVeloc[sidx] * fire->gForce[sidx];
        fire->vdotv += fire->gVeloc[sidx] * fire->gVeloc[sidx];
        fire->fdotf += fire->gForce[sidx] * fire->gForce[sidx];
    }
    
    fire->aveForce = sumForce / update->forceUnits / particle->nAtom;
}
bool cptFire_CheckConverge(Box *box, Particle *particle, Update *update, ConstPressTriFIRE *fire) {
    if (fire->aveForce > __ZeroForce__)
        return false;
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            int sidx = spaceIdx2voigt(idim, jdim);
            if (idim == jdim) {
                double fact = fabs(update->pVirTens[sidx] / fire->Pset - 1.0);
                if (fact > fire->ptol) {
                    return false;
                }
            } else {
                double fact = fabs(update->pVirTens[sidx] / update->pressureUnits);
                if (fact > __ZeroStress__) {
                    return false;
                }
            }
        }
    }
    return true;
}

ConstPressTriFIRE *addConstPressTriFIRE(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getConstPressTriFIRE(update) != NULL) {
        return getConstPressTriFIRE(update);
    }
    
    ConstPressTriFIRE *fire = (ConstPressTriFIRE *)calloc(sizeof(ConstPressTriFIRE), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__ConstPressTriFireRelax__");
    
    cmdArg *cmd = findVariable(var, "cpmin");
    if(!cmd) Abort("--cpmin pTarget or --cpmin");
    if(cmd->cmdArgc == 0){
        fire->isInit = false;
        return fire;
    }
    
    if (cmd->cmdArgc != 1) Abort("--cpmin pTarget or --cpmin");

    fire->Pset = atof(cmd->cmdArgv[0]);
    fire->ptol = __Ptol__;
    if (fire->Pset < 1E-8) {
        Abort("The Pset is too small!");
    }
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return fire;
}
int reInitConstPressTriFIRE(ConstPressTriFIRE *fire, double pTarget){
    memset(fire, '\0', sizeof(ConstPressTriFIRE));
    
    fire->Pset = pTarget;
    fire->ptol = __Ptol__;
    if (fire->Pset < 1E-8) {
        Abort("The Pset is too small!");
    }
    
    fire->dtSet = 5E-3;
    fire->isInit = true;
    return 1;
}
ConstPressTriFIRE *getConstPressTriFIRE(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__ConstPressTriFireRelax__");
    if (whichTool < 0)
        return NULL;
    return (ConstPressTriFIRE *)update->toolkit.toolkit[whichTool];
}
int minConstPressTriFIRE(Box *box, Particle *particle, Update *update) {
    ConstPressTriFIRE *fire = getConstPressTriFIRE(update);
    if (!fire || !fire->isInit) Abort("Call addConstPressTriFireRelax(...) or reInitConstPressTriFireRelax(...)!");

    box->isShapeFixed = false;
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
    double Pset = fire->Pset;
    fire->Pset *= update->pressureUnits;//incorporationg units
    double ePairLastNeg = update->ePair;
    int rVal = 0;
    
    calcForce(box, particle, update);
    cptFire_DotProduct(box, particle, update, fire);
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    uptriMatZeros(fire->gVeloc);
    ePairLastNeg = update->ePair;
    bool isConverge = cptFire_CheckConverge(box, particle, update, fire);
    while (!isConverge) {
        // update x and box
        cptFire_Integrate(box, particle, update, fire);
        
        // check timestep
        calcForce(box, particle, update);
        cptFire_DotProduct(box, particle, update, fire);
        isConverge = cptFire_CheckConverge(box, particle, update, fire);
        currStep++;
        
        // update velocity
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
            for (int sidx = 0; sidx < (DIM * (DIM + 1) / 2); sidx++) {
                fire->gVeloc[sidx] =
                scale1 * fire->gVeloc[sidx] + scale2 * fire->gForce[sidx];
            }
        } else {
            int deltaStep = currStep - last_negative;
            if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 0.9;
                dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * update->timeUnits);
            } else if (deltaStep > 2.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
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
            uptriMatZeros(fire->gVeloc);
        }
    }
    
    fire->Pset = Pset;
    fire->isInit = false;
    return rVal;
}
int delConstPressTriFIRE(Update *update) {
    ConstPressTriFIRE *fire = getConstPressTriFIRE(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__ConstPressTriFireRelax__");
    return 0;
}

#if (DIM == 3 || DIM == 2)
//==============
void cvsFire_DotProduct(Box *box, Particle *particle, Update *update, ConstVolStressFIRE *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    fire->gForce = (update->pVirTens[spaceIdx2voigt(sdim, gdim)] + fire->stressSet) * box->volume;
    double fact = fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0);
    double Qmass = DIM * particle->nAtom * (exp(-fact / fire->stressTol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    double sumForce = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    
    fire->vdotf += fire->gVeloc * fire->gForce;
    fire->vdotv += fire->gVeloc * fire->gVeloc;
    fire->fdotf += fire->gForce * fire->gForce;
    
    fire->aveForce = sumForce / update->forceUnits / particle->nAtom;
}
void cvsFire_Integrate(Box *box, Particle *particle, Update *update, ConstVolStressFIRE *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    double fact = fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0);
    double Qmass = DIM * particle->nAtom * (exp(-fact / fire->stressTol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    double deltaStrain = fire->dt * fire->gVeloc;
    double srate = fabs(deltaStrain) / (fire->dt / update->timeUnits);
    if (srate > maxStrainRate) {
        fire->gVeloc = deltaStrain / fabs(deltaStrain) * maxStrainRate / update->timeUnits;
        deltaStrain = fire->dt * fire->gVeloc;
    }
    fire->deltaStrain += deltaStrain;
    
    box->boxH[spaceIdx2voigt(sdim, gdim)] += deltaStrain * box->boxH[spaceIdx2voigt(gdim, gdim)];
    #if (DIM == 3)
    if (sdim == 0 && gdim == 1) {  // xy
        box->boxH[spaceIdx2voigt(0, 2)] += deltaStrain * box->boxH[spaceIdx2voigt(1, 2)];
    }
    #endif
    setBoxPara(box);
    
    //=============update particle=========
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        pos[sdim] += deltaStrain * pos[gdim];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
}
bool cvsFire_CheckConverge(Box *box, Particle *particle, Update *update, ConstVolStressFIRE *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    if (fabs(fire->deltaStrain) > fire->maxDeltaStrain) {
        fire->rtype = -1;
        return true;
    }
    
    if (fire->aveForce > __ZeroForce__)
        return false;
    if (fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0) > fire->stressTol) {
        return false;
    }
    
    fire->rtype = 0;
    return true;
}

ConstVolStressFIRE *addConstVolStressFIRE(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getConstVolStressFIRE(update))
        Abort("repetitive addConstVolStressFire(...)");
    
    ConstVolStressFIRE *fire = (ConstVolStressFIRE *)calloc(sizeof(ConstVolStressFIRE), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__ConstVolStressFireRelax__");
    
    cmdArg *cmd = findVariable(var, "csmin");
    if(!cmd) Abort("--csmin tStress xy|xz|yz maxDeltaStrain or --csmin");
    if(cmd->cmdArgc == 0){
        fire->isInit = false;
        return fire;
    }
    
    if (cmd->cmdArgc != 3) {
        Abort("--csmin tStress xy|xz|yz maxDeltaStrain");
    }
    
    fire->stressTol = __Stol__;
    fire->dtSet = 5E-3;
    
    if (strcmp(cmd->cmdArgv[1], "xy") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 1;
    }
    #if (DIM == 3)
    else if (strcmp(cmd->cmdArgv[1], "xz") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 2;
    } else if (strcmp(cmd->cmdArgv[1], "yz") == 0) {
        fire->shearDim = 1;
        fire->gradDim = 2;
    }
    #endif
    else {
        Abort("No --csmin tStress xy|xz|yz maxDeltaStrain");
    }
    
    fire->stressSet = atof(cmd->cmdArgv[0]);
    fire->maxDeltaStrain = fabs(atof(cmd->cmdArgv[2]));
    fire->deltaStrain = 0;
    
    fire->isInit = true;
    return fire;
}
int reInitConstVolStressFIRE(ConstVolStressFIRE *fire, double tStress, char *sType, double maxDeltaStrain) {
    memset(fire, '\0', sizeof(ConstVolStressFIRE));
    
    fire->stressTol = __Stol__;
    fire->dtSet = 5E-3;
    
    if (strcmp(sType, "xy") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 1;
    }
#if (DIM == 3)
    else if (strcmp(sType, "xz") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 2;
    } else if (strcmp(sType, "yz") == 0) {
        fire->shearDim = 1;
        fire->gradDim = 2;
    }
#endif
    else {
        safeFprintf(stderr, "--csmin tStress xy|xz|yz maxDeltaStrain");
        return -1;
    }
    fire->stressSet = tStress;
    fire->maxDeltaStrain = maxDeltaStrain;
    fire->deltaStrain = 0;
        
    fire->isInit = true;
    return 1;
}
ConstVolStressFIRE *getConstVolStressFIRE(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__ConstVolStressFireRelax__");
    if (whichTool < 0)
        return NULL;
    return (ConstVolStressFIRE *)update->toolkit.toolkit[whichTool];
}
int minConstVolStressFIRE(Box *box, Particle *particle, Update *update) {
    ConstVolStressFIRE *fire = getConstVolStressFIRE(update);
    if (!fire || !fire->isInit) Abort("call addConstVolStressFire(...)");
    
    box->isShapeFixed = false;
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
    double stressSet = fire->stressSet;
    fire->stressSet *= update->pressureUnits;//incorporationg units
    double ePairLastNeg = 0;
    int rVal = 0;
    
    calcForce(box, particle, update);
    cvsFire_DotProduct(box, particle, update, fire);
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    fire->gVeloc = 0;
    ePairLastNeg = update->ePair;
    
    bool isConverge = cvsFire_CheckConverge(box, particle, update, fire);
    while (!isConverge) {
        // update x
        cvsFire_Integrate(box, particle, update, fire);
        
        // update v
        calcForce(box, particle, update);
        cvsFire_DotProduct(box, particle, update, fire);
        isConverge = cvsFire_CheckConverge(box, particle, update, fire);
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
            fire->gVeloc = fire->gVeloc * scale1 + fire->gForce * scale2;
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
            fire->gVeloc = 0;
        }
        
        if (currStep % 1000 == 0) {
            safeFprintf(stdout, "%d %g %g %g %.15g\n", currStep, fire->deltaStrain, update->ePair / update->energyUnits,
                        update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / update->pressureUnits,
                        fire->aveForce / update->forceUnits);
        }
    }
    
    {
        safeFprintf(stdout, "%d %g %g %g %.15g\n", currStep, fire->deltaStrain, update->ePair / update->energyUnits,
                    update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / update->pressureUnits,
                    fire->aveForce / update->forceUnits);
    }

    
    fire->stressSet = stressSet;
    fire->isInit = false;
    return rVal;
}
int delConstVolStressFIRE(Update *update) {
    ConstVolStressFIRE *fire = getConstVolStressFIRE(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__ConstVolStressFireRelax__");
    return 0;
}
#endif

#endif

#endif
