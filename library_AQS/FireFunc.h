#ifndef __FIREFUNC__
#define __FIREFUNC__

#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

void cbsFire_init(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichVar = findVariable(var, "cvmin");
  if (whichVar < 0) Abort("No --cvmin cmd");
  if (update->cbsFire == NULL) {
    update->cbsFire =
        (minConstBoxShapeFIRE *)calloc(sizeof(minConstBoxShapeFIRE), 1);
  }
  minConstBoxShapeFIRE *fire = update->cbsFire;
  if (fire->isInit) return;

  fire->dtSet = 5E-3;

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  fire->isInit = true;
}
void cbsFire_Integrate(Box *box, Particle *particle, Thermo *thermo,
                       Update *update, Variable *var) {
  minConstBoxShapeFIRE *fire = update->cbsFire;

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 pos = particle->xyz[iatom];
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    vecScaleAdd(veloc, veloc, fire->dt, force);
    vecScaleAdd(pos, pos, fire->dt, veloc);

    particle->xyz[iatom] = pos;
    particle->veloc[iatom] = veloc;
  }

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
}
void cbsFire_DotProduct(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  minConstBoxShapeFIRE *fire = update->cbsFire;

  double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
  double dotPrd = 0.0;
  //__float128 sumForce = 0;
  double sumForce = 0;  //?
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    vecDot(dotPrd, veloc, force);
    vdotf += dotPrd;
    vecDot(dotPrd, veloc, veloc);
    vdotv += dotPrd;
    vecDot(dotPrd, force, force);
    fdotf += dotPrd;
    sumForce += sqrt(dotPrd);
  }
  fire->vdotf = vdotf;
  fire->vdotv = vdotv;
  fire->fdotf = fdotf;
  fire->aveForce = sumForce / thermo->forceUnits / particle->nAtom;
}
void constBoxShapeFireRelax(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  // const-Box-shape FIRE: minimize U(r)
  // The shape of parallelepiped is fixed!
  // The relaxation stops if the system unjammed or the minimum reached.
  // unjamming criteria: pair potential per particle is less than __ZeroEnergy__
  // minimum criteria: average force amplitude is less than __AveZeroForce__

  cbsFire_init(box, particle, thermo, update, var);
  minConstBoxShapeFIRE *fire = update->cbsFire;

  memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
  double dtmax = 0.5 * thermo->timeUnits;
  int currStep = 0, last_negative = 0;
  double alpha = ALPHA0;
  fire->dt = fire->dtSet * thermo->timeUnits;
  double dt0 = fire->dt;

  calcForce(box, particle, thermo, update, var);
  cbsFire_DotProduct(box, particle, thermo, update, var);
  while (!(fire->aveForce <= __AveZeroForce__)) {
    // update x
    cbsFire_Integrate(box, particle, thermo, update, var);

    // update v
    calcForce(box, particle, thermo, update, var);
    cbsFire_DotProduct(box, particle, thermo, update, var);
    currStep++;

    // mixing v and update timestep
    if (fire->vdotf > 0) {
      if (currStep - last_negative > DELAYSTEP) {
        fire->dt *= DT_GROW;
        fire->dt = (fire->dt > dtmax ? dtmax : fire->dt);
        alpha *= ALPHA_SHRINK;
      }

      double scale1 = 1.0 - alpha;
      double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double3 veloc = particle->veloc[iatom];
        double3 force = particle->force[iatom];
        vecScale(veloc, scale1, veloc);
        vecScaleAdd(veloc, veloc, scale2, force);
        particle->veloc[iatom] = veloc;
      }
    } else {
      int deltaStep = currStep - last_negative;
      if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 0.9;
        dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * thermo->timeUnits);
      } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 1.1;
        dtmax = cpuMin(dtmax, 0.5 * thermo->timeUnits);
      }

      fire->dt *= DT_SHRINK;
      dt0 = fire->dt;
      last_negative = currStep;
      alpha = ALPHA0;
      memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
    }

    if (thermo->Epair / thermo->energyUnits < __ZeroEnergy__) {
      return;
    }
  }
}

void cbvFire_init(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichVar = findVariable(var, "cvmin");
  if (whichVar < 0) Abort("No --cvmin cmd");
  if (update->cbvFire == NULL) {
    update->cbvFire =
        (minConstBoxVolFIRE *)calloc(sizeof(minConstBoxVolFIRE), 1);
  }
  minConstBoxVolFIRE *fire = update->cbvFire;
  if (fire->isInit) return;
#ifndef __triBox__
  Abort("BOX error!");
#endif

  fire->dtSet = 2E-3;
  box->isShapeFixed = false;
  particle->isSizeFixed = true;
  fire->isInit = true;
}
void cbvFire_Integrate(Box *box, Particle *particle, Thermo *thermo,
                       Update *update, Variable *var) {
  minConstBoxVolFIRE *fire = update->cbvFire;
  double relativeP, dtFact;
  double3 dtAniso, dtStress;

  //==============box=============
  relativeP = fabs(thermo->ptensor.h0 / thermo->pressure - 1.0);
  dtFact = sqrt(particle->nAtom) * exp(-pow(relativeP / __Ptol__, 2)) + 10.0;
  dtAniso.x = fire->dt / dtFact;
  fire->gVelocAniso.x += dtAniso.x * fire->gForceAniso.x;
  fire->sfactAniso.x = exp(dtAniso.x * fire->gVelocAniso.x);

  relativeP = fabs(thermo->ptensor.h1 / thermo->pressure - 1.0);
  dtFact = sqrt(particle->nAtom) * exp(-pow(relativeP / __Ptol__, 2)) + 10.0;
  dtAniso.y = fire->dt / dtFact;
  fire->gVelocAniso.y += dtAniso.y * fire->gForceAniso.y;
  fire->sfactAniso.y = exp(dtAniso.y * fire->gVelocAniso.y);

  fire->gVelocAniso.z = 0.0;
  fire->sfactAniso.z = 1.0 / (fire->sfactAniso.x * fire->sfactAniso.y);

  //======================
  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h3, 2)) + 10.0;
  dtStress.x = fire->dt / dtFact;
  fire->gVelocStress.x += dtStress.x * fire->gForceStress.x;
  fire->sfactStress.x = dtStress.x * fire->gVelocStress.x;

  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h4, 2)) + 10.0;
  dtStress.y = fire->dt / dtFact;
  fire->gVelocStress.y += dtStress.y * fire->gForceStress.y;
  fire->sfactStress.y = dtStress.y * fire->gVelocStress.y;

  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h5, 2)) + 10.0;
  dtStress.z = fire->dt / dtFact;
  fire->gVelocStress.z += dtStress.z * fire->gForceStress.z;
  fire->sfactStress.z = dtStress.z * fire->gVelocStress.z;
  //=============update box=========
  Hvoigt6 sfact;
  sfact.h0 = fire->sfactAniso.x;
  sfact.h1 = fire->sfactAniso.y;
  sfact.h2 = fire->sfactAniso.z;
  sfact.h3 = fire->sfactStress.x;
  sfact.h4 = fire->sfactStress.y;
  sfact.h5 = fire->sfactStress.z;

  Hvoigt6 boxHvoigt = box->boxHvoigt;
  boxHvoigt.h0 *= sfact.h0;

  boxHvoigt.h5 *= sfact.h0;
  boxHvoigt.h5 += sfact.h5 * boxHvoigt.h1;
  boxHvoigt.h1 *= sfact.h1;

  boxHvoigt.h4 *= sfact.h0;
  boxHvoigt.h4 += sfact.h5 * boxHvoigt.h3 + sfact.h4 * boxHvoigt.h2;
  boxHvoigt.h3 *= sfact.h1;
  boxHvoigt.h3 += sfact.h3 * boxHvoigt.h2;
  boxHvoigt.h2 *= sfact.h2;

  box->boxHvoigt = boxHvoigt;
  thermo->volFrac /= sfact.h0 * sfact.h1 * sfact.h2;
  setBoxPara(box, particle, thermo, update, var);

  setBasicInfo(box, particle, thermo, update, var);

  //=============update particle=========
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 pos = particle->xyz[iatom];
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    pos.x *= sfact.h0;
    pos.x += sfact.h5 * pos.y + sfact.h4 * pos.z;
    pos.y *= sfact.h1;
    pos.y += sfact.h3 * pos.z;
    pos.z *= sfact.h2;

    vecScaleAdd(veloc, veloc, fire->dt, force);
    vecScaleAdd(pos, pos, fire->dt, veloc);

    particle->xyz[iatom] = pos;
    particle->veloc[iatom] = veloc;
  }

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
}
void cbvFire_DotProduct(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  minConstBoxVolFIRE *fire = update->cbvFire;

  fire->gForceAniso.x = (thermo->ptensor.h0 - thermo->pressure) * box->volume;
  fire->gForceAniso.y = (thermo->ptensor.h1 - thermo->pressure) * box->volume;
  fire->gForceAniso.z = 0.0;

  fire->gForceStress.x = thermo->ptensor.h3 * box->volume;
  fire->gForceStress.y = thermo->ptensor.h4 * box->volume;
  fire->gForceStress.z = thermo->ptensor.h5 * box->volume;
  double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
  double dotPrd = 0.0;
  double sumForce = 0.0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    vecDot(dotPrd, veloc, force);
    vdotf += dotPrd;
    vecDot(dotPrd, veloc, veloc);
    vdotv += dotPrd;
    vecDot(dotPrd, force, force);
    fdotf += dotPrd;
    sumForce += sqrt(dotPrd);
  }
  fire->vdotf = vdotf;
  fire->vdotv = vdotv;
  fire->fdotf = fdotf;

  double dotPrdAniso = 0, dotPrdStress = 0;
  vecDot(dotPrdAniso, fire->gVelocAniso, fire->gForceAniso);
  vecDot(dotPrdStress, fire->gVelocStress, fire->gForceStress);
  fire->vdotf += dotPrdAniso + dotPrdStress;

  vecDot(dotPrdAniso, fire->gVelocAniso, fire->gVelocAniso);
  vecDot(dotPrdStress, fire->gVelocStress, fire->gVelocStress);
  fire->vdotv += dotPrdAniso + dotPrdStress;

  vecDot(dotPrdAniso, fire->gForceAniso, fire->gForceAniso);
  vecDot(dotPrdStress, fire->gForceStress, fire->gForceStress);
  fire->fdotf += dotPrdAniso + dotPrdStress;

  fire->aveForce = sumForce / thermo->forceUnits / particle->nAtom;
}
bool cbvFire_CheckConverge(Box *box, Particle *particle, Thermo *thermo,
                           Update *update, Variable *var) {
  minConstBoxVolFIRE *fire = update->cbvFire;
  double3 relativeP;
  relativeP.x = fabs(thermo->ptensor.h0 / thermo->pressure - 1.0);
  relativeP.y = fabs(thermo->ptensor.h1 / thermo->pressure - 1.0);
  relativeP.z = 0.0;
  double3 stress;
  stress.x = fabs(thermo->ptensor.h3 / thermo->pressureUnits);
  stress.y = fabs(thermo->ptensor.h4 / thermo->pressureUnits);
  stress.z = fabs(thermo->ptensor.h5 / thermo->pressureUnits);

  if (relativeP.x <= __Ptol__ && relativeP.y <= __Ptol__ &&
      stress.x <= __StressTol__ && stress.y <= __StressTol__ &&
      stress.z <= __StressTol__ && fire->aveForce <= __AveZeroForce__) {
    return true;
  }

  return false;
}
void constBoxVolFireRelax(Box *box, Particle *particle, Thermo *thermo,
                          Update *update, Variable *var) {
  // const-Box-volume FIRE: minimize U(r)
  // The volume of parallelepiped is fixed!
  // The relaxation stops if the system unjammed or the minimum reached.
  // unjamming criteria: pair potential per particle is less than __ZeroEnergy__
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and the relative difference among (Pxx, Pyy, Pzz) are less than __Ptol__
  // and the absolute value of (Pxy, Pxz, Pxy) are less than __StressTol__
  cbvFire_init(box, particle, thermo, update, var);
  minConstBoxVolFIRE *fire = update->cbvFire;

  memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
  double dtmax = 0.5 * thermo->timeUnits;
  int currStep = 0, last_negative = 0;
  double alpha = ALPHA0;

  fire->dt = fire->dtSet * thermo->timeUnits;
  double dt0 = fire->dt;
  calcForce(box, particle, thermo, update, var);
  cbvFire_DotProduct(box, particle, thermo, update, var);
  bool isConverge = cbvFire_CheckConverge(box, particle, thermo, update, var);
  while (!isConverge) {
    // update x and box
    cbvFire_Integrate(box, particle, thermo, update, var);

    // check timestep
    calcForce(box, particle, thermo, update, var);
    cbvFire_DotProduct(box, particle, thermo, update, var);
    isConverge = cbvFire_CheckConverge(box, particle, thermo, update, var);
    currStep++;

    // update velocity
    if (fire->vdotf > 0) {
      if (currStep - last_negative > DELAYSTEP) {
        fire->dt *= DT_GROW;
        fire->dt = (fire->dt > dtmax ? dtmax : fire->dt);
        alpha *= ALPHA_SHRINK;
      }

      double scale1 = 1.0 - alpha;
      double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double3 veloc = particle->veloc[iatom];
        double3 force = particle->force[iatom];
        vecScale(veloc, scale1, veloc);
        vecScaleAdd(veloc, veloc, scale2, force);
        particle->veloc[iatom] = veloc;
      }
      vecScale(fire->gVelocAniso, scale1, fire->gVelocAniso);
      vecScaleAdd(fire->gVelocAniso, fire->gVelocAniso, scale2,
                  fire->gForceAniso);
      vecScale(fire->gVelocStress, scale1, fire->gVelocStress);
      vecScaleAdd(fire->gVelocStress, fire->gVelocStress, scale2,
                  fire->gForceStress);
    } else {
      int deltaStep = currStep - last_negative;
      if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 0.9;
        dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * thermo->timeUnits);
      } else if (deltaStep > 2.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 1.1;
        dtmax = cpuMin(dtmax, 0.5 * thermo->timeUnits);
      }

      fire->dt *= DT_SHRINK;
      last_negative = currStep;
      alpha = ALPHA0;
      memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
      fire->gVelocAniso = make_double3(0, 0, 0);
      fire->gVelocStress = make_double3(0, 0, 0);
    }

    if (thermo->Epair / thermo->energyUnits < __ZeroEnergy__) {
      return;
    }
  }
}

//=================const press========
void cpiFire_init(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichVar = findVariable(var, "cpmin");
  if (whichVar < 0) Abort("No --cpmin cmd");
  if (update->cpiFire == NULL) {
    update->cpiFire =
        (minConstPressIsoFIRE *)calloc(sizeof(minConstPressIsoFIRE), 1);
  }
  minConstPressIsoFIRE *fire = update->cpiFire;
  if (fire->isInit) return;
  if (var->cmd[whichVar].cmdArgc != 1) {
    Abort("--cpmin pTarget");
  }
  fire->Pset = atof(var->cmd[whichVar].cmdArgv[0]);
  fire->ptol = __Ptol__;
  if (fire->Pset < 1E-8) {
    Abort("The Pset is too small!");
  }

  fire->dtSet = 5E-3;

  box->isShapeFixed = true;
  particle->isSizeFixed = false;
  fire->isInit = true;
}
void cpiFire_Integrate(Box *box, Particle *particle, Thermo *thermo,
                       Update *update, Variable *var) {
  minConstPressIsoFIRE *fire = update->cpiFire;

  //==============volume=============
  double relativeP =
      fabs(thermo->pressure / (fire->Pset * thermo->pressureUnits) - 1.0);
  double dtFact =
      sqrt(particle->nAtom) * exp(-pow(relativeP / fire->ptol, 2)) + 10.0;
  double box_dt = fire->dt / dtFact;
  fire->gVeloc += box_dt * fire->gForce;
  fire->sfact = exp(box_dt * fire->gVeloc);
  double phi = thermo->volFrac / pow(fire->sfact, 3.0);
  if (thermo->volFrac - phi > maxDeltaVF) {
    fire->sfact =
        pow(thermo->volFrac / (thermo->volFrac - maxDeltaVF), 1.0 / 3.0);
    fire->gVeloc = log(fire->sfact) / box_dt;
  } else if (phi - thermo->volFrac > maxDeltaVF) {
    fire->sfact =
        pow(thermo->volFrac / (thermo->volFrac + maxDeltaVF), 1.0 / 3.0);
    fire->gVeloc = log(fire->sfact) / box_dt;
  }
  thermo->volFrac /= pow(fire->sfact, 3);
  particle->meanDiameter /= fire->sfact;
  setBasicInfo(box, particle, thermo, update, var);

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 pos = particle->xyz[iatom];
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    vecScaleAdd(veloc, veloc, fire->dt, force);
    vecScaleAdd(pos, pos, fire->dt, veloc);

    particle->xyz[iatom] = pos;
    particle->veloc[iatom] = veloc;
  }

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
}
void cpiFire_DotProduct(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  minConstPressIsoFIRE *fire = update->cpiFire;

  fire->gForce =
      (thermo->pressure - fire->Pset * thermo->pressureUnits) * box->volume;
  double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
  double dotPrd = 0.0;
  double sumForce = 0.0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];
    vecDot(dotPrd, veloc, force);
    vdotf += dotPrd;
    vecDot(dotPrd, veloc, veloc);
    vdotv += dotPrd;
    vecDot(dotPrd, force, force);
    fdotf += dotPrd;
    sumForce += sqrt(dotPrd);
  }
  fire->vdotf = vdotf + fire->gVeloc * fire->gForce;
  fire->vdotv = vdotv + fire->gVeloc * fire->gVeloc;
  fire->fdotf = fdotf + fire->gForce * fire->gForce;
  fire->aveForce = sumForce / thermo->forceUnits / particle->nAtom;
}
void constPressIsoFireRelax(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  // const-iso-Pressure FIRE: minimize H(r,h) = U(r) + Pset * trace(h).
  // The mean-diameter of particles is adjusted and the box shape is keeped.
  // The relaxation stops if the minimum of H(r,h) reached.
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and the absolute value of (P / Pset - 1.0) is less than __Ptol__

  cpiFire_init(box, particle, thermo, update, var);
  minConstPressIsoFIRE *fire = update->cpiFire;

  memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
  double dtmax = 0.5 * thermo->timeUnits;
  int currStep = 0, last_negative = 0;
  double alpha = ALPHA0;

  fire->dt = fire->dtSet * thermo->timeUnits;
  double dt0 = fire->dt;
  calcForce(box, particle, thermo, update, var);
  cpiFire_DotProduct(box, particle, thermo, update, var);
  double relativeP =
      fabs(thermo->pressure / (fire->Pset * thermo->pressureUnits) - 1.0);
  while (!(relativeP <= fire->ptol && fire->aveForce <= __AveZeroForce__)) {
    // update x and box
    cpiFire_Integrate(box, particle, thermo, update, var);

    // calc force
    calcForce(box, particle, thermo, update, var);
    cpiFire_DotProduct(box, particle, thermo, update, var);
    relativeP =
        fabs(thermo->pressure / (fire->Pset * thermo->pressureUnits) - 1.0);
    currStep++;

    // mixing v and update timestep
    if (fire->vdotf > 0) {
      if (currStep - last_negative > DELAYSTEP) {
        fire->dt *= DT_GROW;
        fire->dt = (fire->dt > dtmax ? dtmax : fire->dt);
        alpha *= ALPHA_SHRINK;
      }

      double scale1 = 1.0 - alpha;
      double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double3 veloc = particle->veloc[iatom];
        double3 force = particle->force[iatom];
        vecScale(veloc, scale1, veloc);
        vecScaleAdd(veloc, veloc, scale2, force);
        particle->veloc[iatom] = veloc;
      }
      fire->gVeloc = fire->gVeloc * scale1 + fire->gForce * scale2;
    } else {
      int deltaStep = currStep - last_negative;
      if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 0.9;
        dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * thermo->timeUnits);
      } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 1.1;
        dtmax = cpuMin(dtmax, 0.5 * thermo->timeUnits);
      }

      fire->dt *= DT_SHRINK;
      last_negative = currStep;
      alpha = ALPHA0;
      memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
      fire->gVeloc = 0;
    }
  }
}

void cptFire_init(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichVar = findVariable(var, "cpmin");
  if (whichVar < 0) Abort("No --cpmin cmd");
  if (update->cptFire == NULL) {
    update->cptFire =
        (minConstPressTriFIRE *)calloc(sizeof(minConstPressTriFIRE), 1);
  }
  minConstPressTriFIRE *fire = update->cptFire;
  if (fire->isInit) return;
#ifndef __triBox__
  Abort("BOX error!");
#endif
  if (var->cmd[whichVar].cmdArgc != 1) {
    Abort("--cpmin pTarget");
  }
  fire->Pset = atof(var->cmd[whichVar].cmdArgv[0]);
  fire->ptol = __Ptol__;
  if (fire->Pset < 1E-8) {
    Abort("The Pset is too small!");
  }

  fire->dtSet = 5E-3;

  box->isShapeFixed = false;
  particle->isSizeFixed = true;
  fire->isInit = true;
}
void cptFire_Integrate(Box *box, Particle *particle, Thermo *thermo,
                       Update *update, Variable *var) {
  minConstPressTriFIRE *fire = update->cptFire;
  double relativeP, dtFact, dtIso;
  double3 dtAniso, dtStress;

  //==============box=============
  relativeP =
      fabs(thermo->ptensor.h0 / (fire->Pset * thermo->pressureUnits) - 1.0);
  dtFact = sqrt(particle->nAtom) * exp(-pow(relativeP / fire->ptol, 2)) + 10.0;
  dtAniso.x = fire->dt / dtFact;
  fire->gVelocAniso.x += dtAniso.x * fire->gForceAniso.x;
  fire->sfactAniso.x = exp(dtAniso.x * fire->gVelocAniso.x);

  relativeP =
      fabs(thermo->ptensor.h1 / (fire->Pset * thermo->pressureUnits) - 1.0);
  dtFact = sqrt(particle->nAtom) * exp(-pow(relativeP / fire->ptol, 2)) + 10.0;
  dtAniso.y = fire->dt / dtFact;
  fire->gVelocAniso.y += dtAniso.y * fire->gForceAniso.y;
  fire->sfactAniso.y = exp(dtAniso.y * fire->gVelocAniso.y);

  relativeP =
      fabs(thermo->ptensor.h2 / (fire->Pset * thermo->pressureUnits) - 1.0);
  dtFact = sqrt(particle->nAtom) * exp(-pow(relativeP / fire->ptol, 2)) + 10.0;
  dtAniso.z = fire->dt / dtFact;
  fire->gVelocAniso.z += dtAniso.z * fire->gForceAniso.z;
  fire->sfactAniso.z = exp(dtAniso.z * fire->gVelocAniso.z);

  // relativeP = fabs(thermo->ptensor.h3) / (fire->Pset *
  // thermo->pressureUnits);
  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h3, 2)) + 10.0;
  dtStress.x = fire->dt / dtFact;
  fire->gVelocStress.x += dtStress.x * fire->gForceStress.x;
  fire->sfactStress.x = dtStress.x * fire->gVelocStress.x;

  // relativeP = fabs(thermo->ptensor.h4) / (fire->Pset *
  // thermo->pressureUnits);
  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h4, 2)) + 10.0;
  dtStress.y = fire->dt / dtFact;
  fire->gVelocStress.y += dtStress.y * fire->gForceStress.y;
  fire->sfactStress.y = dtStress.y * fire->gVelocStress.y;

  // relativeP = fabs(thermo->ptensor.h5) / (fire->Pset *
  // thermo->pressureUnits);
  dtFact = sqrt(particle->nAtom) * exp(-pow(thermo->ptensor.h5, 2)) + 10.0;
  dtStress.z = fire->dt / dtFact;
  fire->gVelocStress.z += dtStress.z * fire->gForceStress.z;
  fire->sfactStress.z = dtStress.z * fire->gVelocStress.z;
  //=============update box=========
  double phi = thermo->volFrac /
               (fire->sfactAniso.x * fire->sfactAniso.y * fire->sfactAniso.z);
  if (thermo->volFrac - phi > maxDeltaVF) {
    double sfact = pow(phi / (thermo->volFrac - maxDeltaVF), 1.0 / 3.0);
    vecScale(fire->sfactAniso, sfact, fire->sfactAniso);
    fire->gVelocAniso.x = log(fire->sfactAniso.x) / dtAniso.x;
    fire->gVelocAniso.y = log(fire->sfactAniso.y) / dtAniso.y;
    fire->gVelocAniso.z = log(fire->sfactAniso.z) / dtAniso.z;
  } else if (phi - thermo->volFrac > maxDeltaVF) {
    double sfact = pow(phi / (thermo->volFrac - maxDeltaVF), 1.0 / 3.0);
    vecScale(fire->sfactAniso, sfact, fire->sfactAniso);
    fire->gVelocAniso.x = log(fire->sfactAniso.x) / dtAniso.x;
    fire->gVelocAniso.y = log(fire->sfactAniso.y) / dtAniso.y;
    fire->gVelocAniso.z = log(fire->sfactAniso.z) / dtAniso.z;
  }

  Hvoigt6 sfact;
  sfact.h0 = fire->sfactAniso.x;
  sfact.h1 = fire->sfactAniso.y;
  sfact.h2 = fire->sfactAniso.z;
  sfact.h3 = fire->sfactStress.x;
  sfact.h4 = fire->sfactStress.y;
  sfact.h5 = fire->sfactStress.z;

  Hvoigt6 boxHvoigt = box->boxHvoigt;
  boxHvoigt.h0 *= sfact.h0;

  boxHvoigt.h5 *= sfact.h0;
  boxHvoigt.h5 += sfact.h5 * boxHvoigt.h1;
  boxHvoigt.h1 *= sfact.h1;

  boxHvoigt.h4 *= sfact.h0;
  boxHvoigt.h4 += sfact.h5 * boxHvoigt.h3 + sfact.h4 * boxHvoigt.h2;
  boxHvoigt.h3 *= sfact.h1;
  boxHvoigt.h3 += sfact.h3 * boxHvoigt.h2;
  boxHvoigt.h2 *= sfact.h2;

  box->boxHvoigt = boxHvoigt;
  thermo->volFrac /= sfact.h0 * sfact.h1 * sfact.h2;
  setBoxPara(box, particle, thermo, update, var);

  setBasicInfo(box, particle, thermo, update, var);

  //=============update particle=========
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 pos = particle->xyz[iatom];
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    pos.x *= sfact.h0;
    pos.x += sfact.h5 * pos.y + sfact.h4 * pos.z;
    pos.y *= sfact.h1;
    pos.y += sfact.h3 * pos.z;
    pos.z *= sfact.h2;

    vecScaleAdd(veloc, veloc, fire->dt, force);
    vecScaleAdd(pos, pos, fire->dt, veloc);

    particle->xyz[iatom] = pos;
    particle->veloc[iatom] = veloc;
  }

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
}
void cptFire_DotProduct(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  minConstPressTriFIRE *fire = update->cptFire;

  fire->gForceAniso.x =
      (thermo->ptensor.h0 - fire->Pset * thermo->pressureUnits) * box->volume;
  fire->gForceAniso.y =
      (thermo->ptensor.h1 - fire->Pset * thermo->pressureUnits) * box->volume;
  fire->gForceAniso.z =
      (thermo->ptensor.h2 - fire->Pset * thermo->pressureUnits) * box->volume;
  fire->gForceStress.x = thermo->ptensor.h3 * box->volume;
  fire->gForceStress.y = thermo->ptensor.h4 * box->volume;
  fire->gForceStress.z = thermo->ptensor.h5 * box->volume;
  double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
  double dotPrd = 0.0;
  double sumForce = 0.0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 veloc = particle->veloc[iatom];
    double3 force = particle->force[iatom];

    vecDot(dotPrd, veloc, force);
    vdotf += dotPrd;
    vecDot(dotPrd, veloc, veloc);
    vdotv += dotPrd;
    vecDot(dotPrd, force, force);
    fdotf += dotPrd;
    sumForce += sqrt(dotPrd);
  }
  fire->vdotf = vdotf;
  fire->vdotv = vdotv;
  fire->fdotf = fdotf;

  double dotPrdAniso = 0, dotPrdStress = 0;
  vecDot(dotPrdAniso, fire->gVelocAniso, fire->gForceAniso);
  vecDot(dotPrdStress, fire->gVelocStress, fire->gForceStress);
  fire->vdotf += dotPrdAniso + dotPrdStress;

  vecDot(dotPrdAniso, fire->gVelocAniso, fire->gVelocAniso);
  vecDot(dotPrdStress, fire->gVelocStress, fire->gVelocStress);
  fire->vdotv += dotPrdAniso + dotPrdStress;

  vecDot(dotPrdAniso, fire->gForceAniso, fire->gForceAniso);
  vecDot(dotPrdStress, fire->gForceStress, fire->gForceStress);
  fire->fdotf += dotPrdAniso + dotPrdStress;

  fire->aveForce = sumForce / thermo->forceUnits / particle->nAtom;
}
bool cptFire_CheckConverge(Box *box, Particle *particle, Thermo *thermo,
                           Update *update, Variable *var) {
  minConstPressTriFIRE *fire = update->cptFire;
  double3 relativeP;
  relativeP.x =
      fabs(thermo->ptensor.h0 / fire->Pset / thermo->pressureUnits - 1.0);
  relativeP.y =
      fabs(thermo->ptensor.h1 / fire->Pset / thermo->pressureUnits - 1.0);
  relativeP.z =
      fabs(thermo->ptensor.h2 / fire->Pset / thermo->pressureUnits - 1.0);
  double3 stress;
  stress.x = fabs(thermo->ptensor.h3 / thermo->pressureUnits);
  stress.y = fabs(thermo->ptensor.h4 / thermo->pressureUnits);
  stress.z = fabs(thermo->ptensor.h5 / thermo->pressureUnits);
  if (relativeP.x <= fire->ptol && relativeP.y <= fire->ptol &&
      relativeP.z <= fire->ptol && stress.x <= __StressTol__ &&
      stress.y <= __StressTol__ && stress.z <= __StressTol__ &&
      fire->aveForce <= __AveZeroForce__) {
    return true;
  }
  return false;
}
void constPressTriFireRelax(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  // const-tri-Pressure-tensor FIRE: minimize H(r,h) = U(r) + Pset * trace(h).
  // All six parameters of box are adjusted.
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and relative difference between Pxx,Pyy,Pzz and Pset are less than __Ptol__
  // and the absolute value of (Pxy, Pxz, Pxy) are less than __StressTol__

  cptFire_init(box, particle, thermo, update, var);
  minConstPressTriFIRE *fire = update->cptFire;

  memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
  double dtmax = 0.5 * thermo->timeUnits;
  int currStep = 0, last_negative = 0;
  double alpha = ALPHA0;

  fire->dt = fire->dtSet * thermo->timeUnits;
  double dt0 = fire->dt;
  calcForce(box, particle, thermo, update, var);
  cptFire_DotProduct(box, particle, thermo, update, var);
  bool isConverge = cptFire_CheckConverge(box, particle, thermo, update, var);
  while (!isConverge) {
    // update x and box
    cptFire_Integrate(box, particle, thermo, update, var);

    // check timestep
    calcForce(box, particle, thermo, update, var);
    cptFire_DotProduct(box, particle, thermo, update, var);
    isConverge = cptFire_CheckConverge(box, particle, thermo, update, var);
    currStep++;

    // update velocity
    if (fire->vdotf > 0) {
      if (currStep - last_negative > DELAYSTEP) {
        fire->dt *= DT_GROW;
        fire->dt = (fire->dt > dtmax ? dtmax : fire->dt);
        alpha *= ALPHA_SHRINK;
      }

      double scale1 = 1.0 - alpha;
      double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double3 veloc = particle->veloc[iatom];
        double3 force = particle->force[iatom];
        vecScale(veloc, scale1, veloc);
        vecScaleAdd(veloc, veloc, scale2, force);
        particle->veloc[iatom] = veloc;
      }
      vecScale(fire->gVelocAniso, scale1, fire->gVelocAniso);
      vecScaleAdd(fire->gVelocAniso, fire->gVelocAniso, scale2,
                  fire->gForceAniso);
      vecScale(fire->gVelocStress, scale1, fire->gVelocStress);
      vecScaleAdd(fire->gVelocStress, fire->gVelocStress, scale2,
                  fire->gForceStress);
    } else {
      int deltaStep = currStep - last_negative;
      if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 0.9;
        dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * thermo->timeUnits);
      } else if (deltaStep > 2.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
        dtmax = dtmax * 1.1;
        dtmax = cpuMin(dtmax, 0.5 * thermo->timeUnits);
      }

      fire->dt *= DT_SHRINK;
      last_negative = currStep;
      alpha = ALPHA0;
      memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
      fire->gVelocAniso = make_double3(0, 0, 0);
      fire->gVelocStress = make_double3(0, 0, 0);
    }
  }
}

#endif