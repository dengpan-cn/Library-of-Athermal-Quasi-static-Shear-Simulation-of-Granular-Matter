#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 constPressIsoShearNL.c -o cpisnl
// run: ./cpisnl --rf rf.bin --shear 1E-4 10.0 --dump --cpmin 1E-4 --sf suffix
void constPressIsoShearInit(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  if (update->sShear) return;

  int whichVar = findVariable(var, "shear");
  if (whichVar < 0) Abort("--shear 1E-4 1.0");
  if (var->cmd[whichVar].cmdArgc != 2) Abort("--shear 1E-4 1.0");
  if (!update->sShear)
    update->sShear = (simpleShear *)calloc(1, sizeof(simpleShear));

  simpleShear *sShear = update->sShear;
  sShear->gamma = 0;
  sShear->deltaGamma = atof(var->cmd[whichVar].cmdArgv[0]);
  sShear->maxGamma = atof(var->cmd[whichVar].cmdArgv[1]);
  if (sShear->maxGamma <= sShear->gamma) Abort("--shear 1E-4 1.0");

  char fname[4096];
  sprintf(fname, "%s/constPressShear_%%13lf_%s.bin", var->cwd, var->sf);
  sShear->shearLog = createReadWriteFile(fname);
  sShear->cinfo = initContactInfo(box, particle, thermo, update, var);
  sShear->dinfo = initDump(box, particle, thermo, update, var);

  constPressIsoFireRelax(box, particle, thermo, update, var);
  computeCoordination(box, particle, thermo, update, var);
  dump(box, particle, thermo, update, var);

  double writeInfo[13];
  writeInfo[0] = sShear->gamma;
  writeInfo[1] = thermo->volFrac;
  writeInfo[2] = thermo->Epair / thermo->energyUnits;
  writeInfo[3] = thermo->ptensor.h0 / thermo->pressureUnits;
  writeInfo[4] = thermo->ptensor.h1 / thermo->pressureUnits;
  writeInfo[5] = thermo->ptensor.h2 / thermo->pressureUnits;
  writeInfo[6] = thermo->ptensor.h3 / thermo->pressureUnits;
  writeInfo[7] = thermo->ptensor.h4 / thermo->pressureUnits;
  writeInfo[8] = thermo->ptensor.h5 / thermo->pressureUnits;
  writeInfo[9] = (double)sShear->cinfo->nRattler / (double)particle->nAtom;
  writeInfo[10] = sShear->cinfo->aveCoordNum;
  writeInfo[11] = sShear->cinfo->aveCoordNumExRattler;
  writeInfo[12] = sShear->cinfo->meanForceExRattler;
  fwrite(writeInfo, sizeof(double), 13, sShear->shearLog);
  fflush(sShear->shearLog);

  sShear->isInit = true;
}
void constPressIsoShearFinalize(Box *box, Particle *particle, Thermo *thermo,
                                Update *update, Variable *var) {
  finalizeDump(box, particle, thermo, update, var);
  finalizeContactInfo(box, particle, thermo, update, var);
  safeCloseFile(update->sShear->shearLog);
  update->sShear->cinfo = NULL;
  update->sShear->dinfo = NULL;
  safeFree(update->sShear);
}
void constPressIsoShear(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  constPressIsoShearInit(box, particle, thermo, update, var);

  simpleShear *sShear = update->sShear;
  int whichStep = 0;
  while (sShear->gamma <= sShear->maxGamma) {
    double deltaStrain = sShear->deltaGamma;
    if (whichStep <= 900) {
      double strain = pow(10.0, -8.0 + 6.0 / 900.0 * whichStep);
      deltaStrain = strain - sShear->gamma;
    }
    instant_simpShearXz(box, particle, thermo, update, var, deltaStrain);
    sShear->gamma += deltaStrain;
    whichStep++;

    constPressIsoFireRelax(box, particle, thermo, update, var);
    computeCoordination(box, particle, thermo, update, var);
    dump(box, particle, thermo, update, var);

    double writeInfo[13];
    writeInfo[0] = sShear->gamma;
    writeInfo[1] = thermo->volFrac;
    writeInfo[2] = thermo->Epair / thermo->energyUnits;
    writeInfo[3] = thermo->ptensor.h0 / thermo->pressureUnits;
    writeInfo[4] = thermo->ptensor.h1 / thermo->pressureUnits;
    writeInfo[5] = thermo->ptensor.h2 / thermo->pressureUnits;
    writeInfo[6] = thermo->ptensor.h3 / thermo->pressureUnits;
    writeInfo[7] = thermo->ptensor.h4 / thermo->pressureUnits;
    writeInfo[8] = thermo->ptensor.h5 / thermo->pressureUnits;
    writeInfo[9] = (double)sShear->cinfo->nRattler / (double)particle->nAtom;
    writeInfo[10] = sShear->cinfo->aveCoordNum;
    writeInfo[11] = sShear->cinfo->aveCoordNumExRattler;
    writeInfo[12] = sShear->cinfo->meanForceExRattler;
    fwrite(writeInfo, sizeof(double), 13, sShear->shearLog);
    fflush(sShear->shearLog);
  }

  constPressIsoShearFinalize(box, particle, thermo, update, var);
}

int main(int argc, char *argv[]) {
  Box *box = (Box *)calloc(1, sizeof(Box));
  Particle *particle = (Particle *)calloc(1, sizeof(Particle));
  Thermo *thermo = (Thermo *)calloc(1, sizeof(Thermo));
  Update *update = (Update *)calloc(1, sizeof(Update));
  Variable *var = (Variable *)calloc(1, sizeof(Variable));

  //======================Get Parameter========================
  parseCmdLine(box, particle, thermo, update, var, argc, argv);
  //===========================================================

  //=============================================================
  double tic = getClock();

  constPressIsoShear(box, particle, thermo, update, var);

  write_data(box, particle, thermo, update, var);
  double toc = getClock();
  //=============================================================

  safeFprintf(stdout, "Total Time: %g (%d:%d), nBuild: %ld, nForce: %ld.\n",
              toc - tic, shiftTimeStamp(tic), shiftTimeStamp(toc),
              update->nebrList.nBuild, update->nebrList.nForce);
}