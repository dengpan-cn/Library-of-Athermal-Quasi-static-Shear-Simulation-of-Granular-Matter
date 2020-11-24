#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 constPressIsoShearNL.c -o cpisnl
// run: ./cpisnl --rf rf.bin --shear 1E-4 10.0 --dump --cpmin 1E-4 --sf suffix
// The format of binary file is described in function "read_data" and
// "write_data"
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

  sShear->cinfo = initContactInfo(box, particle, thermo, update, var);
  sShear->dinfo = initDump(box, particle, thermo, update, var);

  constPressIsoFireRelax(box, particle, thermo, update, var);
  computeCoordination(box, particle, thermo, update, var);
  dump(box, particle, thermo, update, var);

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
