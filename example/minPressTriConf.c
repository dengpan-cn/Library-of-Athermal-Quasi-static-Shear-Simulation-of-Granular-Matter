#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 minPressTriConf.c -o mptc
// run: ./mptc --rf rf.bin --wf wf.bin --cpmin 1E-4 --sf suffix
int main(int argc, char *argv[]) {
  Box *box = (Box *)calloc(1, sizeof(Box));
  Particle *particle = (Particle *)calloc(1, sizeof(Particle));
  Thermo *thermo = (Thermo *)calloc(1, sizeof(Thermo));
  Update *update = (Update *)calloc(1, sizeof(Update));
  Variable *var = (Variable *)calloc(1, sizeof(Variable));

  //======================Get Parameter========================
  parseCmdLine(box, particle, thermo, update, var, argc, argv);
  //===========================================================
  if (findVariable(var, "wf") < 0) {
    Abort("No --wf jp.bin");
  }

  double tic = getClock();

  constPressTriFireRelax(box, particle, thermo, update, var);

  write_data(box, particle, thermo, update, var);

  printf("===========System Info==========\n");
  printf("Number of Particles: %d;\n", particle->nAtom);
  printf("Volume Fraction: %g;\n", thermo->volFrac);
  printf("Min(diameter): %g;\nMax(diameter): %g;\n",
         update->nebrList.minDiameterScale, update->nebrList.maxDiameterScale);
  printf("Simulation Box: \n");
  printf("\t%-8.6e\t%-8.6e\t%-8.6e\n",
         box->boxHvoigt.h0 / thermo->distanceUnits,
         box->boxHvoigt.h5 / thermo->distanceUnits,
         box->boxHvoigt.h4 / thermo->distanceUnits);
  printf("\t%-8.6e\t%-8.6e\t%-8.6e\n", 0.0,
         box->boxHvoigt.h1 / thermo->distanceUnits,
         box->boxHvoigt.h3 / thermo->distanceUnits);
  printf("\t%-8.6e\t%-8.6e\t%-8.6e\n", 0.0, 0.0,
         box->boxHvoigt.h2 / thermo->distanceUnits);
  printf("Pressure tensor: \n");
  printf("\t%-8.6e\t%-8.6e\t%-8.6e\t%-8.6e\t%-8.6e\t%-8.6e\n",
         thermo->ptensor.h0 / thermo->pressureUnits,
         thermo->ptensor.h1 / thermo->pressureUnits,
         thermo->ptensor.h2 / thermo->pressureUnits,
         thermo->ptensor.h3 / thermo->pressureUnits,
         thermo->ptensor.h4 / thermo->pressureUnits,
         thermo->ptensor.h5 / thermo->pressureUnits);
  printf("===========System Info==========\n");

  double toc = getClock();
  safeFprintf(stdout, "Total Time: %g (%d:%d), nBuild: %ld, nForce: %ld.\n",
              toc - tic, shiftTimeStamp(tic), shiftTimeStamp(toc),
              update->nebrList.nBuild, update->nebrList.nForce);
}
