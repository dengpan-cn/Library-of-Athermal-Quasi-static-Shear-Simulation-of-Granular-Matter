#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 minPressTriConf.c -o mptc
// run: ./mptc --rf rf.bin --wf wf.bin --cpmin 1E-4 --sf suffix
// The format of binary file is described in function "read_data" and
// "write_data"
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

  if (findVariable(var, "wf") < 0) Abort("No --wf cmd!");

  contactInfo *cinfo = initContactInfo(box, particle, thermo, update, var);

  constPressTriFireRelax(box, particle, thermo, update, var);
  computeCoordination(box, particle, thermo, update, var);

  printf("Dimension: %d\n", DIM);
  printf("No. of Particles: %d\n", particle->nAtom);
  printf("Mean diameter: %g (unit 1)\n", particle->meanDiameter);
  printf("Min(diameter): %g;\nMax(diameter): %g;\n",
         update->nebrList.minDiameterScale, update->nebrList.maxDiameterScale);
  printf("VolFrac: %g\n", thermo->volFrac);
  printf("Edges of Simulation Box: \n");
  for (int idim = 0; idim < DIM; idim++) {
    printf("\t");
    for (int jdim = 0; jdim < DIM; jdim++) {
      printf("%-8.6e\t", box->boxEdge[idim][jdim] / thermo->distanceUnits);
    }
    printf("\n");
  }
  printf("Epair: %g\n", thermo->Epair / thermo->energyUnits);
  printf("Pressure: %g\n", thermo->pressure / thermo->pressureUnits);
  printf("Ptensor(diagonal):");
  for (int idim = 0; idim < DIM; idim++) {
    printf(" %g",
           thermo->ptensor[spaceIdx(idim, idim)] / thermo->pressureUnits);
  }
  printf("\n");
  printf("Ptensor:");
  for (int idim = 0; idim < DIM; idim++) {
    for (int jdim = idim + 1; jdim < DIM; jdim++) {
      printf(" %g",
             thermo->ptensor[spaceIdx(idim, jdim)] / thermo->pressureUnits);
    }
  }
  printf("\n");
  printf("Number of Rattlers: %d, Z: %g, Z_NR: %g.\n", cinfo->nRattler,
         cinfo->aveCoordNum, cinfo->aveCoordNumExRattler);

  finalizeContactInfo(box, particle, thermo, update, var);
  write_data(box, particle, thermo, update, var);
  double toc = getClock();
  //=============================================================

  safeFprintf(stdout, "Total Time: %g (%d:%d), nBuild: %ld, nForce: %ld.\n",
              toc - tic, shiftTimeStamp(tic), shiftTimeStamp(toc),
              update->nebrList.nBuild, update->nebrList.nForce);
}
