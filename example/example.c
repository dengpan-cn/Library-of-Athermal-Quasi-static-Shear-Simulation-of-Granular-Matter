// header
#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 example.c -o app
// run: ./app --trunc --rf startConf.bin --wf stopConf.bin --cvmin
// --rf startConf.bin: the format of startConf.bin is explained in "read_data"
int main(int argc, char const *argv[]) {
  // allocate memory
  Box *box = (Box *)calloc(1, sizeof(Box));
  Particle *particle = (Particle *)calloc(1, sizeof(Particle));
  Thermo *thermo = (Thermo *)calloc(1, sizeof(Thermo));
  Update *update = (Update *)calloc(1, sizeof(Update));
  Variable *var = (Variable *)calloc(1, sizeof(Variable));

  //======================Get Parameter and read data========================
  parseCmdLine(box, particle, thermo, update, var, argc, argv);
  //=========================================================================

  addVariable(var, "--cvmin");  //"--"
  addVariable(var, "--wf minConstBoxShape.bin");
  constBoxShapeFireRelax(box, particle, thermo, update, var);
  write_data(box, particle, thermo, update, var);
  delVariable(var, "cvmin");  // no "--"
  delVariable(var, "wf");

  addVariable(var, "--cvmin");  //"--"
  addVariable(var, "--wf minConstVol.bin");
  constBoxVolFireRelax(box, particle, thermo, update, var);
  write_data(box, particle, thermo, update, var);
  delVariable(var, "cvmin");  // no "--"
  delVariable(var, "wf");

  addVariable(var, "--cpmin 1E-3");  //"--"
  addVariable(var, "--wf minConstIsoPressure1Ei3.bin");
  constPressIsoFireRelax(box, particle, thermo, update, var);
  write_data(box, particle, thermo, update, var);
  delVariable(var, "cpmin");  // no "--"
  delVariable(var, "wf");

  addVariable(var, "--cpmin 1E-3");  //"--"
  addVariable(var, "--wf minConstPressureTensor1Ei3.bin");
  constPressTriFireRelax(box, particle, thermo, update, var);
  write_data(box, particle, thermo, update, var);
  delVariable(var, "cpmin");  // no "--"
  delVariable(var, "wf");

  return 0;
}
