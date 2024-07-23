#include "FireFunc.h"
#include "SimSubFunc.h"
#include "StructSim.h"
#include "VectorMath.h"

extern FILE *logFile;
extern int truncFileFlag;

// compile: gcc -std=gnu99 -lm -O3 cyclicShear.c -o cycShear
// run: ./cycShear --rf rf.bin --cvmin --sf suffix
// The format of binary file is described in function "read_data" and "write_data"

// To perform 2D simulation, you shuold modify DIM in VectorMath.h.
// To use different criteria, you should modify the parameters in StructSim.h.
// This code is used to perform simulations in article "Jamming is a first-order transition with quenched disorder in amorphous materials sheared by cyclic quasistatic deformations".
// Please contact the authors if you need extra help.

#if (DIM != 2)
#error "DIM should be 2."
#endif

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
    
    printf("Force criterion: mean(|F|) <= %.14lf\n", __ZeroForce__);
    printf("Criterion of unjamming: PE <= %.14lf\n", __ZeroEnergy__);
    
    char fname[4096];
    sprintf(fname, "%s/cvs_%%9lf_%s.bin", var->cwd, var->sf);
    FILE *fout = createReadWriteFile(fname);
    double strain = 0.0;
    for (int icyc = 0; icyc < 100; icyc++) {
        double dstrain = 1e-2;
        //a. shear
        for (int ith = 0; ith < 70; ith++) {
            instant_simpShearXz(box, particle, update, dstrain);
            strain += dstrain;
            
            constBoxShapeFireRelax(box, particle, thermo, update, var);
            
            computeContactInfo(box, particle, update);
            
            {
                double data[9];
                data[0] = icyc;
                data[1] = strain;
                data[2] = update->ePair;
                data[3] = update->pVirTens[spaceIdx2voigt(0, 0)]/update->pressureUnits;
                data[4] = update->pVirTens[spaceIdx2voigt(1, 1)]/update->pressureUnits;
                data[5] = update->pVirTens[spaceIdx2voigt(0, 1)]/update->pressureUnits;
                data[6] = cinfo->aveCoordNum;
                data[7] = cinfo->aveCoordNumExRattler;
                data[8] = cinfo->nRattler;
                safeFwrite(fout, data, sizeof(double), 9);
            }
            
        }
        
        //b. reverse shear
        for (int ith = 0; ith < 70; ith++) {
            instant_simpShearXz(box, particle, update, -dstrain);
            strain -= dstrain;
            
            constBoxShapeFireRelax(box, particle, thermo, update, var);
            
            computeContactInfo(box, particle, update);
            
            {
                double data[9];
                data[0] = icyc;
                data[1] = strain;
                data[2] = update->ePair;
                data[3] = update->pVirTens[spaceIdx2voigt(0, 0)]/update->pressureUnits;
                data[4] = update->pVirTens[spaceIdx2voigt(1, 1)]/update->pressureUnits;
                data[5] = update->pVirTens[spaceIdx2voigt(0, 1)]/update->pressureUnits;
                data[6] = cinfo->aveCoordNum;
                data[7] = cinfo->aveCoordNumExRattler;
                data[8] = cinfo->nRattler;
                safeFwrite(fout, data, sizeof(double), 9);
            }
            
        }
        
    }
    safeCloseFile(fout);
    
    
    double toc = getClock();
    //=============================================================
    
    safeFprintf(stdout, "Total Time: %g (%d:%d), nBuild: %ld, nForce: %ld.\n",
                toc - tic, shiftTimeStamp(tic), shiftTimeStamp(toc),
                update->nebrList.nBuild, update->nebrList.nForce);
}
