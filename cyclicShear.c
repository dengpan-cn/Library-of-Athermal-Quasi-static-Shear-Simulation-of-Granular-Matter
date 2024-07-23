
#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"
#include "FireMinAthermal_ndim.h"

// compile: make
// run: ./a.out --rf rf.bin --cvmin --sf suffix
// The format of binary file is described in function "readConf" and "writeConf".

// To perform 2D simulation, you should modify DIM in Makefile.
// To use different criteria, you should modify the parameters in FireMinAthermal_ndim.h.
// This code is used to perform simulations in manuscript "Jamming is a first-order transition with quenched disorder in amorphous materials sheared by cyclic quasistatic deformations".
// Please contact the authors if you need extra help.

int main(int argc, char const *argv[]) {
    MSG = "Cyclis Shear.";
    
    double tic = getTimeStamp();
    char ticStr[32];
    getTimeString(ticStr);
    //=============================================================
    Variable *var = (Variable *)calloc(1, sizeof(Variable));
    readCmdLineVar(var, argc, argv);
    
    Box *box = (Box *)calloc(1, sizeof(Box));
    Particle *particle = (Particle *)calloc(1, sizeof(Particle));
    Update *update = (Update *)calloc(1, sizeof(Update));
    
    readConf(box, particle, update, var);
    
    contactInfo *cinfo = addContactInfo(box, particle, update);
    ConstBoxShapeFIRE *fire = addConstBoxShapeFIRE(box, particle, update);

    char fname[4096];
    sprintf(fname,"%s/isoComp_%%8lf_%s.bin",var->cwd,var->sf);
    FILE *fout = createFileReadWrite(fname);
    //1. isotropically compressing to target packing fraction;
    while (update->volFrac < 0.8425) {
        instant_inflate(box, particle, update, 1E-4);
        reInitConstBoxShapeFIRE(fire);
        minConstBoxShapeFIRE(box, particle, update);
        computeContactInfo(box, particle, update);
        
        {
            double data[8];
            data[0] = update->volFrac;
            data[1] = update->ePair;
            data[2] = update->pVirTens[spaceIdx2voigt(0, 0)]/update->pressureUnits;
            data[3] = update->pVirTens[spaceIdx2voigt(1, 1)]/update->pressureUnits;
            data[4] = update->pVirTens[spaceIdx2voigt(0, 1)]/update->pressureUnits;
            data[5] = cinfo->aveCoordNum;
            data[6] = cinfo->aveCoordNumExRattler;
            data[7] = cinfo->nRattler;
            safeFwrite(fout, data, sizeof(double), 8);
        }
    }
    instant_inflate(box, particle, update, 0.8425 - update->volFrac);
    reInitConstBoxShapeFIRE(fire);
    minConstBoxShapeFIRE(box, particle, update);
    computeContactInfo(box, particle, update);
    {
        double data[8];
        data[0] = update->volFrac;
        data[1] = update->ePair;
        data[2] = update->pVirTens[spaceIdx2voigt(0, 0)]/update->pressureUnits;
        data[3] = update->pVirTens[spaceIdx2voigt(1, 1)]/update->pressureUnits;
        data[4] = update->pVirTens[spaceIdx2voigt(0, 1)]/update->pressureUnits;
        data[5] = cinfo->aveCoordNum;
        data[6] = cinfo->aveCoordNumExRattler;
        data[7] = cinfo->nRattler;
        safeFwrite(fout, data, sizeof(double), 8);
    }
    safeCloseFile(fout);
    
    writeConf(box, particle, update, var);
    
    //2. cyclic shear and recording data.
    reInitConstBoxShapeFIRE(fire);
    minConstBoxShapeFIRE(box, particle, update);
    computeContactInfo(box, particle, update);
    
    addWriteDumpFile(box, particle, update, var);
    writeDump(box, particle, update);
    sprintf(fname,"%s/cvs_%%9lf_%s.bin",var->cwd,var->sf);
    fout = createFileReadWrite(fname);
    double strain = 0.0;
    {
        double data[9];
        data[0] = 0;
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
    for (int icyc = 0; icyc < 100; icyc++) {
        double dstrain = 0.01;
        //a. shear
        for (int ith = 0; ith < 70; ith++) {
            instant_simpShearXz(box, particle, update, dstrain);
            strain += dstrain;
            reInitConstBoxShapeFIRE(fire);
            minConstBoxShapeFIRE(box, particle, update);
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
            reInitConstBoxShapeFIRE(fire);
            minConstBoxShapeFIRE(box, particle, update);
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
        
        writeDump(box, particle, update);
    }
    safeCloseFile(fout);
    
    //=============================================================
    double toc = getTimeStamp();
    char tocStr[32];
    getTimeString(tocStr);
    safeFprintf(stdout, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    safeFprintf(logFile, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    return EXIT_SUCCESS;
}

