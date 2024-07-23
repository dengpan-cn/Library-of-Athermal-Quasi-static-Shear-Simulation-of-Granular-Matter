#ifndef __SimSubFunc_ndim_dengPan_h__
#define __SimSubFunc_ndim_dengPan_h__

#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*FuncPtrToolkitWriteConf)(FILE *fout, Box *box, Particle *particle, Update *update);

void readCmdLineVar(Variable *var, int argc, char const *argv[]);

//====read and write configurations======
void readConf(Box *box, Particle *particle, Update *update, Variable *var);
void writeConf(Box *box, Particle *particle, Update *update, Variable *var);
int emergWriteConf(Box *box, Particle *particle, Update *update, char *fname);

//====compute force======
#if (DIM == 2 || DIM == 3)
#define image2local(jatom, nlocal, imgageSource) \
    (((jatom) < (nlocal)) ? (jatom) : (imgageSource[(jatom) - (nlocal)]))
#else
#define image2local(jatom, nlocal, imgageSource) (jatom)
#endif

#define Harmonic
//#define Hertzian
void adjustImg(Box *box, Particle *particle);
void calcForce(Box *box, Particle *particle, Update *update);
void calcForceAllPair(Box *box, Particle *particle, Update *update);


bool isNebrListValid(Box *box, Particle *particle, Update *update);
void buildNebrList(Box *box, Particle *particle, Update *update);
void genHalfNebrListPBCimage(Box *box, Particle *particle, Update *update);

int calcExpanCoeff(Box *box, Particle *particle, Update *update, doubleVector *deltaVeloc, double coeff[4]);

//=======box deformation========
void setBoxPara(Box *box);
void calcBoxLenAngle(Box *box, uptriMat para);
void setUnits(Update *update, double distanceUnits);
void calcDistBoxPlane(doubleVector distPlane, doubleVector boxEdge[DIM]);
void reInitSim(Box *box, Particle *particle, Update *update);

void calcTemp(Particle *particle, Update *update);
void calcPressure(Box *box, Particle *particle, Update *update);

void moveBarycenter(Box *box, Particle *particle, Update *update);

void scaleVeloc(Particle *particle, Update *update, double tT);
void genVeloc(Particle *particle, Update *update, double tT);

void instant_inflate(Box *box, Particle *particle, Update *update, double deltaVF);
void instant_simpShearXz(Box *box, Particle *particle, Update *update, double deltaGamma);
#if (DIM == 3)
void instant_simpShearXy(Box *box, Particle *particle, Update *update, double deltaGamma);
void instant_simpShearYz(Box *box, Particle *particle, Update *update, double deltaGamma);
#endif
void instant_pureShearXz(Box *box, Particle *particle, Update *update, double deltaStrain);
void instant_deformation(Box *box, Particle *particle, Update *update, uptriMat transMat);
void normaliseBox(Box *box, Particle *particle, Update *update);

//=============backup and restore simulation===============
void backupSimInfo(Box *box, Particle *particle, Box *bakBox, Particle *bakParticle);
void restoreSimInfo(Box *box, Particle *particle, Update *update, Box *bakBox, Particle *bakParticle);
void freeSimInfo(Box *bakBox, Particle *bakParticle);

//===============================================================================
typedef struct contactInfo {
    int refCnt;

    bool *isSoliton;  // contacting number is 0.
    bool *isRattler;  // Z_i < d + 1;
    bool *isBuckler;  // Z_i = d + 1;

    bool *isNebrContact;
    int allocMem4CoordNum;
    int *nCoordNum;
    int *nCoordNumExRattler;
    double aveCoordNum, aveCoordNumExRattler;
    int nSoliton, nRattler, nBuckler;

    // unjammed: nRattler/nAtom >= 0.5 || Z_NR - 2*DIM - 2*DIM/N_NR < 0
    bool isJammed;
    
    doubleVector *forceExRattler;
    double meanForceExRattler;
} contactInfo;
contactInfo *getContactInfo(Update *update);
contactInfo *addContactInfo(Box *box, Particle *particle, Update *update);
int computeContactInfo(Box *box, Particle *particle, Update *update);
int delContactInfo(Update *update);

//============write and read dump File==========================
typedef struct writeDumpFile {
    int refCnt;
    int revNum;
    FILE *fdump;
} writeDumpFile;
writeDumpFile *getWriteDumpFile(Update *update);
writeDumpFile *addWriteDumpFile(Box *box, Particle *particle, Update *update, Variable *var);
int writeDump(Box *box, Particle *particle, Update *update);
int delWriteDumpFile(Update *update);

typedef struct mmapBinFile {
    int refCnt;
    int fd;               // the descriptor the mapped file
    struct stat binStat;  // the stat of the opened file
    void *dataSection;    // pointer to the mapped memory
    int nStep;            // the number of timestep
    int stepSize;
    int headerSize;
    int revNum;
} mmapBinFile;  // binary dump file
mmapBinFile *openBinFile(char *fname);
int readSimInfo(Box *box, Particle *particle, Update *update, mmapBinFile *binFile);
int readDump(Box *box, Particle *particle, Update *update, mmapBinFile *binFile, int whichStep);
int closeBinFile(mmapBinFile **binFilePtr);

//===============================================================================
typedef struct fullNebrList {
    double rmax;

    int *list;
    int2 *nNebr;
    int maxAllocNebr;

    intVector nBin;
    int *binHead, *binList;
    int allocBinHead;
} fullNebrList;
int buildFullNebrList(Box *box, Particle *particle, double rmax, fullNebrList *fNebrList);
int delFullNebrList(fullNebrList *fNebrList);

void generateUnwarpPos(Box *box, Particle *particle);

//================================
void purgeSim(Box *box, Particle *particle, Update *update);

void showInfo(Box *box, Particle *particle, Update *update);
void showStructInfo(Box *box, Particle *particle, Update *update);

//================================
#if (DIM == 2 || DIM == 3)
void writeLammpsTopo(Box *box, Particle *particle, Update *update, Variable *var);
void writeLammpsTraj(Box *box, Particle *particle, Update *update, Variable *var, int globalTimeStep);
#endif

#ifdef __cplusplus
}
#endif

#endif
