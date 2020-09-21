#ifndef __STRUCTSIM__
#define __STRUCTSIM__

#include "VectorMath.h"

typedef struct Box {
  /*===
  The simulation box is a parallelepiped with edge:
            a = (h0,0,0), b = (h5, h1, 0), c = (h4, h3, h2)
  The edge vectors constitute the matrix:
            h = [[h0, h5, h4];[0, h1, h3];[0, 0, h2]];
  and the h is stroed in boxHvoigt.
  The inverse of matrix h (invH * h = h * invH = I):
            invH = [[ih0, ih5, ih4];[0, ih1, ih3];[0, 0, ih2]];
  is stored in invBoxHvoigt.
  We place the center of the parallelepiped at the origin (0,0,0).
  The front-let-below corner is boxLo and the opposite corner is boxHi.
  The volume of parallelepiped is trace(h) = h0 * h1 * h2;
  ===*/
  double3 boxLo, boxHi;
  Hvoigt6 boxHvoigt, invBoxHvoigt;
  double volume;
  bool isShapeFixed;  // whether the parallelepiped will change or not.
} Box;
typedef struct Particle {
  int nAtom, nAtomType;  // number of particles and particle type
  double3 *xyz;          // position
  double3 *veloc;        // velocity
  int3 *img;             // periodic image
  double3 *force;        // force
  double *diameterScale,
      meanDiameter;  // diameterScale = diameter/mean(diameter);

  //===just for future feature===
  int *type;            // type
  double *mass;         // mass
  double *massPerType;  // mass of each type
  //===just for future feature===

  // The memory are sorted to achieve better cache performance.
  // tag is the global label which is an attribute of particle.
  // id is the index which tells the memory location of particle.
  // id2tag[id] : the particle's label which is stored at location id;
  // tag2id[tag] : the memory location of particle labeled as tag;
  int *id2tag, *tag2id;

  bool isSizeFixed;   // whether the size of particle will change or not.
  bool isForceValid;  // whether the force is valid or not.
} Particle;

typedef struct contactInfo {
  // contact info
  bool *isRattler;
  bool *isNebrContact;
  int allocMem4CoordNum;
  int *nCoordNum;
  int *nCoordNumExRattler;
  double aveCoordNum, aveCoordNumExRattler;
  int nRattler;

  double3 *forceExRattler;
  double lapUexRattler;
  double meanForceExRattler;
} contactInfo;
typedef struct dumpInfo {
  // output configuration
  FILE *fdump;
} dumpInfo;
typedef struct {
  int fd;               // the descriptor the mapped file
  struct stat binStat;  // the stat of the opened file
  void *dataSection;    // pointer to the mapped memory
  int nStep;            // the number of timestep
  int stepSize;         // the bytes of one step
  int headerSize;       // the bytes of the binary dump file
} BinDumpFile;          // binary dump file

typedef struct Thermo {
  double massUnits, energyUnits, distanceUnits;
  double timeUnits, forceUnits, velocityUnits;
  double pressureUnits, volumeUnits;  // units

  double dof;                   // degree of freedom
  double Ekin, Temperature;     // kinetic energy
  Hvoigt6 Kintensor;            // m*v_alpha*v_beta
  double Epair;                 // pair potential
  double virialPair, pressure;  // virial term
  Hvoigt6 ptensor;              // pressure tensor
  double volFrac;               // packing fraction
  double lapU;  // Laplacian of U, units (forceUnits / distanceUnits)

  bool isInit, Edone, Pdone, Tdone;  // flags

  Toolkit toolkit;  // toolkit
} Thermo;

#define __minSkinSet__ 5E-3
typedef struct idPosRadius {
  int id;
  double3 pos;
  double radius;
} idPosRadius;
typedef struct NebrList {
  // neighbour list
  // The neighbour of particle iatom:
  // from list[nNebr[iatom].x] (included) to list[nNebr[iatom].y] (excluded);
  double maxDiameterScale, minDiameterScale;
  double skinSet, rskin;

  double3 binLen;
  int3 nbin;
  int totBin, allocBin;

  double3 *xyzHold;
  double meanDiameterHold;
  Hvoigt6 invBoxHvoigtHold;

  bool isFullStyle;  // default is false.

  //====bin list====
  idPosRadius *ipr;
  int *nAtomPerBin;
  int maxAtomPerBin;

  //====Nebr list====
  int *list;
  int2 *nNebr;
  int maxAllocNebr;

  //===adjacent bin===
  int3 *deltaAdjBin;
  int nAdjBin;
  //====

  int *binHead4sort, *binList4sort;
  int totBin4sort, allocBin4sort;
  int *oid2nid;
  void *buffer;

  long int nBuild, nForce, cntForce;
  bool isInit, isValid, doSort;
  bool compelInit;
} NebrList;

// parameters of FIRE
#define DELAYSTEP 5
#define DT_GROW 1.1
#define ALPHA_SHRINK 0.99
#define DT_SHRINK 0.5
#define ALPHA0 0.1
#define maxDeltaVF 1E-3
#define __Ptol__ 1E-4
#define __StressTol__ 1E-14
#define __AveZeroForce__ 1E-14
#define __ZeroEnergy__ 1E-16

typedef struct minConstBoxShapeFIRE {
  // const-Box-shape FIRE: minimize U(r)
  // The shape of parallelepiped is fixed!
  // The relaxation stops if the system unjammed or the minimum reached.
  // unjamming criteria: pair potential per particle is less than __ZeroEnergy__
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  double dtSet;
  double vdotf, vdotv, fdotf;
  double aveForce;
  double dt;

  bool isInit;
} minConstBoxShapeFIRE;
typedef struct minConstBoxVolFIRE {
  // const-Box-volume FIRE: minimize U(r)
  // The volume of parallelepiped is fixed!
  // The relaxation stops if the system unjammed or the minimum reached.
  // unjamming criteria: pair potential per particle is less than __ZeroEnergy__
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and the relative difference among (Pxx, Pyy, Pzz) are less than __Ptol__
  // and the absolute value of (Pxy, Pxz, Pxy) are less than __StressTol__
  double3 gVelocAniso, sfactAniso, gForceAniso;
  double3 gVelocStress, sfactStress, gForceStress;

  double dtSet;
  double vdotf, vdotv, fdotf;
  double aveForce;
  double dt;

  bool isInit;
} minConstBoxVolFIRE;

typedef struct minConstPressIsoFIRE {
  // const-iso-Pressure FIRE: minimize H(r,h) = U(r) + Pset * trace(h).
  // The mean-diameter of particles is adjusted and the box shape is keeped.
  // The relaxation stops if the minimum of H(r,h) reached.
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and the absolute value of (P / Pset - 1.0) is less than __Ptol__

  double Pset, ptol;
  double gVeloc, sfact, gForce;

  double dtSet;
  double vdotf, vdotv, fdotf;
  double aveForce;
  double dt;

  bool isInit;
} minConstPressIsoFIRE;
typedef struct minConstPressTriFIRE {
  // const-tri-Pressure-tensor FIRE: minimize H(r,h) = U(r) + Pset * trace(h).
  // All six parameters of box are adjusted.
  // minimum criteria: average force amplitude is less than __AveZeroForce__
  // and relative difference between Pxx,Pyy,Pzz and Pset are less than __Ptol__
  // and the absolute value of (Pxy, Pxz, Pxy) are less than __StressTol__

  double Pset, ptol;
  double3 gVelocAniso, sfactAniso, gForceAniso;
  double3 gVelocStress, sfactStress, gForceStress;

  double dtSet;
  double vdotf, vdotv, fdotf;
  double aveForce;
  double dt;

  bool isInit;
} minConstPressTriFIRE;

typedef struct pureShear {
  FILE *shearLog;

  double strain;  // strain = (Lz'-Lz0)/Lz;
  // Lz' = (1+strain)*Lz0; Ly' = Ly0; Lx' = Lx0/(1+strain);
  double deltaStrain, maxStrain;
  double3 refBoxLen;

  bool isInit;
} pureShear;
typedef struct simpleShear {
  double gamma, deltaGamma, maxGamma;
  FILE *shearLog;
  dumpInfo *dinfo;
  contactInfo *cinfo;

  bool isInit;
} simpleShear;

typedef struct Update {
  NebrList nebrList;

  minConstBoxShapeFIRE *cbsFire;
  minConstBoxVolFIRE *cbvFire;
  minConstPressIsoFIRE *cpiFire;
  minConstPressTriFIRE *cptFire;

  pureShear *pShear;
  simpleShear *sShear;
} Update;

#endif