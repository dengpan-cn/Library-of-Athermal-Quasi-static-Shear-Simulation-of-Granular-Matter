#ifndef __VECTORFUNC__
#define __VECTORFUNC__

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>  //-lm
//#include <quadmath.h>  //-lquadmath
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define DIM 3  // dimension of simulation box, currently, it is 3.
#define PI 3.141592653589793238462643  // pi
#ifndef VolUnitSphere
// volume of unit sphere in space with dimension of DIM.
#define VolUnitSphere               \
  (pow(PI, ((double)(DIM)) / 2.0) / \
   exp(lgamma(1 + ((double)(DIM)) / 2.0)))  // volume pre-factor for sphere
#endif

/*===
The simulation box is a parallelepiped with edge:
          a = (h0,0,0), b = (h5, h1, 0), c = (h4, h3, h2)
__orthBox__: h3 = h4 = h5 = 0;
__mono_XZ_Box__: h3 = h5 = 0; shear along XZ direction.
__triBox__: general parallelepiped.
===*/
#define __triBox__
// #define __orthBox__
// #define __mono_XZ_Box__

typedef struct double3 {
  double x, y, z;
} double3;
typedef struct int3 {
  int x, y, z;
} int3;
typedef struct int2 {
  int x, y;
} int2;
// structure for quantities with six elements
typedef struct Hvoigt6 {
  double h0, h1, h2, h3, h4, h5;  // xx, yy, zz, yz, xz, xy
} Hvoigt6;

#ifndef __VEC__
#define __VEC__
// mathematical operation for double3 and int3
#define vecAdd(vSum, vecA, vecB) \
  {                              \
    vSum.x = vecA.x + vecB.x;    \
    vSum.y = vecA.y + vecB.y;    \
    vSum.z = vecA.z + vecB.z;    \
  }
#define vecSub(vSub, vecA, vecB) \
  {                              \
    vSub.x = vecA.x - vecB.x;    \
    vSub.y = vecA.y - vecB.y;    \
    vSub.z = vecA.z - vecB.z;    \
  }
#define vecDot(sDot, vecA, vecB) \
  { sDot = vecA.x * vecB.x + vecA.y * vecB.y + vecA.z * vecB.z; }
#define vecCross(vCross, vecA, vecB)                 \
  {                                                  \
    vCross.x = vecA.y * vecB.z - vecA.z * vecB.y;    \
    vCross.y = -(vecA.x * vecB.z - vecA.z * vecB.x); \
    cross.z = vecA.x * vecB.y - vecA.y * vecB.x;     \
  }
#define vecScale(vScale, scale, vecA) \
  {                                   \
    vScale.x = (scale)*vecA.x;        \
    vScale.y = (scale)*vecA.y;        \
    vScale.z = (scale)*vecA.z;        \
  }
#define vecScaleAdd(vSum, vecA, scale, vecB) \
  {                                          \
    vSum.x = vecA.x + (scale)*vecB.x;        \
    vSum.y = vecA.y + (scale)*vecB.y;        \
    vSum.z = vecA.z + (scale)*vecB.z;        \
  }
#define vecNormP2(sNormP2, vecA) vecDot(sNormP2, vecA, vecA)
#define vecNorm(sNorm, vecA)   \
  {                            \
    vecDot(sNorm, vecA, vecA); \
    sNorm = sqrt(sNorm);       \
  }
#define vecUnit(vUnit, vecA)          \
  {                                   \
    double len = 0.0;                 \
    vecNorm(len, vecA);               \
    vecScale(vUnit, 1.0 / len, vecA); \
  }
#define vecHvoigtMulVec(vHmV, H, vecA)                      \
  {                                                         \
    vHmV.x = H.h0 * vecA.x + H.h5 * vecA.y + H.h4 * vecA.z; \
    vHmV.y = H.h1 * vecA.y + H.h3 * vecA.z;                 \
    vHmV.z = H.h2 * vecA.z;                                 \
  }
#define maxElementVec(sMax, vecA)               \
  {                                             \
    sMax = (vecA.x > vecA.y ? vecA.x : vecA.y); \
    sMax = (sMax > vecA.z ? sMax : vecA.z);     \
  }
#define minElementVec(sMin, vecA)               \
  {                                             \
    sMin = (vecA.x < vecA.y ? vecA.x : vecA.y); \
    sMin = (sMin < vecA.z ? sMin : vecA.z);     \
  }
#endif

// PBC(vector, box): do periodic boundary condition, vector will be changed!
#ifdef __orthBox__
#define PBC(nRij, box)                  \
  {                                     \
    if (nRij.z >= box->boxHi.z) {       \
      nRij.z -= box->boxHvoigt.h2;      \
    } else if (nRij.z < box->boxLo.z) { \
      nRij.z += box->boxHvoigt.h2;      \
    }                                   \
    if (nRij.y >= box->boxHi.y) {       \
      nRij.y -= box->boxHvoigt.h1;      \
    } else if (nRij.y < box->boxLo.y) { \
      nRij.y += box->boxHvoigt.h1;      \
    }                                   \
    if (nRij.x >= box->boxHi.x) {       \
      nRij.x -= box->boxHvoigt.h0;      \
    } else if (nRij.x < box->boxLo.x) { \
      nRij.x += box->boxHvoigt.h0;      \
    }                                   \
  }
#endif

#ifdef __mono_XZ_Box__
#define PBC(nRij, box)                          \
  {                                             \
    if (nRij.z >= box->boxHi.z) {               \
      nRij.z -= box->boxHvoigt.h2;              \
      nRij.x -= box->boxHvoigt.h4;              \
    } else if (nRij.z < box->boxLo.z) {         \
      nRij.z += box->boxHvoigt.h2;              \
      nRij.x += box->boxHvoigt.h4;              \
    }                                           \
    if (nRij.y >= box->boxHi.y)                 \
      nRij.y -= box->boxHvoigt.h1;              \
    else if (nRij.y < box->boxLo.y)             \
      nRij.y += box->boxHvoigt.h1;              \
    if (nRij.x >= 0.5 * box->boxHvoigt.h0)      \
      nRij.x -= box->boxHvoigt.h0;              \
    else if (nRij.x < -0.5 * box->boxHvoigt.h0) \
      nRij.x += box->boxHvoigt.h0;              \
  }
#endif

#ifdef __triBox__
#define PBC(nRij, box)                              \
  {                                                 \
    if (nRij.z >= box->boxHi.z) {                   \
      nRij.z -= box->boxHvoigt.h2;                  \
      nRij.y -= box->boxHvoigt.h3;                  \
      nRij.x -= box->boxHvoigt.h4;                  \
    } else if (nRij.z < box->boxLo.z) {             \
      nRij.z += box->boxHvoigt.h2;                  \
      nRij.y += box->boxHvoigt.h3;                  \
      nRij.x += box->boxHvoigt.h4;                  \
    }                                               \
    if (nRij.y >= 0.5 * box->boxHvoigt.h1) {        \
      nRij.y -= box->boxHvoigt.h1;                  \
      nRij.x -= box->boxHvoigt.h5;                  \
    } else if (nRij.y < -0.5 * box->boxHvoigt.h1) { \
      nRij.y += box->boxHvoigt.h1;                  \
      nRij.x += box->boxHvoigt.h5;                  \
    }                                               \
    if (nRij.x >= 0.5 * box->boxHvoigt.h0) {        \
      nRij.x -= box->boxHvoigt.h0;                  \
    } else if (nRij.x < -0.5 * box->boxHvoigt.h0) { \
      nRij.x += box->boxHvoigt.h0;                  \
    }                                               \
  }
#endif

// global variable
FILE *logFile = NULL;
int truncFileFlag = 0;

#ifndef safeFprintf
#define safeFprintf(fileFp, format, ...)      \
  {                                           \
    if (fileFp != NULL) {                     \
      fprintf(fileFp, format, ##__VA_ARGS__); \
      fflush(fileFp);                         \
    }                                         \
  }
#endif

#ifndef safeCloseFile
#define safeCloseFile(fileFp) \
  {                           \
    if (fileFp != NULL) {     \
      fclose(fileFp);         \
      fileFp = NULL;          \
    }                         \
  }
#endif

#ifndef safeFree
#define safeFree(ptr)  \
  {                    \
    if (ptr != NULL) { \
      free(ptr);       \
      ptr = NULL;      \
    }                  \
  }
#endif

#ifndef Abort
#define Abort(format, ...)                                                     \
  {                                                                            \
    safeFprintf(stderr, "File: %s; Func: %s; Line: %d\n\tError: " format "\n", \
                __FILE__, __func__, __LINE__, ##__VA_ARGS__);                  \
    safeFprintf(logFile,                                                       \
                "File: %s; Func: %s; Line: %d\n\tError: " format "\n",         \
                __FILE__, __func__, __LINE__, ##__VA_ARGS__);                  \
    exit(EXIT_FAILURE);                                                        \
  }
#endif

#ifndef Info
#define Info(format, ...)                                                     \
  {                                                                           \
    safeFprintf(stderr, "File: %s; Func: %s; Line: %d\n\tInfo: " format "\n", \
                __FILE__, __func__, __LINE__, ##__VA_ARGS__);                 \
  }
#endif

#ifndef getClock
#define getClock()                               \
  ({                                             \
    struct timeval sTv;                          \
    gettimeofday(&sTv, NULL);                    \
    double tic = sTv.tv_sec + sTv.tv_usec / 1e6; \
    tic;                                         \
  })
#endif

#ifndef shiftTimeStamp
#define shiftTimeStamp(sec) ((int)(sec - 1.6E9))
#endif

#ifndef getTimeString
#define getTimeString()                                                      \
  ({                                                                         \
    struct timeval sTv;                                                      \
    gettimeofday(&sTv, NULL);                                                \
    struct tm sTm;                                                           \
    localtime_r((time_t *)&sTv.tv_sec, &sTm);                                \
    char str[32];                                                            \
    sprintf(str, "%d-%d-%d-%d-%d", sTm.tm_mon + 1, sTm.tm_mday, sTm.tm_hour, \
            sTm.tm_min, sTm.tm_sec);                                         \
    str;                                                                     \
  })
#endif

#ifndef cpuMin
#define cpuMin(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef cpuMax
#define cpuMax(a, b) ((a) > (b) ? (a) : (b))
#endif

//--trunc or --append
#ifndef createReadWriteFile
#define createReadWriteFile(filePath)                                 \
  ({                                                                  \
    if (filePath == NULL) Abort("No variable!");                      \
    FILE *fp = NULL;                                                  \
    if (access(filePath, F_OK) == 0) {                                \
      if (truncFileFlag == 0) Abort("File %s exist! Exit", filePath); \
      if (truncFileFlag == 1) {                                       \
        Info("File %s exist! Truncated!", filePath);                  \
        fp = fopen(filePath, "wb+");                                  \
      }                                                               \
      if (truncFileFlag == 2) {                                       \
        Info("File %s exist! Appended!", filePath);                   \
        fp = fopen(filePath, "ab+");                                  \
      }                                                               \
    } else {                                                          \
      fp = fopen(filePath, "wb+");                                    \
    }                                                                 \
    if (fp == NULL) Abort("Open File %s Failed!", filePath);          \
    fp;                                                               \
  })
#endif

#ifndef openAppendedFile
#define openAppendedFile(filePath)                           \
  ({                                                         \
    if (filePath == NULL) Abort("No variable!");             \
    FILE *fp = fopen(filePath, "ab+");                       \
    if (fp == NULL) Abort("Open File %s Failed!", filePath); \
    fp;                                                      \
  })
#endif

#ifndef openReadOnlyFile
#define openReadOnlyFile(filePath)                                   \
  ({                                                                 \
    if (filePath == NULL) Abort("No variable!");                     \
    if (access(filePath, F_OK) != 0) Abort("No File %s!", filePath); \
    FILE *fp = fopen(filePath, "r");                                 \
    if (fp == NULL) Abort("Open File %s Failed!", filePath);         \
    fp;                                                              \
  })
#endif

// generate random number with stdandard normal distribution N(0,1).
double rndStdNorm() {
  static int setSeed = 1;
  if (setSeed == 1) {
    srand((unsigned int)time((time_t *)NULL));
    setSeed = 0;
  }

  static int phase = 0;
  static double v1, v2, s;

  double x;
  if (0 == phase) {
    do {
      double u1 = (double)rand() / RAND_MAX;
      double u2 = (double)rand() / RAND_MAX;

      v1 = 2 * u1 - 1;
      v2 = 2 * u2 - 1;
      s = v1 * v1 + v2 * v2;
    } while (1 <= s || 0 == s);
    x = v1 * sqrt(-2 * log(s) / s);
  } else {
    x = v2 * sqrt(-2 * log(s) / s);
  }
  phase = 1 - phase;

  return x;
}
Hvoigt6 make_Hvoigt6(double h0, double h1, double h2, double h3, double h4,
                     double h5) {
  Hvoigt6 tmp;
  tmp.h0 = h0;
  tmp.h1 = h1;
  tmp.h2 = h2;
  tmp.h3 = h3;
  tmp.h4 = h4;
  tmp.h5 = h5;
  return tmp;
}
double3 make_double3(double x, double y, double z) {
  double3 tmp;
  tmp.x = x;
  tmp.y = y;
  tmp.z = z;
  return tmp;
}
int3 make_int3(int x, int y, int z) {
  int3 tmp;
  tmp.x = x;
  tmp.y = y;
  tmp.z = z;
  return tmp;
}
int2 make_int2(int x, int y) {
  int2 tmp;
  tmp.x = x;
  tmp.y = y;
  return tmp;
}
int isEmpty(char *str, int maxlen) {
  if (maxlen <= 0) Abort("Fatal Error!");
  if (strlen(str) >= maxlen) Abort("Fatal Error");
  if (strlen(str) == 0) return 1;

  char *start = str, *stop = str + strlen(str) - 1;
  while (isspace(*stop) && stop >= start) stop--;
  while (isspace(*start) && start <= stop) start++;

  if (start > stop) return 1;

  return 0;
}

typedef struct cmdArg {
  char *cmdType;
  int cmdArgc;
  char **cmdArgv;
} cmdArg;

typedef struct Variable {
  int nVar, maxVar;
  cmdArg *cmd;

  char *cwd, *sf;
} Variable;

int findVariable(Variable *var, char *name) {
  for (int ith = 0; ith < var->nVar; ith++) {
    if (!strcmp(var->cmd[ith].cmdType, name)) return ith;
  }
  return -1;
}
void addVariable(Variable *var, char *inputCmd) {
  if (strlen(inputCmd) < 3) Abort("No cmd!");

  char *tmp = inputCmd;
  while (tmp[0] != '\0' && isspace(tmp[0])) {
    tmp++;
  }
  if (strlen(tmp) <= 2) Abort("Unrecognized cmd %s!", inputCmd);
  if (strncmp(tmp, "--", 2)) Abort("Unrecognized cmd %s!", inputCmd);
  if (isspace(tmp[2])) Abort("Unrecognized cmd %s!", inputCmd);
  tmp += 2;
  char *str = (char *)calloc(strlen(tmp) + 3, sizeof(char));
  sprintf(str, "%s", tmp);

  int narg = 0, maxArg = 8;
  char **arg = (char **)calloc(maxArg, sizeof(char *));
  arg[narg++] = strtok(str, " \t");
  char *ptr = NULL;
  while ((ptr = strtok(NULL, " \t")) != NULL) {
    arg[narg++] = ptr;
    if (narg == maxArg) {
      maxArg += 8;
      arg = (char **)realloc(arg, maxArg * sizeof(char *));
    }
  }

  if (var->nVar == var->maxVar) {
    var->maxVar += 8;
    var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
  }
  if (findVariable(var, str) >= 0) Abort("Repetitive cmd!");
  var->cmd[var->nVar].cmdType = str;
  var->cmd[var->nVar].cmdArgc = narg - 1;
  var->cmd[var->nVar].cmdArgv = NULL;
  if (var->cmd[var->nVar].cmdArgc != 0) {
    var->cmd[var->nVar].cmdArgv =
        (char **)calloc(var->cmd[var->nVar].cmdArgc, sizeof(char *));
  }
  for (int ith = 1; ith < narg; ith++) {
    if (!strncmp(arg[ith], "--", 2)) Abort("Add ONE variable!");
    var->cmd[var->nVar].cmdArgv[ith - 1] = arg[ith];
  }
  var->nVar++;

  free(arg);
}
void delVariable(Variable *var, char *delCmd) {
  int whichVar = findVariable(var, delCmd);
  if (whichVar < 0) return;

  free(var->cmd[whichVar].cmdArgv);
  var->cmd[whichVar].cmdArgv = NULL;
  free(var->cmd[whichVar].cmdType);
  var->cmd[whichVar].cmdType = NULL;
  var->cmd[whichVar].cmdArgc = 0;

  if (whichVar != var->nVar - 1) {
    var->cmd[whichVar].cmdType = var->cmd[var->nVar - 1].cmdType;
    var->cmd[whichVar].cmdArgc = var->cmd[var->nVar - 1].cmdArgc;
    var->cmd[whichVar].cmdArgv = var->cmd[var->nVar - 1].cmdArgv;

    var->cmd[var->nVar - 1].cmdType = NULL;
    var->cmd[var->nVar - 1].cmdArgc = 0;
    var->cmd[var->nVar - 1].cmdArgv = NULL;
  }

  var->nVar--;
}

typedef struct Toolkit {
  int nToolkit, maxToolkit;
  char **toolkitName;
  void **toolkit;
} Toolkit;
int findToolkit(Toolkit *toolkit, char *name) {
  for (int ith = 0; ith < toolkit->nToolkit; ith++) {
    if (!strcmp(name, toolkit->toolkitName[ith])) return ith;
  }
  return -1;
}
void addToolkit(Toolkit *toolkit, void *tool, char *name) {
  if (findToolkit(toolkit, name) >= 0) Abort("Repetitive tool!");

  if (toolkit->maxToolkit <= toolkit->nToolkit) {
    toolkit->maxToolkit += 8;
    toolkit->toolkit = (void **)realloc(toolkit->toolkit,
                                        toolkit->maxToolkit * sizeof(void *));
    toolkit->toolkitName = (char **)realloc(
        toolkit->toolkitName, toolkit->maxToolkit * sizeof(char *));
  }
  toolkit->toolkit[toolkit->nToolkit] = tool;
  toolkit->toolkitName[toolkit->nToolkit] =
      (char *)calloc(strlen(name) + 5, sizeof(char));
  sprintf(toolkit->toolkitName[toolkit->nToolkit], "%s", name);
  toolkit->nToolkit++;
}
void delToolkit(Toolkit *toolkit, char *name) {
  if (findToolkit(toolkit, name) < 0) return;
  int which = findToolkit(toolkit, name);
  free(toolkit->toolkitName[which]);
  free(toolkit->toolkit[which]);

  if (which != toolkit->nToolkit - 1) {
    toolkit->toolkitName[which] = toolkit->toolkitName[toolkit->nToolkit - 1];
    toolkit->toolkit[which] = toolkit->toolkit[toolkit->nToolkit - 1];
  }

  toolkit->nToolkit--;
}

void exchange_double3(double3 *xyz, double3 *buffer, int *oid2nid, int nAtom) {
  for (int oid = 0; oid < nAtom; oid++) {
    int nid = oid2nid[oid];
    buffer[nid] = xyz[oid];
  }
}
void exchange_int3(int3 *img, int3 *buffer, int *oid2nid, int nAtom) {
  for (int oid = 0; oid < nAtom; oid++) {
    int nid = oid2nid[oid];
    buffer[nid] = img[oid];
  }
}
void exchange_double(double *radius, double *buffer, int *oid2nid, int nAtom) {
  for (int oid = 0; oid < nAtom; oid++) {
    int nid = oid2nid[oid];
    buffer[nid] = radius[oid];
  }
}
void exchange_int(int *type, int *buffer, int *oid2nid, int nAtom) {
  for (int oid = 0; oid < nAtom; oid++) {
    int nid = oid2nid[oid];
    buffer[nid] = type[oid];
  }
}

#define heapSortSwap(arr, a, b) \
  {                             \
    int tmp = arr[b];           \
    arr[b] = arr[a];            \
    arr[a] = tmp;               \
  }

void max_heapify(int *arr, int start, int end) {
  int dad = start;
  int son = dad * 2 + 1;
  while (son <= end) {
    if (son + 1 <= end && arr[son] < arr[son + 1]) {
      son++;
    }
    if (arr[dad] > arr[son]) {
      return;
    } else {
      heapSortSwap(arr, dad, son);
      dad = son;
      son = dad * 2 + 1;
    }
  }
}
void heapSortInc(int *arr, int len) {
  // sort array with increasing order
  if (len <= 1) return;

  for (int ith = len / 2 - 1; ith >= 0; ith--) max_heapify(arr, ith, len - 1);
  for (int ith = len - 1; ith > 0; ith--) {
    heapSortSwap(arr, 0, ith);
    max_heapify(arr, 0, ith - 1);
  }
}

#endif