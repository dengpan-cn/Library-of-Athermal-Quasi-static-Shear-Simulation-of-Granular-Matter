#ifndef __VectorMath_ndim_dengPan_h__
#define __VectorMath_ndim_dengPan_h__

#if !defined(DIM) || (defined(__triBox__) == defined(__orthBox__))
#error "Pleade define 1. DIM; 2. __triBox__ or __orthBox__."
#endif

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>  //-lm
// #include <quadmath.h>  //-lquadmath
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

// #define checkMacroFlag

// #define __triBox__
// #define __orthBox__
// #define DIM 3

#define voigtDIM ((DIM * (DIM + 1) / 2))
#define voigtSymMatDIM ((voigtDIM * (voigtDIM + 1) / 2))

#define PI 3.141592653589793238462643

#ifndef VolUnitSphere
// Vsph(radius) = VolUnitSphere * pow(radius,DIM)
#define VolUnitSphere \
(pow(PI, ((double)(DIM)) / 2.0) / exp(lgamma(1 + ((double)(DIM)) / 2.0)))
#endif

// make local variable
#define mkVar(tag, var) var##tag

extern FILE *logFile;
extern int truncFileFlag;         // 0:normal; 1:truncate; 2:restart
extern int screenOutputReadConf;  // 0: not; 1: yes
extern char *MSG;                 // readme of program
#define dumpFileRevNum 1          // 0: SS 3d; 1: SS; 2: HS;

// bit array: char bitArr[BitNslots(numBit)];
typedef char *BitArr;
#define BitMask(whichBit) (1 << ((whichBit) % CHAR_BIT))
#define BitSlot(whichBit) ((whichBit) / CHAR_BIT)
#define BitSet(bitArr, whichBit) ((bitArr)[BitSlot(whichBit)] |= BitMask(whichBit))
#define BitClear(bitArr, whichBit) ((bitArr)[BitSlot(whichBit)] &= ~BitMask(whichBit))
#define BitTest(bitArr, whichBit) ((bitArr)[BitSlot(whichBit)] & BitMask(whichBit))
#define BitNslots(numBit) ((numBit + CHAR_BIT - 1) / CHAR_BIT)
#define BitClearAll(bitArr, numBit)                                          \
{                                                                        \
    size_t mkVar(bitClearAll, nByte) = BitNslots(numBit) * sizeof(char); \
    memset(bitArr, '\0', mkVar(bitClearAll, nByte));                     \
}

#ifndef safeFprintf
#define safeFprintf(fileFp, format, ...)            \
{                                               \
    if (fileFp != NULL) {                       \
        fprintf(fileFp, format, ##__VA_ARGS__); \
        fflush(fileFp);                         \
    }                                           \
}
#endif

#ifndef safeCloseFile
#define safeCloseFile(fileFp) \
{                         \
    if (fileFp != NULL) { \
        fclose(fileFp);   \
        fileFp = NULL;    \
    }                     \
}
#endif

#ifndef safeFwrite
#define safeFwrite(fileFp, ptr, sizeof_byte, nEle)  \
{                                               \
    if (fileFp != NULL) {                       \
        fwrite(ptr, sizeof_byte, nEle, fileFp); \
        fflush(fileFp);                         \
    }                                           \
}
#endif

#ifndef safeFree
#define safeFree(ptr)      \
{                      \
    if (ptr != NULL) { \
        free(ptr);     \
        ptr = NULL;    \
    }                  \
}
#endif

#ifndef Abort
#define Abort(format, ...)                                                          \
{                                                                               \
    safeFprintf(stderr, "File: %s; Line: %d; Func: %s;\n\tError: " format "\n", \
                __FILE__, __LINE__, __func__, ##__VA_ARGS__);                   \
    safeFprintf(logFile,                                                        \
                "File: %s; Line: %d; Func: %s;\n\tError: " format "\n",         \
                __FILE__, __LINE__, __func__, ##__VA_ARGS__);                   \
    exit(EXIT_FAILURE);                                                         \
}
#endif

#ifndef Info
#define Info(format, ...)                                                           \
{                                                                               \
    safeFprintf(stderr, "File: %s; Line: %d;  Func: %s;\n\tInfo: " format "\n", \
                __FILE__, __LINE__, __func__, ##__VA_ARGS__);                   \
    safeFprintf(logFile, "File: %s; Line: %d; Func: %s;\n\tInfo: " format "\n", \
                __FILE__, __LINE__, __func__, ##__VA_ARGS__);                   \
}
#endif

#ifndef getTimeStamp
#define getTimeStamp()                                       \
({                                                       \
    struct timeval sTv;                                  \
    gettimeofday(&sTv, NULL);                            \
    double tic = sTv.tv_sec + sTv.tv_usec / 1e6 - 1.6E9; \
    tic;                                                 \
})
#endif

#ifndef getTimeString
#define getTimeString(str) ({time_t ctime = time(NULL); struct tm *tmptr = gmtime(&ctime); snprintf(str,32,"%d.%d.%d_%d:%d:%d",1900 + tmptr->tm_year, 1 + tmptr->tm_mon, tmptr->tm_mday, 8 + tmptr->tm_hour, tmptr->tm_min, tmptr->tm_sec); })
#endif

#ifndef cpuMin
#define cpuMin(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef cpuMax
#define cpuMax(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef accumTime
#define accumTime(func)              \
({                               \
    double tic = getTimeStamp(); \
    func;                        \
    double toc = getTimeStamp(); \
    toc - tic;                   \
})
#endif

#ifndef isFileExisted
#define isFileExisted(filePath) (access(filePath, F_OK) == 0 ? true : false)
#endif

// create a non-existed file with read and write permission.
// if file is alread existed, it depends on truncFileFlag.
// default: do nothing and Abort the program.
// trunc: truncating the file.
// append: open the file but you can't modify original content.
#ifndef createFileReadWrite
#define createFileReadWrite(filePath)                                                        \
({                                                                                       \
    if (filePath == NULL)                                                                \
        Abort("No variable!");                                                           \
    FILE *fp = NULL;                                                                     \
    if (access(filePath, F_OK) == 0) {                                                   \
        char tstring[32];                                                                \
        getTimeString(tstring);                                                          \
        if (truncFileFlag == 0)                                                          \
            Abort("File \"%s\" exist! Operation Time %s. Exit", filePath, tstring);      \
        if (truncFileFlag == 1) {                                                        \
            Info("File \"%s\" exist! Operation Time %s. Truncated!", filePath, tstring); \
            fp = fopen(filePath, "wb+");                                                 \
        }                                                                                \
        if (truncFileFlag == 2) {                                                        \
            Info("File \"%s\" exist! Operation Time %s. Appended!", filePath, tstring);  \
            fp = fopen(filePath, "ab+");                                                 \
            fseek(fp, 0, SEEK_END);                                                      \
        }                                                                                \
    } else {                                                                             \
        fp = fopen(filePath, "wb+");                                                     \
    }                                                                                    \
    if (fp == NULL)                                                                      \
        Abort("Open File \"%s\" Failed!", filePath);                                     \
    fp;                                                                                  \
})
#endif

// open existed file with read and write permission; moving ReadWrite Position to the end.
#ifndef openExistFileReadWrite
#define openExistFileReadWrite(filePath)                 \
({                                                   \
    if (filePath == NULL)                            \
        Abort("No variable!");                       \
    if (!isFileExisted(filePath))                    \
        Abort("No File %s!", filePath);              \
    FILE *fp = fopen(filePath, "rb+");               \
    fseek(fp, 0, SEEK_END);                          \
    if (fp == NULL)                                  \
        Abort("Open File \"%s\" Failed!", filePath); \
    fp;                                              \
})
#endif

// open existed file with read only permission
#ifndef openExistFileReadOnly
#define openExistFileReadOnly(filePath)              \
({                                               \
    if (filePath == NULL)                        \
        Abort("No variable!");                   \
    if (access(filePath, F_OK) != 0)             \
        Abort("No File %s!", filePath);          \
    FILE *fp = fopen(filePath, "r");             \
    if (fp == NULL)                              \
        Abort("Open File %s Failed!", filePath); \
    fp;                                          \
})
#endif

#if defined(checkMacroFlag)
// do some checks
#define _checkMacroArgs(macroName, arg, var)                                \
{                                                                       \
    if (strstr(#arg, #var)) {                                           \
        fprintf(stderr,                                                 \
                "File: %s; Func: %s; Line: %d\n\tError: There is " #var \
                " in the argument " #arg " of the Macro: " #macroName   \
                ", which may cause side effect, please check "          \
                "carefully.\n",                                         \
                __FILE__, __func__, __LINE__);                          \
        exit(EXIT_FAILURE);                                             \
    }                                                                   \
}
#define checkMacroArgs(macroName, arg, var) _checkMacroArgs(macroName, arg, var)
#else
#define checkMacroArgs(macro, arg, var)
#endif

// column-major linear index
#ifndef mRowCol2lidx
#define mRowCol2lidx(ir, jc, ldd) ((ldd) * (jc) + ir)
#endif

// (i, alpha) to linear_idx
#ifndef aIA2lidx
#define aIA2lidx(ith, alpha) (DIM * (ith) + alpha)
#endif

// (i, alpha; i, beta) to linear_triplet_idx
#ifndef aIAB2diagTriK
#define aIAB2diagTriK(ith, alpha, beta) \
((DIM * DIM) * (ith) + ((alpha) * DIM + (beta)))
#endif

#if !defined(voigtSymRowCol2lidx) && (DIM < 4)

#if (DIM == 3)
static int _voigt2sidx[voigtDIM][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
static int _sidx2voigt[DIM][DIM] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};
#elif (DIM == 2)
static int _voigt2sidx[voigtDIM][2] = {{0, 0}, {1, 1}, {0, 1}};
static int _sidx2voigt[DIM][DIM] = {{0, 2}, {2, 1}};
#endif

#define voigt2SpaceIdx(alpha, beta, voigt)                           \
{                                                                \
    if (voigt >= voigtDIM) {                                     \
        Abort("Dimension of voigt matrix: %d!", voigtDIM);       \
    };                                                           \
    alpha = _voigt2sidx[voigt][0], beta = _voigt2sidx[voigt][1]; \
}
#define spaceIdx2voigt(alpha, beta) (_sidx2voigt[(alpha)][(beta)])

#define voigtSymRowCol2lidx(voigt1, voigt2)                               \
({                                                                    \
    int mkVar(voigtSymRowCol2lidx, min) = cpuMin((voigt1), (voigt2)); \
    int mkVar(voigtSymRowCol2lidx, max) = cpuMax((voigt1), (voigt2)); \
    ((mkVar(voigtSymRowCol2lidx, min)) +                              \
     (mkVar(voigtSymRowCol2lidx, max)) *                              \
     ((mkVar(voigtSymRowCol2lidx, max)) + 1) / 2);                \
})
#endif

typedef double uptriMat[voigtDIM];//up-trianglar part of symmetric matrix.
// uptriMat[spaceIdx2voigt(alpha,beta)]: alpha_beta component of Matrix
typedef double doubleVector[DIM];
typedef doubleVector diagMat;
typedef double *doubleVecPtr;
typedef int intVector[DIM];
typedef int *intVecPtr;
typedef struct int2 {
    int first, second;
} int2;

#if (DIM == 3)

#define uptriMatCpy(dst, src)                                     \
{                                                             \
    for (int mkVar(UMC, _i_) = 0; mkVar(UMC, _i_) < voigtDIM; \
         mkVar(UMC, _i_)++)                                   \
        dst[mkVar(UMC, _i_)] = src[mkVar(UMC, _i_)];          \
}

#define uptriMatZeros(mat)                                        \
{                                                             \
    for (int mkVar(UMZ, _i_) = 0; mkVar(UMZ, _i_) < voigtDIM; \
         mkVar(UMZ, _i_)++)                                   \
        mat[mkVar(UMZ, _i_)] = 0.0;                           \
}
#define uptriMatUnits(mat)                                               \
{                                                                    \
    uptriMatZeros(mat);                                              \
    for (int mkVar(UMU, _i_) = 0; mkVar(UMU, _i_) < DIM;             \
         mkVar(UMU, _i_)++)                                          \
        mat[spaceIdx2voigt(mkVar(UMU, _i_), mkVar(UMU, _i_))] = 1.0; \
}
#define uptriMatScale(matB, s, matA)                              \
{                                                             \
    for (int mkVar(UMS, _i_) = 0; mkVar(UMS, _i_) < voigtDIM; \
         mkVar(UMS, _i_)++)                                   \
        matB[mkVar(UMS, _i_)] = (s) * matA[mkVar(UMS, _i_)];  \
}

#define vCpy(dst, src) \
{ dst[0] = src[0], dst[1] = src[1], dst[2] = src[2]; }
#define vZeros(a) \
{ a[0] = a[1] = a[2] = 0; }
#define vOnes(a) \
{ a[0] = a[1] = a[2] = 1; }
#define vScale(b, s, a)                                          \
{ /* b[i] = s * a[i] */                                      \
    b[0] = (s) * a[0], b[1] = (s) * a[1], b[2] = (s) * a[2]; \
}
#define vAdd(c, a, b)                                               \
{ /* c[i] = a[i] + b[i] */                                      \
    c[0] = a[0] + b[0], c[1] = a[1] + b[1], c[2] = a[2] + b[2]; \
}
#define vShiftAll(a, s)                  \
{ /* a[i] = a[i] + s */              \
    a[0] += s, a[1] += s, a[2] += s; \
}
#define vSub(c, a, b)                                               \
{ /* c[i] = a[i] - b[i] */                                      \
    c[0] = a[0] - b[0], c[1] = a[1] - b[1], c[2] = a[2] - b[2]; \
}
#define sDot(a, b)                              \
({ /* dot(a,b) */                           \
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; \
})
#define vCross(c, a, b)                    \
{ /*c = cross(a,b)*/                   \
    c[0] = a[1] * b[2] - a[2] * b[1];  \
    c[1] = -a[0] * b[2] + a[2] * b[0]; \
    c[2] = a[0] * b[1] - a[1] * b[0];  \
}
#define vFloor(b, a)                                                               \
{ /*b[i] = (int)floor(a[i])*/                                                  \
    b[0] = (int)floor(a[0]), b[1] = (int)floor(a[1]), b[2] = (int)floor(a[2]); \
}
#define vScaleAdd(c, a, s, b)                                                         \
{ /* c[i] = a[i] + s * b[i] */                                                    \
    c[0] = a[0] + (s) * b[0], c[1] = a[1] + (s) * b[1], c[2] = a[2] + (s) * b[2]; \
}
#define sNormP2(a) sDot(a, a) /* np2 = dot(a,a) */
#define sNorm(a)                                      \
({ /* |a| */                                      \
    sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); \
})
#define vUnit(b, a)                             \
{ /* vec = a / |a| */                       \
    double mkVar(unit, _len_) = sNorm(a);   \
    vScale(b, 1.0 / mkVar(unit, _len_), a); \
}
#define sMaxElement(a)                                                         \
({ /* b = max{a[i]} */                                                     \
    double mkVar(maxEle, b) = (a[0] > a[1] ? a[0] : a[1]);                  \
    mkVar(maxEle, b) = (mkVar(maxEle, b) > a[2] ? mkVar(maxEle, b) : a[2]); \
    mkVar(maxEle, b);                                                       \
})
#define sMinElement(a)                                                         \
({ /* b = max{a[i]} */                                                     \
    double mkVar(minEle, b) = (a[0] < a[1] ? a[0] : a[1]);                  \
    mkVar(minEle, b) = (mkVar(minEle, b) < a[2] ? mkVar(minEle, b) : a[2]); \
    mkVar(minEle, b);                                                       \
})
#define unwrapPos(ux, x, img, hMat)                          \
{ /* ux = x + hMat * img */                              \
    ux[0] = x[0] + img[0] * hMat[spaceIdx2voigt(0, 0)] + \
    img[1] * hMat[spaceIdx2voigt(0, 1)] +        \
    img[2] * hMat[spaceIdx2voigt(0, 2)];         \
    ux[1] = x[1] + img[1] * hMat[spaceIdx2voigt(1, 1)] + \
    img[2] * hMat[spaceIdx2voigt(1, 2)];         \
    ux[2] = x[2] + img[2] * hMat[spaceIdx2voigt(2, 2)];  \
}

#define MatMulVec(mVec, hMat, vec)                      \
{ /* mVec = hMat * vec */                           \
    mVec[0] = hMat[spaceIdx2voigt(0, 0)] * vec[0] + \
    hMat[spaceIdx2voigt(0, 1)] * vec[1] + \
    hMat[spaceIdx2voigt(0, 2)] * vec[2];  \
    mVec[1] = hMat[spaceIdx2voigt(1, 1)] * vec[1] + \
    hMat[spaceIdx2voigt(1, 2)] * vec[2];  \
    mVec[2] = hMat[spaceIdx2voigt(2, 2)] * vec[2];  \
}
#define MatMulMat(c, a, b)                                      \
{ /* c = a * b; c!=a, c!=b */                               \
    c[spaceIdx2voigt(0, 0)] =                               \
    a[spaceIdx2voigt(0, 0)] * b[spaceIdx2voigt(0, 0)];  \
    c[spaceIdx2voigt(0, 1)] =                               \
    a[spaceIdx2voigt(0, 0)] * b[spaceIdx2voigt(0, 1)] + \
    a[spaceIdx2voigt(0, 1)] * b[spaceIdx2voigt(1, 1)];  \
    c[spaceIdx2voigt(0, 2)] =                               \
    a[spaceIdx2voigt(0, 0)] * b[spaceIdx2voigt(0, 2)] + \
    a[spaceIdx2voigt(0, 1)] * b[spaceIdx2voigt(1, 2)] + \
    a[spaceIdx2voigt(0, 2)] * b[spaceIdx2voigt(2, 2)];  \
    c[spaceIdx2voigt(1, 1)] =                               \
    a[spaceIdx2voigt(1, 1)] * b[spaceIdx2voigt(1, 1)];  \
    c[spaceIdx2voigt(1, 2)] =                               \
    a[spaceIdx2voigt(1, 1)] * b[spaceIdx2voigt(1, 2)] + \
    a[spaceIdx2voigt(1, 2)] * b[spaceIdx2voigt(2, 2)];  \
    c[spaceIdx2voigt(2, 2)] =                               \
    a[spaceIdx2voigt(2, 2)] * b[spaceIdx2voigt(2, 2)];  \
}
#define diagMulMat(m, diag, a)                                       \
{ /* m = diag(diag) * a; diag: diagMat, m!=a */                  \
    m[spaceIdx2voigt(0, 0)] = diag[0] * a[spaceIdx2voigt(0, 0)]; \
    m[spaceIdx2voigt(0, 1)] = diag[0] * a[spaceIdx2voigt(0, 1)]; \
    m[spaceIdx2voigt(0, 2)] = diag[0] * a[spaceIdx2voigt(0, 2)]; \
    m[spaceIdx2voigt(1, 1)] = diag[1] * a[spaceIdx2voigt(1, 1)]; \
    m[spaceIdx2voigt(1, 2)] = diag[1] * a[spaceIdx2voigt(1, 2)]; \
    m[spaceIdx2voigt(2, 2)] = diag[2] * a[spaceIdx2voigt(2, 2)]; \
}
#define getDiag(diag, m)                   \
{ /* diag = diag(m) */                 \
    diag[0] = m[spaceIdx2voigt(0, 0)]; \
    diag[1] = m[spaceIdx2voigt(1, 1)]; \
    diag[2] = m[spaceIdx2voigt(2, 2)]; \
}
#define invDiag(idiag, diag)         \
{ /* idiag[i] = 1.0 / diag[i] */ \
    idiag[0] = 1.0 / diag[0];    \
    idiag[1] = 1.0 / diag[1];    \
    idiag[2] = 1.0 / diag[2];    \
}

#if defined(__orthBox__)
#define PBC(nRij, box)                                                 \
{ /*adjust the image of starting Point*/                           \
    if (nRij[2] >= 0.5 * box->boxH[spaceIdx2voigt(2, 2)]) {        \
        nRij[2] -= box->boxH[spaceIdx2voigt(2, 2)];                \
    } else if (nRij[2] < -0.5 * box->boxH[spaceIdx2voigt(2, 2)]) { \
        nRij[2] += box->boxH[spaceIdx2voigt(2, 2)];                \
    }                                                              \
    if (nRij[1] >= 0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {        \
        nRij[1] -= box->boxH[spaceIdx2voigt(1, 1)];                \
    } else if (nRij[1] < -0.5 * box->boxH[spaceIdx2voigt(1, 1)]) { \
        nRij[1] += box->boxH[spaceIdx2voigt(1, 1)];                \
    }                                                              \
    if (nRij[0] >= 0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {        \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 0)];                \
    } else if (nRij[0] < -0.5 * box->boxH[spaceIdx2voigt(0, 0)]) { \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 0)];                \
    }                                                              \
}
#elif defined(__triBox__)
#define PBC(nRij, box)                                                   \
{ /*adjust the image of starting Point*/                             \
    if (nRij[2] >= 0.5 * (box)->boxH[spaceIdx2voigt(2, 2)]) {        \
        nRij[2] -= (box)->boxH[spaceIdx2voigt(2, 2)];                \
        nRij[1] -= (box)->boxH[spaceIdx2voigt(1, 2)];                \
        nRij[0] -= (box)->boxH[spaceIdx2voigt(0, 2)];                \
    } else if (nRij[2] < -0.5 * (box)->boxH[spaceIdx2voigt(2, 2)]) { \
        nRij[2] += (box)->boxH[spaceIdx2voigt(2, 2)];                \
        nRij[1] += (box)->boxH[spaceIdx2voigt(1, 2)];                \
        nRij[0] += (box)->boxH[spaceIdx2voigt(0, 2)];                \
    }                                                                \
    if (nRij[1] >= 0.5 * (box)->boxH[spaceIdx2voigt(1, 1)]) {        \
        nRij[1] -= (box)->boxH[spaceIdx2voigt(1, 1)];                \
        nRij[0] -= (box)->boxH[spaceIdx2voigt(0, 1)];                \
    } else if (nRij[1] < -0.5 * (box)->boxH[spaceIdx2voigt(1, 1)]) { \
        nRij[1] += (box)->boxH[spaceIdx2voigt(1, 1)];                \
        nRij[0] += (box)->boxH[spaceIdx2voigt(0, 1)];                \
    }                                                                \
    if (nRij[0] >= 0.5 * (box)->boxH[spaceIdx2voigt(0, 0)]) {        \
        nRij[0] -= (box)->boxH[spaceIdx2voigt(0, 0)];                \
    } else if (nRij[0] < -0.5 * (box)->boxH[spaceIdx2voigt(0, 0)]) { \
        nRij[0] += (box)->boxH[spaceIdx2voigt(0, 0)];                \
    }                                                                \
}
#define imgPBC(nRij, box, img)                                               \
{ /*nRij0 = Ri0 - Rj0, adjust image of j: nRij = Ri0 - (Rj0 + img*box)*/ \
    img[0] = img[1] = img[2] = 0;                                        \
    if (nRij[2] >= 0.5 * box->boxH[spaceIdx2voigt(2, 2)]) {              \
        img[2] = +1;                                                     \
        nRij[2] -= box->boxH[spaceIdx2voigt(2, 2)];                      \
        nRij[1] -= box->boxH[spaceIdx2voigt(1, 2)];                      \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 2)];                      \
    } else if (nRij[2] < -0.5 * box->boxH[spaceIdx2voigt(2, 2)]) {       \
        img[2] = -1;                                                     \
        nRij[2] += box->boxH[spaceIdx2voigt(2, 2)];                      \
        nRij[1] += box->boxH[spaceIdx2voigt(1, 2)];                      \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 2)];                      \
    }                                                                    \
    if (nRij[1] >= 0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {              \
        img[1] = +1;                                                     \
        nRij[1] -= box->boxH[spaceIdx2voigt(1, 1)];                      \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 1)];                      \
    } else if (nRij[1] < -0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {       \
        img[1] = -1;                                                     \
        nRij[1] += box->boxH[spaceIdx2voigt(1, 1)];                      \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 1)];                      \
    }                                                                    \
    if (nRij[0] >= 0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {              \
        img[0] = +1;                                                     \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 0)];                      \
    } else if (nRij[0] < -0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {       \
        img[0] = -1;                                                     \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 0)];                      \
    }                                                                    \
}
#endif

#ifndef box2BoxBounds
#define box2BoxBounds(box, xlo, xhi, ylo, yhi, zlo, zhi, yz, xz, xy)            \
{                                                                           \
    double mkVar(box2BoxBounds, tmp);                                       \
    mkVar(box2BoxBounds, tmp) =                                             \
    box->boxH[spaceIdx2voigt(0, 1)] + box->boxH[spaceIdx2voigt(0, 2)];  \
    mkVar(box2BoxBounds, tmp) =                                             \
    cpuMin(mkVar(box2BoxBounds, tmp), box->boxH[spaceIdx2voigt(0, 1)]); \
    mkVar(box2BoxBounds, tmp) =                                             \
    cpuMin(mkVar(box2BoxBounds, tmp), box->boxH[spaceIdx2voigt(0, 2)]); \
    mkVar(box2BoxBounds, tmp) = cpuMin(mkVar(box2BoxBounds, tmp), 0.0);     \
    xlo = box->cornerLo[0] + mkVar(box2BoxBounds, tmp);                     \
    mkVar(box2BoxBounds, tmp) =                                             \
    box->boxH[spaceIdx2voigt(0, 1)] + box->boxH[spaceIdx2voigt(0, 2)];  \
    mkVar(box2BoxBounds, tmp) =                                             \
    cpuMax(mkVar(box2BoxBounds, tmp), box->boxH[spaceIdx2voigt(0, 1)]); \
    mkVar(box2BoxBounds, tmp) =                                             \
    cpuMax(mkVar(box2BoxBounds, tmp), box->boxH[spaceIdx2voigt(0, 2)]); \
    mkVar(box2BoxBounds, tmp) = cpuMax(mkVar(box2BoxBounds, tmp), 0.0);     \
    xhi = box->cornerLo[0] + box->boxH[spaceIdx2voigt(0, 0)] +              \
    mkVar(box2BoxBounds, tmp);                                        \
    ylo = box->cornerLo[1] + cpuMin(0.0, box->boxH[spaceIdx2voigt(1, 2)]);  \
    yhi = box->cornerLo[1] + box->boxH[spaceIdx2voigt(1, 1)] +              \
    cpuMax(0.0, box->boxH[spaceIdx2voigt(1, 2)]);                     \
    zlo = box->cornerLo[2];                                                 \
    zhi = box->cornerLo[2] + box->boxH[spaceIdx2voigt(2, 2)];               \
    yz = box->boxH[spaceIdx2voigt(1, 2)];                                   \
    xz = box->boxH[spaceIdx2voigt(0, 2)];                                   \
    xy = box->boxH[spaceIdx2voigt(0, 1)];                                   \
}
#endif

#ifndef triBox2LeesEdwards
#define triBox2LeesEdwards(xyz, box)                              \
{                                                             \
    while (xyz[2] >= box->boxH[spaceIdx2voigt(2, 2)] / 2.0) { \
        xyz[2] -= box->boxH[spaceIdx2voigt(2, 2)];            \
        xyz[1] -= box->boxH[spaceIdx2voigt(1, 2)];            \
        xyz[0] -= box->boxH[spaceIdx2voigt(0, 2)];            \
    }                                                         \
    while (xyz[2] < -box->boxH[spaceIdx2voigt(2, 2)] / 2.0) { \
        xyz[2] += box->boxH[spaceIdx2voigt(2, 2)];            \
        xyz[1] += box->boxH[spaceIdx2voigt(1, 2)];            \
        xyz[0] += box->boxH[spaceIdx2voigt(0, 2)];            \
    }                                                         \
    while (xyz[1] >= box->boxH[spaceIdx2voigt(1, 1)] / 2.0) { \
        xyz[1] -= box->boxH[spaceIdx2voigt(1, 1)];            \
        xyz[0] -= box->boxH[spaceIdx2voigt(0, 1)];            \
    }                                                         \
    while (xyz[1] < -box->boxH[spaceIdx2voigt(1, 1)] / 2.0) { \
        xyz[1] += box->boxH[spaceIdx2voigt(1, 1)];            \
        xyz[0] += box->boxH[spaceIdx2voigt(0, 1)];            \
    }                                                         \
    while (xyz[0] >= box->boxH[spaceIdx2voigt(0, 0)] / 2.0) { \
        xyz[0] -= box->boxH[spaceIdx2voigt(0, 0)];            \
    }                                                         \
    while (xyz[0] < -box->boxH[spaceIdx2voigt(0, 0)] / 2.0) { \
        xyz[0] += box->boxH[spaceIdx2voigt(0, 0)];            \
    }                                                         \
}
#endif

#elif (DIM == 2)

#define uptriMatCpy(dst, src)                                     \
{                                                             \
    for (int mkVar(UMC, _i_) = 0; mkVar(UMC, _i_) < voigtDIM; \
         mkVar(UMC, _i_)++)                                   \
        dst[mkVar(UMC, _i_)] = src[mkVar(UMC, _i_)];          \
}

#define uptriMatZeros(mat)                                        \
{                                                             \
    for (int mkVar(UMZ, _i_) = 0; mkVar(UMZ, _i_) < voigtDIM; \
         mkVar(UMZ, _i_)++)                                   \
        mat[mkVar(UMZ, _i_)] = 0.0;                           \
}
#define uptriMatUnits(mat)                                               \
{                                                                    \
    uptriMatZeros(mat);                                              \
    for (int mkVar(UMU, _i_) = 0; mkVar(UMU, _i_) < DIM;             \
         mkVar(UMU, _i_)++)                                          \
        mat[spaceIdx2voigt(mkVar(UMU, _i_), mkVar(UMU, _i_))] = 1.0; \
}
#define uptriMatScale(matB, s, matA)                              \
{                                                             \
    for (int mkVar(UMS, _i_) = 0; mkVar(UMS, _i_) < voigtDIM; \
         mkVar(UMS, _i_)++)                                   \
        matB[mkVar(UMS, _i_)] = (s) * matA[mkVar(UMS, _i_)];  \
}

#define vCpy(dst, src) \
{ dst[0] = src[0], dst[1] = src[1]; }
#define vZeros(a) \
{ a[0] = a[1] = 0; }
#define vOnes(a) \
{ a[0] = a[1] = 1; }
#define vScale(b, s, a)                       \
{ /* b[i] = s * a[i] */                   \
    b[0] = (s) * a[0], b[1] = (s) * a[1]; \
}
#define vAdd(c, a, b)                           \
{ /* c[i] = a[i] + b[i] */                  \
    c[0] = a[0] + b[0], c[1] = a[1] + b[1]; \
}
#define vShiftAll(a, s)       \
{ /* a[i] = a[i] + s */   \
    a[0] += s, a[1] += s; \
}
#define vSub(c, a, b)                           \
{ /* c[i] = a[i] - b[i] */                  \
    c[0] = a[0] - b[0], c[1] = a[1] - b[1]; \
}
#define sDot(a, b)                \
({ /* c = dot(a,b) */         \
    a[0] * b[0] + a[1] * b[1]; \
})
#define vFloor(b, a)                                      \
{ /*b[i] = (int)floor(a[i])*/                         \
    b[0] = (int)floor(a[0]), b[1] = (int)floor(a[1]); \
}
#define vScaleAdd(c, a, s, b)                               \
{ /* c[i] = a[i] + s * b[i] */                          \
    c[0] = a[0] + (s) * b[0], c[1] = a[1] + (s) * b[1]; \
}
#define sNormP2(a) sDot(a, a) /* np2 = dot(a,a) */
#define sNorm(a)                        \
({ /* n = |a| */                    \
    sqrt(a[0] * a[0] + a[1] * a[1]); \
})
#define vUnit(b, a)                             \
{ /* vec = a / |a| */                       \
    double mkVar(unit, _len_) = sNorm(a);   \
    vScale(b, 1.0 / mkVar(unit, _len_), a); \
}
#define sMaxElement(a)              \
({ /* b = max{a[i]} */          \
    (a[0] > a[1] ? a[0] : a[1]); \
})
#define sMinElement(a)              \
({ /* b = min{a[i]} */          \
    (a[0] < a[1] ? a[0] : a[1]); \
})
#define unwrapPos(ux, x, img, hMat)                          \
{ /* ux = x + hMat * img */                              \
    ux[0] = x[0] + img[0] * hMat[spaceIdx2voigt(0, 0)] + \
    img[1] * hMat[spaceIdx2voigt(0, 1)];         \
    ux[1] = x[1] + img[1] * hMat[spaceIdx2voigt(1, 1)];  \
}
#define MatMulVec(mVec, hMat, vec)                      \
{ /* mVec = hMat * vec */                           \
    mVec[0] = hMat[spaceIdx2voigt(0, 0)] * vec[0] + \
    hMat[spaceIdx2voigt(0, 1)] * vec[1];  \
    mVec[1] = hMat[spaceIdx2voigt(1, 1)] * vec[1];  \
}
#define MatMulMat(c, a, b)                                      \
{ /* c = a * b; c!=a, c!=b */                               \
    c[spaceIdx2voigt(0, 0)] =                               \
    a[spaceIdx2voigt(0, 0)] * b[spaceIdx2voigt(0, 0)];  \
    c[spaceIdx2voigt(0, 1)] =                               \
    a[spaceIdx2voigt(0, 0)] * b[spaceIdx2voigt(0, 1)] + \
    a[spaceIdx2voigt(0, 1)] * b[spaceIdx2voigt(1, 1)];  \
    c[spaceIdx2voigt(1, 1)] =                               \
    a[spaceIdx2voigt(1, 1)] * b[spaceIdx2voigt(1, 1)];  \
}
#define diagMulMat(m, diag, a)                                       \
{ /* m = diag(diag) * a; diag: diagMat, m!=a */                  \
    m[spaceIdx2voigt(0, 0)] = diag[0] * a[spaceIdx2voigt(0, 0)]; \
    m[spaceIdx2voigt(0, 1)] = diag[0] * a[spaceIdx2voigt(0, 1)]; \
    m[spaceIdx2voigt(1, 1)] = diag[1] * a[spaceIdx2voigt(1, 1)]; \
}
#define getDiag(diag, m)                   \
{ /* diag = diag(m) */                 \
    diag[0] = m[spaceIdx2voigt(0, 0)]; \
    diag[1] = m[spaceIdx2voigt(1, 1)]; \
}
#define invDiag(idiag, diag)         \
{ /* idiag[i] = 1.0 / diag[i] */ \
    idiag[0] = 1.0 / diag[0];    \
    idiag[1] = 1.0 / diag[1];    \
}

#if defined(__orthBox__)
#define PBC(nRij, box)                                                 \
{ /*adjust the image of starting Point*/                           \
    if (nRij[1] >= 0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {        \
        nRij[1] -= box->boxH[spaceIdx2voigt(1, 1)];                \
    } else if (nRij[1] < -0.5 * box->boxH[spaceIdx2voigt(1, 1)]) { \
        nRij[1] += box->boxH[spaceIdx2voigt(1, 1)];                \
    }                                                              \
    if (nRij[0] >= 0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {        \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 0)];                \
    } else if (nRij[0] < -0.5 * box->boxH[spaceIdx2voigt(0, 0)]) { \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 0)];                \
    }                                                              \
}
#elif defined(__triBox__)
#define PBC(nRij, box)                                                 \
{ /*adjust the image of starting Point*/                           \
    if (nRij[1] >= 0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {        \
        nRij[1] -= box->boxH[spaceIdx2voigt(1, 1)];                \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 1)];                \
    } else if (nRij[1] < -0.5 * box->boxH[spaceIdx2voigt(1, 1)]) { \
        nRij[1] += box->boxH[spaceIdx2voigt(1, 1)];                \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 1)];                \
    }                                                              \
    if (nRij[0] >= 0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {        \
        nRij[0] -= box->boxH[spaceIdx2voigt(0, 0)];                \
    } else if (nRij[0] < -0.5 * box->boxH[spaceIdx2voigt(0, 0)]) { \
        nRij[0] += box->boxH[spaceIdx2voigt(0, 0)];                \
    }                                                              \
}
#endif

#ifndef box2BoxBounds
#define box2BoxBounds(box, xlo, xhi, ylo, yhi, zlo, zhi, yz, xz, xy)              \
{                                                                             \
    double mkVar(box2BoxBounds, tmp);                                         \
    mkVar(box2BoxBounds, tmp) = cpuMin(box->boxH[spaceIdx2voigt(0, 1)], 0.0); \
    xlo = box->cornerLo[0] + mkVar(box2BoxBounds, tmp);                       \
    mkVar(box2BoxBounds, tmp) = cpuMax(box->boxH[spaceIdx2voigt(0, 1)], 0.0); \
    xhi = box->cornerLo[0] + box->boxH[spaceIdx2voigt(0, 0)] +                \
    mkVar(box2BoxBounds, tmp);                                          \
    ylo = box->cornerLo[1];                                                   \
    yhi = box->cornerLo[1] + box->boxH[spaceIdx2voigt(1, 1)];                 \
    zlo = 0;                                                                  \
    zhi = 0;                                                                  \
    yz = 0;                                                                   \
    xz = 0;                                                                   \
    xy = box->boxH[spaceIdx2voigt(0, 1)];                                     \
}
#endif

#ifndef triBox2LeesEdwards
#define triBox2LeesEdwards(xyz, box)                              \
{                                                             \
    while (xyz[1] >= box->boxH[spaceIdx2voigt(1, 1)] / 2.0) { \
        xyz[1] -= box->boxH[spaceIdx2voigt(1, 1)];            \
        xyz[0] -= box->boxH[spaceIdx2voigt(0, 1)];            \
    }                                                         \
    while (xyz[1] < -box->boxH[spaceIdx2voigt(1, 1)] / 2.0) { \
        xyz[1] += box->boxH[spaceIdx2voigt(1, 1)];            \
        xyz[0] += box->boxH[spaceIdx2voigt(0, 1)];            \
    }                                                         \
    while (xyz[0] >= box->boxH[spaceIdx2voigt(0, 0)] / 2.0) { \
        xyz[0] -= box->boxH[spaceIdx2voigt(0, 0)];            \
    }                                                         \
    while (xyz[0] < -box->boxH[spaceIdx2voigt(0, 0)] / 2.0) { \
        xyz[0] += box->boxH[spaceIdx2voigt(0, 0)];            \
    }                                                         \
}
#endif

#else

#if defined(checkMacroFlag)
#define spaceIdx2voigt(alpha, beta)            \
({                                         \
    if ((alpha) < 0 || (alpha) >= DIM)     \
        exit(EXIT_FAILURE);                \
    if ((beta) < 0 || (beta) >= DIM)       \
        exit(EXIT_FAILURE);                \
    if ((alpha) > (beta))                  \
        exit(EXIT_FAILURE);                \
    ((alpha) + (beta) * ((beta) + 1) / 2); \
})
#else
#define spaceIdx2voigt(alpha, beta) ((alpha) + (beta) * ((beta) + 1) / 2)
#endif

#define uptriMatCpy(dst, src)                                     \
{                                                             \
    for (int mkVar(UMC, _i_) = 0; mkVar(UMC, _i_) < voigtDIM; \
         mkVar(UMC, _i_)++)                                   \
        dst[mkVar(UMC, _i_)] = src[mkVar(UMC, _i_)];          \
}

#define uptriMatZeros(mat)                                        \
{                                                             \
    for (int mkVar(UMZ, _i_) = 0; mkVar(UMZ, _i_) < voigtDIM; \
         mkVar(UMZ, _i_)++)                                   \
        mat[mkVar(UMZ, _i_)] = 0.0;                           \
}
#define uptriMatUnits(mat)                                               \
{                                                                    \
    uptriMatZeros(mat);                                              \
    for (int mkVar(UMU, _i_) = 0; mkVar(UMU, _i_) < DIM;             \
         mkVar(UMU, _i_)++)                                          \
        mat[spaceIdx2voigt(mkVar(UMU, _i_), mkVar(UMU, _i_))] = 1.0; \
}
#define uptriMatScale(matB, s, matA)                              \
{                                                             \
    for (int mkVar(UMS, _i_) = 0; mkVar(UMS, _i_) < voigtDIM; \
         mkVar(UMS, _i_)++)                                   \
        matB[mkVar(UMS, _i_)] = (s) * matA[mkVar(UMS, _i_)];  \
}

#define vCpy(dst, src)                                           \
{ /* dst = src */                                            \
    checkMacroArgs(cpyVec, dst, mkVar(cpy, _ijk_));          \
    checkMacroArgs(cpyVec, src, mkVar(cpy, _ijk_));          \
    for (int mkVar(cpy, _ijk_) = 0; mkVar(cpy, _ijk_) < DIM; \
         mkVar(cpy, _ijk_)++)                                \
        dst[mkVar(cpy, _ijk_)] = src[mkVar(cpy, _ijk_)];     \
}
#define vZeros(a)                                                    \
{ /* a[i] = 0 */                                                 \
    checkMacroArgs(zerosVec, a, mkVar(zeros, _ijk_));            \
    for (int mkVar(zeros, _ijk_) = 0; mkVar(zeros, _ijk_) < DIM; \
         mkVar(zeros, _ijk_)++)                                  \
        a[mkVar(zeros, _ijk_)] = 0.0;                            \
}
#define vOnes(a)                                                   \
{ /* a[i] = 1 */                                               \
    checkMacroArgs(onesVec, a, mkVar(ones, _ijk_));            \
    for (int mkVar(ones, _ijk_) = 0; mkVar(ones, _ijk_) < DIM; \
         mkVar(ones, _ijk_)++)                                 \
        a[mkVar(ones, _ijk_)] = 1.0;                           \
}
#define vScale(b, s, a)                                              \
{ /* b[i] = s * a[i] */                                          \
    checkMacroArgs(scaleVec, b, mkVar(scale, _ijk_));            \
    checkMacroArgs(scaleVec, s, mkVar(scale, _ijk_));            \
    checkMacroArgs(scaleVec, a, mkVar(scale, _ijk_));            \
    for (int mkVar(scale, _ijk_) = 0; mkVar(scale, _ijk_) < DIM; \
         mkVar(scale, _ijk_)++)                                  \
        b[mkVar(scale, _ijk_)] = (s) * a[mkVar(scale, _ijk_)];   \
}
#define vAdd(c, a, b)                                                           \
{ /* c[i] = a[i] + b[i] */                                                  \
    checkMacroArgs(vecAdd, c, mkVar(add, _ijk_));                           \
    checkMacroArgs(vecAdd, a, mkVar(add, _ijk_));                           \
    checkMacroArgs(vecAdd, b, mkVar(add, _ijk_));                           \
    for (int mkVar(add, _ijk_) = 0; mkVar(add, _ijk_) < DIM;                \
         mkVar(add, _ijk_)++)                                               \
        c[mkVar(add, _ijk_)] = a[mkVar(add, _ijk_)] + b[mkVar(add, _ijk_)]; \
}
#define vShiftAll(a, s)                                            \
{ /* a[i] = a[i] + s */                                        \
    checkMacroArgs(vShiftAll, a, mkVar(shift, _ijk));          \
    checkMacroArgs(vShiftAll, s, mkVar(shift, _ijk));          \
    for (int mkVar(shift, _ijk) = 0; mkVar(shift, _ijk) < DIM; \
         mkVar(shift, _ijk)++)                                 \
        a[mkVar(shift, _ijk)] += (s);                          \
}
#define vSub(c, a, b)                                                           \
{ /* c[i] = a[i] - b[i] */                                                  \
    checkMacroArgs(vecSub, c, mkVar(sub, _ijk_));                           \
    checkMacroArgs(vecSub, a, mkVar(sub, _ijk_));                           \
    checkMacroArgs(vecSub, b, mkVar(sub, _ijk_));                           \
    for (int mkVar(sub, _ijk_) = 0; mkVar(sub, _ijk_) < DIM;                \
         mkVar(sub, _ijk_)++)                                               \
        c[mkVar(sub, _ijk_)] = a[mkVar(sub, _ijk_)] - b[mkVar(sub, _ijk_)]; \
}
#define sDot(a, b)                                                           \
({ /* c = dot(a,b) */                                                    \
    checkMacroArgs(sDot, a, mkVar(dot, _ijk_));                           \
    checkMacroArgs(sDot, b, mkVar(dot, _ijk_));                           \
    checkMacroArgs(sDot, a, mkVar(dot, _dot_));                           \
    checkMacroArgs(sDot, b, mkVar(dot, _dot_));                           \
    double mkVar(dot, _dot_) = 0.0;                                       \
    for (int mkVar(dot, _ijk_) = 0; mkVar(dot, _ijk_) < DIM;              \
         mkVar(dot, _ijk_)++)                                             \
        mkVar(dot, _dot_) += a[mkVar(dot, _ijk_)] * b[mkVar(dot, _ijk_)]; \
    mkVar(dot, _dot_);                                                    \
})
#define vFloor(b, a)                                                     \
{ /*b[i] = (int)floor(a[i])*/                                        \
    checkMacroArgs(vFloor, a, mkVar(floor, _ijk_));                  \
    checkMacroArgs(vFloor, b, mkVar(floor, _ijk_));                  \
    for (int mkVar(floor, _ijk_) = 0; mkVar(floor, _ijk_) < DIM;     \
         mkVar(floor, _ijk_)++)                                      \
        b[mkVar(floor, _ijk_)] = (int)floor(a[mkVar(floor, _ijk_)]); \
}
#define vScaleAdd(c, a, s, b)                                                \
{ /* c[i] = a[i] + s * b[i] */                                           \
    checkMacroArgs(vecScaleAdd, c, mkVar(scaleAdd, _ijk_));              \
    checkMacroArgs(vecScaleAdd, a, mkVar(scaleAdd, _ijk_));              \
    checkMacroArgs(vecScaleAdd, s, mkVar(scaleAdd, _ijk_));              \
    checkMacroArgs(vecScaleAdd, b, mkVar(scaleAdd, _ijk_));              \
    for (int mkVar(scaleAdd, _ijk_) = 0; mkVar(scaleAdd, _ijk_) < DIM;   \
         mkVar(scaleAdd, _ijk_)++)                                       \
        c[mkVar(scaleAdd, _ijk_)] =                                      \
        a[mkVar(scaleAdd, _ijk_)] + (s) * b[mkVar(scaleAdd, _ijk_)]; \
}
#define sNormP2(a) sDot(a, a) /* np2 = dot(a,a) */
#define sNorm(a)         \
({ /* n = |a| */     \
    sqrt(sNormP2(a)); \
})
#define vUnit(b, a)                             \
{ /* vec = a / |a| */                       \
    double mkVar(unit, _len_) = sNorm(a);   \
    vScale(b, 1.0 / mkVar(unit, _len_), a); \
}
#define sMaxElement(a)                                                            \
({ /* b = max{a[i]} */                                                        \
    checkMacroArgs(maxElement, a, mkVar(maxEle, _max_));                       \
    checkMacroArgs(maxElement, a, mkVar(maxEle, _ijk_));                       \
    double mkVar(maxEle, _max_) = a[0];                                        \
    for (int mkVar(maxEle, _ijk_) = 1; mkVar(maxEle, _ijk_) < DIM;             \
         mkVar(maxEle, _ijk_)++)                                               \
        mkVar(maxEle, _max_) = (mkVar(maxEle, _max_) > a[mkVar(maxEle, _ijk_)] \
                                ? mkVar(maxEle, _max_)                     \
                                : a[mkVar(maxEle, _ijk_)]);                \
    mkVar(maxEle, _max_);                                                      \
})
#define sMinElement(a)                                                            \
({ /* b = min{a[i]} */                                                        \
    checkMacroArgs(minElement, a, mkVar(minEle, _min_));                       \
    checkMacroArgs(minElement, a, mkVar(minEle, _ijk_));                       \
    double mkVar(minEle, _min_) = a[0];                                        \
    for (int mkVar(minEle, _ijk_) = 1; mkVar(minEle, _ijk_) < DIM;             \
         mkVar(minEle, _ijk_)++)                                               \
        mkVar(minEle, _min_) = (mkVar(minEle, _min_) > a[mkVar(minEle, _ijk_)] \
                                ? mkVar(minEle, _min_)                     \
                                : a[mkVar(minEle, _ijk_)]);                \
    mkVar(minEle, _min_);                                                      \
})
#define unwrapPos(ux, x, img, hMat)                                                   \
{ /* ux = x + hMat * img */                                                       \
    checkMacroArgs(unwrapPos, ux, mkVar(unwrap, _ijk_));                          \
    checkMacroArgs(unwrapPos, x, mkVar(unwrap, _ijk_));                           \
    checkMacroArgs(unwrapPos, img, mkVar(unwrap, _ijk_));                         \
    checkMacroArgs(unwrapPos, hMat, mkVar(unwrap, _ijk_));                        \
    checkMacroArgs(unwrapPos, ux, mkVar(unwrap, _kji_));                          \
    checkMacroArgs(unwrapPos, x, mkVar(unwrap, _kji_));                           \
    checkMacroArgs(unwrapPos, img, mkVar(unwrap, _kji_));                         \
    checkMacroArgs(unwrapPos, hMat, mkVar(unwrap, _kji_));                        \
    for (int mkVar(unwrap, _ijk_) = 0; mkVar(unwrap, _ijk_) < DIM;                \
         mkVar(unwrap, _ijk_)++) {                                                \
        ux[mkVar(unwrap, _ijk_)] = x[mkVar(unwrap, _ijk_)];                       \
        for (int mkVar(unwrap, _kji_) = mkVar(unwrap, _ijk_);                     \
             mkVar(unwrap, _kji_) < DIM; mkVar(unwrap, _kji_)++) {                \
            ux[mkVar(unwrap, _ijk_)] +=                                           \
            img[mkVar(unwrap, _kji_)] *                                       \
            hMat[spaceIdx2voigt(mkVar(unwrap, _ijk_), mkVar(unwrap, _kji_))]; \
        }                                                                         \
    }                                                                             \
}
#define MatMulVec(mVec, hMat, vec)                                                   \
{ /* mVec = hMat * vec */                                                        \
    checkMacroArgs(MatMulVec, mVec, mkVar(mmv, _ijk_));                          \
    checkMacroArgs(MatMulVec, vec, mkVar(mmv, _ijk_));                           \
    checkMacroArgs(MatMulVec, hMat, mkVar(mmv, _ijk_));                          \
    checkMacroArgs(MatMulVec, mVec, mkVar(mmv, _kji_));                          \
    checkMacroArgs(MatMulVec, vec, mkVar(mmv, _kji_));                           \
    checkMacroArgs(MatMulVec, hMat, mkVar(mmv, _kji_));                          \
    checkMacroArgs(MatMulVec, mVec, mkVar(mmv, _tmp_));                          \
    checkMacroArgs(MatMulVec, vec, mkVar(mmv, _tmp_));                           \
    checkMacroArgs(MatMulVec, hMat, mkVar(mmv, _tmp_));                          \
    doubleVector mkVar(mmv, _tmp_);                                              \
    vZeros(mkVar(mmv, _tmp_));                                                   \
    for (int mkVar(mmv, _ijk_) = 0; mkVar(mmv, _ijk_) < DIM;                     \
         mkVar(mmv, _ijk_)++) {                                                  \
        for (int mkVar(mmv, _kji_) = mkVar(mmv, _ijk_); mkVar(mmv, _kji_) < DIM; \
             mkVar(mmv, _kji_)++) {                                              \
            mkVar(mmv, _tmp_)[mkVar(mmv, _ijk_)] +=                              \
            hMat[spaceIdx2voigt(mkVar(mmv, _ijk_), mkVar(mmv, _kji_))] *     \
            vec[mkVar(mmv, _kji_)];                                          \
        }                                                                        \
    }                                                                            \
    vCpy(mVec, mkVar(mmv, _tmp_));                                               \
}
#define MatMulMat(c, a, b)                                                           \
{ /* c = a * b; c!=a, c!=b */                                                    \
    checkMacroArgs(MatMulMat, c, mkVar(mmm, _ijk_));                             \
    checkMacroArgs(MatMulMat, a, mkVar(mmm, _ijk_));                             \
    checkMacroArgs(MatMulMat, b, mkVar(mmm, _ijk_));                             \
    checkMacroArgs(MatMulMat, c, mkVar(mmm, _kji_));                             \
    checkMacroArgs(MatMulMat, a, mkVar(mmm, _kji_));                             \
    checkMacroArgs(MatMulMat, b, mkVar(mmm, _kji_));                             \
    checkMacroArgs(MatMulMat, a, mkVar(mmm, _kkk_));                             \
    checkMacroArgs(MatMulMat, b, mkVar(mmm, _kkk_));                             \
    for (int mkVar(mmm, _ijk_) = 0; mkVar(mmm, _ijk_) < DIM;                     \
         mkVar(mmm, _ijk_)++) {                                                  \
        for (int mkVar(mmm, _kji_) = mkVar(mmm, _ijk_); mkVar(mmm, _kji_) < DIM; \
             mkVar(mmm, _kji_)++) {                                              \
            c[spaceIdx2voigt(mkVar(mmm, _ijk_), mkVar(mmm, _kji_))] = 0;         \
            for (int mkVar(mmm, _kkk_) = mkVar(mmm, _ijk_);                      \
                 mkVar(mmm, _kkk_) <= mkVar(mmm, _kji_); mkVar(mmm, _kkk_)++)    \
                c[spaceIdx2voigt(mkVar(mmm, _ijk_), mkVar(mmm, _kji_))] +=       \
                a[spaceIdx2voigt(mkVar(mmm, _ijk_), mkVar(mmm, _kkk_))] *    \
                b[spaceIdx2voigt(mkVar(mmm, _kkk_), mkVar(mmm, _kji_))];     \
        }                                                                        \
    }                                                                            \
}
#define diagMulMat(m, diag, a)                                                       \
{ /* m = diag(diag) * a; diag: diagMat, m!=a */                                  \
    checkMacroArgs(diagMulMat, m, mkVar(dmm, _ijk_));                            \
    checkMacroArgs(diagMulMat, diag, mkVar(dmm, _ijk_));                         \
    checkMacroArgs(diagMulMat, a, mkVar(dmm, _ijk_));                            \
    checkMacroArgs(diagMulMat, m, mkVar(dmm, _kji_));                            \
    checkMacroArgs(diagMulMat, diag, mkVar(dmm, _kji_));                         \
    checkMacroArgs(diagMulMat, a, mkVar(dmm, _kji_));                            \
    checkMacroArgs(diagMulMat, m, a);                                            \
    for (int mkVar(dmm, _ijk_) = 0; mkVar(dmm, _ijk_) < DIM;                     \
         mkVar(dmm, _ijk_)++)                                                    \
        for (int mkVar(dmm, _kji_) = mkVar(dmm, _ijk_); mkVar(dmm, _kji_) < DIM; \
             mkVar(dmm, _kji_)++)                                                \
            m[spaceIdx2voigt(mkVar(dmm, _ijk_), mkVar(dmm, _kji_))] =            \
            diag[mkVar(dmm, _ijk_)] *                                        \
            a[spaceIdx2voigt(mkVar(dmm, _ijk_), mkVar(dmm, _kji_))];         \
}
#define getDiag(diag, m)                                                 \
{ /* diag = diag(m) */                                               \
    checkMacroArgs(getDiag, m, mkVar(gdg, _ijk_));                   \
    checkMacroArgs(getDiag, diag, mkVar(gdg, _ijk_));                \
    for (int mkVar(gdg, _ijk_) = 0; mkVar(gdg, _ijk_) < DIM;         \
         mkVar(gdg, _ijk_)++)                                        \
        diag[mkVar(gdg, _ijk_)] =                                    \
        m[spaceIdx2voigt(mkVar(gdg, _ijk_), mkVar(gdg, _ijk_))]; \
}
#define invDiag(idiag, diag)                                          \
{ /* idiag[i] = 1.0 / diag[i] */                                  \
    checkMacroArgs(invDiag, idiag, mkVar(idg, _ijk_));            \
    checkMacroArgs(invDiag, diag, mkVar(idg, _ijk_));             \
    for (int mkVar(idg, _ijk_) = 0; mkVar(idg, _ijk_) < DIM;      \
         mkVar(idg, _ijk_)++)                                     \
        idiag[mkVar(idg, _ijk_)] = 1.0 / diag[mkVar(idg, _ijk_)]; \
}

#if defined(__triBox__)
#define PBC(nRij, box)                                                        \
{ /*adjust the image of starting Point*/                                  \
    checkMacroArgs(PBC, nRij, mkVar(PBC, _ijk_));                         \
    for (int mkVar(PBC, _ijk_) = DIM - 1; mkVar(PBC, _ijk_) >= 0;         \
         mkVar(PBC, _ijk_)--) {                                           \
        if (nRij[mkVar(PBC, _ijk_)] >=                                    \
            0.5 * box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_),             \
                                           mkVar(PBC, _ijk_))]) {         \
            vSub(nRij, nRij, box->boxEdge[mkVar(PBC, _ijk_)]);            \
        } else if (nRij[mkVar(PBC, _ijk_)] <                              \
                   -0.5 * box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_),     \
                                                   mkVar(PBC, _ijk_))]) { \
            vAdd(nRij, nRij, box->boxEdge[mkVar(PBC, _ijk_)]);            \
        }                                                                 \
    }                                                                     \
}
#elif defined(__orthBox__)
#define PBC(nRij, box)                                                               \
{ /*adjust the image of starting Point*/                                         \
    checkMacroArgs(PBC, nRij, mkVar(PBC, _ijk_));                                \
    for (int mkVar(PBC, _ijk_) = DIM - 1; mkVar(PBC, _ijk_) >= 0;                \
         mkVar(PBC, _ijk_)--) {                                                  \
        if (nRij[mkVar(PBC, _ijk_)] >=                                           \
            0.5 * box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_),                    \
                                           mkVar(PBC, _ijk_))]) {                \
            nRij[mkVar(PBC, _ijk_)] -=                                           \
            box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_), mkVar(PBC, _ijk_))]; \
        } else if (nRij[mkVar(PBC, _ijk_)] <                                     \
                   -0.5 * box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_),            \
                                                   mkVar(PBC, _ijk_))]) {        \
            nRij[mkVar(PBC, _ijk_)] +=                                           \
            box->boxH[spaceIdx2voigt(mkVar(PBC, _ijk_), mkVar(PBC, _ijk_))]; \
        }                                                                        \
    }                                                                            \
}
#endif

#endif

#ifdef __cplusplus
extern "C" {
#endif

double rndUniform(void);  // random number in [0,1];
double rndStdNorm(void);  // random from N(0,1);
bool isEmpty(char *str, int maxlen);

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
cmdArg *findVariable(Variable *var, char *name);  // return its pointer
void addVariable(Variable *var, char *inputCmd);
void delVariable(Variable *var, char *delCmd);

typedef struct Toolkit {
    int nToolkit, maxToolkit;
    char **toolkitName;
    void **toolkit;
    void **funcPtrWriteConf;
} Toolkit;
int findToolkit(Toolkit *toolkit, char *name);
void addToolkit(Toolkit *toolkit, void *tool, void *funcPtr, char *name);  // void *tool is the allocated memory
void delToolkit(Toolkit *toolkit, char *name);                             // allocated memory "void *tool" is freed.

#ifdef __Linux__
void dumpSourceFile(void);
#endif
extern volatile int __nCatchSignal__;  // number of catched signal (delivered by user).
int setSignalHandler(void);

// insitu sort
void heapSortInc(int *arr, int len);                     // sort int
int sortDoubleVector0Inc(doubleVector *base, int nVec);  // sort doubleVector[0]

void exchange_doubleVector(doubleVector *xyz, doubleVector *buffer, int *oid2nid, int nAtom);
void exchange_intVector(intVector *img, intVector *buffer, int *oid2nid, int nAtom);
void exchange_double(double *radius, double *buffer, int *oid2nid, int nAtom);
void exchange_int(int *type, int *buffer, int *oid2nid, int nAtom);

#if (DIM == 2 || DIM == 3)
void toMandelTensor(double *cVoigt, double *cMandel);
#endif

#ifdef __cplusplus
}
#endif

#endif
