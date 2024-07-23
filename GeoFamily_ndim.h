#ifndef __MsgDigAlg5_ndim_dengPan_h__
#define __MsgDigAlg5_ndim_dengPan_h__

#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"

// MD5 Message-Digest Algorithm: adopted from
// https://www.ietf.org/rfc/rfc1321.txt

typedef unsigned char md5FingerPrint[16];

typedef struct GeoFamily {
    //contactInfo *cinfo;
    fullNebrList fNebr;
    int *contactNetwork;
    int maxCNsize, nCNsize;
    
    //int nGeoFamily, maxGeoFamily;
    md5FingerPrint fingerPrint;
    
    //int nConf, maxConf, *geoFamilyIndx;
    bool excludeRatt;
} GeoFamily;

#ifdef __cplusplus
extern "C" {
#endif

GeoFamily *getGeoFamily(Update *update);
GeoFamily *addGeoFamily(Box *box, Particle *particle, Update *update);
int assignGeoFamily(Box *box, Particle *particle, Update *update);
int delGeoFamily(Update *update);

#ifdef __cplusplus
}
#endif

#endif
