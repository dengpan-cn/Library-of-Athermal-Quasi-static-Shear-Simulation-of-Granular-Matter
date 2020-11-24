
#if false
#ifndef __SUBFUNC__
#define __SUBFUNC__

#include "StructSim.h"
#include "VectorMath.h"

//===============================================================================
void calcVolInfo(Box *box, Particle *particle, Thermo *thermo, Update *update,
                 Variable *var) {
  if (thermo->isInit) return;

  double vol = 0;
  double minScale = 1E10, maxScale = 0.0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    vol += VolUnitSphere * pow(particle->diameterScale[iatom], 3);
    double tmp = particle->diameterScale[iatom];
    minScale = (minScale < tmp ? minScale : tmp);
    maxScale = (maxScale > tmp ? maxScale : tmp);
  }
  update->nebrList.minDiameterScale = minScale;
  update->nebrList.maxDiameterScale = maxScale;

  thermo->volFrac = vol * pow(particle->meanDiameter * 0.5, 3) / box->volume;
  thermo->isInit = true;
}
void setBoxPara(Box *box, Particle *particle, Thermo *thermo, Update *update,
                Variable *var) {
  if (box->boxHvoigt.h0 <= 0 || box->boxHvoigt.h1 <= 0 ||
      box->boxHvoigt.h2 <= 0)
    Abort("Not rihgt hand basis!");

  double3 negHalf = make_double3(-0.5, -0.5, -0.5);
  vecHvoigtMulVec(box->boxLo, box->boxHvoigt, negHalf);
  double3 tmp3 =
      make_double3(box->boxHvoigt.h0, box->boxHvoigt.h1, box->boxHvoigt.h2);
  vecAdd(box->boxHi, tmp3, box->boxLo);

  box->volume = box->boxHvoigt.h0 * box->boxHvoigt.h1 * box->boxHvoigt.h2;

  box->invBoxHvoigt.h0 = 1.0 / box->boxHvoigt.h0;
  box->invBoxHvoigt.h1 = 1.0 / box->boxHvoigt.h1;
  box->invBoxHvoigt.h2 = 1.0 / box->boxHvoigt.h2;

  box->invBoxHvoigt.h3 =
      -box->boxHvoigt.h3 / (box->boxHvoigt.h1 * box->boxHvoigt.h2);
  box->invBoxHvoigt.h4 =
      box->boxHvoigt.h3 * box->boxHvoigt.h5 /
          (box->boxHvoigt.h0 * box->boxHvoigt.h1 * box->boxHvoigt.h2) -
      box->boxHvoigt.h4 / (box->boxHvoigt.h0 * box->boxHvoigt.h2);
  box->invBoxHvoigt.h5 =
      -box->boxHvoigt.h5 / (box->boxHvoigt.h0 * box->boxHvoigt.h1);
}
void setBasicInfo(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  calcVolInfo(box, particle, thermo, update, var);

  thermo->massUnits = 1.0;
  thermo->energyUnits = 1.0;
  thermo->distanceUnits = particle->meanDiameter;
  thermo->timeUnits =
      sqrt(thermo->massUnits / thermo->energyUnits) * thermo->distanceUnits;
  thermo->forceUnits = thermo->energyUnits / thermo->distanceUnits;
  thermo->velocityUnits = sqrt(thermo->energyUnits / thermo->massUnits);
  thermo->pressureUnits = thermo->energyUnits / pow(thermo->distanceUnits, DIM);
  thermo->volumeUnits = pow(thermo->distanceUnits, DIM);
}
void adjustImg(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
#ifdef __orthBox__
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 *posXyz = &particle->xyz[iatom];
    int3 *imgXyz = &particle->img[iatom];

    // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5) = invH * (r -
    // boxLo), boxLo = invH * (-0.5, -0.5, -0.5)
    double lamda_x = box->invBoxHvoigt.h0 * posXyz->x + 0.5;
    double lamda_y = box->invBoxHvoigt.h1 * posXyz->y + 0.5;
    double lamda_z = box->invBoxHvoigt.h2 * posXyz->z + 0.5;
    double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
           shift_z = floor(lamda_z);
    imgXyz->x += (int)shift_x;
    imgXyz->y += (int)shift_y;
    imgXyz->z += (int)shift_z;

    double delta_x = box->boxHvoigt.h0 * shift_x;
    double delta_y = box->boxHvoigt.h1 * shift_y;
    double delta_z = box->boxHvoigt.h2 * shift_z;
    posXyz->x -= delta_x;
    posXyz->y -= delta_y;
    posXyz->z -= delta_z;
  }
#endif
#ifdef __mono_XZ_Box__
  Hvoigt6 boxHvoigtHold = box->boxHvoigt;
  double xtiltmax = 0.505 * box->boxHvoigt.h0;
  bool isFlip = false;
  if (box->boxHvoigt.h4 < -xtiltmax) {
    box->boxHvoigt.h4 += box->boxHvoigt.h0;
    isFlip = true;
  } else if (box->boxHvoigt.h4 > xtiltmax) {
    box->boxHvoigt.h4 -= box->boxHvoigt.h0;
    isFlip = true;
  }
  if (isFlip) {
    Info("Flip box: (%g %g %g %g %g %g) ==> (%g %g %g %g %g %g)",
         boxHvoigtHold.h0, boxHvoigtHold.h1, boxHvoigtHold.h2, boxHvoigtHold.h3,
         boxHvoigtHold.h4, boxHvoigtHold.h5, box->boxHvoigt.h0,
         box->boxHvoigt.h1, box->boxHvoigt.h2, box->boxHvoigt.h3,
         box->boxHvoigt.h4, box->boxHvoigt.h5);

    setBoxPara(box, particle, thermo, update, var);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->xyz[iatom].x += boxHvoigtHold.h0 * particle->img[iatom].x +
                                boxHvoigtHold.h4 * particle->img[iatom].z;
      particle->xyz[iatom].y += boxHvoigtHold.h1 * particle->img[iatom].y;
      particle->xyz[iatom].z += boxHvoigtHold.h2 * particle->img[iatom].z;
    }
    memset(particle->img, '\0', particle->nAtom * sizeof(int3));
  }

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 *posXyz = &particle->xyz[iatom];
    int3 *imgXyz = &particle->img[iatom];

    // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5) = invH * (r -
    // boxLo), boxLo = invH * (-0.5, -0.5, -0.5)
    double lamda_x = box->invBoxHvoigt.h0 * posXyz->x +
                     box->invBoxHvoigt.h4 * posXyz->z + 0.5;
    double lamda_y = box->invBoxHvoigt.h1 * posXyz->y + 0.5;
    double lamda_z = box->invBoxHvoigt.h2 * posXyz->z + 0.5;
    double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
           shift_z = floor(lamda_z);
    imgXyz->x += (int)shift_x;
    imgXyz->y += (int)shift_y;
    imgXyz->z += (int)shift_z;

    double delta_x = box->boxHvoigt.h0 * shift_x + box->boxHvoigt.h4 * shift_z;
    double delta_y = box->boxHvoigt.h1 * shift_y;
    double delta_z = box->boxHvoigt.h2 * shift_z;
    posXyz->x -= delta_x;
    posXyz->y -= delta_y;
    posXyz->z -= delta_z;
  }
#endif
#ifdef __triBox__
  Hvoigt6 boxHvoigtHold = box->boxHvoigt;
  double xtiltmax = 0.505 * box->boxHvoigt.h0;
  double ytiltmax = 0.505 * box->boxHvoigt.h1;
  bool isFlip = false;
  if (box->boxHvoigt.h3 < -ytiltmax) {
    box->boxHvoigt.h3 += box->boxHvoigt.h1;
    box->boxHvoigt.h4 += box->boxHvoigt.h5;
    isFlip = true;
  } else if (box->boxHvoigt.h3 >= ytiltmax) {
    box->boxHvoigt.h3 -= box->boxHvoigt.h1;
    box->boxHvoigt.h4 -= box->boxHvoigt.h5;
    isFlip = true;
  }
  if (box->boxHvoigt.h4 < -xtiltmax) {
    box->boxHvoigt.h4 += box->boxHvoigt.h0;
    isFlip = true;
  } else if (box->boxHvoigt.h4 >= xtiltmax) {
    box->boxHvoigt.h4 -= box->boxHvoigt.h0;
    isFlip = true;
  }
  if (box->boxHvoigt.h5 < -xtiltmax) {
    box->boxHvoigt.h5 += box->boxHvoigt.h0;
    isFlip = true;
  } else if (box->boxHvoigt.h5 >= xtiltmax) {
    box->boxHvoigt.h5 -= box->boxHvoigt.h0;
    isFlip = true;
  }
  if (isFlip) {
    Info("Flip box: (%g %g %g %g %g %g) ==> (%g %g %g %g %g %g)",
         boxHvoigtHold.h0, boxHvoigtHold.h1, boxHvoigtHold.h2, boxHvoigtHold.h3,
         boxHvoigtHold.h4, boxHvoigtHold.h5, box->boxHvoigt.h0,
         box->boxHvoigt.h1, box->boxHvoigt.h2, box->boxHvoigt.h3,
         box->boxHvoigt.h4, box->boxHvoigt.h5);

    setBoxPara(box, particle, thermo, update, var);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->xyz[iatom].x += boxHvoigtHold.h0 * particle->img[iatom].x +
                                boxHvoigtHold.h5 * particle->img[iatom].y +
                                boxHvoigtHold.h4 * particle->img[iatom].z;
      particle->xyz[iatom].y += boxHvoigtHold.h1 * particle->img[iatom].y +
                                boxHvoigtHold.h3 * particle->img[iatom].z;
      particle->xyz[iatom].z += boxHvoigtHold.h2 * particle->img[iatom].z;
    }
    memset(particle->img, '\0', particle->nAtom * sizeof(int3));
  }

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 *posXyz = &particle->xyz[iatom];
    int3 *imgXyz = &particle->img[iatom];

    // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5) = invH * (r -
    // boxLo), boxLo = invH * (-0.5, -0.5, -0.5)
    double lamda_x = box->invBoxHvoigt.h0 * posXyz->x +
                     box->invBoxHvoigt.h5 * posXyz->y +
                     box->invBoxHvoigt.h4 * posXyz->z + 0.5;
    double lamda_y = box->invBoxHvoigt.h1 * posXyz->y +
                     box->invBoxHvoigt.h3 * posXyz->z + 0.5;
    double lamda_z = box->invBoxHvoigt.h2 * posXyz->z + 0.5;
    double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
           shift_z = floor(lamda_z);
    imgXyz->x += (int)shift_x;
    imgXyz->y += (int)shift_y;
    imgXyz->z += (int)shift_z;

    double delta_x = box->boxHvoigt.h0 * shift_x + box->boxHvoigt.h5 * shift_y +
                     box->boxHvoigt.h4 * shift_z;
    double delta_y = box->boxHvoigt.h1 * shift_y + box->boxHvoigt.h3 * shift_z;
    double delta_z = box->boxHvoigt.h2 * shift_z;
    posXyz->x -= delta_x;
    posXyz->y -= delta_y;
    posXyz->z -= delta_z;
  }
#endif
}
//===============================================================================
void read_dump(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  int whichVar = findVariable(var, "rd");
  if (whichVar < 0) Abort("--rd dump.bin whichStep");
  if (var->cmd[whichVar].cmdArgc != 2) Abort("--rd dump.bin whichStep");
  int whichStep = (int)atoi(var->cmd[whichVar].cmdArgv[1]);

  int fd = open(var->cmd[whichVar].cmdArgv[0], O_RDONLY);
  struct stat binStat;
  fstat(fd, &binStat);
  void *memPtr = mmap(NULL, binStat.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (memPtr == MAP_FAILED) Abort("Failed to map file!");
  if (particle->xyz != NULL) Abort("--rf or --rd");

  void *data = memPtr;
  particle->nAtom = *((int *)data);
  data = (void *)((char *)data + sizeof(int));
  particle->nAtomType = *((int *)data);
  data = (void *)((char *)data + sizeof(int));

  particle->xyz = (double3 *)calloc(particle->nAtom, sizeof(double3));
  particle->veloc = (double3 *)calloc(particle->nAtom, sizeof(double3));
  particle->force = (double3 *)calloc(particle->nAtom, sizeof(double3));
  particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
  particle->img = (int3 *)calloc(particle->nAtom, sizeof(int3));
  particle->type = (int *)calloc(particle->nAtom, sizeof(int));
  particle->diameterScale = (double *)calloc(particle->nAtom, sizeof(double));
  particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
  particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
  memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->id2tag[iatom] = iatom;
    particle->tag2id[iatom] = iatom;
  }
  data = (void *)((char *)data + sizeof(double) * particle->nAtom);

  particle->massPerType = (double *)calloc(particle->nAtomType, sizeof(double));
  for (int itype = 0; itype < particle->nAtomType; itype++)
    particle->massPerType[itype] = 1.0;
  memcpy(particle->type, data, particle->nAtom * sizeof(int));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
  }

  int headerSize = sizeof(int) + sizeof(int) +
                   particle->nAtom * (sizeof(double) + sizeof(int));
  int stepSize = sizeof(int) + (3 * particle->nAtom + 7) * sizeof(double);
  int nStep = (binStat.st_size - headerSize) / stepSize;
  if (whichStep < 0 || whichStep >= nStep) Abort("Range: [0,%d];", nStep);
  data = ((char *)memPtr + headerSize + whichStep * stepSize + sizeof(int));

  box->boxHvoigt = *((Hvoigt6 *)data);
  data = ((char *)data + sizeof(Hvoigt6));

  particle->meanDiameter = *((double *)data);
  data = ((char *)data + sizeof(double));

  memcpy(particle->xyz, data, particle->nAtom * sizeof(double3));
  memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
  memset(particle->force, '\0', particle->nAtom * sizeof(double3));
  memset(particle->img, '\0', particle->nAtom * sizeof(int3));

  munmap(memPtr, binStat.st_size);
  close(fd);

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}
void read_data(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  int whichVar = findVariable(var, "rf");
  if (whichVar < 0 || var->cmd[whichVar].cmdArgc == 0)
    Abort("read data: --rf lmp.bin");
  FILE *fp = openReadOnlyFile(var->cmd[whichVar].cmdArgv[0]);
  if (particle->xyz != NULL) Abort("--rf or --rd");

  char str[4096];
  fgets(str, 4096, fp);            // the first line
  if (strcmp(str, "binary") == 0)  // binary file
  {
    bool hasMeanDiameter = false;
    while (fread(str, sizeof(char), 32, fp) == 32) {
      if (strcmp(str, "atoms") == 0) {
        fread(&particle->nAtom, sizeof(int), 1, fp);
        if (particle->nAtom <= 0) Abort("Wrong File!");

        particle->xyz = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->veloc = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->force = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
        particle->img = (int3 *)calloc(particle->nAtom, sizeof(int3));
        particle->type = (int *)calloc(particle->nAtom, sizeof(int));
        particle->diameterScale =
            (double *)calloc(particle->nAtom, sizeof(double));
        particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
        particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      } else if (strcmp(str, "atom types") == 0) {
        fread(&particle->nAtomType, sizeof(int), 1, fp);
        if (particle->nAtomType <= 0) Abort("Wrong File!");

        particle->massPerType =
            (double *)calloc(particle->nAtomType, sizeof(double));
        for (int itype = 0; itype < particle->nAtomType; itype++)
          particle->massPerType[itype] = 1.0;
      } else if (strcmp(str, "box Hvoigt") == 0) {
        fread(&box->boxHvoigt.h0, sizeof(double), 1, fp);
        fread(&box->boxHvoigt.h1, sizeof(double), 1, fp);
        fread(&box->boxHvoigt.h2, sizeof(double), 1, fp);
        fread(&box->boxHvoigt.h3, sizeof(double), 1, fp);
        fread(&box->boxHvoigt.h4, sizeof(double), 1, fp);
        fread(&box->boxHvoigt.h5, sizeof(double), 1, fp);
      } else if (strcmp(str, "mean diameter") == 0) {
        fread(&particle->meanDiameter, sizeof(double), 1, fp);
        hasMeanDiameter = true;
      } else if (strcmp(str, "Atoms") == 0) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
          fread(&particle->type[iatom], sizeof(int), 1, fp);
          fread(&particle->diameterScale[iatom], sizeof(double), 1, fp);
          fread(&particle->xyz[iatom].x, sizeof(double), 1, fp);
          fread(&particle->xyz[iatom].y, sizeof(double), 1, fp);
          fread(&particle->xyz[iatom].z, sizeof(double), 1, fp);
          fread(&particle->img[iatom].x, sizeof(int), 1, fp);
          fread(&particle->img[iatom].y, sizeof(int), 1, fp);
          fread(&particle->img[iatom].z, sizeof(int), 1, fp);

          particle->mass[iatom] = particle->massPerType[particle->type[iatom]];

          particle->tag2id[iatom] = iatom;
          particle->id2tag[iatom] = iatom;
        }
      } else if (strcmp(str, "Velocities") == 0) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
          fread(&particle->veloc[iatom].x, sizeof(double), 1, fp);
          fread(&particle->veloc[iatom].y, sizeof(double), 1, fp);
          fread(&particle->veloc[iatom].z, sizeof(double), 1, fp);
        }
      } else
        Abort("Wrong File!");
    }
    if (!hasMeanDiameter) {
      double meanDiameter = 0.0;
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        meanDiameter += particle->diameterScale[iatom];
      }
      meanDiameter = meanDiameter / particle->nAtom;
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->diameterScale[iatom] /= meanDiameter;
      }
      particle->meanDiameter = meanDiameter;
    }

  } else if (strstr(str, "LAMMPS"))  // Lammps style data file
  {
    while (fgets(str, 4096, fp) != NULL) {
      if (strstr(str, "atoms") != NULL) {
        particle->nAtom = atoi(str);
        if (particle->nAtom <= 0) Abort("No atom!");

        particle->xyz = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->veloc = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->force = (double3 *)calloc(particle->nAtom, sizeof(double3));
        particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
        particle->img = (int3 *)calloc(particle->nAtom, sizeof(int3));
        particle->type = (int *)calloc(particle->nAtom, sizeof(int));
        particle->diameterScale =
            (double *)calloc(particle->nAtom, sizeof(double));
        particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
        particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      }
      if (strstr(str, "atom types") != NULL) {
        int type = atoi(str);
        if (type <= 0) Abort("Wrong DATA file!");
        particle->nAtomType = type;
        particle->massPerType = (double *)calloc(type, sizeof(double));
        for (int itype = 0; itype < type; itype++)
          particle->massPerType[itype] = 1.0;
      }

      if (strstr(str, "xlo xhi") != NULL) {
        double xlo, xhi;
        sscanf(str, "%lf %lf", &xlo, &xhi);
        box->boxHvoigt.h0 = xhi - xlo;
      }
      if (strstr(str, "ylo yhi") != NULL) {
        double ylo, yhi;
        sscanf(str, "%lf %lf", &ylo, &yhi);
        box->boxHvoigt.h1 = yhi - ylo;
      }
      if (strstr(str, "zlo zhi") != NULL) {
        double zlo, zhi;
        sscanf(str, "%lf %lf", &zlo, &zhi);
        box->boxHvoigt.h2 = zhi - zlo;
      }
      if (strstr(str, "xy xz yz") != NULL) {
        double xy, xz, yz;
        sscanf(str, "%lf %lf %lf", &xy, &xz, &yz);
        box->boxHvoigt.h3 = yz;
        box->boxHvoigt.h4 = xz;
        box->boxHvoigt.h5 = xy;
      }

      if (strstr(str, "Atoms") != NULL) {
        double meanDiameter = 0.0;
        for (int iatom = 0; iatom < particle->nAtom;) {
          if (feof(fp) && iatom < particle->nAtom) Abort("Wrong dataFile!");
          fgets(str, 4096, fp);
          if (isEmpty(str, 4096)) continue;

          int num, id, type, ix = 0, iy = 0, iz = 0;
          double diam, density, x, y, z;

          num = sscanf(str, "%d %d %lf %lf %lf %lf %lf %d %d %d", &id, &type,
                       &diam, &density, &x, &y, &z, &ix, &iy, &iz);
          if (num != 10) Abort("Wrong format!");
          if (id <= 0 || id > particle->nAtom) Abort("Wrong atom ID.");
          if (type <= 0 || type > particle->nAtomType)
            Abort("Wrong atom type.");

          particle->xyz[id - 1] = make_double3(x, y, z);
          particle->type[id - 1] = type - 1;
          particle->mass[id - 1] =
              particle->massPerType[particle->type[id - 1]];
          particle->diameterScale[id - 1] = diam;
          particle->img[id - 1] = make_int3(ix, iy, iz);
          particle->tag2id[id - 1] = id - 1;
          particle->id2tag[id - 1] = id - 1;
          meanDiameter += diam;

          iatom++;
        }
        meanDiameter = meanDiameter / particle->nAtom;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
          particle->diameterScale[iatom] /= meanDiameter;
        }
        particle->meanDiameter = meanDiameter;
      }
      if (strstr(str, "Velocities") != NULL) {
        for (int iatom = 0; iatom < particle->nAtom;) {
          if (feof(fp) && iatom < particle->nAtom) Abort("Wrong dataFile!");
          fgets(str, 4096, fp);
          if (isEmpty(str, 4096)) continue;

          int num, id;
          double vx, vy, vz;
          num = sscanf(str, "%d %lf %lf %lf", &id, &vx, &vy, &vz);
          particle->veloc[id - 1] = make_double3(vx, vy, vz);
          iatom++;
        }
      }
    }
  }
  safeCloseFile(fp);

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);

  // reset image
  if (truncFileFlag != 2)
    memset(particle->img, '\0', particle->nAtom * sizeof(int3));
}
void readFile(Box *box, Particle *particle, Thermo *thermo, Update *update,
              Variable *var) {
  int cntR = 0;
  if (findVariable(var, "rf") >= 0) {
    read_data(box, particle, thermo, update, var);
    cntR++;
  }
  if (findVariable(var, "rd") >= 0) {
    read_dump(box, particle, thermo, update, var);
    cntR++;
  }
  if (cntR != 1) Abort("--rf or --rd");

#ifdef __orthBox__
  if (fabs(box->boxHvoigt.h3) >= 5E-16 || fabs(box->boxHvoigt.h4) >= 5E-16 ||
      fabs(box->boxHvoigt.h5) >= 5E-16)
    Abort("Fatall Error: Box!");
  box->boxHvoigt.h3 = 0;
  box->boxHvoigt.h4 = 0;
  box->boxHvoigt.h5 = 0;
#endif
#ifdef __mono_XZ_Box__
  if (fabs(box->boxHvoigt.h3) >= 5E-16 || fabs(box->boxHvoigt.h5) >= 5E-16)
    Abort("Fatall Error: Box!");
  box->boxHvoigt.h3 = 0;
  box->boxHvoigt.h5 = 0;
#endif
}
void write_data(Box *box, Particle *particle, Thermo *thermo, Update *update,
                Variable *var) {
  int whichVar = findVariable(var, "wf");
  if (whichVar < 0) return;
  if (var->cmd[whichVar].cmdArgc == 0) Abort("--wf output.bin");

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);

  FILE *fbin = createReadWriteFile(var->cmd[whichVar].cmdArgv[0]);
  {
    char str[4096];
    // header
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "binary");
    str[31] = '\n';
    fwrite(str, sizeof(char), 32, fbin);  // 32 byte;
    // nAtom atoms
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "atoms");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->nAtom, sizeof(int), 1, fbin);
    // nAtomType atom types
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "atom types");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->nAtomType, sizeof(int), 1, fbin);
    // box Hvoigt
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "box Hvoigt");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&box->boxHvoigt.h0, sizeof(double), 1, fbin);
    fwrite(&box->boxHvoigt.h1, sizeof(double), 1, fbin);
    fwrite(&box->boxHvoigt.h2, sizeof(double), 1, fbin);
    fwrite(&box->boxHvoigt.h3, sizeof(double), 1, fbin);
    fwrite(&box->boxHvoigt.h4, sizeof(double), 1, fbin);
    fwrite(&box->boxHvoigt.h5, sizeof(double), 1, fbin);
    // mean diameter
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "mean diameter");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->meanDiameter, sizeof(double), 1, fbin);
    // Atoms Info
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "Atoms");
    fwrite(str, sizeof(char), 32, fbin);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      int idx = particle->tag2id[iatom];

      int type = particle->type[idx];
      double diaScale = particle->diameterScale[idx];
      double3 xyz = particle->xyz[idx];
      int3 img = particle->img[idx];

      fwrite(&type, sizeof(int), 1, fbin);
      fwrite(&diaScale, sizeof(double), 1, fbin);
      fwrite(&xyz.x, sizeof(double), 1, fbin);
      fwrite(&xyz.y, sizeof(double), 1, fbin);
      fwrite(&xyz.z, sizeof(double), 1, fbin);
      fwrite(&img.x, sizeof(int), 1, fbin);
      fwrite(&img.y, sizeof(int), 1, fbin);
      fwrite(&img.z, sizeof(int), 1, fbin);
    }
  }
  safeCloseFile(fbin);
}
//===============================================================================
void parseCmdLine(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var, int argc, char **argv) {
  // find --trunc or --append tag
  for (int iarg = 1; iarg < argc; iarg++) {
    if (!strcmp(argv[iarg], "--trunc")) {
      if (truncFileFlag != 0) Abort("--trunc or --append");
      truncFileFlag = 1;
    } else if (!strcmp(argv[iarg], "--append")) {
      if (truncFileFlag != 0) Abort("--trunc or --append");
      truncFileFlag = 2;
    }
  }

  // parse cmd line
  var->maxVar = 8;
  var->cmd = (cmdArg *)calloc(var->maxVar, sizeof(cmdArg));
  for (int iarg = 1; iarg < argc;) {
    if (!strcmp(argv[iarg], "--log")) {
      if (iarg + 1 >= argc) Abort("Not variable pair!");
      logFile = createReadWriteFile(argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--trunc") ||
               !strcmp(argv[iarg], "--append")) {
      iarg++;
    } else if (!strcmp(argv[iarg], "--cwd")) {
      if (iarg + 1 >= argc) Abort("--cwd .");
      int len = strlen(argv[iarg + 1]) + 2;
      var->cwd = (char *)calloc(len, sizeof(char));
      sprintf(var->cwd, "%s", argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--sf")) {
      if (iarg + 1 >= argc) Abort("--sf id");
      int len = strlen(argv[iarg + 1]) + 2;
      var->sf = (char *)calloc(len, sizeof(char));
      sprintf(var->sf, "%s", argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--skin")) {
      if (iarg + 1 >= argc) Abort("--skin 0.2");
      double skinSet = atof(argv[iarg + 1]);
      if (skinSet <= 0.0) Abort("Wrong skin");
      update->nebrList.skinSet = skinSet;
      iarg += 2;
    } else if (!strncmp(argv[iarg], "--", 2)) {
      if (findVariable(var, argv[iarg] + 2) >= 0)
        Abort("Repetitive cmd: %s", (argv[iarg] + 2));
      if (var->nVar == var->maxVar) {
        var->maxVar += 8;
        var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
      }

      int nsize = strlen(argv[iarg]) - 2 + 1;
      var->cmd[var->nVar].cmdArgc = 0;
      var->cmd[var->nVar].cmdArgv = NULL;
      int cmdArgvStart = iarg + 1;

      iarg++;
      while (iarg < argc && strncmp(argv[iarg], "--", 2)) {
        var->cmd[var->nVar].cmdArgc++;
        nsize += strlen(argv[iarg]) + 1;
        iarg++;
      }
      if (var->cmd[var->nVar].cmdArgc != 0) {
        var->cmd[var->nVar].cmdArgv =
            (char **)calloc(var->cmd[var->nVar].cmdArgc, sizeof(char *));
      }
      char *ptr = (char *)calloc(nsize + 5, sizeof(char));

      var->cmd[var->nVar].cmdType = ptr;
      memcpy(ptr, argv[cmdArgvStart - 1] + 2,
             (strlen(argv[cmdArgvStart - 1]) - 2 + 1) * sizeof(char));
      ptr += strlen(argv[cmdArgvStart - 1]) - 2 + 1;
      for (int ith = cmdArgvStart; ith < iarg; ith++) {
        var->cmd[var->nVar].cmdArgv[ith - cmdArgvStart] = ptr;
        memcpy(ptr, argv[ith], (strlen(argv[ith]) + 1) * sizeof(char));
        ptr += strlen(argv[ith]) + 1;
      }

      var->nVar++;
    } else
      Abort("Unrecognized cmd %s!", argv[iarg]);
  }

  if (var->cwd == NULL) {
    var->cwd = (char *)calloc(3, sizeof(char));
    sprintf(var->cwd, ".");
  }
  if (var->sf == NULL) {
    var->sf = (char *)calloc(9, sizeof(char));
    sprintf(var->sf, "default");
  }
  if (update->nebrList.skinSet <= 0) {
    update->nebrList.skinSet = 0.1 / 1.1;
  }

  readFile(box, particle, thermo, update, var);

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
  printf("===========System Info==========\n");
}
//===============================================================================
void sortParticle(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  if (!update->nebrList.doSort) return;

  NebrList *nebrList = &update->nebrList;
  double3 DistPlane;
  // L(a,b): distance between plane (a,b); {vol / |a cross c|} = boxHvoigt.h2;
  DistPlane.z = box->boxHvoigt.h2;
  // L(c,a): distance between plane (c,a); {vol / |c cross a|}
  double sca = sqrt(pow(box->boxHvoigt.h0 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h0 * box->boxHvoigt.h3, 2));
  DistPlane.y = box->volume / sca;
  // L(b,c): distance between plane (b,c); {vol / |b cross c|}
  double sbc = sqrt(pow(box->boxHvoigt.h1 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h5 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h5 * box->boxHvoigt.h3 -
                            box->boxHvoigt.h1 * box->boxHvoigt.h4,
                        2));
  DistPlane.x = box->volume / sbc;

  double sysRcs = particle->meanDiameter * (1.0 + nebrList->skinSet);
  int3 nbin4sort;
  nbin4sort.x = (int)floor(DistPlane.x / sysRcs);
  nbin4sort.y = (int)floor(DistPlane.y / sysRcs);
  nbin4sort.z = (int)floor(DistPlane.z / sysRcs);
  nbin4sort.x = cpuMax(1, nbin4sort.x);
  nbin4sort.y = cpuMax(1, nbin4sort.y);
  nbin4sort.z = cpuMax(1, nbin4sort.z);
  nebrList->totBin4sort = nbin4sort.x * nbin4sort.y * nbin4sort.z;
  if (nebrList->totBin4sort > nebrList->allocBin4sort) {
    nebrList->binHead4sort = (int *)realloc(
        nebrList->binHead4sort, nebrList->totBin4sort * sizeof(int));
    nebrList->allocBin4sort = nebrList->totBin4sort;
  }

  for (int ibin = 0; ibin < nebrList->totBin4sort; ibin++) {
    nebrList->binHead4sort[ibin] = -1;
  }
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 *posXyz = &particle->xyz[iatom];

    int3 cIdx;
    double lamda_x = box->invBoxHvoigt.h0 * posXyz->x +
                     box->invBoxHvoigt.h5 * posXyz->y +
                     box->invBoxHvoigt.h4 * posXyz->z + 0.5;
    double lamda_y = box->invBoxHvoigt.h1 * posXyz->y +
                     box->invBoxHvoigt.h3 * posXyz->z + 0.5;
    double lamda_z = box->invBoxHvoigt.h2 * posXyz->z + 0.5;
    cIdx.x = (int)floor(lamda_x * nbin4sort.x);
    cIdx.y = (int)floor(lamda_y * nbin4sort.y);
    cIdx.z = (int)floor(lamda_z * nbin4sort.z);

    cIdx.x = (cIdx.x < 0 ? nbin4sort.x - 1 : cIdx.x);
    cIdx.x = (cIdx.x >= nbin4sort.x ? 0 : cIdx.x);
    cIdx.y = (cIdx.y < 0 ? nbin4sort.y - 1 : cIdx.y);
    cIdx.y = (cIdx.y >= nbin4sort.y ? 0 : cIdx.y);
    cIdx.z = (cIdx.z < 0 ? nbin4sort.z - 1 : cIdx.z);
    cIdx.z = (cIdx.z >= nbin4sort.z ? 0 : cIdx.z);

    int binIdx =
        cIdx.z * nbin4sort.y * nbin4sort.x + cIdx.y * nbin4sort.x + cIdx.x;
    nebrList->binList4sort[iatom] = nebrList->binHead4sort[binIdx];
    nebrList->binHead4sort[binIdx] = iatom;
  }

  for (int iz = 0, noid = 0; iz < nbin4sort.z; iz++) {
    for (int iy = 0; iy < nbin4sort.y; iy++) {
      for (int ix = 0; ix < nbin4sort.x; ix++) {
        int cidx = iz * nbin4sort.y * nbin4sort.x + iy * nbin4sort.x + ix;

        int oid = nebrList->binHead4sort[cidx];
        while (oid != -1) {
          nebrList->oid2nid[oid] = noid;
          noid++;
          oid = nebrList->binList4sort[oid];
        }
      }
    }
  }

  exchange_double3(particle->xyz, (double3 *)nebrList->buffer,
                   nebrList->oid2nid, particle->nAtom);
  memcpy(particle->xyz, nebrList->buffer, particle->nAtom * sizeof(double3));

  exchange_double3(particle->veloc, (double3 *)nebrList->buffer,
                   nebrList->oid2nid, particle->nAtom);
  memcpy(particle->veloc, nebrList->buffer, particle->nAtom * sizeof(double3));

  exchange_int3(particle->img, (int3 *)nebrList->buffer, nebrList->oid2nid,
                particle->nAtom);
  memcpy(particle->img, nebrList->buffer, particle->nAtom * sizeof(int3));

  exchange_double(particle->mass, (double *)nebrList->buffer, nebrList->oid2nid,
                  particle->nAtom);
  memcpy(particle->mass, nebrList->buffer, particle->nAtom * sizeof(double));

  exchange_double(particle->diameterScale, (double *)nebrList->buffer,
                  nebrList->oid2nid, particle->nAtom);
  memcpy(particle->diameterScale, nebrList->buffer,
         particle->nAtom * sizeof(double));

  exchange_int(particle->type, (int *)nebrList->buffer, nebrList->oid2nid,
               particle->nAtom);
  memcpy(particle->type, nebrList->buffer, particle->nAtom * sizeof(int));

  exchange_int(particle->id2tag, (int *)nebrList->buffer, nebrList->oid2nid,
               particle->nAtom);
  memcpy(particle->id2tag, nebrList->buffer, particle->nAtom * sizeof(int));

  for (int aid = 0; aid < particle->nAtom; aid++) {
    particle->tag2id[particle->id2tag[aid]] = aid;
  }

  update->nebrList.isValid = false;

  // reset
  nebrList->maxAtomPerBin = 0;
}
void initNebrList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (nebrList->isInit) return;

  if (nebrList->xyzHold == NULL) {
    nebrList->xyzHold = (double3 *)calloc(particle->nAtom, sizeof(double3));
  }
  if (nebrList->nNebr == NULL) {
    nebrList->nNebr = (int2 *)calloc(particle->nAtom, sizeof(int2));
  }

  if (nebrList->binList4sort == NULL) {
    nebrList->binList4sort = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (nebrList->oid2nid == NULL) {
    nebrList->oid2nid = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (nebrList->buffer == NULL) {
    nebrList->buffer = calloc(particle->nAtom, sizeof(double3));
  }

  sortParticle(box, particle, thermo, update, var);

  nebrList->isInit = true;
}
void calcBinParam(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->compelInit && box->isShapeFixed && particle->isSizeFixed &&
      nebrList->deltaAdjBin)
    return;

  double3 DistPlane;
  // L(a,b): distance between plane (a,b); {vol / |a cross c|} = boxHvoigt.h2;
  DistPlane.z = box->boxHvoigt.h2;
  // L(c,a): distance between plane (c,a); {vol / |c cross a|}
  double sca = sqrt(pow(box->boxHvoigt.h0 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h0 * box->boxHvoigt.h3, 2));
  DistPlane.y = box->volume / sca;
  // L(b,c): distance between plane (b,c); {vol / |b cross c|}
  double sbc = sqrt(pow(box->boxHvoigt.h1 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h5 * box->boxHvoigt.h2, 2) +
                    pow(box->boxHvoigt.h5 * box->boxHvoigt.h3 -
                            box->boxHvoigt.h1 * box->boxHvoigt.h4,
                        2));
  DistPlane.x = box->volume / sbc;

  nebrList->rskin =
      nebrList->minDiameterScale * particle->meanDiameter * nebrList->skinSet;
  double sysRcs =
      nebrList->maxDiameterScale * particle->meanDiameter + nebrList->rskin;
  nebrList->nbin.x = (int)floor(DistPlane.x / sysRcs);
  nebrList->nbin.y = (int)floor(DistPlane.y / sysRcs);
  nebrList->nbin.z = (int)floor(DistPlane.z / sysRcs);
  nebrList->nbin.x = cpuMax(nebrList->nbin.x, 1);
  nebrList->nbin.y = cpuMax(nebrList->nbin.y, 1);
  nebrList->nbin.z = cpuMax(nebrList->nbin.z, 1);
  nebrList->totBin = nebrList->nbin.x * nebrList->nbin.y * nebrList->nbin.z;
  nebrList->binLen.x = DistPlane.x / nebrList->nbin.x;
  nebrList->binLen.y = DistPlane.y / nebrList->nbin.y;
  nebrList->binLen.z = DistPlane.z / nebrList->nbin.z;
  if (nebrList->totBin > particle->nAtom) {
    int3 nbin;
    nbin.x = (int)floor(pow(
        particle->nAtom * DistPlane.y * DistPlane.z / DistPlane.x / DistPlane.x,
        1.0 / 3.0));
    nbin.y = (int)floor(DistPlane.y / DistPlane.x * nbin.x);
    nbin.z = (int)floor(DistPlane.z / DistPlane.x * nbin.y);
    double3 binLen;
    binLen.x = DistPlane.x / nbin.x;
    binLen.y = DistPlane.y / nbin.y;
    binLen.z = DistPlane.z / nbin.z;
    minElementVec(sysRcs, binLen);
    double rskin = sysRcs - nebrList->maxDiameterScale * particle->meanDiameter;
    if (rskin >=
        __minSkinSet__ * nebrList->minDiameterScale * particle->meanDiameter) {
      nebrList->nbin = nbin;
      nebrList->binLen = binLen;
      nebrList->rskin = rskin;
    }
  }

  double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
  if (maxRcut >= DistPlane.x || maxRcut >= DistPlane.y ||
      maxRcut >= DistPlane.z)
    Abort("The Box is too samll!");
  int3 xyzNum, bstart, bstop, ndb;
  xyzNum.z = (int)ceil((maxRcut + nebrList->rskin) / nebrList->binLen.z);
  xyzNum.y = (int)ceil((maxRcut + nebrList->rskin) / nebrList->binLen.y);
  xyzNum.x = (int)ceil((maxRcut + nebrList->rskin) / nebrList->binLen.x);

  bstart.z = -xyzNum.z, bstop.z = xyzNum.z, ndb.z = 2 * xyzNum.z + 1;
  bstart.y = -xyzNum.y, bstop.y = xyzNum.y, ndb.y = 2 * xyzNum.y + 1;
  bstart.x = -xyzNum.x, bstop.x = xyzNum.x, ndb.x = 2 * xyzNum.x + 1;
  if (!nebrList->isFullStyle) {  // half style
    int nAdjBin = ndb.z * ndb.y * ndb.x / 2 + 1;
    if (nAdjBin > nebrList->nAdjBin) {
      nebrList->deltaAdjBin =
          (int3 *)realloc(nebrList->deltaAdjBin, nAdjBin * sizeof(int3));
    }
    nebrList->nAdjBin = 0;
    for (int iz = bstart.z; iz <= bstop.z; iz++) {
      for (int iy = bstart.y; iy <= bstop.y; iy++) {
        for (int ix = bstart.x; ix <= bstop.x; ix++) {
          if (iz * ndb.y * ndb.x + iy * ndb.x + ix < 0) continue;
          nebrList->deltaAdjBin[nebrList->nAdjBin++] = make_int3(ix, iy, iz);
        }
      }
    }
  } else {
    int nAdjBin = ndb.z * ndb.y * ndb.x;
    if (nAdjBin > nebrList->nAdjBin) {
      nebrList->deltaAdjBin =
          (int3 *)realloc(nebrList->deltaAdjBin, nAdjBin * sizeof(int3));
    }
    nebrList->nAdjBin = 0;
    for (int iz = bstart.z; iz <= bstop.z; iz++) {
      for (int iy = bstart.y; iy <= bstop.y; iy++) {
        for (int ix = bstart.x; ix <= bstop.x; ix++) {
          nebrList->deltaAdjBin[nebrList->nAdjBin++] = make_int3(ix, iy, iz);
        }
      }
    }
  }

  if (nebrList->totBin > nebrList->allocBin) {
    nebrList->nAtomPerBin =
        (int *)realloc(nebrList->nAtomPerBin, nebrList->totBin * sizeof(int));
    nebrList->allocBin = nebrList->totBin;
  }

  nebrList->compelInit = false;
}
void binParticle(Box *box, Particle *particle, Thermo *thermo, Update *update,
                 Variable *var) {
  NebrList *nebrList = &update->nebrList;

  calcBinParam(box, particle, thermo, update, var);

  bool overFlow = false;
  do {
    overFlow = false;
    memset(nebrList->nAtomPerBin, '\0', nebrList->totBin * sizeof(int));
    int maxPartPerBin = nebrList->maxAtomPerBin;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      idPosRadius data;
      data.id = iatom;
      data.pos = particle->xyz[iatom];
      data.radius =
          0.5 * particle->meanDiameter * particle->diameterScale[iatom];

      int3 cIdx;
      double lamda_x = box->invBoxHvoigt.h0 * data.pos.x +
                       box->invBoxHvoigt.h5 * data.pos.y +
                       box->invBoxHvoigt.h4 * data.pos.z + 0.5;
      double lamda_y = box->invBoxHvoigt.h1 * data.pos.y +
                       box->invBoxHvoigt.h3 * data.pos.z + 0.5;
      double lamda_z = box->invBoxHvoigt.h2 * data.pos.z + 0.5;
      cIdx.x = (int)floor(lamda_x * nebrList->nbin.x);
      cIdx.y = (int)floor(lamda_y * nebrList->nbin.y);
      cIdx.z = (int)floor(lamda_z * nebrList->nbin.z);

      cIdx.x = (cIdx.x < 0 ? nebrList->nbin.x - 1 : cIdx.x);
      cIdx.x = (cIdx.x >= nebrList->nbin.x ? 0 : cIdx.x);
      cIdx.y = (cIdx.y < 0 ? nebrList->nbin.y - 1 : cIdx.y);
      cIdx.y = (cIdx.y >= nebrList->nbin.y ? 0 : cIdx.y);
      cIdx.z = (cIdx.z < 0 ? nebrList->nbin.z - 1 : cIdx.z);
      cIdx.z = (cIdx.z >= nebrList->nbin.z ? 0 : cIdx.z);

      int binIdx = cIdx.z * nebrList->nbin.y * nebrList->nbin.x +
                   cIdx.y * nebrList->nbin.x + cIdx.x;

      if (nebrList->nAtomPerBin[binIdx] < nebrList->maxAtomPerBin) {
        nebrList->ipr[nebrList->maxAtomPerBin * binIdx +
                      nebrList->nAtomPerBin[binIdx]] = data;

        nebrList->nAtomPerBin[binIdx]++;
      } else {
        nebrList->nAtomPerBin[binIdx]++;
        maxPartPerBin = (nebrList->nAtomPerBin[binIdx] > maxPartPerBin)
                            ? nebrList->nAtomPerBin[binIdx]
                            : maxPartPerBin;
        overFlow = true;
      }
    }
    if (overFlow) {
      // realloc
      nebrList->ipr = (idPosRadius *)realloc(
          nebrList->ipr, particle->nAtom * maxPartPerBin * sizeof(idPosRadius));
      nebrList->maxAtomPerBin = maxPartPerBin;
    }
  } while (overFlow);
}
void constructList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  NebrList *nebrList = &update->nebrList;
  memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int2));
  for (int bidx = 0, cntNebr = 0; bidx < nebrList->totBin; bidx++) {
    int zbin = bidx / (nebrList->nbin.y * nebrList->nbin.x);
    int ybin = bidx / nebrList->nbin.x - zbin * nebrList->nbin.y;
    int xbin = bidx % nebrList->nbin.x;
    for (int ith = 0; ith < nebrList->nAtomPerBin[bidx]; ith++) {
      idPosRadius idata = nebrList->ipr[bidx * nebrList->maxAtomPerBin + ith];
      double iRadiusRskin = idata.radius + nebrList->rskin;
      nebrList->nNebr[idata.id].x = cntNebr;
      for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
        int zadj = (zbin + nebrList->deltaAdjBin[adj].z + nebrList->nbin.z) %
                   nebrList->nbin.z;
        int yadj = (ybin + nebrList->deltaAdjBin[adj].y + nebrList->nbin.y) %
                   nebrList->nbin.y;
        int xadj = (xbin + nebrList->deltaAdjBin[adj].x + nebrList->nbin.x) %
                   nebrList->nbin.x;
        int adjBidx = zadj * nebrList->nbin.y * nebrList->nbin.x +
                      yadj * nebrList->nbin.x + xadj;
        for (int jth = 0; jth < nebrList->nAtomPerBin[adjBidx]; jth++) {
          idPosRadius jdata =
              nebrList->ipr[adjBidx * nebrList->maxAtomPerBin + jth];
          if (adj == 0 && jdata.id <= idata.id) continue;

          double3 dRij;
          vecSub(dRij, idata.pos, jdata.pos);
          PBC(dRij, box);
          double rijP2;
          vecNormP2(rijP2, dRij);
          double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
          if (rijP2 > RcutRskinP2) continue;
          if (cntNebr == nebrList->maxAllocNebr) {
            nebrList->maxAllocNebr += 1024;
            int *tmp = (int *)realloc(nebrList->list,
                                      nebrList->maxAllocNebr * sizeof(int));
            if (tmp == NULL) Abort("realloc nebrList failed!");
            nebrList->list = tmp;
          }
          nebrList->list[cntNebr++] = jdata.id;
        }
      }
      nebrList->nNebr[idata.id].y = cntNebr;
    }
  }
}
void buildNebrList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (nebrList->isFullStyle) {
    nebrList->compelInit = true;
  }
  nebrList->isFullStyle = false;

  initNebrList(box, particle, thermo, update, var);

  adjustImg(box, particle, thermo, update, var);
  if ((nebrList->nBuild % 100 == 0)) {
    sortParticle(box, particle, thermo, update, var);
  }
  binParticle(box, particle, thermo, update, var);

  nebrList->meanDiameterHold = particle->meanDiameter;
  nebrList->invBoxHvoigtHold = box->invBoxHvoigt;
  memcpy(nebrList->xyzHold, particle->xyz, particle->nAtom * sizeof(double3));

  constructList(box, particle, thermo, update, var);

  nebrList->nBuild++;
  nebrList->isValid = true;
  nebrList->cntForce = 0;
}
bool isNebrListValid(Box *box, Particle *particle, Thermo *thermo,
                     Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->isInit) {
    nebrList->isValid = false;
    return false;
  }
  if (!nebrList->isValid) {
    nebrList->isValid = false;
    return false;
  }

  if (particle->isSizeFixed && box->isShapeFixed) {
    double maxShiftP2 = 0.25 * nebrList->rskin * nebrList->rskin;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      double3 dRi0;
      vecSub(dRi0, particle->xyz[iatom], nebrList->xyzHold[iatom]);
      double rP2;
      vecNormP2(rP2, dRi0);
      if (rP2 >= maxShiftP2) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else if (!particle->isSizeFixed && box->isShapeFixed) {
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      double3 dRi0;
      vecSub(dRi0, particle->xyz[iatom], nebrList->xyzHold[iatom]);
      double ri0;
      vecNorm(ri0, dRi0);
      double dsig = (particle->meanDiameter - nebrList->meanDiameterHold) *
                    particle->diameterScale[iatom];

      if (ri0 + dsig >= nebrList->rskin * 0.5) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else if (!box->isShapeFixed && particle->isSizeFixed) {
    Hvoigt6 trans;  // trans =  Ht * inv(H0)
    trans.h0 = box->boxHvoigt.h0 * nebrList->invBoxHvoigtHold.h0;
    trans.h1 = box->boxHvoigt.h1 * nebrList->invBoxHvoigtHold.h1;
    trans.h2 = box->boxHvoigt.h2 * nebrList->invBoxHvoigtHold.h2;
    trans.h3 = box->boxHvoigt.h1 * nebrList->invBoxHvoigtHold.h3 +
               box->boxHvoigt.h3 * nebrList->invBoxHvoigtHold.h2;
    trans.h4 = box->boxHvoigt.h0 * nebrList->invBoxHvoigtHold.h4 +
               box->boxHvoigt.h5 * nebrList->invBoxHvoigtHold.h3 +
               box->boxHvoigt.h4 * nebrList->invBoxHvoigtHold.h2;
    trans.h5 = box->boxHvoigt.h0 * nebrList->invBoxHvoigtHold.h5 +
               box->boxHvoigt.h5 * nebrList->invBoxHvoigtHold.h1;
    double sfact = 1;
    sfact = cpuMin(trans.h0, trans.h1);
    sfact = cpuMin(sfact, trans.h2);
    double sYz = trans.h3 * trans.h3;
    sYz = sqrt(1.0 + 0.5 * sYz - sqrt(sYz + 0.25 * sYz * sYz));
    double sXz = trans.h4 * trans.h4;
    sXz = sqrt(1.0 + 0.5 * sXz - sqrt(sXz + 0.25 * sXz * sXz));
    double sXy = trans.h5 * trans.h5;
    sXy = sqrt(1.0 + 0.5 * sXy - sqrt(sXy + 0.25 * sXy * sXy));
    sfact = sfact * sYz * sXz * sXy;
    double rs_eff = nebrList->rskin;
    if (sfact >= 1.0) {
      rs_eff = sfact * rs_eff +
               (sfact - 1.0) *
                   (2.0 * nebrList->minDiameterScale * particle->meanDiameter);
    } else {
      rs_eff = sfact * rs_eff +
               (sfact - 1.0) *
                   (2.0 * nebrList->maxDiameterScale * particle->meanDiameter);
    }
    if (rs_eff <= 0) {
      nebrList->isValid = false;
      return false;
    }
    double maxShiftP2 = 0.25 * rs_eff * rs_eff;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      double3 xyzHold = nebrList->xyzHold[iatom];
      xyzHold.x =
          xyzHold.x * trans.h0 + xyzHold.y * trans.h5 + xyzHold.z * trans.h4;
      xyzHold.y = xyzHold.y * trans.h1 + xyzHold.z * trans.h3;
      xyzHold.z = xyzHold.z * trans.h2;
      double3 dR;
      vecSub(dR, particle->xyz[iatom], xyzHold);
      double rP2;
      vecNormP2(rP2, dR);

      if (rP2 >= maxShiftP2) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else
    Abort("Not Code!");
  nebrList->isValid = true;
  return true;
}

void calcForce(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (particle->isForceValid) return;

  bool isListValid = isNebrListValid(box, particle, thermo, update, var);
  if (!isListValid) {
    if (nebrList->cntForce < 10) {
      nebrList->skinSet *= 1.1;
      nebrList->skinSet = cpuMin(nebrList->skinSet, 0.5);
      nebrList->compelInit = true;
    } else if (nebrList->cntForce > 20) {
      nebrList->skinSet *= 0.9;
      nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__);
      nebrList->compelInit = true;
    }
  } else {
    if (nebrList->skinSet > __minSkinSet__ * (1 + 1E-10) &&
        nebrList->cntForce > 20) {
      nebrList->skinSet *= 0.9;
      nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__);
      nebrList->compelInit = true;
      isListValid = false;
    }
  }
  if (!isListValid) {
    buildNebrList(box, particle, thermo, update, var);
  }

  double3 *force = particle->force;
  memset(force, '\0', particle->nAtom * sizeof(double3));
  thermo->Epair = 0;
  thermo->virialPair = 0;
  thermo->lapU = 0;
  memset(&thermo->ptensor, '\0', sizeof(Hvoigt6));

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 iforce = force[iatom];
    double3 iPos = particle->xyz[iatom];
    double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
    for (int jth = nebrList->nNebr[iatom].x; jth < nebrList->nNebr[iatom].y;
         jth++) {
      int jatom = nebrList->list[jth];
      double3 dRij;
      double sRc =
          iRc + particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
      vecSub(dRij, particle->xyz[iatom], particle->xyz[jatom]);
      PBC(dRij, box);

      double rijP2;
      vecNormP2(rijP2, dRij);
      if (rijP2 >= sRc * sRc) continue;
      double rij = sqrt(rijP2);
      double rdivsig = rij / sRc;
      double fpair = (1.0 - rdivsig) / rij / sRc;
      vecScaleAdd(iforce, iforce, fpair, dRij);
      vecScaleAdd(force[jatom], force[jatom], -fpair, dRij);

      thermo->Epair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
      thermo->lapU += (6.0 * rdivsig * rdivsig - 4.0 * rdivsig) / rijP2;

      thermo->ptensor.h0 += fpair * dRij.x * dRij.x;
      thermo->ptensor.h1 += fpair * dRij.y * dRij.y;
      thermo->ptensor.h2 += fpair * dRij.z * dRij.z;
      thermo->ptensor.h3 += fpair * dRij.y * dRij.z;
      thermo->ptensor.h4 += fpair * dRij.x * dRij.z;
      thermo->ptensor.h5 += fpair * dRij.x * dRij.y;
    }
    force[iatom] = iforce;
  }

  thermo->Epair = thermo->Epair / particle->nAtom;
  thermo->lapU = thermo->lapU / particle->nAtom;
  thermo->virialPair =
      thermo->ptensor.h0 + thermo->ptensor.h1 + thermo->ptensor.h2;

  thermo->ptensor.h0 = thermo->ptensor.h0 / box->volume;
  thermo->ptensor.h1 = thermo->ptensor.h1 / box->volume;
  thermo->ptensor.h2 = thermo->ptensor.h2 / box->volume;
  thermo->ptensor.h3 = thermo->ptensor.h3 / box->volume;
  thermo->ptensor.h4 = thermo->ptensor.h4 / box->volume;
  thermo->ptensor.h5 = thermo->ptensor.h5 / box->volume;
  thermo->pressure =
      (thermo->ptensor.h0 + thermo->ptensor.h1 + thermo->ptensor.h2) / 3.0;

  thermo->Edone = thermo->Pdone = true;

  particle->isForceValid = true;
  nebrList->cntForce++;
  nebrList->nForce++;
}
//===============================================================================
void instant_inflate(Box *box, Particle *particle, Thermo *thermo,
                     Update *update, Variable *var, double deltaVF) {
  // inflate particles to target packing fraction: thermo-volFrac + deltaVF
  calcVolInfo(box, particle, thermo, update, var);
  double volfrac_target = thermo->volFrac + deltaVF;
  double sfact = pow(volfrac_target / thermo->volFrac, 1.0 / 3.0);
  particle->meanDiameter *= sfact;
  thermo->volFrac = sfact * sfact * sfact * thermo->volFrac;

  calcVolInfo(box, particle, thermo, update, var);

  setBasicInfo(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
void instant_simpShearXz(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var, double deltaGamma) {
  // shear box in XZ direction by deltaGamma
  adjustImg(box, particle, thermo, update, var);

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->xyz[iatom].x += particle->xyz[iatom].z * deltaGamma;
  }
  box->boxHvoigt.h4 += deltaGamma * box->boxHvoigt.h2;
  setBoxPara(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
void instant_pureShearXz(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var, double scaleZ) {
  adjustImg(box, particle, thermo, update, var);

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->xyz[iatom].z *= scaleZ;
    particle->xyz[iatom].x /= scaleZ;
  }

  //(h0,0,0)
  box->boxHvoigt.h0 /= scaleZ;
  //(h5,h1,0)
  box->boxHvoigt.h5 /= scaleZ;
  //(h4,h3,h2)
  box->boxHvoigt.h4 /= scaleZ;
  box->boxHvoigt.h2 *= scaleZ;

  setBoxPara(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
void instant_deformation(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var, Hvoigt6 strain) {
  adjustImg(box, particle, thermo, update, var);
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    double3 R0 = particle->xyz[iatom];
    double3 Rt;
    Rt.x = R0.x + R0.x * strain.h0 + R0.y * strain.h5 + R0.z * strain.h4;
    Rt.y = R0.y + R0.y * strain.h1 + R0.z * strain.h3;
    Rt.z = R0.z + R0.z * strain.h2;
    particle->xyz[iatom] = Rt;
  }

  double3 a0 = make_double3(box->boxHvoigt.h0, 0, 0), at;
  at.x = a0.x + a0.x * strain.h0;

  double3 b0 = make_double3(box->boxHvoigt.h5, box->boxHvoigt.h1, 0), bt;
  bt.x = b0.x + b0.x * strain.h0 + b0.y * strain.h5;
  bt.y = b0.y + b0.y * strain.h1;

  double3 c0 = make_double3(box->boxHvoigt.h4, box->boxHvoigt.h3,
                            box->boxHvoigt.h2),
          ct;
  ct.x = c0.x + c0.x * strain.h0 + c0.y * strain.h5 + c0.z * strain.h4;
  ct.y = c0.y + c0.y * strain.h1 + c0.z * strain.h3;
  ct.z = c0.z + c0.z * strain.h2;

  box->boxHvoigt = make_Hvoigt6(at.x, bt.y, ct.z, ct.y, ct.x, bt.x);

  setBoxPara(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
void normaliseBox(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  double sfact = pow(box->volume, 1.0 / 3.0);
  particle->meanDiameter /= sfact;

  box->boxHvoigt.h0 /= sfact;
  box->boxHvoigt.h1 /= sfact;
  box->boxHvoigt.h2 /= sfact;
  box->boxHvoigt.h3 /= sfact;
  box->boxHvoigt.h4 /= sfact;
  box->boxHvoigt.h5 /= sfact;

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    vecScale(particle->xyz[iatom], 1.0 / sfact, particle->xyz[iatom]);
  }

  thermo->isInit = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.isInit = false;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
}
//===============================================================================
void backupSimInfo(Box *box, Particle *particle, Box *bakBox,
                   Particle *bakParticle) {
  if (bakParticle->nAtom != particle->nAtom) {
    bakParticle->nAtom = particle->nAtom;

    bakParticle->xyz = (double3 *)realloc(bakParticle->xyz,
                                          bakParticle->nAtom * sizeof(double3));
    bakParticle->veloc = (double3 *)realloc(
        bakParticle->veloc, bakParticle->nAtom * sizeof(double3));
    bakParticle->mass = (double *)realloc(bakParticle->mass,
                                          bakParticle->nAtom * sizeof(double));
    bakParticle->img =
        (int3 *)realloc(bakParticle->img, bakParticle->nAtom * sizeof(int3));
    bakParticle->type =
        (int *)realloc(bakParticle->type, bakParticle->nAtom * sizeof(int));
    bakParticle->diameterScale = (double *)realloc(
        bakParticle->diameterScale, bakParticle->nAtom * sizeof(double));
    bakParticle->id2tag =
        (int *)realloc(bakParticle->id2tag, bakParticle->nAtom * sizeof(int));
    bakParticle->tag2id =
        (int *)realloc(bakParticle->tag2id, bakParticle->nAtom * sizeof(int));
  }
  if (bakParticle->nAtomType != particle->nAtomType) {
    bakParticle->nAtomType = particle->nAtomType;
    bakParticle->massPerType = (double *)realloc(
        bakParticle->massPerType, bakParticle->nAtomType * sizeof(double));
  }

  memcpy(bakParticle->xyz, particle->xyz, bakParticle->nAtom * sizeof(double3));
  memcpy(bakParticle->veloc, particle->veloc,
         bakParticle->nAtom * sizeof(double3));
  memcpy(bakParticle->mass, particle->mass,
         bakParticle->nAtom * sizeof(double));
  memcpy(bakParticle->img, particle->img, bakParticle->nAtom * sizeof(int3));
  memcpy(bakParticle->type, particle->type, bakParticle->nAtom * sizeof(int));
  memcpy(bakParticle->diameterScale, particle->diameterScale,
         bakParticle->nAtom * sizeof(double));
  memcpy(bakParticle->id2tag, particle->id2tag,
         bakParticle->nAtom * sizeof(int));
  memcpy(bakParticle->tag2id, particle->tag2id,
         bakParticle->nAtom * sizeof(int));

  memcpy(bakParticle->massPerType, particle->massPerType,
         bakParticle->nAtomType * sizeof(double));

  memcpy(&bakBox->boxHvoigt, &box->boxHvoigt, sizeof(Hvoigt6));
}
void restoreSimInfo(Box *box, Particle *particle, Thermo *thermo,
                    Update *update, Variable *var, Box *bakBox,
                    Particle *bakParticle) {
  // restore data
  backupSimInfo(bakBox, bakParticle, box, particle);

  // clear
  box->boxLo = box->boxHi = make_double3(0, 0, 0);
  box->invBoxHvoigt = make_Hvoigt6(0, 0, 0, 0, 0, 0);
  box->volume = 0;

  particle->isForceValid = false;

  thermo->isInit = false;
  thermo->Edone = false;
  thermo->Pdone = false;
  thermo->Tdone = false;

  update->nebrList.isValid = false;
  update->nebrList.isInit = false;
  update->nebrList.compelInit = true;

  // reinit
  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}
void freeSimInfo(Box *bakBox, Particle *bakParticle) {
  safeFree(bakParticle->xyz);
  safeFree(bakParticle->veloc);
  safeFree(bakParticle->mass);
  safeFree(bakParticle->img);
  safeFree(bakParticle->type);
  safeFree(bakParticle->diameterScale);
  safeFree(bakParticle->id2tag);
  safeFree(bakParticle->tag2id);
}
//===============================================================================
void constructList_full(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int2));
  for (int bidx = 0, cntNebr = 0; bidx < nebrList->totBin; bidx++) {
    int zbin = bidx / (nebrList->nbin.y * nebrList->nbin.x);
    int ybin = bidx / nebrList->nbin.x - zbin * nebrList->nbin.y;
    int xbin = bidx % nebrList->nbin.x;
    for (int ith = 0; ith < nebrList->nAtomPerBin[bidx]; ith++) {
      idPosRadius idata = nebrList->ipr[bidx * nebrList->maxAtomPerBin + ith];
      double iRadiusRskin = idata.radius + nebrList->rskin;
      nebrList->nNebr[idata.id].x = cntNebr;
      for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
        int zadj = (zbin + nebrList->deltaAdjBin[adj].z + nebrList->nbin.z) %
                   nebrList->nbin.z;
        int yadj = (ybin + nebrList->deltaAdjBin[adj].y + nebrList->nbin.y) %
                   nebrList->nbin.y;
        int xadj = (xbin + nebrList->deltaAdjBin[adj].x + nebrList->nbin.x) %
                   nebrList->nbin.x;
        int adjBidx = zadj * nebrList->nbin.y * nebrList->nbin.x +
                      yadj * nebrList->nbin.x + xadj;
        for (int jth = 0; jth < nebrList->nAtomPerBin[adjBidx]; jth++) {
          idPosRadius jdata =
              nebrList->ipr[adjBidx * nebrList->maxAtomPerBin + jth];
          if (idata.id == jdata.id) continue;

          double3 dRij;
          vecSub(dRij, idata.pos, jdata.pos);
          PBC(dRij, box);
          double rijP2;
          vecNormP2(rijP2, dRij);
          double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
          if (rijP2 > RcutRskinP2) continue;
          if (cntNebr == nebrList->maxAllocNebr) {
            nebrList->maxAllocNebr += 1024;
            int *tmp = (int *)realloc(nebrList->list,
                                      nebrList->maxAllocNebr * sizeof(int));
            if (tmp == NULL) Abort("realloc nebrList failed!");
            nebrList->list = tmp;
          }
          nebrList->list[cntNebr++] = jdata.id;
        }
      }
      nebrList->nNebr[idata.id].y = cntNebr;
    }
  }
}
void buildNebrList_full(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->isFullStyle) {
    nebrList->compelInit = true;
  }
  nebrList->isFullStyle = true;

  initNebrList(box, particle, thermo, update, var);

  adjustImg(box, particle, thermo, update, var);
  if ((nebrList->nBuild % 100 == 0)) {
    sortParticle(box, particle, thermo, update, var);
  }
  binParticle(box, particle, thermo, update, var);

  nebrList->meanDiameterHold = particle->meanDiameter;
  nebrList->invBoxHvoigtHold = box->invBoxHvoigt;
  memcpy(nebrList->xyzHold, particle->xyz, particle->nAtom * sizeof(double3));

  constructList_full(box, particle, thermo, update, var);

  nebrList->nBuild++;
  nebrList->isValid = true;
  nebrList->cntForce = 0;
}
//===============================================================================
contactInfo *initContactInfo(Box *box, Particle *particle, Thermo *thermo,
                             Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  contactInfo *cinfo = NULL;
  if (whichTool >= 0) {
    cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  } else {
    cinfo = (contactInfo *)calloc(sizeof(contactInfo), 1);
    addToolkit(&thermo->toolkit, (void *)cinfo, "contactInfo");
  }
  if (particle->nAtom <= 0) Abort("No Atom!");

  if (cinfo->isRattler == NULL) {
    cinfo->isRattler = (bool *)calloc(particle->nAtom, sizeof(bool));
  }
  if (cinfo->nCoordNumExRattler == NULL) {
    cinfo->nCoordNumExRattler = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (cinfo->nCoordNum == NULL) {
    cinfo->nCoordNum = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (cinfo->forceExRattler == NULL) {
    cinfo->forceExRattler = (double3 *)calloc(particle->nAtom, sizeof(double3));
  }

  return cinfo;
}
void computeCoordination(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  if (whichTool < 0) return;
  contactInfo *cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  NebrList *nebrList = &update->nebrList;
  if (!isNebrListValid(box, particle, thermo, update, var)) {
    buildNebrList(box, particle, thermo, update, var);
  }
  if (cinfo->allocMem4CoordNum < nebrList->maxAllocNebr) {
    cinfo->allocMem4CoordNum = nebrList->maxAllocNebr;
    cinfo->isNebrContact = (bool *)realloc(
        cinfo->isNebrContact, cinfo->allocMem4CoordNum * sizeof(bool));
  }

  memset(cinfo->isRattler, '\0', particle->nAtom * sizeof(bool));
  memset(cinfo->nCoordNumExRattler, '\0', particle->nAtom * sizeof(int));
  memset(cinfo->nCoordNum, '\0', particle->nAtom * sizeof(int));
  memset(cinfo->isNebrContact, '\0', cinfo->allocMem4CoordNum * sizeof(bool));

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    for (int jth = nebrList->nNebr[iatom].x; jth < nebrList->nNebr[iatom].y;
         jth++) {
      int jatom = nebrList->list[jth];

      double3 dRij;
      double sRc =
          (particle->diameterScale[iatom] + particle->diameterScale[jatom]) *
          particle->meanDiameter * 0.5;
      vecSub(dRij, particle->xyz[iatom], particle->xyz[jatom]);
      PBC(dRij, box);
      double rij;
      vecNorm(rij, dRij);
      if (rij < sRc) {
        cinfo->isNebrContact[jth] = true;
        cinfo->nCoordNum[jatom]++;
        cinfo->nCoordNum[iatom]++;
      }
    }
  }

  int cordNumWR = 0, cordNumNR = 0, nRattler = 0;
  bool done = true;
  memcpy(cinfo->nCoordNumExRattler, cinfo->nCoordNum,
         particle->nAtom * sizeof(int));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->nCoordNumExRattler[iatom] < (DIM + 1)) {
      cinfo->isRattler[iatom] = true;
      done = false;
      nRattler++;
    }
    cordNumWR += cinfo->nCoordNum[iatom];
    cordNumNR += cinfo->nCoordNumExRattler[iatom];
  }
  while (!done) {
    done = true;
    memset(cinfo->nCoordNumExRattler, '\0', particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      if (cinfo->isRattler[iatom]) continue;
      for (int jth = nebrList->nNebr[iatom].x; jth < nebrList->nNebr[iatom].y;
           jth++) {
        int jatom = nebrList->list[jth];
        if (!cinfo->isNebrContact[jth]) continue;
        if (!cinfo->isRattler[jatom]) {
          cinfo->nCoordNumExRattler[iatom]++;
          cinfo->nCoordNumExRattler[jatom]++;
        }
      }
    }

    cordNumNR = 0;
    nRattler = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      if (cinfo->nCoordNumExRattler[iatom] < (DIM + 1)) {
        if (!cinfo->isRattler[iatom]) {
          done = false;
        }
        cinfo->isRattler[iatom] = true;
        nRattler++;
      }
      cordNumNR += cinfo->nCoordNumExRattler[iatom];
    }
  }
  cinfo->nRattler = nRattler;
  cinfo->aveCoordNum = (double)cordNumWR / (double)particle->nAtom;
  if (nRattler != particle->nAtom) {
    cinfo->aveCoordNumExRattler =
        (double)cordNumNR / (double)(particle->nAtom - nRattler);
  } else {
    cinfo->aveCoordNumExRattler = 0;
  }

  if (particle->nAtom == cinfo->nRattler) {
    cinfo->meanForceExRattler = 0;
    cinfo->lapUexRattler = 0;
    return;
  }
  double3 *force = cinfo->forceExRattler;
  double lapU = 0;
  memset(force, '\0', particle->nAtom * sizeof(double3));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->isRattler[iatom]) continue;
    double3 iforce = force[iatom];
    double3 iPos = particle->xyz[iatom];
    double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
    for (int jth = nebrList->nNebr[iatom].x; jth < nebrList->nNebr[iatom].y;
         jth++) {
      int jatom = nebrList->list[jth];
      if (cinfo->isRattler[jatom]) continue;

      double3 dRij;
      double sRc =
          iRc + particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
      vecSub(dRij, particle->xyz[iatom], particle->xyz[jatom]);
      PBC(dRij, box);

      double rijP2;
      vecNormP2(rijP2, dRij);
      if (rijP2 >= sRc * sRc) continue;
      double rij = sqrt(rijP2);
      double rdivsig = rij / sRc;
      double fpair = (1.0 - rdivsig) / rij / sRc;
      vecScaleAdd(iforce, iforce, fpair, dRij);
      vecScaleAdd(force[jatom], force[jatom], -fpair, dRij);

      lapU += (6.0 * rdivsig * rdivsig - 4.0 * rdivsig) / rijP2;
    }
    force[iatom] = iforce;
  }

  double sumForce = 0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->isRattler[iatom]) continue;
    double dotPrd;
    vecDot(dotPrd, force[iatom], force[iatom]);
    sumForce += sqrt(dotPrd);
  }
  cinfo->meanForceExRattler = sumForce / thermo->forceUnits /
                              (double)(particle->nAtom - cinfo->nRattler);

  cinfo->lapUexRattler = lapU / (double)(particle->nAtom - cinfo->nRattler);
}
void finalizeContactInfo(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  if (whichTool < 0) return;
  contactInfo *cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  safeFree(cinfo->isRattler);
  safeFree(cinfo->isNebrContact);
  safeFree(cinfo->nCoordNum);
  safeFree(cinfo->nCoordNumExRattler);
  safeFree(cinfo->forceExRattler);
  delToolkit(&thermo->toolkit, "contactInfo");
}

dumpInfo *initDump(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  int whichVar = findVariable(var, "dump");
  if (whichVar < 0) return (dumpInfo *)NULL;
  dumpInfo *dinfo = (dumpInfo *)calloc(sizeof(dumpInfo), 1);
  addToolkit(&thermo->toolkit, (void *)dinfo, "dumpInfo");

  char fname[4096];
  sprintf(fname, "%s/dump_%s.bin", var->cwd, var->sf);
  dinfo->fdump = createReadWriteFile(fname);

  fwrite(&particle->nAtom, sizeof(int), 1, dinfo->fdump);
  fwrite(&particle->nAtomType, sizeof(int), 1, dinfo->fdump);
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    fwrite(&particle->diameterScale[idx], sizeof(double), 1, dinfo->fdump);
  }
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    fwrite(&particle->type[idx], sizeof(int), 1, dinfo->fdump);
  }
  return dinfo;
}
void dump(Box *box, Particle *particle, Thermo *thermo, Update *update,
          Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "dumpInfo");
  if (whichTool < 0) return;

  FILE *fdump = ((dumpInfo *)thermo->toolkit.toolkit[whichTool])->fdump;
  fwrite(&particle->nAtom, sizeof(int), 1, fdump);
  fwrite(&box->boxHvoigt.h0, sizeof(double), 1, fdump);
  fwrite(&box->boxHvoigt.h1, sizeof(double), 1, fdump);
  fwrite(&box->boxHvoigt.h2, sizeof(double), 1, fdump);
  fwrite(&box->boxHvoigt.h3, sizeof(double), 1, fdump);
  fwrite(&box->boxHvoigt.h4, sizeof(double), 1, fdump);
  fwrite(&box->boxHvoigt.h5, sizeof(double), 1, fdump);
  fwrite(&particle->meanDiameter, sizeof(double), 1, fdump);
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    double3 xyz = make_double3(0, 0, 0);
    vecHvoigtMulVec(xyz, box->boxHvoigt, particle->img[idx]);
    vecAdd(xyz, particle->xyz[idx], xyz);

    fwrite(&xyz.x, sizeof(double), 1, fdump);
    fwrite(&xyz.y, sizeof(double), 1, fdump);
    fwrite(&xyz.z, sizeof(double), 1, fdump);
  }
}
void finalizeDump(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "dumpInfo");
  if (whichTool < 0) return;
  dumpInfo *dinfo = (dumpInfo *)thermo->toolkit.toolkit[whichTool];
  safeCloseFile(dinfo->fdump);
  delToolkit(&thermo->toolkit, "dumpInfo");
}

BinDumpFile *initMapBinFile(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  BinDumpFile *binFile = (BinDumpFile *)calloc(sizeof(BinDumpFile), 1);
  addToolkit(&thermo->toolkit, (void *)binFile, "BinDumpFile");

  int whichVar = findVariable(var, "read_dump");
  if (whichVar < 0 || var->cmd[whichVar].cmdArgc != 1) {
    Abort("--read_dump dump.bin");
  }
  if (binFile->dataSection != NULL) Abort("Fatal Error!");

  binFile->fd = open(var->cmd[whichVar].cmdArgv[0], O_RDONLY);
  fstat(binFile->fd, &binFile->binStat);
  binFile->dataSection = mmap(NULL, binFile->binStat.st_size, PROT_READ,
                              MAP_PRIVATE, binFile->fd, 0);
  if (binFile->dataSection == MAP_FAILED) {
    Abort("Map Binary file Failed !\n");
  }

  void *data = binFile->dataSection;
  int nElement = *((int *)data);
  binFile->headerSize =
      sizeof(int) * 2 + nElement * (sizeof(double) + sizeof(int));
  binFile->stepSize = sizeof(int) + 6 * sizeof(double) + sizeof(double) +
                      nElement * 3 * sizeof(double);
  binFile->nStep =
      (binFile->binStat.st_size - binFile->headerSize) / binFile->stepSize;

  return binFile;
}
void finalizeMapBinFile(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "BinDumpFile");
  if (whichTool < 0) return;
  BinDumpFile *binFile = (BinDumpFile *)thermo->toolkit.toolkit[whichTool];
  munmap(binFile->dataSection, binFile->binStat.st_size);
  close(binFile->fd);
  binFile->fd = -1;
  binFile->dataSection = NULL;
  delToolkit(&thermo->toolkit, "BinDumpFile");
}
void getBinDumpData(Box *box, Particle *particle, Thermo *thermo,
                    Update *update, Variable *var, int whichStep) {
  int whichTool = findToolkit(&thermo->toolkit, "BinDumpFile");
  if (whichTool < 0) return;
  BinDumpFile *binFile = (BinDumpFile *)thermo->toolkit.toolkit[whichTool];
  if (whichStep >= binFile->nStep || whichStep < 0) {
    Abort("Wrong Step!");
  }

  {
    void *data = binFile->dataSection;
    particle->nAtom = *((int *)data);
    data = (void *)((char *)data + sizeof(int));
    particle->nAtomType = *((int *)data);
    data = (void *)((char *)data + sizeof(int));

    if (particle->xyz == NULL) {
      particle->xyz = (double3 *)calloc(particle->nAtom, sizeof(double3));
      particle->veloc = (double3 *)calloc(particle->nAtom, sizeof(double3));
      particle->force = (double3 *)calloc(particle->nAtom, sizeof(double3));
      particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
      particle->img = (int3 *)calloc(particle->nAtom, sizeof(int3));
      particle->type = (int *)calloc(particle->nAtom, sizeof(int));
      particle->diameterScale =
          (double *)calloc(particle->nAtom, sizeof(double));
      particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
      particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      particle->massPerType =
          (double *)calloc(particle->nAtomType, sizeof(double));
      for (int itype = 0; itype < particle->nAtomType; itype++)
        particle->massPerType[itype] = 1.0;
    }

    memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
    data = (void *)((char *)data + sizeof(double) * particle->nAtom);
    memcpy(particle->type, data, particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->id2tag[iatom] = iatom;
      particle->tag2id[iatom] = iatom;
    }

    data = ((char *)binFile->dataSection + binFile->headerSize +
            whichStep * binFile->stepSize + sizeof(int));

    box->boxHvoigt = *((Hvoigt6 *)data);
    data = ((char *)data + sizeof(Hvoigt6));

    particle->meanDiameter = *((double *)data);
    data = ((char *)data + sizeof(double));

    memcpy(particle->xyz, data, particle->nAtom * sizeof(double3));

    memset(particle->veloc, '\0', particle->nAtom * sizeof(double3));
    memset(particle->force, '\0', particle->nAtom * sizeof(double3));
    memset(particle->img, '\0', particle->nAtom * sizeof(int3));
  }

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;
  update->nebrList.skinSet = __minSkinSet__;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}

#endif

#else

#ifndef __SUBFUNC__
#define __SUBFUNC__

#include "StructSim.h"
#include "VectorMath.h"

//===============================================================================
void calcVolInfo(Box *box, Particle *particle, Thermo *thermo, Update *update,
                 Variable *var) {
  if (thermo->isInit) return;

  double vol = 0;
  double minScale = 1E10, maxScale = 0.0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    vol += VolUnitSphere * pow(particle->diameterScale[iatom], DIM);
    double tmp = particle->diameterScale[iatom];
    minScale = (minScale < tmp ? minScale : tmp);
    maxScale = (maxScale > tmp ? maxScale : tmp);
  }
  update->nebrList.minDiameterScale = minScale;
  update->nebrList.maxDiameterScale = maxScale;

  thermo->volFrac = vol * pow(particle->meanDiameter * 0.5, DIM) / box->volume;
  thermo->isInit = true;
}
void setBoxPara(Box *box, Particle *particle, Thermo *thermo, Update *update,
                Variable *var) {
  box->dim = DIM;
  box->volume = 1.0;
  for (int idim = 0; idim < DIM; idim++) {
    if (box->boxH[spaceIdx(idim, idim)] <= 0.0) {
      Abort("Not rihgt hand basis!");
    }
    box->volume *= box->boxH[spaceIdx(idim, idim)];

    zerosVec(box->boxEdge[idim]);
    for (int jdim = 0; jdim <= idim; jdim++) {
      box->boxEdge[idim][jdim] = box->boxH[spaceIdx(jdim, idim)];
    }
  }

  zerosVec(box->cornerLo);
  for (int idim = 0; idim < DIM; idim++) {
    vecAdd(box->cornerLo, box->cornerLo, box->boxEdge[idim]);
  }
  scaleVec(box->cornerLo, -0.5, box->cornerLo);
  scaleVec(box->cornerHi, -1.0, box->cornerLo);

  zerosUptriMat(box->invBoxH);
  for (int idim = DIM - 1; idim >= 0; idim--) {
    box->invBoxH[spaceIdx(idim, idim)] = 1.0 / box->boxH[spaceIdx(idim, idim)];
    for (int jdim = idim + 1; jdim < DIM; jdim++) {
      double cij = 0.0;
      for (int k = idim + 1; k <= jdim; k++) {
        cij += box->boxH[spaceIdx(idim, k)] * box->invBoxH[spaceIdx(k, jdim)];
      }
      box->invBoxH[spaceIdx(idim, jdim)] =
          -cij / box->boxH[spaceIdx(idim, idim)];
    }
  }
}
void setBasicInfo(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  calcVolInfo(box, particle, thermo, update, var);

  thermo->massUnits = 1.0;
  thermo->energyUnits = 1.0;
  thermo->distanceUnits = particle->meanDiameter;
  thermo->timeUnits =
      sqrt(thermo->massUnits / thermo->energyUnits) * thermo->distanceUnits;
  thermo->forceUnits = thermo->energyUnits / thermo->distanceUnits;
  thermo->velocityUnits = sqrt(thermo->energyUnits / thermo->massUnits);
  thermo->pressureUnits = thermo->energyUnits / pow(thermo->distanceUnits, DIM);
  thermo->volumeUnits = pow(thermo->distanceUnits, DIM);
}
void adjustImg(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
#ifdef __orthBox__
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    doubleVecPtr posPtr = particle->pos[iatom];
    intVecPtr imgPtr = particle->img[iatom];

    // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
    doubleVector lamda;
    MatMulVec(lamda, box->invBoxH, posPtr);
    vecShiftAll(lamda, 0.5);

    intVector shiftImg;
    for (int idim = 0; idim < DIM; idim++) {
      shiftImg[idim] = (int)floor(lamda[idim]);
    }
    vecAdd(imgPtr, imgPtr, shiftImg);

    doubleVector shiftPos;
    MatMulVec(shiftPos, box->boxH, shiftImg);
    vecSub(posPtr, posPtr, shiftPos);
  }
#endif
#ifdef __triBox__
  bool isFlip = false;
  doubleVector maxTilt;
  for (int idim = 0; idim < DIM; idim++) {
    maxTilt[idim] = 0.505 * box->boxH[spaceIdx(idim, idim)];
  }
  doubleVector boxEdge[DIM];
  memcpy(boxEdge, box->boxEdge, DIM * sizeof(doubleVector));
  for (int idim = DIM - 1; idim >= 0; idim--) {
    for (int tilt = idim - 1; tilt >= 0; tilt--) {
      if (boxEdge[idim][tilt] >= maxTilt[tilt]) {
        vecSub(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
        isFlip = true;
      } else if (boxEdge[idim][tilt] < -maxTilt[tilt]) {
        vecAdd(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
        isFlip = true;
      }
    }
  }
  if (isFlip) {
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      unwrapPos(particle->pos[iatom], particle->pos[iatom],
                particle->img[iatom], box->boxH);
    }
    memset(particle->img, '\0', particle->nAtom * sizeof(intVector));

    for (int idim = 0; idim < DIM; idim++) {
      for (int jdim = idim; jdim < DIM; jdim++) {
        box->boxH[spaceIdx(idim, jdim)] = boxEdge[jdim][idim];
      }
    }
    setBoxPara(box, particle, thermo, update, var);
  }

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    doubleVecPtr posPtr = particle->pos[iatom];
    intVecPtr imgPtr = particle->img[iatom];

    // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
    doubleVector lamda;
    MatMulVec(lamda, box->invBoxH, posPtr);
    vecShiftAll(lamda, 0.5);

    intVector shiftImg;
    for (int idim = 0; idim < DIM; idim++) {
      shiftImg[idim] = (int)floor(lamda[idim]);
    }
    vecAdd(imgPtr, imgPtr, shiftImg);

    doubleVector shiftPos;
    MatMulVec(shiftPos, box->boxH, shiftImg);
    vecSub(posPtr, posPtr, shiftPos);
  }
#endif
}
//===============================================================================
void read_dump(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  int whichVar = findVariable(var, "rd");
  if (whichVar < 0) Abort("--rd dump.bin whichStep");
  if (var->cmd[whichVar].cmdArgc != 2) Abort("--rd dump.bin whichStep");
  int whichStep = (int)atoi(var->cmd[whichVar].cmdArgv[1]);

  int fd = open(var->cmd[whichVar].cmdArgv[0], O_RDONLY);
  struct stat binStat;
  fstat(fd, &binStat);
  void *memPtr = mmap(NULL, binStat.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (memPtr == MAP_FAILED) Abort("Failed to map file!");
  if (particle->pos != NULL) Abort("--rf or --rd");

  // header
  char *header = (char *)memPtr;
  void *data = memPtr;
  int dumpRevNum = 0;
  if (strcmp(header, "Revised Binary File") == 0) {
    int *revNum = (int *)(header + 32);
    dumpRevNum = revNum[0];
    int *dim = (int *)(revNum + 1);
    if (dim[0] != DIM) {
      Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0],
            DIM);
    }
    data = (void *)(dim + 1);
  } else {
    if (DIM != 3) {
      Abort("The dumpfile is for d = 3, while the code is for d = %d!", DIM);
    }
  }
  box->dim = DIM;

  particle->nAtom = *((int *)data);
  data = (void *)((char *)data + sizeof(int));
  particle->nAtomType = *((int *)data);
  data = (void *)((char *)data + sizeof(int));

  particle->pos = (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
  particle->veloc =
      (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
  particle->force =
      (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
  particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
  particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
  particle->type = (int *)calloc(particle->nAtom, sizeof(int));
  particle->diameterScale = (double *)calloc(particle->nAtom, sizeof(double));
  particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
  particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
  memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->id2tag[iatom] = iatom;
    particle->tag2id[iatom] = iatom;
  }
  data = (void *)((char *)data + sizeof(double) * particle->nAtom);

  particle->massPerType = (double *)calloc(particle->nAtomType, sizeof(double));
  for (int itype = 0; itype < particle->nAtomType; itype++)
    particle->massPerType[itype] = 1.0;
  memcpy(particle->type, data, particle->nAtom * sizeof(int));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
  }

  int headerSize = sizeof(int) + sizeof(int) +
                   particle->nAtom * (sizeof(double) + sizeof(int));
  if (dumpRevNum == 1) {
    headerSize += 32 * sizeof(char) + 2 * sizeof(int);
  }

  int stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) +
                 particle->nAtom * sizeof(doubleVector);

  int nStep = (binStat.st_size - headerSize) / stepSize;
  if (whichStep < 0 || whichStep >= nStep) Abort("Range: [0,%d];", nStep);
  data = ((char *)memPtr + headerSize + whichStep * stepSize + sizeof(int));

  if (dumpRevNum == 0) {
    double *tmp = (double *)data;
    box->boxH[spaceIdx(0, 0)] = tmp[0];
    box->boxH[spaceIdx(1, 1)] = tmp[1];
    box->boxH[spaceIdx(2, 2)] = tmp[2];
    box->boxH[spaceIdx(1, 2)] = tmp[3];
    box->boxH[spaceIdx(0, 2)] = tmp[4];
    box->boxH[spaceIdx(0, 1)] = tmp[5];
    data = ((char *)data + 6 * sizeof(double));
  } else {
    memcpy(box->boxH, data, sizeof(uptriMat));
    data = ((char *)data + sizeof(uptriMat));
  }

  particle->meanDiameter = *((double *)data);
  data = ((char *)data + sizeof(double));

  memcpy(particle->pos, data, particle->nAtom * sizeof(doubleVector));
  memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
  memset(particle->force, '\0', particle->nAtom * sizeof(doubleVector));
  memset(particle->img, '\0', particle->nAtom * sizeof(intVector));

  munmap(memPtr, binStat.st_size);
  close(fd);

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}
void read_data(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  int whichVar = findVariable(var, "rf");
  if (whichVar < 0 || var->cmd[whichVar].cmdArgc == 0)
    Abort("read data: --rf lmp.bin");
  FILE *fp = openReadOnlyFile(var->cmd[whichVar].cmdArgv[0]);
  if (particle->pos != NULL) Abort("--rf or --rd");

  char str[4096];
  fgets(str, 4096, fp);
  if (strcmp(str, "binary") == 0) {
    bool hasMeanDiameter = false, hasDimension = false;
    while (fread(str, sizeof(char), 32, fp) == 32) {
      if (strcmp(str, "dimension") == 0) {
        fread(&box->dim, sizeof(int), 1, fp);
        if (box->dim != DIM)
          Abort("The file is for d = %d, while the code is for d = %d!",
                box->dim, DIM);
        hasDimension = true;
      } else if (strcmp(str, "atoms") == 0) {
        fread(&particle->nAtom, sizeof(int), 1, fp);
        if (particle->nAtom <= 0) Abort("Wrong File!");

        particle->pos =
            (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
        particle->veloc =
            (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
        particle->force =
            (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
        particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
        particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
        particle->type = (int *)calloc(particle->nAtom, sizeof(int));
        particle->diameterScale =
            (double *)calloc(particle->nAtom, sizeof(double));
        particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
        particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      } else if (strcmp(str, "atom types") == 0) {
        fread(&particle->nAtomType, sizeof(int), 1, fp);
        if (particle->nAtomType <= 0) Abort("Wrong File!");

        particle->massPerType =
            (double *)calloc(particle->nAtomType, sizeof(double));
        for (int itype = 0; itype < particle->nAtomType; itype++)
          particle->massPerType[itype] = 1.0;
      } else if (strcmp(str, "box Hvoigt") == 0) {
        if (!hasDimension) {
          box->dim = 3;
          if (box->dim != DIM) Abort("Dimension is not consistent!");

          fread(&box->boxH[spaceIdx(0, 0)], sizeof(double), 1, fp);  // xx
          fread(&box->boxH[spaceIdx(1, 1)], sizeof(double), 1, fp);  // yy
          fread(&box->boxH[spaceIdx(2, 2)], sizeof(double), 1, fp);  // zz
          fread(&box->boxH[spaceIdx(1, 2)], sizeof(double), 1, fp);  // yz
          fread(&box->boxH[spaceIdx(0, 2)], sizeof(double), 1, fp);  // xz
          fread(&box->boxH[spaceIdx(0, 1)], sizeof(double), 1, fp);  // xy
        } else {
          fread(&box->boxH, sizeof(uptriMat), 1, fp);
        }
      } else if (strcmp(str, "mean diameter") == 0) {
        fread(&particle->meanDiameter, sizeof(double), 1, fp);
        hasMeanDiameter = true;
      } else if (strcmp(str, "Atoms") == 0) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
          fread(&particle->type[iatom], sizeof(int), 1, fp);
          fread(&particle->diameterScale[iatom], sizeof(double), 1, fp);
          fread(&particle->pos[iatom], sizeof(doubleVector), 1, fp);
          fread(&particle->img[iatom], sizeof(intVector), 1, fp);

          particle->mass[iatom] = particle->massPerType[particle->type[iatom]];

          particle->tag2id[iatom] = iatom;
          particle->id2tag[iatom] = iatom;
        }
      } else if (strcmp(str, "Velocities") == 0) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
          fread(&particle->veloc[iatom], sizeof(doubleVector), 1, fp);
        }
      } else
        Abort("Wrong File!");
    }
    if (!hasMeanDiameter) {
      double meanDiameter = 0.0;
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        meanDiameter += particle->diameterScale[iatom];
      }
      meanDiameter = meanDiameter / particle->nAtom;
      for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->diameterScale[iatom] /= meanDiameter;
      }
      particle->meanDiameter = meanDiameter;
    }
  } else
    Abort("Wrong File!");
  safeCloseFile(fp);

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);

  // reset image
  if (truncFileFlag != 2)
    memset(particle->img, '\0', particle->nAtom * sizeof(intVector));
}
void readFile(Box *box, Particle *particle, Thermo *thermo, Update *update,
              Variable *var) {
  int cntR = 0;
  if (findVariable(var, "rf") >= 0) {
    read_data(box, particle, thermo, update, var);
    cntR++;
  }
  if (findVariable(var, "rd") >= 0) {
    read_dump(box, particle, thermo, update, var);
    cntR++;
  }
  if (cntR != 1) Abort("--rf or --rd");

#ifdef __orthBox__
  for (int idim = 0; idim < DIM; idim++) {
    for (int jdim = idim + 1; jdim < DIM; jdim++) {
      if (fabs(box->boxH[spaceIdx(idim, jdim)]) >= 5E-16)
        Abort("Not orthogonal Box!");
      box->boxH[spaceIdx(idim, jdim)] = 0;
    }
  }
#endif
}
void write_data(Box *box, Particle *particle, Thermo *thermo, Update *update,
                Variable *var) {
  int whichVar = findVariable(var, "wf");
  if (whichVar < 0) return;
  if (var->cmd[whichVar].cmdArgc == 0) Abort("--wf output.bin");

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);

  FILE *fbin = createReadWriteFile(var->cmd[whichVar].cmdArgv[0]);
  {
    char str[4096];
    // header
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "binary");
    str[31] = '\n';
    fwrite(str, sizeof(char), 32, fbin);  // 32 byte;
    // dimension
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "dimension");
    fwrite(str, sizeof(char), 32, fbin);  // 32 byte;
    fwrite(&box->dim, sizeof(int), 1, fbin);
    // nAtom atoms
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "atoms");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->nAtom, sizeof(int), 1, fbin);
    // nAtomType atom types
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "atom types");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->nAtomType, sizeof(int), 1, fbin);
    // box Hvoigt
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "box Hvoigt");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&box->boxH, sizeof(uptriMat), 1, fbin);
    // mean diameter
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "mean diameter");
    fwrite(str, sizeof(char), 32, fbin);
    fwrite(&particle->meanDiameter, sizeof(double), 1, fbin);
    // Atoms Info
    memset(str, '\0', 4096 * sizeof(char));
    sprintf(str, "Atoms");
    fwrite(str, sizeof(char), 32, fbin);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      int idx = particle->tag2id[iatom];

      int type = particle->type[idx];
      double diaScale = particle->diameterScale[idx];
      doubleVecPtr posPtr = particle->pos[idx];
      intVecPtr imgPtr = particle->img[idx];

      fwrite(&type, sizeof(int), 1, fbin);
      fwrite(&diaScale, sizeof(double), 1, fbin);
      fwrite(posPtr, sizeof(doubleVector), 1, fbin);
      fwrite(imgPtr, sizeof(intVector), 1, fbin);
    }
  }
  safeCloseFile(fbin);
}
//===============================================================================
void parseCmdLine(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var, int argc, char **argv) {
  // find --trunc or --append tag
  for (int iarg = 1; iarg < argc; iarg++) {
    if (!strcmp(argv[iarg], "--trunc")) {
      if (truncFileFlag != 0) Abort("--trunc or --append");
      truncFileFlag = 1;
    } else if (!strcmp(argv[iarg], "--append")) {
      if (truncFileFlag != 0) Abort("--trunc or --append");
      truncFileFlag = 2;
    }
  }

  // parse cmd line
  var->maxVar = 8;
  var->cmd = (cmdArg *)calloc(var->maxVar, sizeof(cmdArg));
  for (int iarg = 1; iarg < argc;) {
    if (!strcmp(argv[iarg], "--log")) {
      if (iarg + 1 >= argc) Abort("Not variable pair!");
      logFile = createReadWriteFile(argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--trunc") ||
               !strcmp(argv[iarg], "--append")) {
      iarg++;
    } else if (!strcmp(argv[iarg], "--cwd")) {
      if (iarg + 1 >= argc) Abort("--cwd .");
      int len = strlen(argv[iarg + 1]) + 2;
      var->cwd = (char *)calloc(len, sizeof(char));
      sprintf(var->cwd, "%s", argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--sf")) {
      if (iarg + 1 >= argc) Abort("--sf id");
      int len = strlen(argv[iarg + 1]) + 2;
      var->sf = (char *)calloc(len, sizeof(char));
      sprintf(var->sf, "%s", argv[iarg + 1]);
      iarg += 2;
    } else if (!strcmp(argv[iarg], "--skin")) {
      if (iarg + 1 >= argc) Abort("--skin 0.2");
      double skinSet = atof(argv[iarg + 1]);
      if (skinSet <= 0.0) Abort("Wrong skin");
      update->nebrList.skinSet = skinSet;
      iarg += 2;
    } else if (!strncmp(argv[iarg], "--", 2)) {
      if (findVariable(var, argv[iarg] + 2) >= 0)
        Abort("Repetitive cmd: %s", (argv[iarg] + 2));
      if (var->nVar == var->maxVar) {
        var->maxVar += 8;
        var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
      }

      int nsize = strlen(argv[iarg]) - 2 + 1;
      var->cmd[var->nVar].cmdArgc = 0;
      var->cmd[var->nVar].cmdArgv = NULL;
      int cmdArgvStart = iarg + 1;

      iarg++;
      while (iarg < argc && strncmp(argv[iarg], "--", 2)) {
        var->cmd[var->nVar].cmdArgc++;
        nsize += strlen(argv[iarg]) + 1;
        iarg++;
      }
      if (var->cmd[var->nVar].cmdArgc != 0) {
        var->cmd[var->nVar].cmdArgv =
            (char **)calloc(var->cmd[var->nVar].cmdArgc, sizeof(char *));
      }
      char *ptr = (char *)calloc(nsize + 5, sizeof(char));

      var->cmd[var->nVar].cmdType = ptr;
      memcpy(ptr, argv[cmdArgvStart - 1] + 2,
             (strlen(argv[cmdArgvStart - 1]) - 2 + 1) * sizeof(char));
      ptr += strlen(argv[cmdArgvStart - 1]) - 2 + 1;
      for (int ith = cmdArgvStart; ith < iarg; ith++) {
        var->cmd[var->nVar].cmdArgv[ith - cmdArgvStart] = ptr;
        memcpy(ptr, argv[ith], (strlen(argv[ith]) + 1) * sizeof(char));
        ptr += strlen(argv[ith]) + 1;
      }

      var->nVar++;
    } else
      Abort("Unrecognized cmd %s!", argv[iarg]);
  }

  if (var->cwd == NULL) {
    var->cwd = (char *)calloc(3, sizeof(char));
    sprintf(var->cwd, ".");
  }
  if (var->sf == NULL) {
    var->sf = (char *)calloc(9, sizeof(char));
    sprintf(var->sf, "default");
  }
  if (update->nebrList.skinSet <= 0) {
    update->nebrList.skinSet = 0.1 / 1.1;
  }
  if (findVariable(var, "read_dump") >= 0) return;

  readFile(box, particle, thermo, update, var);

  printf("===========System Info==========\n");
  printf("Dimension: %d\n", DIM);
  printf("Number of Particles: %d;\n", particle->nAtom);
  printf("Volume Fraction: %g;\n", thermo->volFrac);
  printf("Min(diameter): %g;\nMax(diameter): %g;\n",
         update->nebrList.minDiameterScale, update->nebrList.maxDiameterScale);
  printf("Edges of Simulation Box: \n");
  for (int idim = 0; idim < DIM; idim++) {
    for (int jdim = 0; jdim < DIM; jdim++) {
      printf("\t");
      printf("%-8.6e\t", box->boxEdge[idim][jdim] / thermo->distanceUnits);
    }
    printf("\n");
  }
  printf("===========System Info==========\n");
}
//===============================================================================
void calcDistBoxPlane(doubleVector distPlane, doubleVector boxEdge[DIM]) {
  doubleVector edge[DIM];

#if false
  for (int idim = 0; idim < DIM; idim++) {
    for (int ith = 0, ord = 0; ith < DIM; ith++) {
      if (ith != idim) {
        cpyVec(edge[ord], boxEdge[ith]);
        ord++;
      }
    }
    cpyVec(edge[DIM - 1], boxEdge[idim]);

    for (int ith = 0; ith < DIM; ith++) {
      vecUnit(edge[ith], edge[ith]);
      for (int jth = ith + 1; jth < DIM; jth++) {
        double dotProd = 0.0;
        vecDot(dotProd, edge[ith], edge[jth]);
        vecScaleAdd(edge[jth], edge[jth], -dotProd, edge[ith]);
      }
    }

    vecDot(distPlane[idim], boxEdge[idim], edge[DIM - 1]);
    distPlane[idim] = fabs(distPlane[idim]);
  }
#endif

  for (int idim = 0; idim < DIM; idim++) {
    zerosVec(edge[idim]);
    edge[idim][idim] = 1.0;
  }
  for (int idim = DIM - 1; idim >= 0; idim--) {
    for (int ith = idim; ith < DIM - 1; ith++) {
      cpyVec(edge[ith], boxEdge[ith + 1]);
      memset(edge[ith], '\0', idim * sizeof(double));
    }
    cpyVec(edge[DIM - 1], boxEdge[idim]);
    memset(edge[DIM - 1], '\0', idim * sizeof(double));

    for (int ith = idim; ith < DIM; ith++) {
      vecUnit(edge[ith], edge[ith]);
      for (int jth = ith + 1; jth < DIM; jth++) {
        double dotProd = 0.0;
        vecDot(dotProd, edge[ith], edge[jth]);
        vecScaleAdd(edge[jth], edge[jth], -dotProd, edge[ith]);
      }
    }

    vecDot(distPlane[idim], boxEdge[idim], edge[DIM - 1]);
    distPlane[idim] = fabs(distPlane[idim]);
  }
}

void sortParticle(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  if (!update->nebrList.doSort) return;

  NebrList *nebrList = &update->nebrList;
  doubleVector distPlane;
  calcDistBoxPlane(distPlane, box->boxEdge);

  double sysRcs = particle->meanDiameter * (1.0 + nebrList->skinSet);
  nebrList->totBin4sort = 1;
  intVector nbin4sort;
  for (int idim = 0; idim < DIM; idim++) {
    nbin4sort[idim] = (int)floor(distPlane[idim] / sysRcs);
    nbin4sort[idim] = cpuMax(nbin4sort[idim], 1);
    nebrList->totBin4sort *= nbin4sort[idim];
  }
  if (nebrList->totBin4sort > nebrList->allocBin4sort) {
    nebrList->binHead4sort = (int *)realloc(
        nebrList->binHead4sort, nebrList->totBin4sort * sizeof(int));
    nebrList->allocBin4sort = nebrList->totBin4sort;
  }

  for (int ibin = 0; ibin < nebrList->totBin4sort; ibin++) {
    nebrList->binHead4sort[ibin] = -1;
  }
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    doubleVecPtr pos = particle->pos[iatom];

    intVector cIdx;
    doubleVector lamda;
    MatMulVec(lamda, box->invBoxH, pos);
    vecShiftAll(lamda, 0.5);
    int binIdx = 0;
    for (int idim = DIM - 1; idim >= 0; idim--) {
      cIdx[idim] = (int)floor(lamda[idim] * nbin4sort[idim]);
      cIdx[idim] = (cIdx[idim] < 0 ? nbin4sort[idim] - 1 : cIdx[idim]);
      cIdx[idim] = (cIdx[idim] >= nbin4sort[idim] ? 0 : cIdx[idim]);
      binIdx = binIdx * nbin4sort[idim] + cIdx[idim];
    }
    nebrList->binList4sort[iatom] = nebrList->binHead4sort[binIdx];
    nebrList->binHead4sort[binIdx] = iatom;
  }

  for (int idx = 0, noid = 0; idx < nebrList->totBin4sort; idx++) {
    int oid = nebrList->binHead4sort[idx];
    while (oid != -1) {
      nebrList->oid2nid[oid] = noid;
      noid++;
      oid = nebrList->binList4sort[oid];
    }
  }

  exchange_doubleVector(particle->pos, (doubleVector *)nebrList->buffer,
                        nebrList->oid2nid, particle->nAtom);
  memcpy(particle->pos, nebrList->buffer,
         particle->nAtom * sizeof(doubleVector));

  exchange_doubleVector(particle->veloc, (doubleVector *)nebrList->buffer,
                        nebrList->oid2nid, particle->nAtom);
  memcpy(particle->veloc, nebrList->buffer,
         particle->nAtom * sizeof(exchange_doubleVector));

  exchange_intVector(particle->img, (intVector *)nebrList->buffer,
                     nebrList->oid2nid, particle->nAtom);
  memcpy(particle->img, nebrList->buffer, particle->nAtom * sizeof(intVector));

  exchange_double(particle->mass, (double *)nebrList->buffer, nebrList->oid2nid,
                  particle->nAtom);
  memcpy(particle->mass, nebrList->buffer, particle->nAtom * sizeof(double));

  exchange_double(particle->diameterScale, (double *)nebrList->buffer,
                  nebrList->oid2nid, particle->nAtom);
  memcpy(particle->diameterScale, nebrList->buffer,
         particle->nAtom * sizeof(double));

  exchange_int(particle->type, (int *)nebrList->buffer, nebrList->oid2nid,
               particle->nAtom);
  memcpy(particle->type, nebrList->buffer, particle->nAtom * sizeof(int));

  exchange_int(particle->id2tag, (int *)nebrList->buffer, nebrList->oid2nid,
               particle->nAtom);
  memcpy(particle->id2tag, nebrList->buffer, particle->nAtom * sizeof(int));

  for (int aid = 0; aid < particle->nAtom; aid++) {
    particle->tag2id[particle->id2tag[aid]] = aid;
  }

  update->nebrList.isValid = false;

  // reset
  nebrList->maxAtomPerBin = 0;
}
void initNebrList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (nebrList->isInit) return;

  if (nebrList->xyzHold == NULL) {
    nebrList->xyzHold =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
  }
  if (nebrList->nNebr == NULL) {
    nebrList->nNebr = (int2 *)calloc(particle->nAtom, sizeof(int2));
  }

  if (nebrList->binList4sort == NULL) {
    nebrList->binList4sort = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (nebrList->oid2nid == NULL) {
    nebrList->oid2nid = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (nebrList->buffer == NULL) {
    nebrList->buffer = calloc(particle->nAtom, sizeof(doubleVector));
  }

  nebrList->isInit = true;
}

void calcBinParam(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->compelInit && box->isShapeFixed && particle->isSizeFixed &&
      nebrList->deltaAdjBin)
    return;

  doubleVector distPlane;
  calcDistBoxPlane(distPlane, box->boxEdge);

  nebrList->rskin =
      nebrList->minDiameterScale * particle->meanDiameter * nebrList->skinSet;
  double sysRcs =
      nebrList->maxDiameterScale * particle->meanDiameter + nebrList->rskin;
  double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
  nebrList->totBin = 1;
  for (int idim = 0; idim < DIM; idim++) {
    nebrList->nbin[idim] = (int)floor(distPlane[idim] / sysRcs);
    nebrList->nbin[idim] = cpuMax(nebrList->nbin[idim], 1);
    nebrList->totBin *= nebrList->nbin[idim];
    nebrList->binLen[idim] = distPlane[idim] / nebrList->nbin[idim];
    if (maxRcut >= distPlane[idim]) Abort("The Box is too samll!");
  }

  intVector xyzNum, bstart, bstop, ndb;
  int nAdjBin = 1;
  for (int idim = 0; idim < DIM; idim++) {
    xyzNum[idim] =
        (int)ceil((maxRcut + nebrList->rskin) / nebrList->binLen[idim]);
    bstart[idim] = -xyzNum[idim];
    bstop[idim] = xyzNum[idim];
    ndb[idim] = 2 * xyzNum[idim] + 1;
    nAdjBin *= ndb[idim];
  }

  if (!nebrList->isFullStyle) {
    int nRecord = nAdjBin / 2 + 1;
    if (nRecord > nebrList->nAdjBin)
      nebrList->deltaAdjBin = (intVector *)realloc(nebrList->deltaAdjBin,
                                                   nRecord * sizeof(intVector));

    nebrList->nAdjBin = 0;
    // for (int idx = nAdjBin - 1; idx >= 0; idx--) {
    for (int idx = 0; idx < nAdjBin; idx++) {
      intVector delta;
      int In = idx;
      for (int idim = 0; idim < DIM; idim++) {
        delta[idim] = In % ndb[idim] + bstart[idim];
        In = In / ndb[idim];
      }
      In = 0;
      for (int idim = DIM - 1; idim >= 0; idim--) {
        In = In * ndb[idim] + delta[idim];
      }
      if (In < 0) continue;
      cpyVec(nebrList->deltaAdjBin[nebrList->nAdjBin], delta);
      nebrList->nAdjBin++;
    }
  } else {
    if (nAdjBin > nebrList->nAdjBin)
      nebrList->deltaAdjBin = (intVector *)realloc(nebrList->deltaAdjBin,
                                                   nAdjBin * sizeof(intVector));
    nebrList->nAdjBin = 0;
    for (int idx = nAdjBin - 1; idx >= 0; idx--) {
      intVector delta;
      int In = idx;
      for (int idim = 0; idim < DIM; idim++) {
        delta[idim] = In % ndb[idim] + bstart[idim];
        In = In / ndb[idim];
      }
      cpyVec(nebrList->deltaAdjBin[nebrList->nAdjBin], delta);
      nebrList->nAdjBin++;
    }
  }

  if (nebrList->totBin > nebrList->allocBin) {
    nebrList->nAtomPerBin =
        (int *)realloc(nebrList->nAtomPerBin, nebrList->totBin * sizeof(int));
    nebrList->allocBin = nebrList->totBin;
  }

  nebrList->compelInit = false;
}
void binParticle(Box *box, Particle *particle, Thermo *thermo, Update *update,
                 Variable *var) {
  NebrList *nebrList = &update->nebrList;

  calcBinParam(box, particle, thermo, update, var);

  bool overFlow = false;
  do {
    overFlow = false;
    memset(nebrList->nAtomPerBin, '\0', nebrList->totBin * sizeof(int));
    int maxPartPerBin = nebrList->maxAtomPerBin;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      idPosRadius data;
      data.id = iatom;
      cpyVec(data.pos, particle->pos[iatom]);
      data.radius =
          0.5 * particle->meanDiameter * particle->diameterScale[iatom];

      intVector cIdx;
      doubleVector lamda;
      MatMulVec(lamda, box->invBoxH, data.pos);
      vecShiftAll(lamda, 0.5);
      int binIdx = 0;
      for (int idim = DIM - 1; idim >= 0; idim--) {
        cIdx[idim] = (int)floor(lamda[idim] * nebrList->nbin[idim]);
        cIdx[idim] = (cIdx[idim] < 0 ? nebrList->nbin[idim] - 1 : cIdx[idim]);
        cIdx[idim] = (cIdx[idim] >= nebrList->nbin[idim] ? 0 : cIdx[idim]);
        binIdx = binIdx * nebrList->nbin[idim] + cIdx[idim];
      }

      if (nebrList->nAtomPerBin[binIdx] < nebrList->maxAtomPerBin) {
        nebrList->ipr[nebrList->maxAtomPerBin * binIdx +
                      nebrList->nAtomPerBin[binIdx]] = data;

        nebrList->nAtomPerBin[binIdx]++;
      } else {
        nebrList->nAtomPerBin[binIdx]++;
        maxPartPerBin = (nebrList->nAtomPerBin[binIdx] > maxPartPerBin)
                            ? nebrList->nAtomPerBin[binIdx]
                            : maxPartPerBin;
        overFlow = true;
      }
    }
    if (overFlow) {
      // realloc
      nebrList->ipr = (idPosRadius *)realloc(
          nebrList->ipr, particle->nAtom * maxPartPerBin * sizeof(idPosRadius));
      nebrList->maxAtomPerBin = maxPartPerBin;
    }
  } while (overFlow);
}
void constructList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  NebrList *nebrList = &update->nebrList;
  memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int2));
  for (int bidx = 0, cntNebr = 0; bidx < nebrList->totBin; bidx++) {
    intVector binIdxVec;
    int In = bidx;
    for (int idim = 0; idim < DIM; idim++) {
      binIdxVec[idim] = In % nebrList->nbin[idim];
      In = In / nebrList->nbin[idim];
    }
    for (int ith = 0; ith < nebrList->nAtomPerBin[bidx]; ith++) {
      idPosRadius idata = nebrList->ipr[bidx * nebrList->maxAtomPerBin + ith];
      double iRadiusRskin = idata.radius + nebrList->rskin;
      nebrList->nNebr[idata.id].first = cntNebr;
      for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
        intVector binAdjVec;
        for (int idim = 0; idim < DIM; idim++) {
          binAdjVec[idim] =
              (binIdxVec[idim] + nebrList->deltaAdjBin[adj][idim] +
               nebrList->nbin[idim]) %
              nebrList->nbin[idim];
        }
        int adjBidx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
          adjBidx = adjBidx * nebrList->nbin[idim] + binAdjVec[idim];
        }
        for (int jth = 0; jth < nebrList->nAtomPerBin[adjBidx]; jth++) {
          idPosRadius jdata =
              nebrList->ipr[adjBidx * nebrList->maxAtomPerBin + jth];
          if (adjBidx == bidx && jdata.id <= idata.id) continue;

          doubleVector dRij;
          vecSub(dRij, idata.pos, jdata.pos);
          PBC(dRij, box);
          double rijP2;
          vecNormP2(rijP2, dRij);
          double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
          if (rijP2 > RcutRskinP2) continue;
          if (cntNebr == nebrList->maxAllocNebr) {
            nebrList->maxAllocNebr += 1024;
            int *tmp = (int *)realloc(nebrList->list,
                                      nebrList->maxAllocNebr * sizeof(int));
            if (tmp == NULL) Abort("realloc nebrList failed!");
            nebrList->list = tmp;
          }
          nebrList->list[cntNebr++] = jdata.id;
        }
      }
      nebrList->nNebr[idata.id].second = cntNebr;
    }
  }
}
void buildNebrList(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (nebrList->isFullStyle) {
    nebrList->compelInit = true;
  }
  nebrList->isFullStyle = false;

  initNebrList(box, particle, thermo, update, var);

  adjustImg(box, particle, thermo, update, var);
  if ((nebrList->nBuild % 100 == 0)) {
    sortParticle(box, particle, thermo, update, var);
  }
  binParticle(box, particle, thermo, update, var);

  nebrList->meanDiameterHold = particle->meanDiameter;
  cpyUptriMat(nebrList->invBoxHold, box->invBoxH);
  memcpy(nebrList->xyzHold, particle->pos,
         particle->nAtom * sizeof(doubleVector));

  constructList(box, particle, thermo, update, var);

  nebrList->nBuild++;
  nebrList->isValid = true;
  nebrList->cntForce = 0;
}
bool isNebrListValid(Box *box, Particle *particle, Thermo *thermo,
                     Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->isInit) {
    nebrList->isValid = false;
    return false;
  }
  if (!nebrList->isValid) {
    nebrList->isValid = false;
    return false;
  }

  if (particle->isSizeFixed && box->isShapeFixed) {
    double maxShiftP2 = 0.25 * nebrList->rskin * nebrList->rskin;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      doubleVector dRt0;
      vecSub(dRt0, particle->pos[iatom], nebrList->xyzHold[iatom]);
      double rP2;
      vecNormP2(rP2, dRt0);
      if (rP2 >= maxShiftP2) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else if (!particle->isSizeFixed && box->isShapeFixed) {
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      doubleVector dRt0;
      vecSub(dRt0, particle->pos[iatom], nebrList->xyzHold[iatom]);
      double ri0;
      vecNorm(ri0, dRt0);
      double dsig = (particle->meanDiameter - nebrList->meanDiameterHold) *
                    particle->diameterScale[iatom];

      if (ri0 + dsig >= nebrList->rskin * 0.5) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else if (!box->isShapeFixed && particle->isSizeFixed) {
    uptriMat trans;  // trans =  Ht * inv(H0)
    MatMulMat(trans, box->boxH, nebrList->invBoxHold);
    diagMat volTrans;
    getDiag(volTrans, trans);

    double sfact = 1;
    minElement(sfact, volTrans);

    invDiag(volTrans, volTrans);
    uptriMat shapeTrans;
    diagMulMat(shapeTrans, volTrans, trans);

    for (int idim = 0; idim < DIM; idim++) {
      for (int jdim = idim + 1; jdim < DIM; jdim++) {
        double eps2 = pow(shapeTrans[spaceIdx(idim, jdim)], 2);
        sfact *= sqrt(1.0 + eps2 / 2.0 - sqrt(eps2 + eps2 * eps2 / 4.0));
      }
    }

    double rs_eff = nebrList->rskin;
    if (sfact >= 1.0) {
      rs_eff = sfact * rs_eff +
               (sfact - 1.0) *
                   (2.0 * nebrList->minDiameterScale * particle->meanDiameter);
    } else {
      rs_eff = sfact * rs_eff +
               (sfact - 1.0) *
                   (2.0 * nebrList->maxDiameterScale * particle->meanDiameter);
    }
    if (rs_eff <= 0) {
      nebrList->isValid = false;
      return false;
    }
    double maxShiftP2 = 0.25 * rs_eff * rs_eff;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      doubleVector posHold;
      cpyVec(posHold, nebrList->xyzHold[iatom]);
      MatMulVec(posHold, trans, posHold);

      doubleVector dR;
      vecSub(dR, particle->pos[iatom], posHold);
      double rP2;
      vecNormP2(rP2, dR);

      if (rP2 >= maxShiftP2) {
        nebrList->isValid = false;
        return false;
      }
    }
  } else
    Abort("Not Code!");
  nebrList->isValid = true;
  return true;
}

void calcForce(Box *box, Particle *particle, Thermo *thermo, Update *update,
               Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (particle->isForceValid) return;

  bool isListValid = isNebrListValid(box, particle, thermo, update, var);
  if (nebrList->isFullStyle) {
    nebrList->isValid = false;
    isListValid = false;
  }
  if (!isListValid) {
    if (nebrList->cntForce < 10) {
      nebrList->skinSet *= 1.1;
      nebrList->skinSet = cpuMin(nebrList->skinSet, 0.5);
      nebrList->compelInit = true;
    } else if (nebrList->cntForce > 20) {
      nebrList->skinSet *= 0.9;
      nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__);
      nebrList->compelInit = true;
    }
  } else {
    if (nebrList->skinSet > __minSkinSet__ * (1 + 1E-10) &&
        nebrList->cntForce > 20) {
      nebrList->skinSet *= 0.9;
      nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__);
      nebrList->compelInit = true;
      isListValid = false;
    }
  }
  if (!isListValid) {
    buildNebrList(box, particle, thermo, update, var);
  }

  doubleVector *force = particle->force;
  memset(force, '\0', particle->nAtom * sizeof(doubleVector));
  thermo->Epair = thermo->virialPair = thermo->lapU = 0;
  zerosUptriMat(thermo->ptensor);

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    doubleVector iforce;
    cpyVec(iforce, force[iatom]);
    doubleVector iPos;
    cpyVec(iPos, particle->pos[iatom]);
    double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
    for (int jth = nebrList->nNebr[iatom].first;
         jth < nebrList->nNebr[iatom].second; jth++) {
      int jatom = nebrList->list[jth];
      double sRc =
          iRc + particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
      doubleVector dRij;
      vecSub(dRij, iPos, particle->pos[jatom]);
      PBC(dRij, box);

      double rijP2;
      vecNormP2(rijP2, dRij);
      if (rijP2 >= sRc * sRc) continue;
      double rij = sqrt(rijP2);
      double rdivsig = rij / sRc;
      double fpair = (1.0 - rdivsig) / rij / sRc;
      vecScaleAdd(iforce, iforce, fpair, dRij);
      vecScaleAdd(force[jatom], force[jatom], -fpair, dRij);

      thermo->Epair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
      // thermo->lapU += (6.0 * rdivsig * rdivsig - 4.0 * rdivsig) / rijP2;

      for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
          thermo->ptensor[spaceIdx(idim, jdim)] +=
              fpair * dRij[idim] * dRij[jdim];
        }
      }
    }
    cpyVec(force[iatom], iforce);
  }

  thermo->Epair = thermo->Epair / particle->nAtom;
  thermo->lapU = thermo->lapU / particle->nAtom;

  thermo->virialPair = 0.0;
  thermo->pressure = 0.0;
  for (int idim = 0; idim < DIM; idim++) {
    thermo->virialPair += thermo->ptensor[spaceIdx(idim, idim)];
    for (int jdim = idim; jdim < DIM; jdim++) {
      thermo->ptensor[spaceIdx(idim, jdim)] /= box->volume;
    }
    thermo->pressure += thermo->ptensor[spaceIdx(idim, idim)];
  }
  thermo->pressure /= DIM;

  thermo->Edone = thermo->Pdone = true;

  particle->isForceValid = true;
  nebrList->cntForce++;
  nebrList->nForce++;
}

//===============================================================================
void instant_inflate(Box *box, Particle *particle, Thermo *thermo,
                     Update *update, Variable *var, double deltaVF) {
  calcVolInfo(box, particle, thermo, update, var);
  double volfrac_target = thermo->volFrac + deltaVF;
  double sfact = pow(volfrac_target / thermo->volFrac, 1.0 / DIM);
  particle->meanDiameter *= sfact;
  thermo->volFrac = pow(sfact, DIM) * thermo->volFrac;

  calcVolInfo(box, particle, thermo, update, var);

  setBasicInfo(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
#ifdef __triBox__
void instant_simpShearXz(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var, double deltaGamma) {
  adjustImg(box, particle, thermo, update, var);

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    particle->pos[iatom][0] += particle->pos[iatom][DIM - 1] * deltaGamma;
  }

  box->boxH[spaceIdx(0, DIM - 1)] +=
      deltaGamma * box->boxH[spaceIdx(DIM - 1, DIM - 1)];

  setBoxPara(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
void instant_deformation(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var, uptriMat transMartix) {
  // Rt = R0 + transMartix * R0;
  adjustImg(box, particle, thermo, update, var);

  uptriMat trans;
  cpyUptriMat(trans, transMartix);
  for (int idim = 0; idim < DIM; idim++) {
    trans[spaceIdx(idim, idim)] += 1.0;
  }

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    MatMulVec(particle->pos[iatom], trans, particle->pos[iatom]);
  }

  // update box
  for (int iedge = 0; iedge < DIM; iedge++) {
    for (int idim = 0; idim < DIM; idim++) {
      box->boxEdge[iedge][idim] *= trans[spaceIdx(idim, idim)];
      for (int jdim = idim + 1; jdim < DIM; jdim++) {
        box->boxEdge[iedge][idim] +=
            trans[spaceIdx(idim, jdim)] * box->boxEdge[iedge][jdim];
      }
    }
  }

  for (int idim = 0; idim < DIM; idim++) {
    for (int jdim = idim; jdim < DIM; jdim++) {
      box->boxH[spaceIdx(idim, jdim)] = box->boxEdge[jdim][idim];
    }
  }
  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);

  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
}
#endif

void normaliseBox(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  double sfact = pow(box->volume, 1.0 / DIM);
  particle->meanDiameter /= sfact;

  for (int idim = 0; idim < DIM; idim++) {
    for (int jdim = idim; jdim < DIM; jdim++) {
      box->boxH[spaceIdx(idim, jdim)] /= sfact;
    }
  }

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    scaleVec(particle->pos[iatom], 1.0 / sfact, particle->pos[iatom]);
  }

  thermo->isInit = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  particle->isForceValid = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.isInit = false;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
}

//===============================================================================
void constructList_full(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int2));
  for (int bidx = 0, cntNebr = 0; bidx < nebrList->totBin; bidx++) {
    intVector binIdxVec;
    int In = bidx;
    for (int idim = 0; idim < DIM; idim++) {
      binIdxVec[idim] = In % nebrList->nbin[idim];
      In = In / nebrList->nbin[idim];
    }
    for (int ith = 0; ith < nebrList->nAtomPerBin[bidx]; ith++) {
      idPosRadius idata = nebrList->ipr[bidx * nebrList->maxAtomPerBin + ith];
      double iRadiusRskin = idata.radius + nebrList->rskin;
      nebrList->nNebr[idata.id].first = cntNebr;
      for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
        intVector binAdjVec;
        for (int idim = 0; idim < DIM; idim++) {
          binAdjVec[idim] =
              (binIdxVec[idim] + nebrList->deltaAdjBin[adj][idim] +
               nebrList->nbin[idim]) %
              nebrList->nbin[idim];
        }
        int adjBidx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
          adjBidx = adjBidx * nebrList->nbin[idim] + binAdjVec[idim];
        }
        for (int jth = 0; jth < nebrList->nAtomPerBin[adjBidx]; jth++) {
          idPosRadius jdata =
              nebrList->ipr[adjBidx * nebrList->maxAtomPerBin + jth];
          if (idata.id == jdata.id) continue;

          doubleVector dRij;
          vecSub(dRij, idata.pos, jdata.pos);
          PBC(dRij, box);
          double rijP2;
          vecNormP2(rijP2, dRij);
          double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
          if (rijP2 > RcutRskinP2) continue;
          if (cntNebr == nebrList->maxAllocNebr) {
            nebrList->maxAllocNebr += 1024;
            int *tmp = (int *)realloc(nebrList->list,
                                      nebrList->maxAllocNebr * sizeof(int));
            if (tmp == NULL) Abort("realloc nebrList failed!");
            nebrList->list = tmp;
          }
          nebrList->list[cntNebr++] = jdata.id;
        }
      }
      nebrList->nNebr[idata.id].second = cntNebr;
    }
  }
}
void buildNebrList_full(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  NebrList *nebrList = &update->nebrList;
  if (!nebrList->isFullStyle) {
    nebrList->compelInit = true;
  }
  nebrList->isFullStyle = true;

  initNebrList(box, particle, thermo, update, var);

  adjustImg(box, particle, thermo, update, var);
  if ((nebrList->nBuild % 100 == 0)) {
    sortParticle(box, particle, thermo, update, var);
  }
  binParticle(box, particle, thermo, update, var);

  nebrList->meanDiameterHold = particle->meanDiameter;
  cpyUptriMat(nebrList->invBoxHold, box->invBoxH);
  memcpy(nebrList->xyzHold, particle->pos,
         particle->nAtom * sizeof(doubleVector));

  constructList_full(box, particle, thermo, update, var);

  nebrList->nBuild++;
  nebrList->isValid = true;
  nebrList->cntForce = 0;
}

//===============================================================================
void backupSimInfo(Box *box, Particle *particle, Box *bakBox,
                   Particle *bakParticle) {
  if (bakParticle->nAtom != particle->nAtom) {
    bakParticle->nAtom = particle->nAtom;

    bakParticle->pos = (doubleVector *)realloc(
        bakParticle->pos, bakParticle->nAtom * sizeof(doubleVector));
    bakParticle->veloc = (doubleVector *)realloc(
        bakParticle->veloc, bakParticle->nAtom * sizeof(doubleVector));
    bakParticle->mass = (double *)realloc(bakParticle->mass,
                                          bakParticle->nAtom * sizeof(double));
    bakParticle->img = (intVector *)realloc(
        bakParticle->img, bakParticle->nAtom * sizeof(intVector));
    bakParticle->type =
        (int *)realloc(bakParticle->type, bakParticle->nAtom * sizeof(int));
    bakParticle->diameterScale = (double *)realloc(
        bakParticle->diameterScale, bakParticle->nAtom * sizeof(double));
    bakParticle->id2tag =
        (int *)realloc(bakParticle->id2tag, bakParticle->nAtom * sizeof(int));
    bakParticle->tag2id =
        (int *)realloc(bakParticle->tag2id, bakParticle->nAtom * sizeof(int));
  }
  if (bakParticle->nAtomType != particle->nAtomType) {
    bakParticle->nAtomType = particle->nAtomType;
    bakParticle->massPerType = (double *)realloc(
        bakParticle->massPerType, bakParticle->nAtomType * sizeof(double));
  }

  memcpy(bakParticle->pos, particle->pos,
         bakParticle->nAtom * sizeof(doubleVector));
  memcpy(bakParticle->veloc, particle->veloc,
         bakParticle->nAtom * sizeof(doubleVector));
  memcpy(bakParticle->mass, particle->mass,
         bakParticle->nAtom * sizeof(double));
  memcpy(bakParticle->img, particle->img,
         bakParticle->nAtom * sizeof(intVector));
  memcpy(bakParticle->type, particle->type, bakParticle->nAtom * sizeof(int));
  memcpy(bakParticle->diameterScale, particle->diameterScale,
         bakParticle->nAtom * sizeof(double));
  memcpy(bakParticle->id2tag, particle->id2tag,
         bakParticle->nAtom * sizeof(int));
  memcpy(bakParticle->tag2id, particle->tag2id,
         bakParticle->nAtom * sizeof(int));

  memcpy(bakParticle->massPerType, particle->massPerType,
         bakParticle->nAtomType * sizeof(double));

  memcpy(&bakBox->boxH, &box->boxH, sizeof(uptriMat));
}
void restoreSimInfo(Box *box, Particle *particle, Thermo *thermo,
                    Update *update, Variable *var, Box *bakBox,
                    Particle *bakParticle) {
  // restore data
  backupSimInfo(bakBox, bakParticle, box, particle);

  // clear
  particle->isForceValid = false;

  thermo->isInit = false;
  thermo->Edone = false;
  thermo->Pdone = false;
  thermo->Tdone = false;

  update->nebrList.isValid = false;
  update->nebrList.isInit = false;
  update->nebrList.compelInit = true;

  // reinit
  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}
void freeSimInfo(Box *bakBox, Particle *bakParticle) {
  safeFree(bakParticle->pos);
  safeFree(bakParticle->veloc);
  safeFree(bakParticle->mass);
  safeFree(bakParticle->img);
  safeFree(bakParticle->type);
  safeFree(bakParticle->diameterScale);
  safeFree(bakParticle->id2tag);
  safeFree(bakParticle->tag2id);
}

//===============================================================================
contactInfo *initContactInfo(Box *box, Particle *particle, Thermo *thermo,
                             Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  contactInfo *cinfo = NULL;
  if (whichTool >= 0) {
    cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  } else {
    cinfo = (contactInfo *)calloc(sizeof(contactInfo), 1);
    addToolkit(&thermo->toolkit, (void *)cinfo, "contactInfo");
  }
  if (particle->nAtom <= 0) Abort("No Atom!");

  if (cinfo->isRattler == NULL) {
    cinfo->isRattler = (bool *)calloc(particle->nAtom, sizeof(bool));
  }
  if (cinfo->nCoordNumExRattler == NULL) {
    cinfo->nCoordNumExRattler = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (cinfo->nCoordNum == NULL) {
    cinfo->nCoordNum = (int *)calloc(particle->nAtom, sizeof(int));
  }
  if (cinfo->forceExRattler == NULL) {
    cinfo->forceExRattler =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
  }

  return cinfo;
}
void computeCoordination(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  if (whichTool < 0) return;
  contactInfo *cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  NebrList *nebrList = &update->nebrList;
  if (!isNebrListValid(box, particle, thermo, update, var)) {
    buildNebrList(box, particle, thermo, update, var);
  }
  if (cinfo->allocMem4CoordNum < nebrList->maxAllocNebr) {
    cinfo->allocMem4CoordNum = nebrList->maxAllocNebr;
    cinfo->isNebrContact = (bool *)realloc(
        cinfo->isNebrContact, cinfo->allocMem4CoordNum * sizeof(bool));
  }

  memset(cinfo->isRattler, '\0', particle->nAtom * sizeof(bool));
  memset(cinfo->nCoordNumExRattler, '\0', particle->nAtom * sizeof(int));
  memset(cinfo->nCoordNum, '\0', particle->nAtom * sizeof(int));
  memset(cinfo->isNebrContact, '\0', cinfo->allocMem4CoordNum * sizeof(bool));

  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    for (int jth = nebrList->nNebr[iatom].first;
         jth < nebrList->nNebr[iatom].second; jth++) {
      int jatom = nebrList->list[jth];

      doubleVector dRij;
      double sRc =
          (particle->diameterScale[iatom] + particle->diameterScale[jatom]) *
          particle->meanDiameter * 0.5;
      vecSub(dRij, particle->pos[iatom], particle->pos[jatom]);
      PBC(dRij, box);
      double rij;
      vecNorm(rij, dRij);
      if (rij < sRc) {
        cinfo->isNebrContact[jth] = true;
        cinfo->nCoordNum[jatom]++;
        cinfo->nCoordNum[iatom]++;
      }
    }
  }

  int cordNumWR = 0, cordNumNR = 0, nRattler = 0;
  bool done = true;
  memcpy(cinfo->nCoordNumExRattler, cinfo->nCoordNum,
         particle->nAtom * sizeof(int));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->nCoordNumExRattler[iatom] < (DIM + 1)) {
      cinfo->isRattler[iatom] = true;
      done = false;
      nRattler++;
    }
    cordNumWR += cinfo->nCoordNum[iatom];
    cordNumNR += cinfo->nCoordNumExRattler[iatom];
  }
  while (!done) {
    done = true;
    memset(cinfo->nCoordNumExRattler, '\0', particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      if (cinfo->isRattler[iatom]) continue;
      for (int jth = nebrList->nNebr[iatom].first;
           jth < nebrList->nNebr[iatom].second; jth++) {
        int jatom = nebrList->list[jth];
        if (!cinfo->isNebrContact[jth]) continue;
        if (!cinfo->isRattler[jatom]) {
          cinfo->nCoordNumExRattler[iatom]++;
          cinfo->nCoordNumExRattler[jatom]++;
        }
      }
    }

    cordNumNR = 0;
    nRattler = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      if (cinfo->nCoordNumExRattler[iatom] < (DIM + 1)) {
        if (!cinfo->isRattler[iatom]) {
          done = false;
        }
        cinfo->isRattler[iatom] = true;
        nRattler++;
      }
      cordNumNR += cinfo->nCoordNumExRattler[iatom];
    }
  }
  cinfo->nRattler = nRattler;
  cinfo->aveCoordNum = (double)cordNumWR / (double)particle->nAtom;
  if (nRattler != particle->nAtom) {
    cinfo->aveCoordNumExRattler =
        (double)cordNumNR / (double)(particle->nAtom - nRattler);
  } else {
    cinfo->aveCoordNumExRattler = 0;
  }

  if (particle->nAtom == cinfo->nRattler) {
    cinfo->meanForceExRattler = 0;
    cinfo->lapUexRattler = 0;
    return;
  }
  doubleVector *force = cinfo->forceExRattler;
  // double lapU = 0;
  memset(force, '\0', particle->nAtom * sizeof(doubleVector));
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->isRattler[iatom]) continue;
    doubleVector iforce;
    cpyVec(iforce, force[iatom]);
    doubleVector iPos;
    cpyVec(iPos, particle->pos[iatom]);
    double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
    for (int jth = nebrList->nNebr[iatom].first;
         jth < nebrList->nNebr[iatom].second; jth++) {
      int jatom = nebrList->list[jth];
      if (cinfo->isRattler[jatom]) continue;

      doubleVector dRij;
      double sRc =
          iRc + particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
      vecSub(dRij, iPos, particle->pos[jatom]);
      PBC(dRij, box);

      double rijP2;
      vecNormP2(rijP2, dRij);
      if (rijP2 >= sRc * sRc) continue;
      double rij = sqrt(rijP2);
      double rdivsig = rij / sRc;
      double fpair = (1.0 - rdivsig) / rij / sRc;
      vecScaleAdd(iforce, iforce, fpair, dRij);
      vecScaleAdd(force[jatom], force[jatom], -fpair, dRij);

      // lapU += (6.0 * rdivsig * rdivsig - 4.0 * rdivsig) / rijP2;
    }
    cpyVec(force[iatom], iforce);
  }

  double sumForce = 0;
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    if (cinfo->isRattler[iatom]) continue;
    double dotPrd;
    vecDot(dotPrd, force[iatom], force[iatom]);
    sumForce += sqrt(dotPrd);
  }
  cinfo->meanForceExRattler = sumForce / thermo->forceUnits /
                              (double)(particle->nAtom - cinfo->nRattler);

  // cinfo->lapUexRattler = lapU / (double)(particle->nAtom - cinfo->nRattler);
}
void finalizeContactInfo(Box *box, Particle *particle, Thermo *thermo,
                         Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "contactInfo");
  if (whichTool < 0) return;
  contactInfo *cinfo = (contactInfo *)thermo->toolkit.toolkit[whichTool];
  safeFree(cinfo->isRattler);
  safeFree(cinfo->isNebrContact);
  safeFree(cinfo->nCoordNum);
  safeFree(cinfo->nCoordNumExRattler);
  safeFree(cinfo->forceExRattler);
  delToolkit(&thermo->toolkit, "contactInfo");
}

dumpInfo *initDump(Box *box, Particle *particle, Thermo *thermo, Update *update,
                   Variable *var) {
  int whichVar = findVariable(var, "dump");
  if (whichVar < 0) return (dumpInfo *)NULL;
  dumpInfo *dinfo = (dumpInfo *)calloc(sizeof(dumpInfo), 1);
  addToolkit(&thermo->toolkit, (void *)dinfo, "dumpInfo");

  char fname[4096];
  sprintf(fname, "%s/dump_%s.bin", var->cwd, var->sf);
  dinfo->fdump = createReadWriteFile(fname);

  // header
  char str[32];
  memset(str, '\0', 32 * sizeof(char));
  sprintf(str, "Revised Binary File");
  str[31] = '\n';
  fwrite(str, sizeof(char), 32, dinfo->fdump);
  dinfo->revNum = dumpFileRevNum;
  fwrite(&dinfo->revNum, sizeof(int), 1, dinfo->fdump);

  fwrite(&box->dim, sizeof(int), 1, dinfo->fdump);
  fwrite(&particle->nAtom, sizeof(int), 1, dinfo->fdump);
  fwrite(&particle->nAtomType, sizeof(int), 1, dinfo->fdump);
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    fwrite(&particle->diameterScale[idx], sizeof(double), 1, dinfo->fdump);
  }
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    fwrite(&particle->type[idx], sizeof(int), 1, dinfo->fdump);
  }
  return dinfo;
}
void dump(Box *box, Particle *particle, Thermo *thermo, Update *update,
          Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "dumpInfo");
  if (whichTool < 0) return;

  FILE *fdump = ((dumpInfo *)thermo->toolkit.toolkit[whichTool])->fdump;
  fwrite(&particle->nAtom, sizeof(int), 1, fdump);
  fwrite(&box->boxH, sizeof(uptriMat), 1, fdump);
  fwrite(&particle->meanDiameter, sizeof(double), 1, fdump);
  for (int iatom = 0; iatom < particle->nAtom; iatom++) {
    int idx = particle->tag2id[iatom];
    doubleVector uxyz;
    unwrapPos(uxyz, particle->pos[idx], particle->img[idx], box->boxH);

    fwrite(&uxyz, sizeof(doubleVector), 1, fdump);
  }
}
void finalizeDump(Box *box, Particle *particle, Thermo *thermo, Update *update,
                  Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "dumpInfo");
  if (whichTool < 0) return;
  dumpInfo *dinfo = (dumpInfo *)thermo->toolkit.toolkit[whichTool];
  safeCloseFile(dinfo->fdump);
  delToolkit(&thermo->toolkit, "dumpInfo");
}

BinDumpFile *initMapBinFile(Box *box, Particle *particle, Thermo *thermo,
                            Update *update, Variable *var) {
  BinDumpFile *binFile = (BinDumpFile *)calloc(sizeof(BinDumpFile), 1);
  addToolkit(&thermo->toolkit, (void *)binFile, "BinDumpFile");

  int whichVar = findVariable(var, "read_dump");
  if (whichVar < 0 || var->cmd[whichVar].cmdArgc != 1) {
    Abort("--read_dump dump.bin");
  }
  if (binFile->dataSection != NULL) Abort("Fatal Error!");

  binFile->fd = open(var->cmd[whichVar].cmdArgv[0], O_RDONLY);
  fstat(binFile->fd, &binFile->binStat);
  binFile->dataSection = mmap(NULL, binFile->binStat.st_size, PROT_READ,
                              MAP_PRIVATE, binFile->fd, 0);
  if (binFile->dataSection == MAP_FAILED) {
    Abort("Map Binary file Failed !\n");
  }

  // header
  char *header = (char *)binFile->dataSection;
  void *data = binFile->dataSection;
  int *nElement = (int *)data;
  if (strcmp(header, "Revised Binary File") == 0) {
    int *revNum = (int *)(header + 32);
    binFile->revNum = revNum[0];
    if (binFile->revNum == 0) Abort("Wrong Binary File!");

    int *dim = (int *)(revNum + 1);
    if (dim[0] != DIM) {
      Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0],
            DIM);
    }
    nElement = (int *)(dim + 1);
  } else {
    if (DIM != 3) {
      Abort("The dumpfile is for d = 3, while the code is for d = %d!", DIM);
    }
    binFile->revNum = 0;
  }
  box->dim = DIM;

  binFile->headerSize =
      sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
  if (binFile->revNum == 1) {
    binFile->headerSize += 32 * sizeof(char) + 2 * sizeof(int);
  }
  binFile->stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) +
                      nElement[0] * sizeof(doubleVector);
  binFile->nStep =
      (binFile->binStat.st_size - binFile->headerSize) / binFile->stepSize;
  // read first step
  {
    void *data = (void *)(binFile->dataSection);
    if (binFile->revNum == 1) {
      data = data + 32 * sizeof(char) + 2 * sizeof(int);
    }

    particle->nAtom = *((int *)data);
    data = (void *)((char *)data + sizeof(int));
    particle->nAtomType = *((int *)data);
    data = (void *)((char *)data + sizeof(int));

    if (particle->pos == NULL) {
      particle->pos =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->veloc =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->force =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
      particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
      particle->type = (int *)calloc(particle->nAtom, sizeof(int));
      particle->diameterScale =
          (double *)calloc(particle->nAtom, sizeof(double));
      particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
      particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      particle->massPerType =
          (double *)calloc(particle->nAtomType, sizeof(double));
      for (int itype = 0; itype < particle->nAtomType; itype++)
        particle->massPerType[itype] = 1.0;
    }

    memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
    data = (void *)((char *)data + sizeof(double) * particle->nAtom);
    memcpy(particle->type, data, particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->id2tag[iatom] = iatom;
      particle->tag2id[iatom] = iatom;
    }

    data = ((char *)binFile->dataSection + binFile->headerSize + sizeof(int));

    if (binFile->revNum == 0) {
      double *tmp = (double *)data;
      box->boxH[spaceIdx(0, 0)] = tmp[0];
      box->boxH[spaceIdx(1, 1)] = tmp[1];
      box->boxH[spaceIdx(2, 2)] = tmp[2];
      box->boxH[spaceIdx(1, 2)] = tmp[3];
      box->boxH[spaceIdx(0, 2)] = tmp[4];
      box->boxH[spaceIdx(0, 1)] = tmp[5];
      data = ((char *)data + 6 * sizeof(double));
    } else {
      memcpy(box->boxH, data, sizeof(uptriMat));
      data = ((char *)data + sizeof(uptriMat));
    }

    particle->meanDiameter = *((double *)data);
    data = ((char *)data + sizeof(double));

    memcpy(particle->pos, data, particle->nAtom * sizeof(doubleVector));

    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    memset(particle->force, '\0', particle->nAtom * sizeof(doubleVector));
    memset(particle->img, '\0', particle->nAtom * sizeof(intVector));

    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    particle->isForceValid = false;
    thermo->Edone = thermo->Pdone = thermo->Tdone = false;
    thermo->isInit = false;
    update->nebrList.isFullStyle = false;
    update->nebrList.isInit = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.nForce = 0;
    update->nebrList.cntForce = 0;
    update->nebrList.doSort = true;
    update->nebrList.skinSet = __minSkinSet__;

    setBoxPara(box, particle, thermo, update, var);
    setBasicInfo(box, particle, thermo, update, var);
    adjustImg(box, particle, thermo, update, var);
  }

  return binFile;
}
void finalizeMapBinFile(Box *box, Particle *particle, Thermo *thermo,
                        Update *update, Variable *var) {
  int whichTool = findToolkit(&thermo->toolkit, "BinDumpFile");
  if (whichTool < 0) return;
  BinDumpFile *binFile = (BinDumpFile *)thermo->toolkit.toolkit[whichTool];
  munmap(binFile->dataSection, binFile->binStat.st_size);
  close(binFile->fd);
  binFile->fd = -1;
  binFile->dataSection = NULL;
  delToolkit(&thermo->toolkit, "BinDumpFile");
}
void getBinDumpData(Box *box, Particle *particle, Thermo *thermo,
                    Update *update, Variable *var, int whichStep) {
  int whichTool = findToolkit(&thermo->toolkit, "BinDumpFile");
  if (whichTool < 0) return;
  BinDumpFile *binFile = (BinDumpFile *)thermo->toolkit.toolkit[whichTool];
  if (whichStep >= binFile->nStep || whichStep < 0) {
    Abort("Wrong Step!");
  }

  {
    void *data = (void *)(binFile->dataSection);
    if (binFile->revNum == 1) {
      data = data + 32 * sizeof(char) + 2 * sizeof(int);
    }
    particle->nAtom = *((int *)data);
    data = (void *)((char *)data + sizeof(int));
    particle->nAtomType = *((int *)data);
    data = (void *)((char *)data + sizeof(int));

    if (particle->pos == NULL) {
      particle->pos =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->veloc =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->force =
          (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
      particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
      particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
      particle->type = (int *)calloc(particle->nAtom, sizeof(int));
      particle->diameterScale =
          (double *)calloc(particle->nAtom, sizeof(double));
      particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
      particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
      particle->massPerType =
          (double *)calloc(particle->nAtomType, sizeof(double));
      for (int itype = 0; itype < particle->nAtomType; itype++)
        particle->massPerType[itype] = 1.0;
    }

    memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
    data = (void *)((char *)data + sizeof(double) * particle->nAtom);
    memcpy(particle->type, data, particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
      particle->id2tag[iatom] = iatom;
      particle->tag2id[iatom] = iatom;
    }

    data = ((char *)binFile->dataSection + binFile->headerSize +
            whichStep * binFile->stepSize + sizeof(int));

    if (binFile->revNum == 0) {
      double *tmp = (double *)data;
      box->boxH[spaceIdx(0, 0)] = tmp[0];
      box->boxH[spaceIdx(1, 1)] = tmp[1];
      box->boxH[spaceIdx(2, 2)] = tmp[2];
      box->boxH[spaceIdx(1, 2)] = tmp[3];
      box->boxH[spaceIdx(0, 2)] = tmp[4];
      box->boxH[spaceIdx(0, 1)] = tmp[5];
      data = ((char *)data + 6 * sizeof(double));
    } else {
      memcpy(box->boxH, data, sizeof(uptriMat));
      data = ((char *)data + sizeof(uptriMat));
    }

    particle->meanDiameter = *((double *)data);
    data = ((char *)data + sizeof(double));

    memcpy(particle->pos, data, particle->nAtom * sizeof(doubleVector));

    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    memset(particle->force, '\0', particle->nAtom * sizeof(doubleVector));
    memset(particle->img, '\0', particle->nAtom * sizeof(intVector));
  }

  box->isShapeFixed = true;
  particle->isSizeFixed = true;
  particle->isForceValid = false;
  thermo->Edone = thermo->Pdone = thermo->Tdone = false;
  thermo->isInit = false;
  update->nebrList.isFullStyle = false;
  update->nebrList.isInit = false;
  update->nebrList.isValid = false;
  update->nebrList.compelInit = true;
  update->nebrList.nForce = 0;
  update->nebrList.cntForce = 0;
  update->nebrList.doSort = true;
  update->nebrList.skinSet = __minSkinSet__;

  setBoxPara(box, particle, thermo, update, var);
  setBasicInfo(box, particle, thermo, update, var);
  adjustImg(box, particle, thermo, update, var);
}
#endif

#endif