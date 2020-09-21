#include "VectorMath.h"

// convert between lammps data file and binary file
int main(int argc, char *argv[]) {
  if (argc != 2) {
    Abort("./app datafile");
  }

  int nAtom = 0, nAtomType = 0;
  Hvoigt6 boxHvoigt = make_Hvoigt6(0, 0, 0, 0, 0, 0);
  double3 *xyz = NULL, *veloc = NULL;
  int *type = NULL;
  double *diameterScale = NULL;
  int3 *img = NULL;
  double meanDiameter = 1.0;

  FILE *fp = openReadOnlyFile(argv[1]);
  char str[4096];
  fgets(str, 4096, fp);
  if (strcmp(str, "binary") == 0)  // binary file
  {
    bool hasMeanDiameter = false;
    while (fread(str, sizeof(char), 32, fp) == 32) {
      if (strcmp(str, "atoms") == 0) {
        fread(&nAtom, sizeof(int), 1, fp);
        if (nAtom <= 0) Abort("Wrong File!");

        xyz = (double3 *)calloc(nAtom, sizeof(double3));
        img = (int3 *)calloc(nAtom, sizeof(int3));
        type = (int *)calloc(nAtom, sizeof(int));
        diameterScale = (double *)calloc(nAtom, sizeof(double));
      } else if (strcmp(str, "atom types") == 0) {
        fread(&nAtomType, sizeof(int), 1, fp);
        if (nAtomType <= 0) Abort("Wrong File!");
      } else if (strcmp(str, "box Hvoigt") == 0) {
        fread(&boxHvoigt.h0, sizeof(double), 1, fp);
        fread(&boxHvoigt.h1, sizeof(double), 1, fp);
        fread(&boxHvoigt.h2, sizeof(double), 1, fp);
        fread(&boxHvoigt.h3, sizeof(double), 1, fp);
        fread(&boxHvoigt.h4, sizeof(double), 1, fp);
        fread(&boxHvoigt.h5, sizeof(double), 1, fp);
      } else if (strcmp(str, "mean diameter") == 0) {
        fread(&meanDiameter, sizeof(double), 1, fp);
        hasMeanDiameter = true;
      } else if (strcmp(str, "Atoms") == 0) {
        for (int iatom = 0; iatom < nAtom; iatom++) {
          fread(&type[iatom], sizeof(int), 1, fp);
          fread(&diameterScale[iatom], sizeof(double), 1, fp);
          fread(&xyz[iatom].x, sizeof(double), 1, fp);
          fread(&xyz[iatom].y, sizeof(double), 1, fp);
          fread(&xyz[iatom].z, sizeof(double), 1, fp);
          fread(&img[iatom].x, sizeof(int), 1, fp);
          fread(&img[iatom].y, sizeof(int), 1, fp);
          fread(&img[iatom].z, sizeof(int), 1, fp);
        }
      } else if (strcmp(str, "Velocities") == 0) {
        veloc = (double3 *)calloc(nAtom, sizeof(double3));
        for (int iatom = 0; iatom < nAtom; iatom++) {
          fread(&veloc[iatom].x, sizeof(double), 1, fp);
          fread(&veloc[iatom].y, sizeof(double), 1, fp);
          fread(&veloc[iatom].z, sizeof(double), 1, fp);
        }
      } else
        Abort("Wrong File!");
    }

    {
      char str[4096];
      sprintf(str, "%s.data", argv[1]);
      FILE *fout = createReadWriteFile(str);

      safeFprintf(fout,
                  "LAMMPS compatible data file. atom_style: sphere (id type "
                  "diameter density x y z ix iy iz).\n\n");
      safeFprintf(fout, "\t%d atoms\n", nAtom);
      safeFprintf(fout, "\t%d atom types\n", nAtomType);

      safeFprintf(
          fout,
          "\n\t%.16g %.16g xlo xhi\n\t%.16g %.16g ylo yhi\n\t%.16g %.16g "
          "zlo zhi\n\t%.16g %.16g %.16g xy xz yz\n",
          -boxHvoigt.h0 / 2.0, boxHvoigt.h0 / 2.0, -boxHvoigt.h1 / 2.0,
          boxHvoigt.h1 / 2.0, -boxHvoigt.h2 / 2.0, boxHvoigt.h2 / 2.0,
          boxHvoigt.h5, boxHvoigt.h4, boxHvoigt.h3);

      safeFprintf(fout, "\nAtoms\n\n");
      for (int ith = 0; ith < nAtom; ith++) {
        safeFprintf(fout, "%d %d %.16g 1 %.16g %.16g %.16g %d %d %d\n", ith + 1,
                    type[ith] + 1, diameterScale[ith] * meanDiameter,
                    xyz[ith].x, xyz[ith].y, xyz[ith].z, img[ith].x, img[ith].y,
                    img[ith].z);
      }

      if (veloc != NULL) {
        safeFprintf(fout, "\nVelocities\n\n");
        for (int ith = 0; ith < nAtom; ith++) {
          safeFprintf(fout, "%d %.16g %.16g %.16g\n", ith + 1, veloc[ith].x,
                      veloc[ith].y, veloc[ith].z);
        }
      }

      safeCloseFile(fout);
    }

  } else if (strstr(str, "LAMMPS"))  // Lammps style data file
  {
    while (fgets(str, 4096, fp) != NULL) {
      if (strstr(str, "atoms") != NULL) {
        nAtom = atoi(str);
        if (nAtom <= 0) Abort("No atom!");
        xyz = (double3 *)calloc(nAtom, sizeof(double3));

        img = (int3 *)calloc(nAtom, sizeof(int3));
        type = (int *)calloc(nAtom, sizeof(int));
        diameterScale = (double *)calloc(nAtom, sizeof(double));
      }
      if (strstr(str, "atom types") != NULL) {
        int type = atoi(str);
        if (type <= 0) Abort("Wrong DATA file!");
        nAtomType = type;
      }

      if (strstr(str, "xlo xhi") != NULL) {
        double xlo, xhi;
        sscanf(str, "%lf %lf", &xlo, &xhi);
        boxHvoigt.h0 = xhi - xlo;
      }
      if (strstr(str, "ylo yhi") != NULL) {
        double ylo, yhi;
        sscanf(str, "%lf %lf", &ylo, &yhi);
        boxHvoigt.h1 = yhi - ylo;
      }
      if (strstr(str, "zlo zhi") != NULL) {
        double zlo, zhi;
        sscanf(str, "%lf %lf", &zlo, &zhi);
        boxHvoigt.h2 = zhi - zlo;
      }
      if (strstr(str, "xy xz yz") != NULL) {
        double xy, xz, yz;
        sscanf(str, "%lf %lf %lf", &xy, &xz, &yz);
        boxHvoigt.h3 = yz;
        boxHvoigt.h4 = xz;
        boxHvoigt.h5 = xy;
      }

      if (strstr(str, "Atoms") != NULL) {
        meanDiameter = 0.0;
        for (int iatom = 0; iatom < nAtom;) {
          if (feof(fp) && iatom < nAtom) Abort("Wrong dataFile!");
          fgets(str, 4096, fp);
          if (isEmpty(str, 4096)) continue;

          int num, id, itype, ix = 0, iy = 0, iz = 0;
          double diam, density, x, y, z;

          num = sscanf(str, "%d %d %lf %lf %lf %lf %lf %d %d %d", &id, &itype,
                       &diam, &density, &x, &y, &z, &ix, &iy, &iz);
          if (num != 10) Abort("Wrong format!");
          if (id <= 0 || id > nAtom) Abort("Wrong atom ID.");
          if (itype <= 0 || itype > nAtomType) Abort("Wrong atom type.");

          xyz[id - 1] = make_double3(x, y, z);
          type[id - 1] = itype - 1;
          diameterScale[id - 1] = diam;
          img[id - 1] = make_int3(ix, iy, iz);
          meanDiameter += diam;

          iatom++;
        }
        meanDiameter = meanDiameter / nAtom;
        for (int iatom = 0; iatom < nAtom; iatom++) {
          diameterScale[iatom] /= meanDiameter;
        }
      }
      if (strstr(str, "Velocities") != NULL) {
        veloc = (double3 *)calloc(nAtom, sizeof(double3));
        for (int iatom = 0; iatom < nAtom;) {
          if (feof(fp) && iatom < nAtom) Abort("Wrong dataFile!");
          fgets(str, 4096, fp);
          if (isEmpty(str, 4096)) continue;

          int num, id;
          double vx, vy, vz;
          num = sscanf(str, "%d %lf %lf %lf", &id, &vx, &vy, &vz);
          veloc[id - 1] = make_double3(vx, vy, vz);
          iatom++;
        }
      }
    }

    {
      char str[4096];
      sprintf(str, "%s.bin", argv[1]);
      FILE *fout = createReadWriteFile(str);

      // header
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "binary");
      str[31] = '\n';
      fwrite(str, sizeof(char), 32, fout);  // 32 byte;
      // nAtom atoms
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "atoms");
      fwrite(str, sizeof(char), 32, fout);
      fwrite(&nAtom, sizeof(int), 1, fout);
      // nAtomType atom types
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "atom types");
      fwrite(str, sizeof(char), 32, fout);
      fwrite(&nAtomType, sizeof(int), 1, fout);
      // box Hvoigt
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "box Hvoigt");
      fwrite(str, sizeof(char), 32, fout);
      fwrite(&boxHvoigt, sizeof(double), 6, fout);
      // mean diameter
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "mean diameter");
      fwrite(str, sizeof(char), 32, fout);
      fwrite(&meanDiameter, sizeof(double), 1, fout);
      // Atoms Info
      memset(str, '\0', 4096 * sizeof(char));
      sprintf(str, "Atoms");
      fwrite(str, sizeof(char), 32, fout);
      for (int iatom = 0; iatom < nAtom; iatom++) {
        fwrite(&type[iatom], sizeof(int), 1, fout);
        fwrite(&diameterScale[iatom], sizeof(double), 1, fout);
        fwrite(&xyz[iatom].x, sizeof(double), 1, fout);
        fwrite(&xyz[iatom].y, sizeof(double), 1, fout);
        fwrite(&xyz[iatom].z, sizeof(double), 1, fout);
        fwrite(&img[iatom].x, sizeof(int), 1, fout);
        fwrite(&img[iatom].y, sizeof(int), 1, fout);
        fwrite(&img[iatom].z, sizeof(int), 1, fout);
      }
      if (veloc != NULL) {  // velocities Info
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "Velocities");
        fwrite(str, sizeof(char), 32, fout);
        for (int iatom = 0; iatom < nAtom; iatom++) {
          fwrite(&veloc[iatom].x, sizeof(double), 1, fout);
          fwrite(&veloc[iatom].y, sizeof(double), 1, fout);
          fwrite(&veloc[iatom].z, sizeof(double), 1, fout);
        }
      }

      safeCloseFile(fout);
    }
  }

  fclose(fp);
}
