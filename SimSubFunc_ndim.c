#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"

#ifdef __cplusplus
#error "This file must NOT be compiled by C++ compliler!"
#endif

void readCmdLineVar(Variable *var, int argc, char const *argv[]) {
    if (setSignalHandler() < 0)
        Abort("Failed when call setSignalHandler()!");
    
    if (screenOutputReadConf)
        safeFprintf(stderr, "Reading cmdline args ...\n");
    // find --trunc or --append tag
    bool nocite = true;
    for (int iarg = 1; iarg < argc; iarg++) {
        if (!strcmp(argv[iarg], "--trunc")) {
            if (truncFileFlag != 0)
                Abort("--trunc or --append");
            truncFileFlag = 1;
        } else if (!strcmp(argv[iarg], "--append")) {
            if (truncFileFlag != 0)
                Abort("--trunc or --append");
            truncFileFlag = 2;
        } else if (!strcmp(argv[iarg], "--cite")) {
            nocite = false;
            screenOutputReadConf += 10000;
        } else if (!strcmp(argv[iarg], "--dumpSourceFile")) {
#ifdef __Linux__
            dumpSourceFile();
#endif
            exit(EXIT_SUCCESS);
        } else if (!strcmp(argv[iarg], "--help")) {
            logFile = stdout;
            safeFprintf(logFile, "Info at compiling time:\n");
            safeFprintf(logFile, "\tCompile time: %s %s\n", __DATE__, __TIME__);
#ifdef __CompilePath__
            safeFprintf(logFile, "\tCompile path: %s\n", __CompilePath__);
#endif
#ifdef __SourceFileName__
            safeFprintf(logFile, "\tSource Files: %s\n", __SourceFileName__);
#endif
#if defined(__triBox__)
            safeFprintf(logFile, "\tBox: triclinic\n");
#elif defined(__orthBox__)
            safeFprintf(logFile, "\tBox: orthogonal\n");
#endif
            safeFprintf(logFile, "\tDIM: %d\n", DIM);
            safeFprintf(logFile, "\tdumpFileRevNum: %d\n", dumpFileRevNum);
            
            char *ptr = getcwd(NULL, 0);
            safeFprintf(logFile, "\tRunning path: %s\n", ptr);
            safeFree(ptr);
            
            exit(EXIT_SUCCESS);
        }
    }
    
    // parse cmd line
    var->maxVar = 8;
    var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
    char *msg = NULL;
    // variable start with "--" and the next character is not "-".
    for (int iarg = 1; iarg < argc;) {
        if (!strcmp(argv[iarg], "--log")) {
            iarg += 1;
            if (iarg < argc && strncmp(argv[iarg], "--", 2)) {
                msg = (char *)argv[iarg];
                iarg += 1;
            }
        } else if (!strcmp(argv[iarg], "--trunc") ||
                   !strcmp(argv[iarg], "--append")) {
            iarg++;
        } else if (!strcmp(argv[iarg], "--cite")) {
            iarg++;
        } else if (!strcmp(argv[iarg], "--cwd")) {
            if (iarg + 1 >= argc)
                Abort("--cwd .");
            int len = (int)strlen(argv[iarg + 1]) + 2;
            var->cwd = (char *)calloc(len, sizeof(char));
            sprintf(var->cwd, "%s", argv[iarg + 1]);
            iarg += 2;
        } else if (!strcmp(argv[iarg], "--sf")) {
            if (iarg + 1 >= argc)
                Abort("--sf id");
            int len = (int)strlen(argv[iarg + 1]) + 2;
            var->sf = (char *)calloc(len, sizeof(char));
            sprintf(var->sf, "%s", argv[iarg + 1]);
            iarg += 2;
        } else if (!strncmp(argv[iarg], "--", 2)) {
            if (strlen(argv[iarg]) == 2) {
                Abort("--varName varArg ...");
            } else if (argv[iarg][2] == '-') {
                Abort("--varName varArg ...");
            }
            
            if (findVariable(var, (char *)argv[iarg]))
                Abort("Repetitive cmd: %s", argv[iarg]);
            if (var->nVar == var->maxVar) {
                var->maxVar += 8;
                var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
            }
            
            int nsize = (int)strlen(argv[iarg]) + 1;
            var->cmd[var->nVar].cmdArgc = 0;
            var->cmd[var->nVar].cmdArgv = NULL;
            int cmdArgvStart = iarg + 1;
            
            iarg++;
            while (iarg < argc) {
                bool isVar = true;
                if (strncmp(argv[iarg], "--", 2))
                    isVar = false;
                if (isVar) {
                    if (strlen(argv[iarg]) <= 2)
                        isVar = false;
                    else if (argv[iarg][2] == '-')
                        isVar = false;
                }
                if (isVar)
                    break;
                
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
            memcpy(ptr, argv[cmdArgvStart - 1], (strlen(argv[cmdArgvStart - 1]) + 1) * sizeof(char));
            ptr += strlen(argv[cmdArgvStart - 1]) + 1;
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
    
    // logfile
    if (!nocite) {
        char timeString[32];
        getTimeString(timeString);
        
        char fname[40960] = {0};
        snprintf(fname, 40960, "%s/logFile_%s_%d_%s.dat", var->cwd, timeString, getpid(), var->sf);
        
        logFile = createFileReadWrite(fname);
        safeFprintf(logFile, "Info at compiling time:\n");
        safeFprintf(logFile, "\tCompile time: %s %s\n", __DATE__, __TIME__);
#ifdef __CompilePath__
        safeFprintf(logFile, "\tCompile path: %s\n", __CompilePath__);
#endif
#ifdef __SourceFileName__
        safeFprintf(logFile, "\tSource Files: %s\n", __SourceFileName__);
#endif
#if defined(__triBox__)
        safeFprintf(logFile, "\tBox: triclinic\n");
#elif defined(__orthBox__)
        safeFprintf(logFile, "\tBox: orthogonal\n");
#endif
        safeFprintf(logFile, "\tDIM: %d\n", DIM);
        safeFprintf(logFile, "\tdumpFileRevNum: %d\n", dumpFileRevNum);
        
        safeFprintf(logFile, "Info at running time:\n");
        safeFprintf(logFile, "\tProg. Name: %s\n", argv[0]);
        safeFprintf(logFile, "\tRunning time: %s\n", timeString);
        safeFprintf(logFile, "\tCmdLine arguments: ");
        for (int iarg = 1; iarg < argc; iarg++) {
            safeFprintf(logFile, "%s ", argv[iarg]);
        }
        safeFprintf(logFile, "\n");
        
        char *ptr = getcwd(NULL, 0);
        safeFprintf(logFile, "\tRunning path: %s\n", ptr);
        safeFree(ptr);
        
        if (MSG)
            safeFprintf(logFile, "Readme of Program:\n\t%s\n", MSG);
        if (msg)
            safeFprintf(logFile, "Extra info from cmdline: \n\t%s\n", msg);
    }
    
    // screen output
    if (screenOutputReadConf) {
        safeFprintf(stderr, "Prog. Name: \n\t%s\n", argv[0]);
        safeFprintf(stderr, "CmdLine arguments: \n\t");
        for (int iarg = 1; iarg < argc; iarg++) {
            safeFprintf(stderr, "%s ", argv[iarg]);
        }
        safeFprintf(stderr, "\n");
        
        char *ptr = getcwd(NULL, 0);
        safeFprintf(stderr, "Running path: \n\t%s\n", ptr);
        safeFree(ptr);
        if (MSG)
            safeFprintf(stderr, "Readme of Program:\n\t%s\n", MSG);
        if (msg)
            safeFprintf(stderr, "Extra info from cmdline: \n\t%s\n", msg);
    }
}

//===============================================================================
void initConfInfo(Box *box, Particle *particle, Update *update) {
    double vol = 0;
    double minScale = particle->diameterScale[0],
    maxScale = particle->diameterScale[0];
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vol += VolUnitSphere * pow(particle->diameterScale[iatom], DIM);
        double tmp = particle->diameterScale[iatom];
        minScale = (minScale < tmp ? minScale : tmp);
        maxScale = (maxScale > tmp ? maxScale : tmp);
    }
    update->nebrList.minDiameterScale = minScale;
    update->nebrList.maxDiameterScale = maxScale;
    
    update->volFrac = vol * pow(particle->meanDiameter * 0.5, DIM) / box->volume;
    update->dof = particle->nAtom * DIM;
}

void readConf_dump(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--rd");
    if (!cmd || cmd->cmdArgc != 2)
        Abort("--rd dump.bin whichStep");
    int whichStep = (int)atoi(cmd->cmdArgv[1]);
    
    mmapBinFile *binFile = openBinFile(cmd->cmdArgv[0]);
    if (whichStep < 0)
        whichStep = binFile->nStep - 1;
    if (whichStep >= binFile->nStep) {
        Abort("Step %d is out of Range: [-1,%d];", whichStep, binFile->nStep - 1);
    }
    
    readDump(box, particle, update, binFile, whichStep);
    
    closeBinFile(&binFile);
}

void readConf_data(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--rf");
    if (!cmd || cmd->cmdArgc == 0)
        Abort("read data: --rf lmp.bin");
    FILE *fp = openExistFileReadOnly(cmd->cmdArgv[0]);
    if (particle->pos != NULL)
        Abort("--rf or --rd");
    
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
                if (particle->nAtom <= 0)
                    Abort("Wrong File!");
                
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
                if (particle->nAtomType <= 0)
                    Abort("Wrong File!");
                
                particle->massPerType =
                (double *)calloc(particle->nAtomType, sizeof(double));
                for (int itype = 0; itype < particle->nAtomType; itype++)
                    particle->massPerType[itype] = 1.0;
            } else if (strcmp(str, "box Hvoigt") == 0) {
                if (!hasDimension) {
                    box->dim = 3;
                    if (box->dim != DIM)
                        Abort("Dimension is not consistent!");
                    
                    fread(&box->boxH[spaceIdx2voigt(0, 0)], sizeof(double), 1, fp);  // xx
                    fread(&box->boxH[spaceIdx2voigt(1, 1)], sizeof(double), 1, fp);  // yy
                    fread(&box->boxH[spaceIdx2voigt(2, 2)], sizeof(double), 1, fp);  // zz
                    fread(&box->boxH[spaceIdx2voigt(1, 2)], sizeof(double), 1, fp);  // yz
                    fread(&box->boxH[spaceIdx2voigt(0, 2)], sizeof(double), 1, fp);  // xz
                    fread(&box->boxH[spaceIdx2voigt(0, 1)], sizeof(double), 1, fp);  // xy
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
            } else {
                int sizeByte = 0;
                fread(&sizeByte, sizeof(int), 1, fp);
                void *ptr = calloc(sizeByte, 1);
                fread(ptr, 1, sizeByte, fp);
                if (findToolkit(&update->toolkit, str) >= 0)
                    Abort("Wrong binary file!");
                addToolkit(&update->toolkit, ptr, NULL, str);
            }
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
    } else if (strstr(str, "LAMMPS compatible data file.")) {
        // Lammps style data file
        box->dim = DIM;
        while (fgets(str, 4096, fp) != NULL) {
            if (strstr(str, "atoms") != NULL) {
                particle->nAtom = atoi(str);
                if (particle->nAtom <= 0)
                    Abort("No atom!");
                
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
            }
            if (strstr(str, "atom types") != NULL) {
                particle->nAtomType = atoi(str);
                if (particle->nAtomType <= 0)
                    Abort("Wrong DATA file!");
                
                particle->massPerType =
                (double *)calloc(particle->nAtomType, sizeof(double));
                for (int itype = 0; itype < particle->nAtomType; itype++)
                    particle->massPerType[itype] = 1.0;
            }
            
            if (strstr(str, "xlo xhi") != NULL) {
                double xlo, xhi;
                sscanf(str, "%lf %lf", &xlo, &xhi);
                box->boxH[spaceIdx2voigt(0, 0)] = xhi - xlo;
            }
            if (strstr(str, "ylo yhi") != NULL) {
                double ylo, yhi;
                sscanf(str, "%lf %lf", &ylo, &yhi);
                box->boxH[spaceIdx2voigt(1, 1)] = yhi - ylo;
            }
#if (DIM == 3)
            if (strstr(str, "zlo zhi") != NULL) {
                double zlo, zhi;
                sscanf(str, "%lf %lf", &zlo, &zhi);
                box->boxH[spaceIdx2voigt(2, 2)] = zhi - zlo;
            }
#endif
            if (strstr(str, "xy xz yz") != NULL) {
                double xy, xz, yz;
                sscanf(str, "%lf %lf %lf", &xy, &xz, &yz);
#if (DIM == 3)
                box->boxH[spaceIdx2voigt(1, 2)] = yz;
                box->boxH[spaceIdx2voigt(0, 2)] = xz;
                box->boxH[spaceIdx2voigt(0, 1)] = xy;
#elif (DIM == 2)
                box->boxH[spaceIdx2voigt(0, 1)] = xy;
#else
                Abort("Only 2D and 3D are vailid!");
#endif
            }
            
            if (strstr(str, "Atoms") != NULL) {
                double meanDiameter = 0.0;
                for (int iatom = 0; iatom < particle->nAtom;) {
                    if (feof(fp) && iatom < particle->nAtom)
                        Abort("Wrong dataFile!");
                    fgets(str, 4096, fp);
                    if (isEmpty(str, 4096))
                        continue;
                    
                    int num, id, type, ix = 0, iy = 0, iz = 0;
                    double diam, density, x, y, z;
                    
                    num = sscanf(str, "%d %d %lf %lf %lf %lf %lf %d %d %d", &id, &type,
                                 &diam, &density, &x, &y, &z, &ix, &iy, &iz);
                    if (num != 10)
                        Abort("Wrong format!");
                    if (id <= 0 || id > particle->nAtom)
                        Abort("Wrong atom ID.");
                    if (type <= 0 || type > particle->nAtomType)
                        Abort("Wrong atom type.");
                    
#if (DIM == 3)
                    particle->pos[id - 1][0] = x;
                    particle->pos[id - 1][1] = y;
                    particle->pos[id - 1][2] = z;
#elif (DIM == 2)
                    particle->pos[id - 1][0] = x;
                    particle->pos[id - 1][1] = y;
#endif
                    
                    particle->type[id - 1] = type - 1;
                    particle->mass[id - 1] =
                    particle->massPerType[particle->type[id - 1]];
                    particle->diameterScale[id - 1] = diam;
#if (DIM == 3)
                    particle->img[id - 1][0] = ix;
                    particle->img[id - 1][1] = iy;
                    particle->img[id - 1][2] = iz;
#elif (DIM == 2)
                    particle->img[id - 1][0] = ix;
                    particle->img[id - 1][1] = iy;
#endif
                    
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
                    if (feof(fp) && iatom < particle->nAtom)
                        Abort("Wrong dataFile!");
                    fgets(str, 4096, fp);
                    if (isEmpty(str, 4096))
                        continue;
                    
                    int num, id;
                    double vx, vy, vz;
                    num = sscanf(str, "%d %lf %lf %lf", &id, &vx, &vy, &vz);
                    
#if (DIM == 3)
                    particle->veloc[id - 1][0] = vx;
                    particle->veloc[id - 1][1] = vy;
                    particle->veloc[id - 1][2] = vz;
#elif (DIM == 2)
                    particle->veloc[id - 1][0] = vx;
                    particle->veloc[id - 1][1] = vy;
#endif
                    
                    iatom++;
                }
            }
        }
    } else
        Abort("Wrong File!");
    safeCloseFile(fp);
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    particle->isForceValid = false;
    update->Edone = update->Pdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.nForce = 0;
    update->nebrList.cntForce = 0;
    update->nebrList.doSort = true;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
    adjustImg(box, particle);
}

void readConf(Box *box, Particle *particle, Update *update, Variable *var) {
    update->nebrList.skinSet = -1.0;
    update->nebrList.nDelay = -1;
    
    int cntR = 0;
    cmdArg *cmd = findVariable(var, "--rf");
    if (cmd) {
        readConf_data(box, particle, update, var);
        cntR++;
    }
    cmd = findVariable(var, "--rd");
    if (cmd && cmd->cmdArgc == 2) {
        if (cntR == 1)
            Abort("--rf or --rd");
        readConf_dump(box, particle, update, var);
        cntR++;
    }
    if (cntR != 1)
        Abort("--rd dump.bin step or --rf conf.bin");  //--rd dump.bin or no "--rd"
    // and "--rf"
    
    cmd = findVariable(var, "--skin");
    if (cmd) {
        if (cmd->cmdArgc <= 0)
            Abort("--skin 0.1");
        update->nebrList.skinSet = atof(cmd->cmdArgv[0]) / 1.1;
    } else {
        update->nebrList.skinSet = 0.1 / 1.1;
    }
    
#ifdef __orthBox__
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            if (fabs(box->boxH[spaceIdx2voigt(idim, jdim)]) >= 5E-16)
                Abort("Not orthogonal Box!");
            box->boxH[spaceIdx2voigt(idim, jdim)] = 0;
        }
    }
#endif
    
    if (screenOutputReadConf) {
        printf("===========System Info==========\n");
        printf("Dimension: %d\n", DIM);
        printf("Number of Particles: %d;\n", particle->nAtom);
        printf("Volume Fraction: %g;\n", update->volFrac);
        printf("Min(diameter): %g;\nMax(diameter): %g;\n",
               update->nebrList.minDiameterScale,
               update->nebrList.maxDiameterScale);
        printf("Edges of Simulation Box: \n");
        for (int iedge = 0; iedge < DIM; iedge++) {
            printf("\t");
            for (int jdim = 0; jdim < DIM; jdim++) {
                printf("%-8.6e\t", box->boxEdge[iedge][jdim] / update->distanceUnits);
            }
            printf("\n");
        }
#if defined(Harmonic)
        printf("Interaction between particles: Harmonic\n");
#elif defined(Hertzian)
        printf("Interaction between particles: Hertzian\n");
#else
#error "No Interaction is defined!"
#endif
        printf("===========System Info==========\n");
        
        if (screenOutputReadConf >= 10000)
            screenOutputReadConf -= 10000;
    }
}
void writeConf(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--wf");
    if (!cmd)
        return;
    if (cmd->cmdArgc == 0)
        Abort("--wf output.bin");
    
    FILE *fbin = NULL;
    fbin = createFileReadWrite(cmd->cmdArgv[0]);
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
        if (update->isThermalRun) {
            // Atoms Info
            memset(str, '\0', 4096 * sizeof(char));
            sprintf(str, "Velocities");
            fwrite(str, sizeof(char), 32, fbin);
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                int idx = particle->tag2id[iatom];
                
                doubleVecPtr velocPtr = particle->veloc[idx];
                fwrite(velocPtr, sizeof(doubleVector), 1, fbin);
            }
        }
        // toolkit
        for (int itool = 0; itool < update->toolkit.nToolkit; itool++) {
            FuncPtrToolkitWriteConf ptr = (FuncPtrToolkitWriteConf)update->toolkit.funcPtrWriteConf[itool];
            if (ptr) {
                ptr(fbin, box, particle, update);
            }
        }
    }
    safeCloseFile(fbin);
}
int emergWriteConf(Box *box, Particle *particle, Update *update, char *fname) {
    if (access(fname, F_OK) == 0) {
        char tsring[32];
        getTimeString(tsring);
        Info("File \"%s\" will be truncated! Operation Time: %s", fname, tsring);
    }
    FILE *fbin = fopen(fname, "wb+");
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
        // toolkit
        for (int itool = 0; itool < update->toolkit.nToolkit; itool++) {
            FuncPtrToolkitWriteConf ptr = (FuncPtrToolkitWriteConf)update->toolkit.funcPtrWriteConf[itool];
            if (ptr) {
                ptr(fbin, box, particle, update);
            }
        }
        fflush(fbin);
    }
    safeCloseFile(fbin);
    return 0;
}

//===============================================================================
void zeroMomentum(Box *box, Particle *particle, Update *update) {
    doubleVector sV;
    vZeros(sV);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vAdd(sV, sV, particle->veloc[iatom]);
    }
    vScale(sV, 1.0 / particle->nAtom, sV);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vSub(particle->veloc[iatom], particle->veloc[iatom], sV);
    }
    
    update->Edone = update->Pdone = update->Tdone = update->Wdone = false;
    particle->isForceValid = false;
}
void moveBarycenter(Box *box, Particle *particle, Update *update) {
    doubleVector sCent;
    vZeros(sCent);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector imgPos;
        unwrapPos(imgPos, particle->pos[iatom], particle->img[iatom], box->boxH);
        vAdd(sCent, imgPos, sCent);
    }
    vScale(sCent, 1.0 / particle->nAtom, sCent);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vSub(particle->pos[iatom], particle->pos[iatom], sCent);
    }
    
    update->Edone = update->Pdone = update->Tdone = update->Wdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
}

void calcDistBoxPlane(doubleVector distPlane, doubleVector boxEdge[DIM]) {
    doubleVector edge[DIM];
    
    for (int idim = 0; idim < DIM; idim++) {
        vZeros(edge[idim]);
        edge[idim][idim] = 1.0;
    }
    for (int idim = DIM - 1; idim >= 0; idim--) {
        for (int ith = idim; ith < DIM - 1; ith++) {
            vCpy(edge[ith], boxEdge[ith + 1]);
            memset(edge[ith], '\0', idim * sizeof(double));
        }
        vCpy(edge[DIM - 1], boxEdge[idim]);
        memset(edge[DIM - 1], '\0', idim * sizeof(double));
        
        for (int ith = idim; ith < DIM; ith++) {
            vUnit(edge[ith], edge[ith]);
            for (int jth = ith + 1; jth < DIM; jth++) {
                double dotProd = sDot(edge[ith], edge[jth]);
                vScaleAdd(edge[jth], edge[jth], -dotProd, edge[ith]);
            }
        }
        
        distPlane[idim] = sDot(boxEdge[idim], edge[DIM - 1]);
        distPlane[idim] = fabs(distPlane[idim]);
    }
}

void calcTemp(Particle *particle, Update *update) {
    if (!update->isThermalRun) {
        update->tempKin = 0;
        uptriMatZeros(update->kinTens);
        update->Tdone = true;
        return;
    }
    
    double totKin = 0.0;
    uptriMat kinTens;
    uptriMatZeros(kinTens);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        totKin += sNormP2(particle->veloc[iatom]);
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim; jdim < DIM; jdim++) {
                kinTens[spaceIdx2voigt(idim, jdim)] += particle->veloc[iatom][idim] * particle->veloc[iatom][jdim];
            }
        }
    }
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            kinTens[spaceIdx2voigt(idim, jdim)] /= update->dof;
        }
    }
    uptriMatCpy(update->kinTens, kinTens);
    
    update->tempKin = totKin / update->dof;
    update->Tdone = true;
}
void calcPressure(Box *box, Particle *particle, Update *update) {
    if (!update->Wdone) {
        Info("virial is not computed!");
        update->pVirKin = 0;
        return;
    }
    if (update->Pdone)
        return;
    
    if (!update->Tdone) {
        calcTemp(particle, update);
    }
    
    update->pVirKin = 0.0;
    update->pVir = 0.0;
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            update->pVirKinTens[spaceIdx2voigt(idim, jdim)] = update->pVirTens[spaceIdx2voigt(idim, jdim)] + particle->nAtom / box->volume * update->kinTens[spaceIdx2voigt(idim, jdim)];
        }
        update->pVirKin += update->pVirKinTens[spaceIdx2voigt(idim, idim)];
        update->pVir += update->pVirTens[spaceIdx2voigt(idim, idim)];
    }
    update->pVirKin /= DIM;
    update->pVir /= DIM;
    
    update->Pdone = true;
}
void scaleVeloc(Particle *particle, Update *update, double tT) {
    if (!update->Tdone) {
        double totKin = 0.0;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            totKin += sNormP2(particle->veloc[iatom]);
        }
        update->tempKin = 0.5 * totKin / update->dof;
        update->Tdone = true;
    }
    
    double Tcurr = update->tempKin;
    if (Tcurr < 1E-12) {
        genVeloc(particle, update, tT);
    } else {
        double sfact = sqrt(tT / Tcurr);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            vScale(particle->veloc[iatom], sfact, particle->veloc[iatom]);
        }
        update->tempKin = Tcurr * sfact * sfact;
        update->Tdone = true;
        update->Pdone = false;
    }
    
    return;
}
void genVeloc(Particle *particle, Update *update, double tT) {
    doubleVector sVeloc;
    vZeros(sVeloc);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        for (int idim = 0; idim < DIM; idim++) {
            particle->veloc[iatom][idim] = rndStdNorm();
        }
        vAdd(sVeloc, sVeloc, particle->veloc[iatom]);
    }
    vScale(sVeloc, 1.0 / particle->nAtom, sVeloc);
    double totKin = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vSub(particle->veloc[iatom], particle->veloc[iatom], sVeloc);
        totKin += sNormP2(particle->veloc[iatom]);
    }
    double Tcurr = totKin / update->dof;
    double sfact = sqrt(tT / Tcurr);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScale(particle->veloc[iatom], sfact, particle->veloc[iatom]);
    }
    update->tempKin = Tcurr * sfact * sfact;
    update->Tdone = true;
    update->Pdone = false;
}

#define RcsBinLenRatio 1
void calcBinParam(Box *box, Particle *particle, Update *update) {
    NebrList *nebrList = &update->nebrList;
    if (!nebrList->compelInit && box->isShapeFixed && particle->isSizeFixed &&
        nebrList->binList)
        return;
    // nvt
    
    nebrList->maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
    // checking Rskin
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
    for (int idim = 1; idim < DIM; idim++) {
        minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
    }
    double maxRskin = (0.5 * minAxByCz - nebrList->maxRcut) * 0.99;
    double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
    if (maxRskinSet < __minSkinSet__) {
        double minRskin = __minSkinSet__ * nebrList->minDiameterScale * particle->meanDiameter;
        double minLen = (nebrList->maxRcut + minRskin) * 2.0;
        double minBoxVol = pow(minLen, DIM);
        double maxPhi = update->volFrac * box->volume / minBoxVol;
        
        Abort(
              "The box is too small, the PBC method breaks down! Please use large N, "
              "or rewrite codes! The 0.5*min(box(.,.)) is %g, max(Rcut) is %g, "
              "min(Rskin) is %g. The max(Rcut)+min(Rskin) should be smaller "
              "than 0.5*min(ax,by,cz). The estimated max(volFrac) is %g, currently, "
              "it is %g.",
              0.5 * minAxByCz / update->distanceUnits,
              nebrList->maxRcut / update->distanceUnits,
              minRskin / update->distanceUnits, maxPhi, update->volFrac);
    }
    
    // calculate Rskin
    nebrList->rskin = nebrList->minDiameterScale * particle->meanDiameter * nebrList->skinSet;
    nebrList->rskin = cpuMin(nebrList->rskin, maxRskin);
    double sysRcs = (nebrList->maxRcut + nebrList->rskin) / RcsBinLenRatio;
    
#if (DIM == 3 || DIM == 2)
    //==========================================
    nebrList->totBin = 1;
    nebrList->totBinExt = 1;
    intVector bstart, bstop, ndb;
    intVector nExtBin;
    int nAdjBin = 1;
    for (int idim = 0; idim < DIM; idim++) {
        int inbin = (int)floor(distPlane[idim] / sysRcs);
        // maxRskin require that inbin >= 2;
        inbin = cpuMax(inbin, 1);
        double ibinLen = distPlane[idim] / inbin;
        int ixyzNum = (int)ceil((nebrList->maxRcut + nebrList->rskin) / ibinLen);
        int ibstart = -ixyzNum;
        int ibstop = ixyzNum;
        int indb = 2 * ixyzNum + 1;
        if (indb > inbin) {
            inbin = 2 * (inbin / 2) + 1;
            // in this situation inbin >= 3
            ibinLen = distPlane[idim] / inbin;
            ixyzNum = (inbin / 2);
            ibstart = -ixyzNum;
            ibstop = ixyzNum;
            indb = 2 * ixyzNum + 1;
        }
        
        nExtBin[idim] = ixyzNum;
        nebrList->nbin[idim] = inbin;
        nebrList->nbinExt[idim] = inbin + 2 * ixyzNum;
        nebrList->binLen[idim] = ibinLen;
        nebrList->totBin *= nebrList->nbin[idim];
        nebrList->totBinExt *= nebrList->nbinExt[idim];
        nebrList->nStencil[idim] = ixyzNum;
        
        bstart[idim] = ibstart;
        bstop[idim] = ibstop;
        ndb[idim] = indb;
        nAdjBin *= ndb[idim];
    }
    if (nebrList->totBinExt > nebrList->allocBinExt) {
        nebrList->allocBinExt = nebrList->totBinExt;
        nebrList->isSourceBin = (bool *)realloc(
                                                nebrList->isSourceBin, nebrList->allocBinExt * sizeof(bool));
        nebrList->stencilSource = (int2 *)realloc(
                                                  nebrList->stencilSource,
                                                  (nebrList->totBinExt - nebrList->totBin) * sizeof(int2));
    }
    if (nAdjBin / 2 > nebrList->allocAdjBin) {
        nebrList->allocAdjBin = nAdjBin / 2;
        nebrList->adjBinList = (int *)realloc(nebrList->adjBinList,
                                              nebrList->allocAdjBin * sizeof(int));
    }
    nebrList->nAdjBin = 0;
    for (int idx = 0; idx < nAdjBin; idx++) {
        intVector delta;
        int In = idx;
        for (int idim = 0; idim < DIM; idim++) {
            delta[idim] = (In % ndb[idim]) + bstart[idim];
            In = In / ndb[idim];
        }
        In = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            In = In * ndb[idim] + delta[idim];
        }
        if (In <= 0)
            continue;
        
        In = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            In = In * nebrList->nbinExt[idim] + delta[idim];
        }
        nebrList->adjBinList[nebrList->nAdjBin++] = In;
    }
    
    // check very bin
    
#if (DIM == 3)
    for (int bz = 0; bz < nebrList->nbinExt[2]; bz++) {
        for (int by = 0; by < nebrList->nbinExt[1]; by++) {
            for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
                int ibin = (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                bool isSource = false;
                if (bx >= nebrList->nStencil[0] &&
                    bx < nebrList->nStencil[0] + nebrList->nbin[0] &&
                    by >= nebrList->nStencil[1] &&
                    by < nebrList->nStencil[1] + nebrList->nbin[1] &&
                    bz >= nebrList->nStencil[2] &&
                    bz < nebrList->nStencil[2] + nebrList->nbin[2])
                    isSource = true;
                
                nebrList->isSourceBin[ibin] = isSource;
            }
        }
    }
    
    // stencil-xlo
    int cnt = 0;
    for (int idx = 0; idx < nebrList->nStencil[0]; idx++) {
        int bx = idx + nebrList->nbin[0];
        for (int by = nebrList->nStencil[1];
             by < nebrList->nbin[1] + nebrList->nStencil[1]; by++) {
            for (int bz = nebrList->nStencil[2];
                 bz < nebrList->nbin[2] + nebrList->nStencil[2]; bz++) {
                nebrList->stencilSource[cnt].first =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + idx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[0] = cnt;
    // stencil-xhi
    for (int idx = nebrList->nbinExt[0] - nebrList->nStencil[0];
         idx < nebrList->nbinExt[0]; idx++) {
        int bx = idx - nebrList->nbin[0];
        for (int by = nebrList->nStencil[1];
             by < nebrList->nbin[1] + nebrList->nStencil[1]; by++) {
            for (int bz = nebrList->nStencil[2];
                 bz < nebrList->nbin[2] + nebrList->nStencil[2]; bz++) {
                nebrList->stencilSource[cnt].first =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + idx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[1] = cnt;
    // stencil-ylo
    for (int idy = 0; idy < nebrList->nStencil[1]; idy++) {
        int by = idy + nebrList->nbin[1];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            for (int bz = nebrList->nStencil[2];
                 bz < nebrList->nbin[2] + nebrList->nStencil[2]; bz++) {
                nebrList->stencilSource[cnt].first =
                (bz * nebrList->nbinExt[1] + idy) * nebrList->nbinExt[0] + bx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[2] = cnt;
    // stencil-yhi
    for (int idy = nebrList->nbinExt[1] - nebrList->nStencil[1];
         idy < nebrList->nbinExt[1]; idy++) {
        int by = idy - nebrList->nbin[1];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            for (int bz = nebrList->nStencil[2];
                 bz < nebrList->nbin[2] + nebrList->nStencil[2]; bz++) {
                nebrList->stencilSource[cnt].first =
                (bz * nebrList->nbinExt[1] + idy) * nebrList->nbinExt[0] + bx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[3] = cnt;
    
    // stencil-zlo
    for (int idz = 0; idz < nebrList->nStencil[2]; idz++) {
        int bz = idz + nebrList->nbin[2];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            for (int by = 0; by < nebrList->nbinExt[1]; by++) {
                nebrList->stencilSource[cnt].first =
                (idz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[4] = cnt;
    // stencil-zhi
    for (int idz = nebrList->nbinExt[2] - nebrList->nStencil[2];
         idz < nebrList->nbinExt[2]; idz++) {
        int bz = idz - nebrList->nbin[2];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            for (int by = 0; by < nebrList->nbinExt[1]; by++) {
                nebrList->stencilSource[cnt].first =
                (idz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                nebrList->stencilSource[cnt].second =
                (bz * nebrList->nbinExt[1] + by) * nebrList->nbinExt[0] + bx;
                cnt++;
            }
        }
    }
    nebrList->exLayerNum[5] = cnt;
    
#else
    for (int by = 0; by < nebrList->nbinExt[1]; by++) {
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            int ibin = by * nebrList->nbinExt[0] + bx;
            bool isSource = false;
            if (bx >= nebrList->nStencil[0] &&
                bx < nebrList->nStencil[0] + nebrList->nbin[0] &&
                by >= nebrList->nStencil[1] &&
                by < nebrList->nStencil[1] + nebrList->nbin[1])
                isSource = true;
            
            nebrList->isSourceBin[ibin] = isSource;
        }
    }
    
    // stencil-xlo
    int cnt = 0;
    for (int idx = 0; idx < nebrList->nStencil[0]; idx++) {
        int bx = idx + nebrList->nbin[0];
        for (int by = nebrList->nStencil[1];
             by < nebrList->nbin[1] + nebrList->nStencil[1]; by++) {
            nebrList->stencilSource[cnt].first = (by)*nebrList->nbinExt[0] + idx;
            nebrList->stencilSource[cnt].second = (by)*nebrList->nbinExt[0] + bx;
            cnt++;
        }
    }
    nebrList->exLayerNum[0] = cnt;
    // stencil-xhi
    for (int idx = nebrList->nbinExt[0] - nebrList->nStencil[0];
         idx < nebrList->nbinExt[0]; idx++) {
        int bx = idx - nebrList->nbin[0];
        for (int by = nebrList->nStencil[1];
             by < nebrList->nbin[1] + nebrList->nStencil[1]; by++) {
            nebrList->stencilSource[cnt].first = (by)*nebrList->nbinExt[0] + idx;
            nebrList->stencilSource[cnt].second = (by)*nebrList->nbinExt[0] + bx;
            cnt++;
        }
    }
    nebrList->exLayerNum[1] = cnt;
    // stencil-ylo
    for (int idy = 0; idy < nebrList->nStencil[1]; idy++) {
        int by = idy + nebrList->nbin[1];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            nebrList->stencilSource[cnt].first = (idy)*nebrList->nbinExt[0] + bx;
            nebrList->stencilSource[cnt].second = (by)*nebrList->nbinExt[0] + bx;
            cnt++;
        }
    }
    nebrList->exLayerNum[2] = cnt;
    // stencil-yhi
    for (int idy = nebrList->nbinExt[1] - nebrList->nStencil[1];
         idy < nebrList->nbinExt[1]; idy++) {
        int by = idy - nebrList->nbin[1];
        for (int bx = 0; bx < nebrList->nbinExt[0]; bx++) {
            nebrList->stencilSource[cnt].first = (idy)*nebrList->nbinExt[0] + bx;
            nebrList->stencilSource[cnt].second = (by)*nebrList->nbinExt[0] + bx;
            cnt++;
        }
    }
    nebrList->exLayerNum[3] = cnt;
    
#endif
    
    if (nebrList->totBinExt > nebrList->maxBinHead) {
        nebrList->binHead = (idPosRadius *)realloc(
                                                   nebrList->binHead, nebrList->totBinExt * sizeof(idPosRadius));
        nebrList->maxBinHead = nebrList->totBinExt;
    }
    if (nebrList->binList == NULL) {
        nebrList->binList =
        (idPosRadius *)calloc(particle->nAtom, sizeof(idPosRadius));
    }
    
#if (DIM == 3)
    // calc norm vector of box-plane: DIM == 3
    for (int idim = 0; idim < DIM; idim++) {
        doubleVecPtr v1 = NULL, v2 = NULL;
        switch (idim) {
        case 0:
            v1 = box->boxEdge[1];
            v2 = box->boxEdge[2];
            break;
        case 1:
            v1 = box->boxEdge[2];
            v2 = box->boxEdge[0];
            break;
        case 2:
            v1 = box->boxEdge[0];
            v2 = box->boxEdge[1];
            break;
        }
        vCross(nebrList->normPlane[idim], v1, v2);
        vUnit(nebrList->normPlane[idim], nebrList->normPlane[idim]);
    }
#else
    // calc norm vector of box-line: DIM == 2
    for (int idim = 0; idim < DIM; idim++) {
        if (idim == 0) {
            nebrList->normPlane[idim][0] = box->boxEdge[1][1];
            nebrList->normPlane[idim][1] = -box->boxEdge[1][0];
        } else {
            nebrList->normPlane[idim][0] = box->boxEdge[0][1];
            nebrList->normPlane[idim][1] = -box->boxEdge[0][0];
        }
        vUnit(nebrList->normPlane[idim], nebrList->normPlane[idim]);
    }
#endif
#else
    nebrList->totBin = 1;
    nebrList->totBinExt = 1;
    intVector bstart, bstop, ndb;
    int nAdjBin = 1;
    for (int idim = 0; idim < DIM; idim++) {
        int inbin = (int)floor(distPlane[idim] / sysRcs);
        inbin = cpuMax(inbin, 1);
        double ibinLen = distPlane[idim] / inbin;
        int ixyzNum = (int)ceil((nebrList->maxRcut + nebrList->rskin) / ibinLen);
        int ibstart = -ixyzNum;
        int ibstop = ixyzNum;
        int indb = 2 * ixyzNum + 1;
        if (indb > inbin) {
            inbin = 2 * (inbin / 2) + 1;
            ibinLen = distPlane[idim] / inbin;
            ixyzNum = (inbin / 2);
            ibstart = -ixyzNum;
            ibstop = ixyzNum;
            indb = 2 * ixyzNum + 1;
        }
        
        nebrList->nbin[idim] = inbin;
        nebrList->nbinExt[idim] = inbin;
        nebrList->binLen[idim] = ibinLen;
        nebrList->totBin *= nebrList->nbin[idim];
        nebrList->totBinExt *= nebrList->nbinExt[idim];
        nebrList->nStencil[idim] = 0;
        
        bstart[idim] = ibstart;
        bstop[idim] = ibstop;
        ndb[idim] = indb;
        nAdjBin *= ndb[idim];
    }
    
    if (nebrList->totBinExt > nebrList->allocBinExt) {
        nebrList->allocBinExt = nebrList->totBinExt;
        nebrList->isSourceBin = NULL;
        nebrList->stencilSource = NULL;
    }
    if (nAdjBin / 2 > nebrList->allocAdjBin) {
        nebrList->allocAdjBin = nAdjBin / 2;
        nebrList->adjBinList = (intVector *)realloc(nebrList->adjBinList, nebrList->allocAdjBin * sizeof(intVector));
    }
    
    nebrList->nAdjBin = 0;
    for (int idx = 0; idx < nAdjBin; idx++) {
        intVector delta;
        int In = idx;
        for (int idim = 0; idim < DIM; idim++) {
            delta[idim] = (In % ndb[idim]) + bstart[idim];
            In = In / ndb[idim];
        }
        In = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            In = In * ndb[idim] + delta[idim];
        }
        if (In <= 0)
            continue;
        
        vCpy(nebrList->adjBinList[nebrList->nAdjBin], delta);
        nebrList->nAdjBin++;
    }
    
    if (nebrList->totBinExt > nebrList->maxBinHead) {
        nebrList->binHead = (idPosRadius *)realloc(
                                                   nebrList->binHead, nebrList->totBinExt * sizeof(idPosRadius));
        nebrList->maxBinHead = nebrList->totBinExt;
    }
    if (nebrList->binList == NULL) {
        nebrList->binList =
        (idPosRadius *)calloc(particle->nAtom, sizeof(idPosRadius));
    }
    
#endif
    
    nebrList->compelInit = false;
}

void binParticle(Box *box, Particle *particle, Update *update) {
    NebrList *nebrList = &update->nebrList;
    calcBinParam(box, particle, update);
    
    for (int ibin = 0; ibin < nebrList->totBinExt; ibin++) {
        nebrList->binHead[ibin].id = -1;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        idPosRadius data;
        data.id = iatom;
        vCpy(data.pos, particle->pos[iatom]);
        data.radius = 0.5 * particle->meanDiameter * particle->diameterScale[iatom];
        
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, data.pos);
        vShiftAll(lamda, 0.5);
        int binIdx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            // local box
            int cIdx = (int)floor(lamda[idim] * nebrList->nbin[idim]);
            cIdx = (cIdx < 0 ? nebrList->nbin[idim] - 1 : cIdx);
            cIdx = (cIdx >= nebrList->nbin[idim] ? 0 : cIdx);
            // extra stencil layer
            cIdx += nebrList->nStencil[idim];
            binIdx = binIdx * nebrList->nbinExt[idim] + cIdx;
        }
        
        nebrList->binList[data.id] = nebrList->binHead[binIdx];
        nebrList->binHead[binIdx] = data;
    }
    
    // build ghost/image atoms
#if (DIM == 3 || DIM == 2)
    // stencil-xlo
    int ith = 0;
    nebrList->nImage = 0;
    double mRcsP2 = pow(nebrList->maxRcut + nebrList->rskin, 2);
    for (; ith < nebrList->exLayerNum[0]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerHi);
            double pPlane = sDot(nebrList->normPlane[0], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] -= box->boxH[spaceIdx2voigt(0, 0)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[0] = nebrList->nImage;
    // stencil-xhi
    for (; ith < nebrList->exLayerNum[1]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerLo);
            double pPlane = sDot(nebrList->normPlane[0], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] += box->boxH[spaceIdx2voigt(0, 0)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[1] = nebrList->nImage;
    // stencil-ylo
    for (; ith < nebrList->exLayerNum[2]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerHi);
            double pPlane = sDot(nebrList->normPlane[1], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] -= box->boxH[spaceIdx2voigt(0, 1)];
            jdata.pos[1] -= box->boxH[spaceIdx2voigt(1, 1)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[2] = nebrList->nImage;
    // stencil-yhi
    for (; ith < nebrList->exLayerNum[3]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerLo);
            double pPlane = sDot(nebrList->normPlane[1], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] += box->boxH[spaceIdx2voigt(0, 1)];
            jdata.pos[1] += box->boxH[spaceIdx2voigt(1, 1)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[3] = nebrList->nImage;
    
#if (DIM == 3)
    // stencil-zlo
    for (; ith < nebrList->exLayerNum[4]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerHi);
            double pPlane = sDot(nebrList->normPlane[2], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] -= box->boxH[spaceIdx2voigt(0, 2)];
            jdata.pos[1] -= box->boxH[spaceIdx2voigt(1, 2)];
            jdata.pos[2] -= box->boxH[spaceIdx2voigt(2, 2)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[4] = nebrList->nImage;
    // stencil-zhi
    for (; ith < nebrList->exLayerNum[5]; ith++) {
        int tobin = nebrList->stencilSource[ith].first;
        int frombin = nebrList->stencilSource[ith].second;
        for (idPosRadius jdata = nebrList->binHead[frombin]; jdata.id >= 0;
             jdata = nebrList->binList[jdata.id]) {
            int jnext = jdata.id;  // store
            
            doubleVector boxCornP;
            vSub(boxCornP, jdata.pos, box->cornerLo);
            double pPlane = sDot(nebrList->normPlane[2], boxCornP);
            if (pPlane * pPlane > mRcsP2)
                continue;  // exclude those images
            
            jdata.id = particle->nAtom + nebrList->nImage;
            jdata.pos[0] += box->boxH[spaceIdx2voigt(0, 2)];
            jdata.pos[1] += box->boxH[spaceIdx2voigt(1, 2)];
            jdata.pos[2] += box->boxH[spaceIdx2voigt(2, 2)];
            if (nebrList->nImage >= nebrList->maxImage) {
                nebrList->maxImage += 1024;
                int nAllmax = particle->nAtom + nebrList->maxImage;
                particle->pos = (doubleVector *)realloc(particle->pos,
                                                        nAllmax * sizeof(doubleVector));
                nebrList->binList = (idPosRadius *)realloc(
                                                           nebrList->binList, nAllmax * sizeof(idPosRadius));
                nebrList->imageParent = (int *)realloc(
                                                       nebrList->imageParent, nebrList->maxImage * sizeof(int));
                nebrList->imgageSource = (int *)realloc(
                                                        nebrList->imgageSource, nebrList->maxImage * sizeof(int));
                if (!(particle->pos || nebrList->binList || nebrList->imageParent ||
                      nebrList->imgageSource))
                    Abort("Realloc Error!");
            }
            nebrList->imageParent[nebrList->nImage++] = jnext;
            
            nebrList->binList[jdata.id] = nebrList->binHead[tobin];
            nebrList->binHead[tobin] = jdata;
            
            jdata.id = jnext;  // restore
        }
    }
    nebrList->imagePartCount[5] = nebrList->nImage;
#endif
    
    // build source
    for (int iiamge = 0; iiamge < nebrList->nImage; iiamge++) {
        int jsource = nebrList->imageParent[iiamge];
        while (jsource >= particle->nAtom) {
            jsource = nebrList->imageParent[jsource - particle->nAtom];
        }
        nebrList->imgageSource[iiamge] = jsource;
    }
    
#endif
}
void constructList(Box *box, Particle *particle, Update *update) {
    NebrList *nebrList = &update->nebrList;
    
#if (DIM == 3 || DIM == 2)
    for (int gbin = 0, cntNebr = 0; gbin < nebrList->totBinExt; gbin++) {
        if (!nebrList->isSourceBin[gbin])
            continue;
        for (idPosRadius idata = nebrList->binHead[gbin]; idata.id >= 0;
             idata = nebrList->binList[idata.id]) {
            double iRadiusRskin = idata.radius + nebrList->rskin;
            nebrList->nNebr[idata.id].first = cntNebr;
            
            // loop over rest of atoms in i's bin
            for (idPosRadius jdata = nebrList->binList[idata.id]; jdata.id >= 0;
                 jdata = nebrList->binList[jdata.id]) {
                doubleVector dRij;
                vSub(dRij, idata.pos, jdata.pos);
                double rijP2 = sNormP2(dRij);
                double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                if (rijP2 > RcutRskinP2)
                    continue;
                
                if (cntNebr == nebrList->maxAllocNebr) {
                    nebrList->maxAllocNebr += 10240;
                    nebrList->list = (int *)realloc(nebrList->list,
                                                    nebrList->maxAllocNebr * sizeof(int));
                    if (nebrList->list == NULL)
                        Abort("realloc nebrList failed!");
                }
                nebrList->list[cntNebr++] = jdata.id;
            }
            
            // loop over all atoms in other bins in stencil
            for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                int jbin = gbin + nebrList->adjBinList[adj];
                for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                     jdata = nebrList->binList[jdata.id]) {
                    doubleVector dRij;
                    vSub(dRij, idata.pos, jdata.pos);
                    double rijP2 = sNormP2(dRij);
                    double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                    if (rijP2 > RcutRskinP2)
                        continue;
                    
                    if (cntNebr == nebrList->maxAllocNebr) {
                        nebrList->maxAllocNebr += 10240;
                        nebrList->list = (int *)realloc(
                                                        nebrList->list, nebrList->maxAllocNebr * sizeof(int));
                        if (nebrList->list == NULL)
                            Abort("realloc nebrList failed!");
                    }
                    nebrList->list[cntNebr++] = jdata.id;
                }
            }
            
            nebrList->nNebr[idata.id].second = cntNebr;
        }
    }
#else
    for (int gbin = 0, cntNebr = 0; gbin < nebrList->totBinExt; gbin++) {
        intVector gbinVec;
        for (int idim = 0, itmp = gbin; idim < DIM; idim++) {
            gbinVec[idim] = itmp % nebrList->nbinExt[idim];
            itmp = itmp / nebrList->nbinExt[idim];
        }
        
        for (idPosRadius idata = nebrList->binHead[gbin]; idata.id >= 0;
             idata = nebrList->binList[idata.id]) {
            double iRadiusRskin = idata.radius + nebrList->rskin;
            nebrList->nNebr[idata.id].first = cntNebr;
            
            // loop over rest of atoms in i's bin
            for (idPosRadius jdata = nebrList->binList[idata.id]; jdata.id >= 0;
                 jdata = nebrList->binList[jdata.id]) {
                doubleVector dRij;
                vSub(dRij, idata.pos, jdata.pos);
                // apply PBC for DIM >= 4;
                PBC(dRij, box);
                
                double rijP2 = sNormP2(dRij);
                double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                if (rijP2 > RcutRskinP2)
                    continue;
                
                if (cntNebr == nebrList->maxAllocNebr) {
                    nebrList->maxAllocNebr += 10240;
                    nebrList->list = (int *)realloc(nebrList->list,
                                                    nebrList->maxAllocNebr * sizeof(int));
                    if (nebrList->list == NULL)
                        Abort("realloc nebrList failed!");
                }
                nebrList->list[cntNebr++] = jdata.id;
            }
            
            // loop over all atoms in other bins in stencil
            for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                int jbin = 0;
                for (int idim = DIM - 1; idim >= 0; idim--) {
                    jbin = jbin * nebrList->nbinExt[idim] +
                    (gbinVec[idim] + nebrList->adjBinList[adj][idim] +
                     nebrList->nbinExt[idim]) %
                    nebrList->nbinExt[idim];
                }
                
                for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                     jdata = nebrList->binList[jdata.id]) {
                    doubleVector dRij;
                    vSub(dRij, idata.pos, jdata.pos);
                    // apply PBC for DIM >= 4;
                    PBC(dRij, box);
                    
                    double rijP2 = sNormP2(dRij);
                    double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                    if (rijP2 > RcutRskinP2)
                        continue;
                    
                    if (cntNebr == nebrList->maxAllocNebr) {
                        nebrList->maxAllocNebr += 10240;
                        nebrList->list = (int *)realloc(
                                                        nebrList->list, nebrList->maxAllocNebr * sizeof(int));
                        if (nebrList->list == NULL)
                            Abort("realloc nebrList failed!");
                    }
                    nebrList->list[cntNebr++] = jdata.id;
                }
            }
            
            nebrList->nNebr[idata.id].second = cntNebr;
        }
    }
#endif
}

bool isNebrListValid(Box *box, Particle *particle, Update *update) {
    NebrList *nebrList = &update->nebrList;
    if (nebrList->xyzHold == NULL) {
        nebrList->isValid = false;
        return false;
    }
    if (!nebrList->isValid) {
        nebrList->isValid = false;
        return false;
    }
    
    if (nebrList->cntForce < nebrList->nDelay)
        return true;
    
    if (particle->isSizeFixed && box->isShapeFixed) {
        double maxShiftP2 = 0.25 * nebrList->rskin * nebrList->rskin;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVector dRt0;
            vSub(dRt0, particle->pos[iatom], nebrList->xyzHold[iatom]);
            double rP2 = sNormP2(dRt0);
            if (rP2 >= maxShiftP2) {
                if (nebrList->cntForce == nebrList->nDelay)
                    nebrList->nDangerous++;
                nebrList->isValid = false;
                return false;
            }
        }
    } else if (!particle->isSizeFixed && box->isShapeFixed) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVector dRt0;
            vSub(dRt0, particle->pos[iatom], nebrList->xyzHold[iatom]);
            double ri0 = sNorm(dRt0);
            double dsig = (particle->meanDiameter - nebrList->meanDiameterHold) *
            particle->diameterScale[iatom];
            
            if (ri0 + dsig >= nebrList->rskin * 0.5) {
                if (nebrList->cntForce == nebrList->nDelay)
                    nebrList->nDangerous++;
                nebrList->isValid = false;
                return false;
            }
        }
    } else if (!box->isShapeFixed && particle->isSizeFixed) {
        uptriMat trans;  // trans =  Ht * inv(H0)
        MatMulMat(trans, box->boxH, nebrList->invBoxHold);
        diagMat volTrans;
        getDiag(volTrans, trans);
        
        double sfact = sMinElement(volTrans);
        
        invDiag(volTrans, volTrans);
        uptriMat shapeTrans;
        diagMulMat(shapeTrans, volTrans, trans);
        
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim + 1; jdim < DIM; jdim++) {
                double eps2 = pow(shapeTrans[spaceIdx2voigt(idim, jdim)], 2);
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
            if (nebrList->cntForce == nebrList->nDelay)
                nebrList->nDangerous++;
            nebrList->isValid = false;
            return false;
        }
        double maxShiftP2 = 0.25 * rs_eff * rs_eff;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVector posHold;
            vCpy(posHold, nebrList->xyzHold[iatom]);
            MatMulVec(posHold, trans, posHold);
            
            doubleVector dR;
            vSub(dR, particle->pos[iatom], posHold);
            double rP2 = sNormP2(dR);
            
            if (rP2 >= maxShiftP2) {
                if (nebrList->cntForce == nebrList->nDelay)
                    nebrList->nDangerous++;
                nebrList->isValid = false;
                return false;
            }
        }
    } else
        Abort("Not Code!");
    nebrList->isValid = true;
    return true;
}

void adjustImg(Box *box, Particle *particle) {
#ifdef __orthBox__
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr posPtr = particle->pos[iatom];
        intVecPtr imgPtr = particle->img[iatom];
        
        // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, posPtr);
        vShiftAll(lamda, 0.5);
        
        intVector shiftImg;
        vFloor(shiftImg, lamda);
        vAdd(imgPtr, imgPtr, shiftImg);
        
        doubleVector shiftPos;
        MatMulVec(shiftPos, box->boxH, shiftImg);
        vSub(posPtr, posPtr, shiftPos);
    }
#endif
#ifdef __triBox__
    bool isFlip = false;
    doubleVector maxTilt;
    for (int idim = 0; idim < DIM; idim++) {
        maxTilt[idim] = 0.505 * box->boxH[spaceIdx2voigt(idim, idim)];
    }
    doubleVector boxEdge[DIM];
    memcpy(boxEdge, box->boxEdge, DIM * sizeof(doubleVector));
    for (int idim = DIM - 1; idim >= 0; idim--) {
        for (int tilt = idim - 1; tilt >= 0; tilt--) {
            if (boxEdge[idim][tilt] >= maxTilt[tilt]) {
                vSub(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
                isFlip = true;
            } else if (boxEdge[idim][tilt] < -maxTilt[tilt]) {
                vAdd(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
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
                box->boxH[spaceIdx2voigt(idim, jdim)] = boxEdge[jdim][idim];
            }
        }
        setBoxPara(box);
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr posPtr = particle->pos[iatom];
        intVecPtr imgPtr = particle->img[iatom];
        
        // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, posPtr);
        vShiftAll(lamda, 0.5);
        
        intVector shiftImg;
        vFloor(shiftImg, lamda);
        vAdd(imgPtr, imgPtr, shiftImg);
        
        doubleVector shiftPos;
        MatMulVec(shiftPos, box->boxH, shiftImg);
        vSub(posPtr, posPtr, shiftPos);
    }
#endif
}
void sortParticle(Box *box, Particle *particle, Update *update) {
    if (particle->isSortForbidden)
        return;
    if (!update->nebrList.doSort)
        return;
    
    particle->sortFlag++;
    
    adjustImg(box, particle);
    
    NebrList *nebrList = &update->nebrList;
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    
    double sysRcs =
    particle->meanDiameter * (1.0 + nebrList->skinSet) / RcsBinLenRatio;
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
        vShiftAll(lamda, 0.5);
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
           particle->nAtom * sizeof(doubleVector));
    
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
    
    nebrList->isValid = false;
}
void buildNebrList(Box *box, Particle *particle, Update *update) {
    NebrList *nebrList = &update->nebrList;
    if (nebrList->xyzHold == NULL) {
        int nAtom = particle->nAtom;
        nebrList->xyzHold = (doubleVector *)calloc(nAtom, sizeof(doubleVector));
        nebrList->nNebr = (int2 *)calloc(nAtom, sizeof(int2));
        nebrList->binList4sort = (int *)calloc(nAtom, sizeof(int));
        nebrList->oid2nid = (int *)calloc(nAtom, sizeof(int));
        nebrList->buffer = calloc(nAtom, sizeof(doubleVector));
    }
    if (nebrList->boxPtrSave == NULL) {
        nebrList->boxPtrSave = box;
        nebrList->particlePtrSave = particle;
    } else {
        if (nebrList->boxPtrSave != box || nebrList->particlePtrSave != particle) {
            Abort("Inconsistency between NebrList and its host!");
        }
    }
    
    if (update->isThermalRun) {
        if ((nebrList->nBuild % 50 == 0)) {
            zeroMomentum(box, particle, update);
        }
        if ((nebrList->nBuild % 100 == 0)) {
            moveBarycenter(box, particle, update);
        }
    }
    if ((nebrList->nBuild % 200 == 0)) {
        sortParticle(box, particle, update);
    }
    
    adjustImg(box, particle);
    binParticle(box, particle, update);
    
    nebrList->meanDiameterHold = particle->meanDiameter;
    uptriMatCpy(nebrList->invBoxHold, box->invBoxH);
    memcpy(nebrList->xyzHold, particle->pos,
           particle->nAtom * sizeof(doubleVector));
    
    constructList(box, particle, update);
    
    nebrList->nBuild++;
    nebrList->isValid = true;
    nebrList->cntForce = 0;
}

void calcForceAllPair(Box *box, Particle *particle, Update *update) {
    if (particle->isForceValid)
        return;
    
    doubleVector *force = particle->force;
    memset(force, '\0', particle->nAtom * sizeof(doubleVector));
    double ePair = 0, maxOvlp = 0;
    uptriMat pVirTens, fabricTensor;
    uptriMatZeros(pVirTens);
    uptriMatZeros(fabricTensor);
    int nContact = 0;
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector iforce;
        vCpy(iforce, force[iatom]);
        doubleVector iPos;
        vCpy(iPos, particle->pos[iatom]);
        double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
        for (int jatom = iatom + 1; jatom < particle->nAtom; jatom++) {
            double jRc = particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
            double sRc = iRc + jRc;
            
            doubleVector dRij;
            vSub(dRij, iPos, particle->pos[jatom]);
            PBC(dRij, box);
            
            double rijP2 = sNormP2(dRij);
            if (rijP2 >= sRc * sRc)
                continue;
            double rij = sqrt(rijP2);
            maxOvlp = cpuMax(maxOvlp, 1.0 - rij / sRc);
            
#if defined(Harmonic)
            double rdivsig = rij / sRc;
            double fpair = (1.0 - rdivsig) / rij / sRc;
            ePair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
#elif defined(Hertzian)
            double delta = sRc - rij;
            // lammps style hertzian
            //  double fpair = sqrt(iRc * jRc / sRc * delta) * delta / rij;
            //  ePair += 0.4 * sqrt(iRc * jRc / sRc * delta) * delta * delta;
            double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
            ePair += (2.0 / 5.0) * pow(delta / sRc, 5.0 / 2.0);
#endif
            
            vScaleAdd(iforce, iforce, fpair, dRij);
            vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
            
            nContact++;
            for (int idim = 0; idim < DIM; idim++) {
                for (int jdim = idim; jdim < DIM; jdim++) {
                    pVirTens[spaceIdx2voigt(idim, jdim)] +=
                    fpair * dRij[idim] * dRij[jdim];
                    fabricTensor[spaceIdx2voigt(idim, jdim)] +=
                    dRij[idim] * dRij[jdim] / rijP2;
                }
            }
        }
        vCpy(force[iatom], iforce);
    }
    
    update->ePair = ePair / particle->nAtom;
    update->maxOvlp = maxOvlp;
    uptriMatCpy(update->pVirTens, pVirTens);
    uptriMatCpy(update->fabricTensor, fabricTensor);
    update->nContact = nContact;
    update->pVir = 0.0;
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            update->pVirTens[spaceIdx2voigt(idim, jdim)] /= box->volume;
            update->fabricTensor[spaceIdx2voigt(idim, jdim)] *= (-(double)DIM / (double)nContact);
        }
        update->pVir += update->pVirTens[spaceIdx2voigt(idim, idim)];
    }
    update->pVir /= DIM;
    
    update->Edone = update->Pdone = true;
    particle->isForceValid = true;
}

//===============================================================================
#if (DIM == 2 || DIM == 3)
void generateImageAtom(Box *box, Particle *particle, NebrList *nebrList) {
    // stencil-xlo
    int ith = 0;
    for (; ith < nebrList->imagePartCount[0]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] -= box->boxH[spaceIdx2voigt(0, 0)];
    }
    // stencil-xhi
    for (; ith < nebrList->imagePartCount[1]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] += box->boxH[spaceIdx2voigt(0, 0)];
    }
    // stencil-ylo
    for (; ith < nebrList->imagePartCount[2]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] -= box->boxH[spaceIdx2voigt(0, 1)];
        particle->pos[particle->nAtom + ith][1] -= box->boxH[spaceIdx2voigt(1, 1)];
    }
    // stencil-yhi
    for (; ith < nebrList->imagePartCount[3]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] += box->boxH[spaceIdx2voigt(0, 1)];
        particle->pos[particle->nAtom + ith][1] += box->boxH[spaceIdx2voigt(1, 1)];
    }
#if (DIM == 3)
    // stencil-zlo
    for (; ith < nebrList->imagePartCount[4]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] -= box->boxH[spaceIdx2voigt(0, 2)];
        particle->pos[particle->nAtom + ith][1] -= box->boxH[spaceIdx2voigt(1, 2)];
        particle->pos[particle->nAtom + ith][2] -= box->boxH[spaceIdx2voigt(2, 2)];
    }
    // stencil-zhi
    for (; ith < nebrList->imagePartCount[5]; ith++) {
        int iorig = nebrList->imageParent[ith];
        vCpy(particle->pos[particle->nAtom + ith], particle->pos[iorig]);
        particle->pos[particle->nAtom + ith][0] += box->boxH[spaceIdx2voigt(0, 2)];
        particle->pos[particle->nAtom + ith][1] += box->boxH[spaceIdx2voigt(1, 2)];
        particle->pos[particle->nAtom + ith][2] += box->boxH[spaceIdx2voigt(2, 2)];
    }
#endif
}
#else
#define generateImageAtom(box, particle, nebrList) ;
#endif
/*=======
 For DIM == 2 or DIM == 3, the PBC is treated during building nebrList by
 including image particles in the Particle structure.
 ========*/
void calcForce(Box *box, Particle *particle, Update *update) {
    // half-style NebrList
    NebrList *nebrList = &update->nebrList;
    if (particle->isForceValid)
        return;
    
#if false
    {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVecPtr posPtr = particle->pos[iatom];
            intVecPtr imgPtr = particle->img[iatom];
            
            // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
            doubleVector lamda;
            MatMulVec(lamda, box->invBoxH, posPtr);
            vShiftAll(lamda, 0.5);
            
            intVector shiftImg;
            vFloor(shiftImg, lamda);
            vAdd(imgPtr, imgPtr, shiftImg);
            
            doubleVector shiftPos;
            MatMulVec(shiftPos, box->boxH, shiftImg);
            vSub(posPtr, posPtr, shiftPos);
        }
        
        
        
        doubleVector *force = particle->force;
        memset(force, '\0', particle->nAtom * sizeof(doubleVector));
        double ePair = 0, maxOvlp = 0;
        uptriMat pVirTens, fabricTensor;
        uptriMatZeros(pVirTens);
        uptriMatZeros(fabricTensor);
        int nContact = 0;
        
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVector iforce;
            vCpy(iforce, force[iatom]);
            doubleVector iPos;
            vCpy(iPos, particle->pos[iatom]);
            double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
            for(int jatom = iatom + 1; jatom < particle->nAtom;jatom++){
                double jRc = particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
                double sRc = iRc + jRc;
                
                doubleVector dRij;
                vSub(dRij, iPos, particle->pos[jatom]);
                
                { /*adjust the image of starting Point*/
                    while (dRij[2] >= 0.5 * box->boxH[spaceIdx2voigt(2, 2)]) {
                        dRij[2] -= box->boxH[spaceIdx2voigt(2, 2)];
                        dRij[1] -= box->boxH[spaceIdx2voigt(1, 2)];
                        dRij[0] -= box->boxH[spaceIdx2voigt(0, 2)];
                    }
                    while (dRij[2] < -0.5 * box->boxH[spaceIdx2voigt(2, 2)]) {
                        dRij[2] += box->boxH[spaceIdx2voigt(2, 2)];
                        dRij[1] += box->boxH[spaceIdx2voigt(1, 2)];
                        dRij[0] += box->boxH[spaceIdx2voigt(0, 2)];
                    }
                    while (dRij[1] >= 0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {
                        dRij[1] -= box->boxH[spaceIdx2voigt(1, 1)];
                        dRij[0] -= box->boxH[spaceIdx2voigt(0, 1)];
                    }
                    while (dRij[1] < -0.5 * box->boxH[spaceIdx2voigt(1, 1)]) {
                        dRij[1] += box->boxH[spaceIdx2voigt(1, 1)];
                        dRij[0] += box->boxH[spaceIdx2voigt(0, 1)];
                    }
                    while (dRij[0] >= 0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {
                        dRij[0] -= box->boxH[spaceIdx2voigt(0, 0)];
                    }
                    while (dRij[0] < -0.5 * box->boxH[spaceIdx2voigt(0, 0)]) {
                        dRij[0] += box->boxH[spaceIdx2voigt(0, 0)];
                    }
                }
                
                double rijP2 = sNormP2(dRij);
                if (rijP2 >= sRc * sRc)
                    continue;
                double rij = sqrt(rijP2);
                maxOvlp = cpuMax(maxOvlp, 1.0 - rij / sRc);
                
#if defined(Harmonic)
                double rdivsig = rij / sRc;
                double fpair = (1.0 - rdivsig) / rij / sRc;
                ePair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
#elif defined(Hertzian)
                double delta = sRc - rij;
                // lammps style hertzian
                //  double fpair = sqrt(iRc * jRc / sRc * delta) * delta / rij;
                //  ePair += 0.4 * sqrt(iRc * jRc / sRc * delta) * delta * delta;
                double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
                ePair += (2.0 / 5.0) * pow(delta / sRc, 5.0 / 2.0);
#endif
                
                vScaleAdd(iforce, iforce, fpair, dRij);
                vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
                
                nContact++;
                for (int idim = 0; idim < DIM; idim++) {
                    for (int jdim = idim; jdim < DIM; jdim++) {
                        pVirTens[spaceIdx2voigt(idim, jdim)] +=
                        fpair * dRij[idim] * dRij[jdim];
                        fabricTensor[spaceIdx2voigt(idim, jdim)] +=
                        dRij[idim] * dRij[jdim] / rijP2;
                    }
                }
            }
            vCpy(force[iatom], iforce);
        }
        
        update->ePair = ePair / particle->nAtom;
        update->maxOvlp = maxOvlp;
        uptriMatCpy(update->pVirTens, pVirTens);
        uptriMatCpy(update->fabricTensor, fabricTensor);
        update->nContact = nContact;
        update->pVir = 0.0;
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim; jdim < DIM; jdim++) {
                update->pVirTens[spaceIdx2voigt(idim, jdim)] /= box->volume;
                update->fabricTensor[spaceIdx2voigt(idim, jdim)] *=
                (-(double)DIM / (double)nContact);
            }
            update->pVir += update->pVirTens[spaceIdx2voigt(idim, idim)];
        }
        update->pVir /= DIM;
        
        update->Edone = update->Pdone = true;
        particle->isForceValid = true;
        nebrList->cntForce++;
        nebrList->nForce++;
    }
#endif
    
    bool isListValid = isNebrListValid(box, particle, update);
    if (update->isThermalRun) {
        if (!isListValid) {
            if (nebrList->cntForce < nebrList->nRebuildMax) {
                nebrList->skinSet *= 1.1;
                nebrList->skinSet = cpuMin(nebrList->skinSet, __maxSkinSet__);
                nebrList->compelInit = true;
            } else if (nebrList->cntForce > nebrList->nRebuildMax * 2) {
                nebrList->skinSet *= 0.9;
                nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__ * 10.0);
                nebrList->compelInit = true;
            }
        } else {
            if (nebrList->skinSet > 10.0 * __minSkinSet__ * (1 + 1E-10) &&
                nebrList->cntForce > nebrList->nRebuildMax * 2) {
                nebrList->skinSet *= 0.9;
                nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__ * 10.0);
                nebrList->compelInit = true;
                isListValid = false;
            }
        }
    } else {
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
    }
    if (!isListValid) {
        buildNebrList(box, particle, update);
    }
    
    doubleVector *force = particle->force;
    memset(force, '\0', particle->nAtom * sizeof(doubleVector));
    double ePair = 0, maxOvlp = 0;
    uptriMat pVirTens, fabricTensor;
    uptriMatZeros(pVirTens);
    uptriMatZeros(fabricTensor);
    int nContact = 0;
    
    generateImageAtom(box, particle, nebrList);
    
    if (update->isThermalRun) {
        switch (update->evFlag) {
        case 2:
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                doubleVector iforce;
                vCpy(iforce, force[iatom]);
                doubleVector iPos;
                vCpy(iPos, particle->pos[iatom]);
                double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
                for (int jth = nebrList->nNebr[iatom].first;
                     jth < nebrList->nNebr[iatom].second; jth++) {
                    int jimage = nebrList->list[jth];
                    int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
                    double jRc = particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
                    double sRc = iRc + jRc;
                    
                    doubleVector dRij;
                    vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
                    PBC(dRij, box);
#endif
                    double rijP2 = sNormP2(dRij);
                    if (rijP2 >= sRc * sRc)
                        continue;
                    double rij = sqrt(rijP2);
                    
#if defined(Harmonic)
                    double rdivsig = rij / sRc;
                    double fpair = (1.0 - rdivsig) / rij / sRc;
                    ePair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
#elif defined(Hertzian)
                    double delta = sRc - rij;
                    double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
                    ePair += (2.0 / 5.0) * pow(delta / sRc, 5.0 / 2.0);
#endif
                    
                    vScaleAdd(iforce, iforce, fpair, dRij);
                    vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
                    
                    for (int idim = 0; idim < DIM; idim++) {
                        for (int jdim = idim; jdim < DIM; jdim++) {
                            pVirTens[spaceIdx2voigt(idim, jdim)] += fpair * dRij[idim] * dRij[jdim];
                        }
                    }
                }
                vCpy(force[iatom], iforce);
            }
            update->ePair = ePair / particle->nAtom;
            uptriMatCpy(update->pVirTens, pVirTens);
            for (int idim = 0; idim < DIM; idim++) {
                for (int jdim = idim; jdim < DIM; jdim++) {
                    update->pVirTens[spaceIdx2voigt(idim, jdim)] /= box->volume;
                }
            }
            update->Edone = update->Wdone = true;
            break;
        default:
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                doubleVector iforce;
                vCpy(iforce, force[iatom]);
                doubleVector iPos;
                vCpy(iPos, particle->pos[iatom]);
                double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
                for (int jth = nebrList->nNebr[iatom].first;
                     jth < nebrList->nNebr[iatom].second; jth++) {
                    int jimage = nebrList->list[jth];
                    int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
                    double jRc =
                    particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
                    double sRc = iRc + jRc;
                    
                    doubleVector dRij;
                    vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
                    PBC(dRij, box);
#endif
                    double rijP2 = sNormP2(dRij);
                    if (rijP2 >= sRc * sRc)
                        continue;
                    double rij = sqrt(rijP2);
                    
#if defined(Harmonic)
                    double rdivsig = rij / sRc;
                    double fpair = (1.0 - rdivsig) / rij / sRc;
                    // ePair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
#elif defined(Hertzian)
                    double delta = sRc - rij;
                    double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
                    // ePair += (2.0 / 5.0) * pow(delta / sRc, 5.0 / 2.0);
#endif
                    
                    vScaleAdd(iforce, iforce, fpair, dRij);
                    vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
                }
                vCpy(force[iatom], iforce);
            }
            update->ePair = 0;
            uptriMatCpy(update->pVirTens, pVirTens);
            update->Edone = update->Wdone = false;
            break;
        }
        
        particle->isForceValid = true;
        nebrList->cntForce++;
        nebrList->nForce++;
    } else {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVector iforce;
            vCpy(iforce, force[iatom]);
            doubleVector iPos;
            vCpy(iPos, particle->pos[iatom]);
            double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
            for (int jth = nebrList->nNebr[iatom].first;
                 jth < nebrList->nNebr[iatom].second; jth++) {
                int jimage = nebrList->list[jth];
                int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
                double jRc =
                particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
                double sRc = iRc + jRc;
                
                doubleVector dRij;
                vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
                PBC(dRij, box);
#endif
                double rijP2 = sNormP2(dRij);
                if (rijP2 >= sRc * sRc)
                    continue;
                double rij = sqrt(rijP2);
                maxOvlp = cpuMax(maxOvlp, 1.0 - rij / sRc);
                
#if defined(Harmonic)
                double rdivsig = rij / sRc;
                double fpair = (1.0 - rdivsig) / rij / sRc;
                ePair += 0.5 * (1.0 - rdivsig) * (1.0 - rdivsig);
#elif defined(Hertzian)
                double delta = sRc - rij;
                // lammps style hertzian
                //  double fpair = sqrt(iRc * jRc / sRc * delta) * delta / rij;
                //  ePair += 0.4 * sqrt(iRc * jRc / sRc * delta) * delta * delta;
                double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
                ePair += (2.0 / 5.0) * pow(delta / sRc, 5.0 / 2.0);
#endif
                
                vScaleAdd(iforce, iforce, fpair, dRij);
                vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
                
                nContact++;
                for (int idim = 0; idim < DIM; idim++) {
                    for (int jdim = idim; jdim < DIM; jdim++) {
                        pVirTens[spaceIdx2voigt(idim, jdim)] +=
                        fpair * dRij[idim] * dRij[jdim];
                        fabricTensor[spaceIdx2voigt(idim, jdim)] +=
                        dRij[idim] * dRij[jdim] / rijP2;
                    }
                }
            }
            vCpy(force[iatom], iforce);
        }
        
        update->ePair = ePair / particle->nAtom;
        update->maxOvlp = maxOvlp;
        uptriMatCpy(update->pVirTens, pVirTens);
        uptriMatCpy(update->fabricTensor, fabricTensor);
        update->nContact = nContact;
        update->pVir = 0.0;
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim; jdim < DIM; jdim++) {
                update->pVirTens[spaceIdx2voigt(idim, jdim)] /= box->volume;
                update->fabricTensor[spaceIdx2voigt(idim, jdim)] *=
                (-(double)DIM / (double)nContact);
            }
            update->pVir += update->pVirTens[spaceIdx2voigt(idim, idim)];
        }
        update->pVir /= DIM;
        
        update->Edone = update->Pdone = true;
        particle->isForceValid = true;
        nebrList->cntForce++;
        nebrList->nForce++;
    }
}

void genHalfNebrListPBCimage(Box *box, Particle *particle, Update *update){
    NebrList *nebrList = &update->nebrList;
    
    bool isListValid = isNebrListValid(box, particle, update);
    if (update->isThermalRun) {
        if (!isListValid) {
            if (nebrList->cntForce < nebrList->nRebuildMax) {
                nebrList->skinSet *= 1.1;
                nebrList->skinSet = cpuMin(nebrList->skinSet, __maxSkinSet__);
                nebrList->compelInit = true;
            } else if (nebrList->cntForce > nebrList->nRebuildMax * 2) {
                nebrList->skinSet *= 0.9;
                nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__ * 10.0);
                nebrList->compelInit = true;
            }
        } else {
            if (nebrList->skinSet > 10.0 * __minSkinSet__ * (1 + 1E-10) &&
                nebrList->cntForce > nebrList->nRebuildMax * 2) {
                nebrList->skinSet *= 0.9;
                nebrList->skinSet = cpuMax(nebrList->skinSet, __minSkinSet__ * 10.0);
                nebrList->compelInit = true;
                isListValid = false;
            }
        }
    } else {
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
    }
    if (!isListValid) {
        buildNebrList(box, particle, update);
    }
        
    generateImageAtom(box, particle, nebrList);
}


// compute energy expansion coefficients.
// input: |deltaVeloc| = 1.0;
// dE = dot(grad E, dVel) * x + 1/2 * contract(Hess, dVel) * x^2 + 1/6 * contract(M3, dVel) * x^3 + 1/24 * contract(M4, dVel) * x^4 +  o(x^5);
// c1 = dot(grad E, dVel); c2 = contract(Hess, dVel); c3 = contract(M3, dVel), c4 = contract(M4, dVel);
int calcExpanCoeff(Box *box, Particle *particle, Update *update, doubleVector *deltaVeloc, double coeff[4]) {// output unit of coeff[..] is energy.
    //build Half Neighbour List
    calcForce(box, particle, update);
    
    double lenP2 = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        lenP2 += sNormP2(deltaVeloc[iatom]);
    }
    if (lenP2 < 1E-24) {
        Info("Invalid deltaVeloc");
        return -1;
    }
    
    NebrList *nebrList = &update->nebrList;
    double y1 = 0, y2 = 0, y3 = 0, y4 = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector iPos;
        vCpy(iPos, particle->pos[iatom]);
        double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
        for (int jth = nebrList->nNebr[iatom].first; jth < nebrList->nNebr[iatom].second; jth++) {
            int jimage = nebrList->list[jth];
            int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
            double jRc = particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
            double sRc = iRc + jRc;
            doubleVector dRij;
            vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
            PBC(dRij, box);
#endif
            double rijP2 = sNormP2(dRij);
            if (rijP2 >= sRc * sRc)
                continue;
            double rij = sqrt(rijP2);
            // Hij * dVij
            doubleVector uRij;  // unit vecter along dRij
            vScale(uRij, 1.0 / rij, dRij);
            doubleVector dVij;
            vSub(dVij, deltaVeloc[iatom], deltaVeloc[jatom]);
            
#if defined(Harmonic)
            double rdivsig = rij / sRc;
            double fpair = (1.0 - rdivsig) / rij / sRc;
            
            double tij = -(1.0 - rdivsig) / sRc;  // first derivative
            double sij = 1 / sRc / sRc;           // second derivative
            double cij = 0.0;                     // third derivative
            double dij = 0.0;                     // fourth derivative
#elif defined(Hertzian)
            double delta = sRc - rij;
            double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
            
            double tij = -sqrt(delta / sRc) * delta / sRc / sRc;
            double sij = 3.0 / 2.0 * sqrt(delta / sRc) / sRc / sRc;
            double cij = -3.0 / 4.0 / sqrt(delta / sRc) / sRc / sRc / sRc;
            double dij = -3.0 / 8.0 * pow(delta / sRc, -1.5) * pow(sRc, -4.0);
#endif
            // First order expansion
            y1 += tij * sDot(uRij, dVij);
            
            // Second order expansion
            {
                double afact = sij - tij / rij;
                double bfact = tij / rij;
                y2 += afact * pow(sDot(uRij, dVij), 2.0) + bfact * sNormP2(dVij);
            }
            // Third order expansion
            {
                double afact = cij - 3.0 * sij / rij + 3.0 * tij / rij / rij;
                double bfact = sij / rij - tij / rij / rij;
                y3 += afact * pow(sDot(uRij, dVij), 3.0) + 3.0 * bfact * sDot(uRij, dVij) * sNormP2(dVij);
            }
            // Fourth order expansion
            {
                double afact = dij - 6.0 * cij / rij + 15.0 * sij / rij / rij - 15.0 * tij / rij / rij / rij;
                double bfact = 6.0 * cij / rij - 18.0 * sij / rij / rij + 18.0 * tij / rij / rij / rij;
                double cfact = 3.0 * sij / rij / rij - 3.0 * tij / rij / rij / rij;
                y4 += afact * pow(sDot(uRij, dVij), 4.0) + bfact * pow(sDot(uRij, dVij), 2.0) * sNormP2(dVij) + cfact * pow(sNormP2(dVij), 2.0);
            }
        }
    }
    
    coeff[0] = y1 / sqrt(lenP2) * update->distanceUnits;
    coeff[1] = y2 / lenP2 * pow(update->distanceUnits, 2.0);
    coeff[2] = y3 / pow(lenP2, 1.5) * pow(update->distanceUnits, 3.0);
    coeff[3] = y4 / lenP2 / lenP2  * pow(update->distanceUnits, 4.0);
    
    return 0;
}


//===============================================================================
void setBoxPara(Box *box) {
    box->dim = DIM;
    box->volume = 1.0;
    for (int idim = 0; idim < DIM; idim++) {
        if (box->boxH[spaceIdx2voigt(idim, idim)] <= 0.0) {
            Abort("Not rihgt hand basis!");
        }
        box->volume *= box->boxH[spaceIdx2voigt(idim, idim)];
        
        vZeros(box->boxEdge[idim]);
        for (int jdim = 0; jdim <= idim; jdim++) {
            box->boxEdge[idim][jdim] = box->boxH[spaceIdx2voigt(jdim, idim)];
        }
    }
    
    vZeros(box->cornerLo);
    for (int idim = 0; idim < DIM; idim++) {
        vAdd(box->cornerLo, box->cornerLo, box->boxEdge[idim]);
    }
    vScale(box->cornerLo, -0.5, box->cornerLo);
    vScale(box->cornerHi, -1.0, box->cornerLo);
    
    uptriMatZeros(box->invBoxH);
    for (int idim = DIM - 1; idim >= 0; idim--) {
        box->invBoxH[spaceIdx2voigt(idim, idim)] =
        1.0 / box->boxH[spaceIdx2voigt(idim, idim)];
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            double cij = 0.0;
            for (int k = idim + 1; k <= jdim; k++) {
                cij += box->boxH[spaceIdx2voigt(idim, k)] *
                box->invBoxH[spaceIdx2voigt(k, jdim)];
            }
            box->invBoxH[spaceIdx2voigt(idim, jdim)] =
            -cij / box->boxH[spaceIdx2voigt(idim, idim)];
        }
    }
}
void calcBoxLenAngle(Box *box, uptriMat para) {
    for (int idim = 0; idim < DIM; idim++) {
        para[spaceIdx2voigt(idim, idim)] = box->boxH[spaceIdx2voigt(idim, idim)];
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            double len1 = sNorm(box->boxEdge[idim]);
            double len2 = sNorm(box->boxEdge[jdim]);
            double c12 = sDot(box->boxEdge[idim], box->boxEdge[jdim]);
            para[spaceIdx2voigt(idim, jdim)] = acos(c12 / (len1 * len2)) / PI * 180.0;
        }
    }
}
void setUnits(Update *update, double distanceUnits) {
    update->massUnits = 1.0;
    update->energyUnits = 1.0;
    update->distanceUnits = distanceUnits;
    
    update->timeUnits =
    sqrt(update->massUnits / update->energyUnits) * update->distanceUnits;
    update->forceUnits = update->energyUnits / update->distanceUnits;
    update->velocityUnits = sqrt(update->energyUnits / update->massUnits);
    update->pressureUnits = update->energyUnits / pow(update->distanceUnits, DIM);
    update->volumeUnits = pow(update->distanceUnits, DIM);
}

void instant_inflate(Box *box, Particle *particle, Update *update, double deltaVF) {
    double volfrac_target = update->volFrac + deltaVF;
    double sfact = pow(volfrac_target / update->volFrac, 1.0 / DIM);
    particle->meanDiameter *= sfact;
    update->volFrac = pow(sfact, DIM) * update->volFrac;
    setUnits(update, particle->meanDiameter);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}

#if defined(__triBox__)
void instant_simpShearXz(Box *box, Particle *particle, Update *update, double deltaGamma) {
    adjustImg(box, particle);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->pos[iatom][0] += particle->pos[iatom][DIM - 1] * deltaGamma;
    }
    
    box->boxH[spaceIdx2voigt(0, DIM - 1)] += deltaGamma * box->boxH[spaceIdx2voigt(DIM - 1, DIM - 1)];
    
    setBoxPara(box);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}
#if (DIM == 3)
void instant_simpShearXy(Box *box, Particle *particle, Update *update, double deltaGamma) {
    adjustImg(box, particle);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->pos[iatom][0] += particle->pos[iatom][1] * deltaGamma;
    }
    box->boxH[spaceIdx2voigt(0, 2)] += box->boxH[spaceIdx2voigt(1, 2)] * deltaGamma;
    box->boxH[spaceIdx2voigt(0, 1)] += box->boxH[spaceIdx2voigt(1, 1)] * deltaGamma;
    
    setBoxPara(box);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}
void instant_simpShearYz(Box *box, Particle *particle, Update *update, double deltaGamma) {
    adjustImg(box, particle);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->pos[iatom][1] += particle->pos[iatom][2] * deltaGamma;
    }
    
    box->boxH[spaceIdx2voigt(1, 2)] += box->boxH[spaceIdx2voigt(2, 2)] * deltaGamma;
    setBoxPara(box);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}
#endif
#endif

void instant_pureShearXz(Box *box, Particle *particle, Update *update, double deltaStrain) {
    adjustImg(box, particle);
    
    double scaleX = exp(0.5 * deltaStrain);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        particle->pos[iatom][DIM - 1] /= scaleX;
        particle->pos[iatom][0] *= scaleX;
    }
    
    for (int idim = 0; idim < DIM - 1; idim++) {
        box->boxH[spaceIdx2voigt(0, idim)] *= scaleX;
    }
    box->boxH[spaceIdx2voigt(DIM - 1, DIM - 1)] /= scaleX;
    
    setBoxPara(box);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}
void instant_deformation(Box *box, Particle *particle, Update *update, uptriMat transMat) {
    // Rt = R0 + transMat * R0;
    adjustImg(box, particle);
    
    uptriMat trans;
    uptriMatCpy(trans, transMat);
    
    for (int idim = 0; idim < DIM; idim++) {
        trans[spaceIdx2voigt(idim, idim)] += 1.0;
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        MatMulVec(particle->pos[iatom], trans, particle->pos[iatom]);
    }
    
    // update edge
    for (int iedge = 0; iedge < DIM; iedge++) {
        for (int idim = 0; idim < DIM; idim++) {
            box->boxEdge[iedge][idim] *= trans[spaceIdx2voigt(idim, idim)];
            for (int jdim = idim + 1; jdim < DIM; jdim++) {
                box->boxEdge[iedge][idim] +=
                trans[spaceIdx2voigt(idim, jdim)] * box->boxEdge[iedge][jdim];
            }
        }
    }
    // construct boxH from edge
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            box->boxH[spaceIdx2voigt(idim, jdim)] = box->boxEdge[jdim][idim];
        }
    }
    setBoxPara(box);
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
}
void normaliseBox(Box *box, Particle *particle, Update *update) {
    double sfact = pow(box->volume, 1.0 / DIM);
    particle->meanDiameter /= sfact;
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            box->boxH[spaceIdx2voigt(idim, jdim)] /= sfact;
        }
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScale(particle->pos[iatom], 1.0 / sfact, particle->pos[iatom]);
    }
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    
    setBoxPara(box);
    setUnits(update, particle->meanDiameter);
}

//===============================================================================
void backupSimInfo(Box *box, Particle *particle, Box *bakBox, Particle *bakParticle) {
    if (bakParticle->nAtom != particle->nAtom) {
        bakParticle->nAtom = particle->nAtom;
        
        bakParticle->pos = (doubleVector *)realloc(
                                                   bakParticle->pos, bakParticle->nAtom * sizeof(doubleVector));
        bakParticle->veloc = (doubleVector *)realloc(
                                                     bakParticle->veloc, bakParticle->nAtom * sizeof(doubleVector));
        bakParticle->force =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
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
    memcpy(bakParticle->force, particle->force,
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
    
    bakParticle->meanDiameter = particle->meanDiameter;
    bakParticle->isSizeFixed = particle->isSizeFixed;
    bakParticle->isForceValid = particle->isForceValid;
    
    memcpy(bakBox, box, sizeof(Box));
}
void restoreSimInfo(Box *box, Particle *particle, Update *update, Box *bakBox, Particle *bakParticle) {
    // restore data
    backupSimInfo(bakBox, bakParticle, box, particle);
    
    // clear
    particle->isForceValid = false;
    
    update->Edone = false;
    update->Pdone = false;
    update->Tdone = false;
    
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    
    // reinit
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
}
void freeSimInfo(Box *bakBox, Particle *bakParticle) {
    safeFree(bakParticle->pos);
    safeFree(bakParticle->uPos);
    safeFree(bakParticle->veloc);
    safeFree(bakParticle->mass);
    safeFree(bakParticle->img);
    safeFree(bakParticle->type);
    safeFree(bakParticle->diameterScale);
    safeFree(bakParticle->id2tag);
    safeFree(bakParticle->tag2id);
}
void reInitSim(Box *box, Particle *particle, Update *update) {
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    particle->isForceValid = false;
    update->Edone = update->Pdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.nForce = 0;
    update->nebrList.cntForce = 0;
    update->nebrList.doSort = true;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
}

void showInfo(Box *box, Particle *particle, Update *update) {
    calcForce(box, particle, update);
    
    printf("===========System Info==========\n");
    printf("Dimension: %d\n", DIM);
    printf("No. of Particles: %d\n", particle->nAtom);
    printf("Mean diameter: %g (unit 1)\n", particle->meanDiameter);
    printf("Min(diameter): %g;\nMax(diameter): %g;\n",
           update->nebrList.minDiameterScale, update->nebrList.maxDiameterScale);
    printf("VolFrac: %g\n", update->volFrac);
    printf("Edges of Simulation Box: \n");
    for (int iedge = 0; iedge < DIM; iedge++) {
        printf("\t");
        for (int jdim = 0; jdim < DIM; jdim++) {
            printf("%-8.6e\t", box->boxEdge[iedge][jdim] / update->distanceUnits);
        }
        printf("\n");
    }
    printf("Epair: %g\n", update->ePair / update->energyUnits);
    {
        double sumForce = 0;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            doubleVecPtr force = particle->force[iatom];
            sumForce += sNorm(force);
        }
        printf("AveForce: %g\n", sumForce / update->forceUnits / particle->nAtom);
    }
    printf("Pressure: %g\n", update->pVir / update->pressureUnits);
    printf("Normal Ptensor:\n\t");
    for (int idim = 0; idim < DIM; idim++) {
        printf("%g ", update->pVirTens[spaceIdx2voigt(idim, idim)] /
               update->pressureUnits);
    }
    printf("\nTangential Ptensor:\n\t");
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            printf("%g ", update->pVirTens[spaceIdx2voigt(idim, jdim)] /
                   update->pressureUnits);
        }
    }
    printf("\n");
    printf("Normal fabricTensor:\n\t");
    for (int idim = 0; idim < DIM; idim++) {
        printf("%g ", update->fabricTensor[spaceIdx2voigt(idim, idim)]);
    }
    printf("\nTangential fabricTensor:\n\t");
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            printf("%g ", update->fabricTensor[spaceIdx2voigt(idim, jdim)]);
        }
    }
    printf("\n");
    
    double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
    for (int idim = 0; idim < DIM; idim++) {
        minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
    }
    double maxRcut = update->nebrList.maxDiameterScale * particle->meanDiameter;
    double minRskin = __minSkinSet__ * update->nebrList.minDiameterScale *
    particle->meanDiameter;
    
    double minLen = (maxRcut + minRskin) * 2.0;
    double minBoxVol = pow(minLen, DIM);
    double maxPhi = update->volFrac * box->volume / minBoxVol;
    printf(
           "The 0.5*min(box(.,.)) is %g, max(Rcut) is %g, "
           "min(Rskin) is %g. The max(volFrac) is ~ %g.\n",
           0.5 * minAxByCz / update->distanceUnits, maxRcut / update->distanceUnits,
           minRskin / update->distanceUnits, maxPhi);
    if (maxRcut + minRskin >= 0.5 * minAxByCz) {
        printf("The PBC method breaks down!\n");
    }
    printf("===========System Info==========\n");
}

void showStructInfo(Box *box, Particle *particle, Update *update) {
    printf("===========Structure Info==========\n");
    contactInfo *cinfo = addContactInfo(box, particle, update);
    computeContactInfo(box, particle, update);
    
    printf("Solition: %d, Rattlers: %d\n", cinfo->nSoliton, cinfo->nRattler);
    printf("Z: %g; Z_NR: %g; Ratio_R: %g\n", cinfo->aveCoordNum,
           cinfo->aveCoordNumExRattler,
           (double)cinfo->nRattler / (double)particle->nAtom);
    printf("AveForceExRattler: %g\n", cinfo->meanForceExRattler);
    printf("===========Structure Info==========\n");
    
    delContactInfo(update);
}

//===============================================================================
contactInfo *getContactInfo(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "contactInfo");
    if (whichTool < 0)
        return (contactInfo *)NULL;
    return (contactInfo *)update->toolkit.toolkit[whichTool];
}
contactInfo *addContactInfo(Box *box, Particle *particle, Update *update) {
    int whichTool = findToolkit(&update->toolkit, "contactInfo");
    contactInfo *cinfo = NULL;
    if (whichTool >= 0) {
        cinfo = (contactInfo *)update->toolkit.toolkit[whichTool];
    } else {
        cinfo = (contactInfo *)calloc(sizeof(contactInfo), 1);
        addToolkit(&update->toolkit, (void *)cinfo, NULL, "contactInfo");
    }
    if (particle->nAtom <= 0)
        Abort("No Atom!");
    cinfo->refCnt++;
    
    if (cinfo->isRattler == NULL) {
        cinfo->isSoliton = (bool *)calloc(particle->nAtom, sizeof(bool));
        cinfo->isRattler = (bool *)calloc(particle->nAtom, sizeof(bool));
        cinfo->isBuckler = (bool *)calloc(particle->nAtom, sizeof(bool));
        cinfo->nCoordNumExRattler = (int *)calloc(particle->nAtom, sizeof(int));
        cinfo->nCoordNum = (int *)calloc(particle->nAtom, sizeof(int));
        cinfo->forceExRattler =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
    }
    
    return cinfo;
}
int computeContactInfo(Box *box, Particle *particle, Update *update) {
    contactInfo *cinfo = getContactInfo(update);
    if (cinfo == NULL)
        return -1;
    NebrList *nebrList = &update->nebrList;
    
    if (!isNebrListValid(box, particle, update)) {
        buildNebrList(box, particle, update);
    }
    if (cinfo->allocMem4CoordNum < nebrList->maxAllocNebr) {
        cinfo->allocMem4CoordNum = nebrList->maxAllocNebr;
        cinfo->isNebrContact = (bool *)realloc(cinfo->isNebrContact, cinfo->allocMem4CoordNum * sizeof(bool));
    }
    
    memset(cinfo->isSoliton, '\0', particle->nAtom * sizeof(bool));
    memset(cinfo->isRattler, '\0', particle->nAtom * sizeof(bool));
    memset(cinfo->isBuckler, '\0', particle->nAtom * sizeof(bool));
    memset(cinfo->nCoordNumExRattler, '\0', particle->nAtom * sizeof(int));
    memset(cinfo->nCoordNum, '\0', particle->nAtom * sizeof(int));
    memset(cinfo->isNebrContact, '\0', cinfo->allocMem4CoordNum * sizeof(bool));
    cinfo->isJammed = false;
    
    generateImageAtom(box, particle, nebrList);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector iPos;
        vCpy(iPos, particle->pos[iatom]);
        double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
        for (int jth = nebrList->nNebr[iatom].first;
             jth < nebrList->nNebr[iatom].second; jth++) {
            int jimage = nebrList->list[jth];
            int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
            double jRc =
            particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
            double sRc = iRc + jRc;
            
            doubleVector dRij;
            vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
            PBC(dRij, box);
#endif
            double rijP2 = sNormP2(dRij);
            if (rijP2 >= sRc * sRc) {
                continue;
            } else {
                cinfo->isNebrContact[jth] = true;
                cinfo->nCoordNum[jatom]++;
                cinfo->nCoordNum[iatom]++;
            }
        }
    }
    
    int cordNumWR = 0, cordNumNR = 0, nRattler = 0, nSoliton = 0;
    bool done = true;
    memcpy(cinfo->nCoordNumExRattler, cinfo->nCoordNum, particle->nAtom * sizeof(int));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        if (cinfo->nCoordNum[iatom] == 0) {
            cinfo->isSoliton[iatom] = true;
            nSoliton++;
        }
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
            if (cinfo->isRattler[iatom])
                continue;
            for (int jth = nebrList->nNebr[iatom].first;
                 jth < nebrList->nNebr[iatom].second; jth++) {
                int jimage = nebrList->list[jth];
                int jatom =
                image2local(jimage, particle->nAtom, nebrList->imgageSource);
                
                if (!cinfo->isNebrContact[jth])
                    continue;
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
    cinfo->nSoliton = nSoliton;
    cinfo->aveCoordNum = (double)cordNumWR / (double)particle->nAtom;
    if (nRattler != particle->nAtom) {
        cinfo->aveCoordNumExRattler = (double)cordNumNR / (double)(particle->nAtom - nRattler);
        
        if (nRattler > 0.5 * particle->nAtom) {
            cinfo->isJammed = false;
        } else {
            // with compression boundary freedom. Ref.: PRE, 90, 022138.
            double dZ = cinfo->aveCoordNumExRattler - 2 * DIM + (2 * DIM - 2) / (particle->nAtom - cinfo->nRattler);
            if (dZ < -1.0 / (particle->nAtom - cinfo->nRattler))
                cinfo->isJammed = false;
            else
                cinfo->isJammed = true;
        }
        
    } else {
        cinfo->isJammed = false;
        cinfo->aveCoordNumExRattler = 0;
    }
    
    cinfo->nBuckler = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        if (cinfo->nCoordNumExRattler[iatom] == DIM + 1) {
            cinfo->isBuckler[iatom] = true;
            cinfo->nBuckler++;
        } else
            cinfo->isBuckler[iatom] = false;
    }
    
    if (particle->nAtom == cinfo->nRattler) {
        cinfo->meanForceExRattler = 0;
        return 0;
    }
    doubleVector *force = cinfo->forceExRattler;
    memset(force, '\0', particle->nAtom * sizeof(doubleVector));
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        if (cinfo->isRattler[iatom])
            continue;
        doubleVector iforce;
        vCpy(iforce, force[iatom]);
        doubleVector iPos;
        vCpy(iPos, particle->pos[iatom]);
        double iRc = particle->diameterScale[iatom] * particle->meanDiameter * 0.5;
        for (int jth = nebrList->nNebr[iatom].first;
             jth < nebrList->nNebr[iatom].second; jth++) {
            int jimage = nebrList->list[jth];
            int jatom = image2local(jimage, particle->nAtom, nebrList->imgageSource);
            if (cinfo->isRattler[jatom])
                continue;
            
            doubleVector dRij;
            double jRc =
            particle->diameterScale[jatom] * particle->meanDiameter * 0.5;
            double sRc = iRc + jRc;
            vSub(dRij, iPos, particle->pos[jimage]);
#if !(DIM == 2 || DIM == 3)
            PBC(dRij, box);
#endif
            double rijP2 = sNormP2(dRij);
            if (rijP2 >= sRc * sRc)
                continue;
            double rij = sqrt(rijP2);
            
#if defined(Harmonic)
            double rdivsig = rij / sRc;
            double fpair = (1.0 - rdivsig) / rij / sRc;
#elif defined(Hertzian)
            double delta = sRc - rij;
            // lammps style hertzian
            //  double fpair = sqrt(iRc * jRc / sRc * delta) * delta / rij;
            //  ePair += 0.4 * sqrt(iRc * jRc / sRc * delta) * delta * delta;
            double fpair = sqrt(delta / sRc) * delta / sRc / sRc / rij;
#endif
            
            vScaleAdd(iforce, iforce, fpair, dRij);
            vScaleAdd(force[jatom], force[jatom], -fpair, dRij);
        }
        vCpy(force[iatom], iforce);
    }
    
    double sumForce = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        if (cinfo->isRattler[iatom])
            continue;
        sumForce += sNorm(force[iatom]);
    }
    cinfo->meanForceExRattler = sumForce / update->forceUnits /
    (double)(particle->nAtom - cinfo->nRattler);
    
    return 0;
}
int delContactInfo(Update *update) {
    contactInfo *cinfo = getContactInfo(update);
    if (!cinfo)
        return -1;
    
    cinfo->refCnt--;
    if (cinfo->refCnt <= 0) {
        safeFree(cinfo->isSoliton);
        safeFree(cinfo->isRattler);
        safeFree(cinfo->isNebrContact);
        safeFree(cinfo->nCoordNum);
        safeFree(cinfo->nCoordNumExRattler);
        safeFree(cinfo->forceExRattler);
        delToolkit(&update->toolkit, "contactInfo");
    }
    return 0;
}

//===============================================================================
writeDumpFile *getWriteDumpFile(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "writeDumpFile");
    if (whichTool < 0)
        return (writeDumpFile *)NULL;
    return (writeDumpFile *)update->toolkit.toolkit[whichTool];
}
writeDumpFile *addWriteDumpFile(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--dump");
    if (!cmd) {
        return (writeDumpFile *)NULL;
    } else if (getWriteDumpFile(update) != NULL) {
        Info("repetitive initWriteDumpFile!");
        return getWriteDumpFile(update);
    }
    writeDumpFile *dinfo = (writeDumpFile *)calloc(sizeof(writeDumpFile), 1);
    addToolkit(&update->toolkit, (void *)dinfo, NULL, "writeDumpFile");
    
    char fname[4096];
    if (cmd->cmdArgc == 1) {
        sprintf(fname, "%s/dump_%s_%s.bin", var->cwd, cmd->cmdArgv[0],
                var->sf);
    } else {
        sprintf(fname, "%s/dump_%s.bin", var->cwd, var->sf);
    }
    
    bool isFileExist = isFileExisted(fname);
    dinfo->fdump = createFileReadWrite(fname);
    if ((!isFileExist) || (isFileExist && (truncFileFlag != 2))) {
        // write header
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
        fflush(dinfo->fdump);
    } else {
        {//do mininor checking
            fseek(dinfo->fdump, 0, SEEK_SET);
            
            char str[32];
            fread(str, sizeof(char), 32, dinfo->fdump);
            if (strcmp(str, "Revised Binary File") == 0) {
                int dumpRevNum = -1;
                fread(&dumpRevNum, sizeof(int), 1, dinfo->fdump);
                if (dumpRevNum != dumpFileRevNum)
                    Abort("The dump file is for Rev. %d, code is for %d.", dumpRevNum,
                          dumpFileRevNum);
                int dim = -1;
                fread(&dim, sizeof(int), 1, dinfo->fdump);
                if (dim != box->dim)
                    Abort("The dumpfile is for d = %d, dim of initial file is %d.", dim,
                          box->dim);
            } else {
                fseek(dinfo->fdump, 0, SEEK_SET);
                if (DIM != 3)
                    Abort("The dumpfile is for d = 3, the code is for d = %d.", DIM);
            }
            
            int nAtom = -1;
            fread(&nAtom, sizeof(int), 1, dinfo->fdump);
            if (particle->nAtom != nAtom) Abort("Not Consistent!");
            
            int nAtomType = -1;
            fread(&nAtomType, sizeof(int), 1, dinfo->fdump);
            if (particle->nAtomType != nAtomType) Abort("Not Consistent!");
        }
        fseek(dinfo->fdump, 0, SEEK_END);
    }
    
    return dinfo;
}
int writeDump(Box *box, Particle *particle, Update *update) {
    writeDumpFile *wdf = getWriteDumpFile(update);
    if (wdf == NULL)
        return -1;
    
    safeFwrite(wdf->fdump, &particle->nAtom, sizeof(int), 1);
    safeFwrite(wdf->fdump, &box->boxH, sizeof(uptriMat), 1);
    safeFwrite(wdf->fdump, &particle->meanDiameter, sizeof(double), 1);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        doubleVector uxyz;
        unwrapPos(uxyz, particle->pos[idx], particle->img[idx], box->boxH);
        
        safeFwrite(wdf->fdump, uxyz, sizeof(doubleVector), 1);
    }
    return 0;
}
int delWriteDumpFile(Update *update) {
    writeDumpFile *wdf = getWriteDumpFile(update);
    if (!wdf)
        return -1;
    safeCloseFile(wdf->fdump);
    wdf->refCnt--;
    delToolkit(&update->toolkit, "writeDumpFile");
    return 0;
}

mmapBinFile *openBinFile(char *fname) {
    mmapBinFile *binFile = (mmapBinFile *)calloc(sizeof(mmapBinFile), 1);
    
    binFile->fd = open(fname, O_RDONLY);
    fstat(binFile->fd, &binFile->binStat);
    binFile->dataSection = mmap(NULL, binFile->binStat.st_size, PROT_READ,
                                MAP_PRIVATE, binFile->fd, 0);
    if (binFile->dataSection == MAP_FAILED) {
        Abort("Map Binary file:\"%s\" Failed !\n", fname);
        return NULL;
    }
    
    // header
    char *header = (char *)binFile->dataSection;
    void *data = binFile->dataSection;
    int *nElement = (int *)data;
    if (strcmp(header, "Revised Binary File") == 0) {
        int *revNum = (int *)(header + 32);
        binFile->revNum = revNum[0];
        if (binFile->revNum != 1)
            Abort("Wrong Binary File!");
        
        int *dim = (int *)(revNum + 1);
        if (dim[0] != DIM)
            Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0], DIM);
        nElement = (int *)(dim + 1);
        
        Info("revNum of %s: %d. (0: SS 3d; 1: SS; 2: HS;)", fname, binFile->revNum);
    } else if (strcmp(header, "Revised (HS) Binary File") == 0) {
        int *revNum = (int *)(header + 32);
        binFile->revNum = revNum[0];
        if (binFile->revNum != 2)
            Abort("Wrong Binary File!");
        
        int *dim = (int *)(revNum + 1);
        if (dim[0] != DIM)
            Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0], DIM);
        nElement = (int *)(dim + 1);
        
        Info("revNum of %s: %d. (0: SS 3d; 1: SS; 2: HS;)", fname, binFile->revNum);
    } else {
        if (strcmp(header, "binary") == 0 || strcmp(header, "LAMMPS compatible data file.") == 0) {
            Abort("The file \"%s\" is not dump file.", fname);
        }
        binFile->revNum = 0;
        if (DIM != 3)
            Abort("The dumpfile is for d = 3, while the code is for d = %d!", DIM);
        
        Info("We assume the format of DUMP file \"%s\": 0. (0: SS 3d; 1: SS; 2: HS;).\n\tPlease check carefully!", fname);
    }
    
    switch (binFile->revNum) {
    case 0:
        binFile->headerSize = sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    case 1:
        binFile->headerSize = 32 * sizeof(char) + 2 * sizeof(int) + sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    case 2:
        binFile->headerSize = 32 * sizeof(char) + 2 * sizeof(int) + sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    default:
        Abort("Not code!");
        break;
    }
    
    // We assume the stepSize < INT_MAX.
    switch (binFile->revNum) {
    case 0:
        binFile->stepSize = sizeof(int) + sizeof(double) * 6 + sizeof(double) + nElement[0] * sizeof(doubleVector);
        break;
    case 1:
        binFile->stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) + nElement[0] * sizeof(doubleVector);
        break;
    case 2:
        binFile->stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) * 2 + nElement[0] * sizeof(doubleVector);
        break;
    default:
        Abort("Not code!");
        break;
    }
    binFile->nStep = (int)((binFile->binStat.st_size - binFile->headerSize) / binFile->stepSize);
    
    return binFile;
}
int readSimInfo(Box *box, Particle *particle, Update *update, mmapBinFile *binFile) {
    return readDump(box, particle, update, binFile, 0);
}
int readDump(Box *box, Particle *particle, Update *update, mmapBinFile *binFile, int whichStep) {
    if (whichStep >= binFile->nStep || whichStep < 0) {
        Abort("Step %d is out of Range: [0,%d];", whichStep, binFile->nStep - 1);
    }
    
    {
        void *data = (void *)(binFile->dataSection);
        if (binFile->revNum != 0) {
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
        
        // increase step by step!
        data = (char *)binFile->dataSection + binFile->headerSize;
        for (int ith = 0; ith < whichStep; ith++) {
            data = (char *)data + binFile->stepSize;
        }
        data = (char *)data + sizeof(int);
        
        if (binFile->revNum == 0) {
            double *tmp = (double *)data;
            box->boxH[spaceIdx2voigt(0, 0)] = tmp[0];
            box->boxH[spaceIdx2voigt(1, 1)] = tmp[1];
            box->boxH[spaceIdx2voigt(2, 2)] = tmp[2];
            box->boxH[spaceIdx2voigt(1, 2)] = tmp[3];
            box->boxH[spaceIdx2voigt(0, 2)] = tmp[4];
            box->boxH[spaceIdx2voigt(0, 1)] = tmp[5];
            data = ((char *)data + 6 * sizeof(double));
        } else {
            memcpy(box->boxH, data, sizeof(uptriMat));
            data = ((char *)data + sizeof(uptriMat));
        }
        
        particle->meanDiameter = *((double *)data);
        data = ((char *)data + sizeof(double));
        if (binFile->revNum == 2) {
            // double runtimeReal = *((double *)data);
            data = ((char *)data + sizeof(double));
        }
        
        memcpy(particle->pos, data, particle->nAtom * sizeof(doubleVector));
        
        memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
        memset(particle->force, '\0', particle->nAtom * sizeof(doubleVector));
        memset(particle->img, '\0', particle->nAtom * sizeof(intVector));
    }
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    particle->isForceValid = false;
    update->Edone = update->Pdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.doSort = true;
    update->nebrList.skinSet = __minSkinSet__;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
    reInitSim(box, particle, update);
    adjustImg(box, particle);
    
    generateUnwarpPos(box, particle);
    
    return 0;
}
int closeBinFile(mmapBinFile **binFilePtr) {
    if (binFilePtr[0] == NULL)
        return 0;
    munmap(binFilePtr[0]->dataSection, binFilePtr[0]->binStat.st_size);
    close(binFilePtr[0]->fd);
    safeFree(binFilePtr[0]);
    return 0;
}

//===============================================================================
int buildFullNebrList(Box *box, Particle *particle, double rmax, fullNebrList *fNebrList) {
    adjustImg(box, particle);
    
    // bin List
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    double mindistPlane = sMinElement(distPlane);
    if (rmax < 0) {
        fNebrList->rmax = mindistPlane / 2.0 * 0.995;
    } else if (rmax > mindistPlane / 2.0 * 0.995) {
        Info("reset rmax to %g.", mindistPlane / 2.0 * 0.995);
        fNebrList->rmax = mindistPlane / 2.0 * 0.995;
    } else
        fNebrList->rmax = rmax;
    
    // build full-style NebrList
    if (fNebrList->nNebr == NULL) {
        fNebrList->nNebr = (int2 *)calloc(particle->nAtom, sizeof(int2));
    }
    if (fNebrList->binList == NULL) {
        fNebrList->binList = (int *)calloc(particle->nAtom, sizeof(int));
    }
    
    int totBin = 1;
    for (int idim = 0; idim < DIM; idim++) {
        fNebrList->nBin[idim] = (int)floor(distPlane[idim] / fNebrList->rmax);
        fNebrList->nBin[idim] = cpuMax(fNebrList->nBin[idim], 3);
        totBin *= fNebrList->nBin[idim];
    }
    if (totBin > fNebrList->allocBinHead) {
        fNebrList->allocBinHead = totBin;
        fNebrList->binHead =
        (int *)realloc(fNebrList->binHead, totBin * sizeof(int));
    }
    
    intVector nbin;
    vCpy(nbin, fNebrList->nBin);
    int *binHead = fNebrList->binHead;
    int *binList = fNebrList->binList;
    for (int ibin = 0; ibin < totBin; ibin++) {
        binHead[ibin] = -1;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        
        intVector cIdx;
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, pos);
        vShiftAll(lamda, 0.5);
        int binIdx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            cIdx[idim] = (int)floor(lamda[idim] * nbin[idim]);
            cIdx[idim] = (cIdx[idim] < 0 ? nbin[idim] - 1 : cIdx[idim]);
            cIdx[idim] = (cIdx[idim] >= nbin[idim] ? 0 : cIdx[idim]);
            binIdx = binIdx * nbin[idim] + cIdx[idim];
        }
        binList[iatom] = binHead[binIdx];
        binHead[binIdx] = iatom;
    }
    
#if (DIM == 3)
    for (int ibz = 0, cntNebr = 0; ibz < nbin[2]; ibz++) {
        for (int iby = 0; iby < nbin[1]; iby++) {
            for (int ibx = 0; ibx < nbin[0]; ibx++) {
                int ibin = ibz * nbin[1] * nbin[0] + iby * nbin[0] + ibx;
                for (int iatom = binHead[ibin]; iatom >= 0; iatom = binList[iatom]) {
                    fNebrList->nNebr[iatom].first = cntNebr;
                    for (int dz = -1; dz <= 1; dz++) {
                        int jbz = (ibz + dz + nbin[2]) % nbin[2];
                        for (int dy = -1; dy <= 1; dy++) {
                            int jby = (iby + dy + nbin[1]) % nbin[1];
                            for (int dx = -1; dx <= 1; dx++) {
                                int jbx = (ibx + dx + nbin[0]) % nbin[0];
                                int jbin = jbz * nbin[1] * nbin[0] + jby * nbin[0] + jbx;
                                for (int jatom = binHead[jbin]; jatom >= 0;
                                     jatom = binList[jatom]) {
                                    if (iatom == jatom)
                                        continue;
                                    
                                    doubleVector drij;
                                    vSub(drij, particle->pos[iatom], particle->pos[jatom]);
                                    PBC(drij, box);
                                    double dr = sNorm(drij);
                                    if (dr > fNebrList->rmax)
                                        continue;
                                    
                                    if (cntNebr == fNebrList->maxAllocNebr) {
                                        fNebrList->maxAllocNebr += 10240;
                                        fNebrList->list = (int *)realloc(
                                                                         fNebrList->list, fNebrList->maxAllocNebr * sizeof(int));
                                        if (fNebrList->list == NULL)
                                            Abort("realloc nebrList failed!");
                                    }
                                    fNebrList->list[cntNebr++] = jatom;
                                }
                            }
                        }
                    }
                    fNebrList->nNebr[iatom].second = cntNebr;
                }
            }
        }
    }
#elif (DIM == 2)
    for (int iby = 0, cntNebr = 0; iby < nbin[1]; iby++) {
        for (int ibx = 0; ibx < nbin[0]; ibx++) {
            int ibin = iby * nbin[0] + ibx;
            for (int iatom = binHead[ibin]; iatom >= 0; iatom = binList[iatom]) {
                fNebrList->nNebr[iatom].first = cntNebr;
                for (int dy = -1; dy <= 1; dy++) {
                    int jby = (iby + dy + nbin[1]) % nbin[1];
                    for (int dx = -1; dx <= 1; dx++) {
                        int jbx = (ibx + dx + nbin[0]) % nbin[0];
                        int jbin = jby * nbin[0] + jbx;
                        for (int jatom = binHead[jbin]; jatom >= 0;
                             jatom = binList[jatom]) {
                            if (iatom == jatom)
                                continue;
                            doubleVector drij;
                            vSub(drij, particle->pos[iatom], particle->pos[jatom]);
                            PBC(drij, box);
                            double dr = sNorm(drij);
                            if (dr > fNebrList->rmax)
                                continue;
                            
                            if (cntNebr == fNebrList->maxAllocNebr) {
                                fNebrList->maxAllocNebr += 10240;
                                fNebrList->list = (int *)realloc(
                                                                 fNebrList->list, fNebrList->maxAllocNebr * sizeof(int));
                                if (fNebrList->list == NULL)
                                    Abort("realloc nebrList failed!");
                            }
                            fNebrList->list[cntNebr++] = jatom;
                        }
                    }
                }
                
                fNebrList->nNebr[iatom].second = cntNebr;
            }
        }
    }
#else
    
    int nAdjBin = pow(3, DIM);
    intVector *adjBinVec = (intVector *)calloc(nAdjBin, sizeof(intVector));
    for (int idx = 0; idx < nAdjBin; idx++) {
        intVector delta;
        for (int idim = 0, itmp = idx; idim < DIM; idim++) {
            delta[idim] = itmp % 3 - 1;
            itmp = itmp / 3;
        }
        vCpy(adjBinVec[idx], delta);
    }
    
    for (int gbin = 0, cntNebr = 0; gbin < totBin; gbin++) {
        intVector gbinVec;
        for (int idim = 0, itmp = gbin; idim < DIM; idim++) {
            gbinVec[idim] = itmp % nbin[idim];
            itmp = itmp / nbin[idim];
        }
        
        for (int iatom = binHead[gbin]; iatom >= 0; iatom = binList[iatom]) {
            fNebrList->nNebr[iatom].first = cntNebr;
            
            // loop over all atoms in adjacent bins
            for (int adj = 0; adj < nAdjBin; adj++) {
                int jbin = 0;
                for (int idim = DIM - 1; idim >= 0; idim--) {
                    jbin =
                    jbin * nbin[idim] +
                    (gbinVec[idim] + adjBinVec[adj][idim] + nbin[idim]) % nbin[idim];
                }
                
                for (int jatom = binHead[jbin]; jatom >= 0; jatom = binList[jatom]) {
                    if (iatom == jatom)
                        continue;
                    
                    doubleVector dRij;
                    vSub(dRij, particle->pos[iatom], particle->pos[jatom]);
                    PBC(dRij, box);
                    
                    double dr = sNorm(dRij);
                    if (dr > fNebrList->rmax)
                        continue;
                    
                    if (cntNebr == fNebrList->maxAllocNebr) {
                        fNebrList->maxAllocNebr += 10240;
                        fNebrList->list = (int *)realloc(
                                                         fNebrList->list, fNebrList->maxAllocNebr * sizeof(int));
                        if (fNebrList->list == NULL)
                            Abort("realloc nebrList failed!");
                    }
                    fNebrList->list[cntNebr++] = jatom;
                }
            }
            
            fNebrList->nNebr[iatom].second = cntNebr;
        }
    }
    
    safeFree(adjBinVec);
#endif
    
    return 0;
}
int delFullNebrList(fullNebrList *fNebrList) {
    safeFree(fNebrList->list);
    safeFree(fNebrList->nNebr);
    safeFree(fNebrList->binHead);
    safeFree(fNebrList->binList);
    return 0;
}

void generateUnwarpPos(Box *box, Particle *particle) {
    if (particle->pos == NULL)
        return;
    
    if (particle->uPos == NULL) {
        particle->uPos =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        unwrapPos(particle->uPos[iatom], particle->pos[iatom], particle->img[iatom],
                  box->boxH);
    }
}

//===============================================================================
void purgeSim(Box *box, Particle *particle, Update *update) {
    safeFprintf(stdout, "Purging ...\n");
    
    memset(box, '\0', sizeof(Box));
    
    safeFree(particle->pos);
    safeFree(particle->veloc);
    safeFree(particle->img);
    safeFree(particle->force);
    safeFree(particle->diameterScale);
    safeFree(particle->type);
    safeFree(particle->massPerType);
    safeFree(particle->mass);
    safeFree(particle->id2tag);
    safeFree(particle->tag2id);
    memset(particle, '\0', sizeof(Particle));
    
    safeFree(update->nebrList.imageParent);
    safeFree(update->nebrList.imgageSource);
    safeFree(update->nebrList.isSourceBin);
    safeFree(update->nebrList.stencilSource);
    safeFree(update->nebrList.xyzHold);
    safeFree(update->nebrList.binHead);
    safeFree(update->nebrList.binList);
    safeFree(update->nebrList.adjBinList);
    safeFree(update->nebrList.list);
    safeFree(update->nebrList.nNebr);
    
    safeFree(update->nebrList.binHead4sort);
    safeFree(update->nebrList.binList4sort);
    safeFree(update->nebrList.oid2nid);
    safeFree(update->nebrList.buffer);
    
    if (update->toolkit.nToolkit != 0) {
        Info("Warning!");
    }
    while (update->toolkit.nToolkit != 0) {
        safeFprintf(stdout, "\t\"%s\" is not released!\n",
                    update->toolkit.toolkitName[0]);
        delToolkit(&update->toolkit, update->toolkit.toolkitName[0]);
    }
    safeFree(update->toolkit.toolkit);
    safeFree(update->toolkit.toolkitName);
    safeFree(update->toolkit.funcPtrWriteConf);
    memset(update, '\0', sizeof(Update));
    
    safeFprintf(stdout, "Done!\n");
}

#if (DIM == 2 || DIM == 3)
void writeLammpsTopo(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--topo");
    if (!cmd)
        return;
    if (cmd->cmdArgc != 1)
        Abort("./app --topo topo.data");
    
    FILE *fout = createFileReadWrite(cmd->cmdArgv[0]);
    
    safeFprintf(fout,
                "LAMMPS compatible data file. atom_style: sphere (id type "
                "diameter density x y z ix iy iz).\n\n");
    safeFprintf(fout, "\t%d atoms\n", particle->nAtom);
    safeFprintf(fout, "\t%d atom types\n", particle->nAtomType);
    
    double distanceUnits = update->distanceUnits;
#if (DIM == 3)
    safeFprintf(fout,
                "\n\t%g %g xlo xhi\n\t%g %g ylo yhi\n\t%g %g "
                "zlo zhi\n\t%g %g %g xy xz yz\n",
                box->cornerLo[0] / distanceUnits,
                box->cornerLo[0] / distanceUnits +
                box->boxH[spaceIdx2voigt(0, 0)] / distanceUnits,
                box->cornerLo[1] / distanceUnits,
                box->cornerLo[1] / distanceUnits +
                box->boxH[spaceIdx2voigt(1, 1)] / distanceUnits,
                box->cornerLo[2] / distanceUnits,
                box->cornerLo[2] / distanceUnits +
                box->boxH[spaceIdx2voigt(2, 2)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 1)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 2)] / distanceUnits,
                box->boxH[spaceIdx2voigt(1, 2)] / distanceUnits);
#elif (DIM == 2)
    safeFprintf(fout,
                "\n\t%g %g xlo xhi\n\t%g %g ylo yhi\n\t0 0 "
                "zlo zhi\n\t%g 0 0 xy xz yz\n",
                box->cornerLo[0] / distanceUnits,
                box->cornerLo[0] / distanceUnits +
                box->boxH[spaceIdx2voigt(0, 0)] / distanceUnits,
                box->cornerLo[1] / distanceUnits,
                box->cornerLo[1] / distanceUnits +
                box->boxH[spaceIdx2voigt(1, 1)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 1)] / distanceUnits);
#endif
    
    safeFprintf(fout, "\nAtoms\n\n");
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        
        int type = particle->type[idx];
        double diameter =
        particle->diameterScale[idx] * particle->meanDiameter / distanceUnits;
        doubleVecPtr posPtr = particle->pos[idx];
        intVecPtr imgPtr = particle->img[idx];
        
#if (DIM == 3)
        safeFprintf(fout, "%d %d %g 1 %g %g %g %d %d %d\n", iatom + 1, type + 1,
                    diameter, posPtr[0] / distanceUnits, posPtr[1] / distanceUnits,
                    posPtr[2] / distanceUnits, imgPtr[0], imgPtr[1], imgPtr[2]);
#elif (DIM == 2)
        safeFprintf(fout, "%d %d %g 1 %g %g 0 %d %d 0\n", iatom + 1, type + 1,
                    diameter, posPtr[0] / distanceUnits, posPtr[1] / distanceUnits,
                    imgPtr[0], imgPtr[1]);
#endif
    }
    
    safeCloseFile(fout);
}
void writeLammpsTraj(Box *box, Particle *particle, Update *update, Variable *var, int globalTimeStep) {
    cmdArg *cmd = findVariable(var, "--traj");
    if (!cmd)
        return;
    if (cmd->cmdArgc != 1)
        Abort("./app --traj lmpTraj.dat");
    
    FILE *fout = NULL;
    if (isFileExisted(cmd->cmdArgv[0])) {
        fout = openExistFileReadWrite(cmd->cmdArgv[0]);
        fseek(fout, 0, SEEK_END);
    } else {
        fout = createFileReadWrite(cmd->cmdArgv[0]);
    }
    
    fprintf(fout, "ITEM: TIMESTEP\n");
    fprintf(fout, "%d\n", globalTimeStep);
    
    double distanceUnits = update->distanceUnits;
    fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fout, "%d\n", particle->nAtom);
    fprintf(fout, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
    double xlo, xhi, ylo, yhi, zlo, zhi, yz, xz, xy;
    box2BoxBounds(box, xlo, xhi, ylo, yhi, zlo, zhi, yz, xz, xy);
    fprintf(fout, "%g %g %g\n", xlo / distanceUnits, xhi / distanceUnits,
            xy / distanceUnits);
    fprintf(fout, "%g %g %g\n", ylo / distanceUnits, yhi / distanceUnits,
            xz / distanceUnits);
    fprintf(fout, "%g %g %g\n", zlo / distanceUnits, zhi / distanceUnits,
            yz / distanceUnits);
    fprintf(fout, "ITEM: ATOMS xu yu zu c_LEx c_LEy c_LEz\n");
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        doubleVector uxyz;
        unwrapPos(uxyz, particle->pos[idx], particle->img[idx], box->boxH);
        doubleVector uxyz_LE;
        vCpy(uxyz_LE, uxyz);
        triBox2LeesEdwards(uxyz_LE, box);
        
#if (DIM == 3)
        fprintf(fout, "%g %g %g %g %g %g\n", uxyz[0] / distanceUnits,
                uxyz[1] / distanceUnits, uxyz[2] / distanceUnits,
                uxyz_LE[0] / distanceUnits, uxyz_LE[1] / distanceUnits,
                uxyz_LE[2] / distanceUnits);
#elif (DIM == 2)
        fprintf(fout, "%g %g 0 %g %g 0\n", uxyz[0] / distanceUnits,
                uxyz[1] / distanceUnits, uxyz_LE[0] / distanceUnits,
                uxyz_LE[1] / distanceUnits);
#endif
    }
    safeCloseFile(fout);
}
#endif
