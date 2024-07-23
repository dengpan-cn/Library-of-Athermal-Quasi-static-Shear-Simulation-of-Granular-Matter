#include "VectorMath_ndim.h"

#ifdef __cplusplus
#error "This file must NOT be compiled by C++ compliler!"
#endif

FILE *logFile = NULL;
int truncFileFlag = 0;         // 0:normal;1:truncate;2:restart
int screenOutputReadConf = 0;  // 0: not; >= 10000: internal state;
char *MSG = NULL;

double rndUniform(void) { // random number belong to  U[0,1]
    static int setSeed = 1;

    static double *u;
    static int i97, j97;
    static double c, cd, cm;
    if (setSeed == 1) {
        srand((unsigned int)time((time_t *)NULL) + (unsigned int)getpid());
        for (int ith = 0; ith < (unsigned int)getpid() % 10000; ith++) {
            rand();
        }

        int seed = rand() % 900000000;
        while (seed <= 0 || seed > 900000000) {
            seed = rand() % 900000000;
        }
        u = (double *)calloc(97 + 1, sizeof(double));
        int ij = (seed - 1) / 30082;
        int kl = (seed - 1) - 30082 * ij;
        int i = (ij / 177) % 177 + 2;
        int j = ij % 177 + 2;
        int k = (kl / 169) % 178 + 1;
        int l = kl % 169;
        for (int ii = 1; ii <= 97; ii++) {
            double s = 0.0;
            double t = 0.5;
            for (int jj = 1; jj <= 24; jj++) {
                int m = ((i * j) % 179) * k % 179;
                i = j;
                j = k;
                k = m;
                l = (53 * l + 1) % 169;
                if ((l * m) % 64 >= 32) s = s + t;
                t = 0.5 * t;
            }
            u[ii] = s;
        }

        c = 362436.0 / 16777216.0;
        cd = 7654321.0 / 16777216.0;
        cm = 16777213.0 / 16777216.0;
        i97 = 97;
        j97 = 33;

        setSeed = 0;

        rndUniform();
    }

    // u,i97,j97,c,cd,cm
    double uni = u[i97] - u[j97];
    if (uni < 0.0) uni += 1.0;
    u[i97] = uni;
    i97--;
    if (i97 == 0) i97 = 97;
    j97--;
    if (j97 == 0) j97 = 97;
    c -= cd;
    if (c < 0.0) c += cm;
    uni -= c;
    if (uni < 0.0) uni += 1.0;

    return uni;
}
double rndStdNorm(void) {  // random number belong to N(0,1)
    //====
    static bool genStdNorm = true;
    static double otherRndStd = 0;

    double rndStdNormNum = 0.0;
    if (genStdNorm) {  // generate two rndStdNorm number.
        double v1, v2, s;
        do {
            v1 = 2.0 * rndUniform() - 1.0;
            v2 = 2.0 * rndUniform() - 1.0;
            s = v1 * v1 + v2 * v2;
        } while ((s >= 1.0) || (s == 0.0));
        double fac = sqrt(-2.0 * log(s) / s);

        otherRndStd = v1 * fac;    // save another rndStdNorm number
        rndStdNormNum = v2 * fac;  // return this
    } else {
        rndStdNormNum = otherRndStd;
    }
    genStdNorm = !(genStdNorm);

    return rndStdNormNum;
}

bool isEmpty(char *str, int maxlen) {
    if (maxlen <= 0) Abort("Fatal Error!");
    if (strlen(str) >= maxlen) Abort("Fatal Error");
    if (strlen(str) == 0) return true;

    char *start = str, *stop = str + strlen(str) - 1;
    while (isspace(*stop) && stop >= start) stop--;
    while (isspace(*start) && start <= stop) start++;

    if (start > stop) return true;
    return false;
}

cmdArg *findVariable(Variable *var, char *name) {
    char *ptr = name;
    if (strlen(ptr) == 0) Abort("Not legal variable name!");
    if (ptr[0] == '-') {
        if (strlen(ptr) < 3 || (ptr[2] == '-')) Abort("%s is not legal variable name!", ptr);
        ptr += 2;
    }

    for (int ith = 0; ith < var->nVar; ith++) {
        if (!strcmp(var->cmd[ith].cmdType + 2, ptr)) return (var->cmd + ith);
    }
    return NULL;
}
void addVariable(Variable *var, char *inputCmd) {
    if (var == NULL) Abort("var is NULL");
    if (strlen(inputCmd) < 3) Abort("No cmd!");

    char *tmp = inputCmd;
    while (tmp[0] != '\0' && isspace(tmp[0])) {
        tmp++;
    }
    if (strlen(tmp) <= 2) Abort("Unrecognized cmd %s!", inputCmd);
    if (strncmp(tmp, "--", 2)) Abort("Unrecognized cmd %s!", inputCmd);
    if (isspace(tmp[2]) || (tmp[2] == '-')) Abort("Unrecognized cmd %s!", inputCmd);
    char *str = (char *)calloc(strlen(tmp) + 5, sizeof(char));
    sprintf(str, "%s", tmp);
    tmp += 2;

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
    if (findVariable(var, str))
        Abort("Repetitive cmd! Please delete one before add another one!");
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
    cmdArg *cmd = findVariable(var, delCmd);
    if (cmd == NULL) return;

    safeFree(cmd->cmdArgv);
    safeFree(cmd->cmdType);
    cmd->cmdArgc = 0;

    int whichVar = (int)(cmd - var->cmd);
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

int findToolkit(Toolkit *toolkit, char *name) {
    for (int ith = 0; ith < toolkit->nToolkit; ith++) {
        if (!strcmp(name, toolkit->toolkitName[ith])) return ith;
    }
    return -1;
}
void addToolkit(Toolkit *toolkit, void *tool, void *funcPtr, char *name) {
    if (strlen(name) >= 32) Abort("\"%s\" is too long!", name);
    if (findToolkit(toolkit, name) >= 0)
        Abort("Repetitive tool! Please delete one before add another one!");

    if (toolkit->maxToolkit <= toolkit->nToolkit) {
        toolkit->maxToolkit += 8;
        toolkit->toolkit = (void **)realloc(toolkit->toolkit,
                                            toolkit->maxToolkit * sizeof(void *));
        toolkit->funcPtrWriteConf = (void **)realloc(toolkit->funcPtrWriteConf,
                                                     toolkit->maxToolkit * sizeof(void *));
        toolkit->toolkitName = (char **)realloc(
            toolkit->toolkitName, toolkit->maxToolkit * sizeof(char *));
    }
    toolkit->toolkit[toolkit->nToolkit] = tool;
    toolkit->funcPtrWriteConf[toolkit->nToolkit] = funcPtr;
    toolkit->toolkitName[toolkit->nToolkit] =
        (char *)calloc(strlen(name) + 5, sizeof(char));
    sprintf(toolkit->toolkitName[toolkit->nToolkit], "%s", name);
    toolkit->nToolkit++;
}
void delToolkit(Toolkit *toolkit, char *name) {
    if (findToolkit(toolkit, name) < 0) return;
    int which = findToolkit(toolkit, name);
    safeFree(toolkit->toolkitName[which]);
    safeFree(toolkit->toolkit[which]);

    if (which != toolkit->nToolkit - 1) {
        toolkit->toolkitName[which] = toolkit->toolkitName[toolkit->nToolkit - 1];
        toolkit->toolkit[which] = toolkit->toolkit[toolkit->nToolkit - 1];
        toolkit->funcPtrWriteConf[which] = toolkit->funcPtrWriteConf[toolkit->nToolkit - 1];
    }

    toolkit->nToolkit--;
}

#ifdef __Linux__
// this function is used to dump all source files, compiled by the following Make cmd.
/*
$(objSourceFile): $(SourceFileName)
    echo "//=====" > build/sourceFile
    for fname in $(SourceFileName); do echo //$$fname >> build/sourceFile; cat $$fname >> build/sourceFile; done
    objcopy -I binary -B i386 -O elf64-x86-64 build/sourceFile $(objSourceFile)
    rm build/sourceFile
*/
extern void *_binary_build_sourceFile_start;
void dumpSourceFile(void) {
    printf("%s", (char *)&_binary_build_sourceFile_start);
}
#endif

volatile int __nCatchSignal__ = 0;
void incCatchCount(int sig) {
    __nCatchSignal__++;
    fprintf(stderr, "Number of catched signal: %d.\n", __nCatchSignal__);
    fflush(stderr);
}
int setSignalHandler(void) {
    if(screenOutputReadConf) safeFprintf(stderr, "SetSignalHandler ...\n");
    
    sigset_t sigSet;
    sigemptyset(&sigSet);
    sigprocmask(SIG_SETMASK, &sigSet, NULL);  // deliver all signals to process.

    struct sigaction action;
    action.sa_handler = incCatchCount;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_RESTART;

    // catching the following signals
    if (sigaction(SIGHUP, &action, NULL) == -1) return -1;   // 1
    if (sigaction(SIGINT, &action, NULL) == -1) return -1;   // 2
    if (sigaction(SIGQUIT, &action, NULL) == -1) return -1;  // 3

    return 0;
}

#define heapSortSwap(arr, a, b) \
    {                           \
        int tmp = arr[b];       \
        arr[b] = arr[a];        \
        arr[a] = tmp;           \
    }
void max_heapify(int *arr, int start, int end) {
    // 建立父节点指标和子节点指标
    int dad = start;
    int son = dad * 2 + 1;
    while (son <= end) {  // 若子节点指标在范围内才做比较
        if (son + 1 <= end && arr[son] < arr[son + 1]) {
            // 先比较两个子节点大小，选择最大的
            son++;
        }
        if (arr[dad] > arr[son]) {  // 如果父节点大於子节点代表调整完毕，直接跳出函数
            return;
        } else {  // 否则交换父子内容再继续子节点和孙节点比较
            heapSortSwap(arr, dad, son);
            dad = son;
            son = dad * 2 + 1;
        }
    }
}
void heapSortInc(int *arr, int len) {
    if (len <= 1) return;

    // 将数组初始化为大顶堆
    for (int ith = len / 2 - 1; ith >= 0; ith--) max_heapify(arr, ith, len - 1);
    // 将第一个元素和已排好元素前一位做交换，再重新调整，直到排序完毕
    for (int ith = len - 1; ith > 0; ith--) {
        heapSortSwap(arr, 0, ith);
        max_heapify(arr, 0, ith - 1);
    }
}
#undef heapSortSwap

int compareDoubleVector(const void *a, const void*b){
    doubleVector* aPtr = (doubleVector*)a;
    doubleVector* bPtr = (doubleVector*)b;
    
    if(aPtr[0][0] >= bPtr[0][0]) return +1;
    else return -1;
}
int sortDoubleVector0Inc(doubleVector *base, int nVec){
    qsort(base,nVec,sizeof(doubleVector),compareDoubleVector);
    return 1;
}

void exchange_doubleVector(doubleVector *xyz, doubleVector *buffer, int *oid2nid, int nAtom) {
    for (int oid = 0; oid < nAtom; oid++) {
        int nid = oid2nid[oid];
        vCpy(buffer[nid], xyz[oid]);
    }
}
void exchange_intVector(intVector *img, intVector *buffer, int *oid2nid, int nAtom) {
    for (int oid = 0; oid < nAtom; oid++) {
        int nid = oid2nid[oid];
        vCpy(buffer[nid], img[oid]);
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


#if (DIM == 3)
void toMandelTensor(double *cVoigt, double *cMandel) {
    for (int voigt1 = 0; voigt1 < voigtDIM; voigt1++) {
        for (int voigt2 = voigt1; voigt2 < voigtDIM; voigt2++) {
            cMandel[voigtSymRowCol2lidx(voigt1, voigt2)] =
                cVoigt[voigtSymRowCol2lidx(voigt1, voigt2)];
            if (voigt1 >= (voigtDIM / 2) && voigt2 > (voigtDIM / 2)) {
                cMandel[voigtSymRowCol2lidx(voigt1, voigt2)] *= 2.0;
            } else if (voigt1 < (voigtDIM / 2) || voigt2 >= (voigtDIM / 2)) {
                cMandel[voigtSymRowCol2lidx(voigt1, voigt2)] *= sqrt(2.0);
            }
        }
    }
}
#elif (DIM == 2)
void toMandelTensor(double *cVoigt, double *cMandel) {
    // just copying
    for (int voigt1 = 0; voigt1 < voigtDIM; voigt1++) {
        for (int voigt2 = voigt1; voigt2 < voigtDIM; voigt2++) {
            cMandel[voigtSymRowCol2lidx(voigt1, voigt2)] =
                cVoigt[voigtSymRowCol2lidx(voigt1, voigt2)];
        }
    }
}
#endif
