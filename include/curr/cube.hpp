# 1 "cube.h"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "cube.h"



# 1 "wavelet.h" 1



# 1 "configs.h" 1





# 1 "const.h" 1
# 7 "configs.h" 2

class Configs
{
 char *group[128];
 char *name[128];
 char *value[128];
 int elements;



public:
 Configs(void);
 virtual ~Configs(void);

 int ReadConfigFile(char *file);

 char* GetParam(char* group, char* param_name);
 float GetParamDouble(char* group, char* param_name);
 int GetParamInt(char* group, char* param_name);
};
# 5 "wavelet.h" 2

class Wavelet
{
public:

 int* points;


 float* values;


 int wavelet_used_values;


 int max_values;

 Wavelet(Configs* config);
 virtual ~Wavelet(void);

 int ReadFromFile(char* file);
 double PointValue(int point);
 double IndexValue(int index);
};
# 5 "cube.h" 2






class Cube
{
public:
 long x_max;
 long y_max;
 long z_max;

 int size;

 int type;

 char *name;

 float* values;

 double* sum_x;
 double* sum_x2;
 double* r;
 double GlobalSum_x;
 double GlobalSum_x2;

 int layers_num;
 double* layersCorr;
 double corrGlobal;
 double corrAvg;

 Cube(long x, long y, long z, int layers_num, int type);
 virtual ~Cube(void);

 int ReadFromFile(char* file);
 int ReadFromFileFast(char* file);
 int SaveToFile(char* file, Cube *mask, double null_value, int root);

 void Init();
 void Init(double);

 void SetCell(int x, int y, int z, double value);
 double GetCell(int x, int y, int z);

 void DeleteValues();

 Cube *MakeReflectionCoefsCube();
 Cube *MakeSyntheticCube(Wavelet *);
 void LoadNewLayers(unsigned short int );
 void GenerateLayersCorr(struct layers_t*, Cube*);
 void CalculateCorrData(struct layers_t*);
 void CalculateStaticCorrSums();
 double GetLayerCorrAverage();
 double GetGlobalCorrAverage();
 double* copyLayersCorr();
 char* CorrAverageDump();
 void DumpCRs();
 long Coords(int x, int y, int z);
};
