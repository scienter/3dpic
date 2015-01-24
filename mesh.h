#include "particle.h"
#include "laser.h"
#define FIRST 1
#define SECOND 2
#define THIRD 3


typedef struct _Domain 
{
   int dimension;

   int fieldType;
   int currentType;
   int interpolationType;

   int maxStep;
   int saveStep;
   int saveStart;
   int dumpStart;
   int dumpSave;
   int fieldSave;
   int particleSave;
   int rhoSave;

   int nx;             //Total domain
   int ny;             //Total domain
   int nz;             //Total domain
   int nxSub;          //Each core has sub domain        
   int nySub;          //Each core has sub domain        
   int nzSub;          //Each core has sub domain        
   int istart;
   int iend;
   int jstart;
   int jend;
   int kstart;
   int kend;
   int minXSub;        //Each core has start mesh point in total domain
   int maxXSub;
   int minYSub;        //Each core has start mesh point in total domain
   int maxYSub;
   int minZSub;        //Each core has start mesh point in total domain
   int maxZSub;
   int numberInCell;
   int moving;         //Moving domain option. 1:on
         
   float lambda;
   float omega;
   float divisionLambda;
   float dt;
   float dx;
   float dy;
   float dz;

   //MPI parameter
   int L;
   int M;
   int N;
   int nextCore;
   int prevCore;
   
   double *YplusJ;
   double *YminusJ;
   double *ZplusJ;
   double *ZminusJ;
   float *minusYC;
   float *plusYC;
   float *minusZC;
   float *plusZC;
   float *minusY;
   float *plusY;
   float *minusZ;
   float *plusZ;

   float ***Ex;    
   float ***Bx;    
   float ***Pr;    
   float ***Pl;    
   float ***Sr;    
   float ***Sl;    
   float ***ExC;    
   float ***BxC;    
   float ***PrC;    
   float ***PlC;    
   float ***SrC;    
   float ***SlC;    
   float ***Jx;    
   float ***Jy;    
   float ***Jz;    
   float ***JxOld;    
   float ***JyOld;    
   float ***JzOld;    

   //sharing mesh
   float *Ufield;
   float *Dfield;
   float *Rfield;
   float *Lfield;
   int numShareR;
   int numShareL;
   int numShareU;
   int numShareD;
   
   struct _UpPML **UpPML;    
   struct _DnPML **DnPML;    
   struct _Particle ***particle;    
   struct _Boost **boost;    

   //Plasma load
   struct _LoadList *loadList;
   int nSpecies;

   //Plasma load
   struct _LaserList *laserList;
   int nLaser;

   //Boost
   int boostOn;
   float gamma;
   float beta;
   int minT;	//boost frame's step
   int maxT;	//boost frame's step
   int boostSaveStep;	//lab frame's step
   
   //Probe
   int probeNum;
   int *probeX;
   int *probeY;
   int *probeZ;
   struct _Probe **probe;

   //PML
   int pmlOn;
   int pmlCell;   
}  Domain; 

typedef struct _Boost
{
   float x;
   float y;
   float E1;
   float B1;
   float Pr;
   float Pl;
   float Sr;
   float Sl;   
}  Boost;

typedef struct _UpPML 
{
   float Dx;
   float Dy;
   float Dz;
   float Ex;
   float Ey;
   float Ez;
   float Hx;
   float Hy;
   float Hz;
   float Bx;
   float By;
   float Bz;
}  UpPML;

typedef struct _DnPML 
{
   float Bzx;
   float Bzy;
}  DnPML;

typedef struct _Particle 
{
   float rho;
   // Particle List Header
   ptclHead **head;            
}  Particle;

typedef struct _External 
{
   float Ex;
   float Ey;
   float Ez;
   float Bx;
   float By;
   float Bz;
}  External;

typedef struct _Probe
{
   float E1;
   float Pr;
   float Pl;
   float B1;
   float Sr;
   float Sl;
}  Probe;
