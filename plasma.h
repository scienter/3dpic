#define Electron 	1
#define HPlus0 	 	100
#define HPlus1 	 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0          600
#define CPlus1          601
#define CPlus2          602
#define CPlus3          603
#define CPlus4          604
#define CPlus5          605
#define CPlus6          606
#define userDefined   	9999999

#define Polygon    1
#define Circle    2
#define Exp    3

typedef struct _LoadList  {
   int type;
   int species;
   float superP;
   float density;
   float numberInCell;
   float criticalDensity;
   float index;  
   float num;      //exceeded number of particle which is less than 1
   int xnodes;     //longitudinal point number
   float *xn;      //longitudinal density (1 is P->density)
   float *xpoint;    //longitudinal point
   int ynodes;     //transverse point number
   float *yn;      //transverse density (1 is P->density)
   float *ypoint;    //transverse point
   int znodes;     //transverse point number
   float *zn;      //transverse density (1 is P->density)
   float *zpoint;    //transverse point
   float cx;       //circular plasma's center X
   float cy;       //circular plasma's center Y


   int pointPosition;   
   float p1;
   float p2;
   float p3;

   float beta;
   float gamma;
   float mass;
   int charge;
   
   float temperature;
   int withNextSpcs;
   int withPrevSpcs;

   struct _LoadList *next;
} LoadList;
