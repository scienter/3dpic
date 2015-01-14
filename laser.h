typedef struct _LaserList  {
   int polarity;
   float lambda;
   float omega;
   float amplitude;   //unit is a0.
   float rU;
   float rD;
   float flat;
   int loadPointX; 
   int loadPointY; 

   float rayleighLength;
   float beamWaist;
   float focus;
   int direction;

   struct _LaserList *next;
} LaserList;
