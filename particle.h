

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    float x; 
    float oldX;   
    float y; 
    float oldY;   
    float z; 
    float oldZ;   
    float p1;    //momentum  
    float p2;
    float p3;
    float Ex;    
    float Ey;    
    float Ez;    
    float Bx;    
    float By;    
    float Bz;    
    float index; 
    struct _ptclList *next;
} ptclList;

