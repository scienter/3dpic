#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void rearrangeParticles3D(Domain *D)
{
    Particle ***particle;
    particle=D->particle;

    int i,j,k,s,intX=0,intY=0,intZ=0,cnt,deleteFlag=0;
    int istart,iend,jstart,jend,kstart,kend;
    float x,y,z;
    ptclList *p,*New,*prev,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    if(D->fieldType==1)
    {
      for(i=istart; i<=iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
            for(s=0; s<D->nSpecies; s++)
            {
              cnt=1;
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                if(cnt==1)
                  prev=p;
                deleteFlag=0;
              
                x=p->x;   y=p->y;   z=p->z;
                if(x>=1.0)  {
                  intX=(int)x;
                  x-=intX;
                  deleteFlag=1;
                }
                else if(x<0) {              
                  intX=(int)(x-1);
                  x-=intX;
                  deleteFlag=1;
                } 
                else   intX=0;
                if(y>=1.0)  {
                   intY=(int)y;
                   y-=intY;
                   deleteFlag=1;
                }
                else if(y<0) {
                  intY=(int)(y-1);
                  y-=intY;
                  deleteFlag=1;
                } 
                else   intY=0;       
                if(z>=1.0)  {
                   intZ=(int)z;
                   z-=intZ;
                   deleteFlag=1;
                }
                else if(z<0) {
                  intZ=(int)(z-1);
                  z-=intZ;
                  deleteFlag=1;
                } 
                else   intZ=0;       
                if(deleteFlag==1)
                {
                  if(cnt==1)
                  {
                    particle[i][j][k].head[s]->pt = p->next;
                    p->x=x;   p->y=y;   p->z=z;   
                    p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                    particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                    p=particle[i][j][k].head[s]->pt;
                    cnt=1;
                  }
                  else
                  {
                    prev->next = p->next;
                    p->x=x;   p->y=y;   p->z=z;  
                    p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                    particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                    p=prev->next;
                  }
                }		//End of if(deleteFlag==1)
                else
                {
                  prev=p;
                  p=p->next;
                  cnt++;
                }              
              }	//End if while(p)
            }			//End of for(s)
    }          	//End of fieldType=1
}

