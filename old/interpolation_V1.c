#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
/*
void interpolation2D_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,ii,jj,i1,j1,istart,iend,jstart,jend,s;
   float x,y,Pr,Pl,Sr,Sl,E1,B1;
   float totalPr,totalPl,totalSr,totalSl,totalE1,totalB1;
   float Wx[3],Wy[3];
   float extPr,extPl,extSr,extSl,extE1,extB1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;
   extB1=Ext->B1;
   extPr=Ext->Pr;
   extPl=Ext->Pl;
   extSr=Ext->Sr;
   extSl=Ext->Sl;
 
   if(D->fieldType==1)
   {   
     FieldDSX **field;
     field=D->fieldDSX;

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
       {
         for(s=0; s<D->nSpecies; s++)
         {
           p=particle[i][j].head[s]->pt;
           while(p)
           {
             x=p->x;  y=p->y;
             Wx[0]=0.5*(1-x)*(1-x);
             Wx[1]=0.75-(0.5-x)*(0.5-x);
             Wx[2]=0.5*x*x;
             Wy[0]=0.5*(1-y)*(1-y);
             Wy[1]=0.75-(0.5-y)*(0.5-y);
             Wy[2]=0.5*y*y;

             totalPr=totalPl=totalSr=totalSl=0;
             for(jj=0; jj<3; jj++)
             {
               Pr=Pl=Sr=Sl=0;
               for(ii=0; ii<3; ii++)
               {
                 Pr+=field[i-1+ii][j-1+jj].Pr*Wx[ii];
                 Pl+=field[i-1+ii][j-1+jj].Pl*Wx[ii];
                 Sr+=field[i-1+ii][j-1+jj].Sr*Wx[ii];
                 Sl+=field[i-1+ii][j-1+jj].Sl*Wx[ii];
               }
               totalPr+=Pr*Wy[jj];
               totalPl+=Pl*Wy[jj];
               totalSr+=Sr*Wy[jj];
               totalSl+=Sl*Wy[jj];
             }

             i1=(int)(i+x+0.5);
             j1=(int)(j+y+0.5);
             x=i+x-i1;
             y=j+y-j1;
             Wx[0]=0.5*(0.5-x)*(0.5-x);
             Wx[1]=0.75-x*x;
             Wx[2]=0.5*(x+0.5)*(x+0.5);
             Wy[0]=0.5*(0.5-y)*(0.5-y);
             Wy[1]=0.75-y*y;
             Wy[2]=0.5*(y+0.5)*(y+0.5);

             totalE1=totalB1=0;
             for(jj=0; jj<3; jj++)
             {
               E1=B1=0;
               for(ii=0; ii<3; ii++)
               {
                 E1+=field[i1-1+ii][j1-1+jj].E1*Wx[ii];
                 B1+=field[i1-1+ii][j1-1+jj].B1*Wx[ii];
               }
               totalE1+=E1*Wy[jj];
               totalB1+=B1*Wy[jj];
             }


             p->E1=totalE1+extE1; p->Pr=totalPr+extPr; p->Pl=totalPl+extPl;
             p->B1=totalB1+extB1; p->Sr=totalSr+extSr; p->Sl=totalSl+extSl;
         
             p=p->next;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1

}
*/

void interpolation3D_1st(Domain *D,External *Ext)
{
   int i,j,k,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s,cnt;
   float Ex,Pr,Pl,Bx,Sr,Sl,extEx,extEy,extEz,extBx,extBy,extBz,x,y,z,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;

   Particle ***particle;
   particle=D->particle;

   extEx=Ext->Ex;
   extEy=Ext->Ey;
   extEz=Ext->Ez;
   extBx=Ext->Bx;
   extBy=Ext->By;
   extBz=Ext->Bz;
   if(D->fieldType==1)
   {   
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
         for(k=kstart; k<kend; k++)
         {
           for(s=0; s<D->nSpecies; s++)
           {
             cnt=0;
             p=particle[i][j][k].head[s]->pt;
             while(p)
             {
               x=p->x;  y=p->y;  z=p->z;
               i1=(int)(i+x+0.5);
               j1=(int)(j+y+0.5);
               k1=(int)(k+z+0.5);
               x1=x+0.5-((int)(x+0.5));
               y1=y+0.5-((int)(y+0.5));
               z1=z+0.5-((int)(z+0.5));

               Bx=(1-x)*(1-y1)*(1-z1)*D->Bx[i][j1-1][k1-1]
                 +    x*(1-y1)*(1-z1)*D->Bx[i+1][j1-1][k1-1]
                 +(1-x)*    y1*(1-z1)*D->Bx[i][j1][k1-1]
                 +    x*    y1*(1-z1)*D->Bx[i+1][j1][k1-1]
                 +(1-x)*(1-y1)*    z1*D->Bx[i][j1-1][k1]
                 +    x*(1-y1)*    z1*D->Bx[i+1][j1-1][k1]
                 +(1-x)*    y1*    z1*D->Bx[i][j1][k1]
                 +    x*    y1*    z1*D->Bx[i+1][j1][k1];
               Ex=(1-x)*(1-y)*(1-z)*D->Ex[i][j][k]
                 +    x*(1-y)*(1-z)*D->Ex[i+1][j][k]
                 +(1-x)*    y*(1-z)*D->Ex[i][j+1][k]
                 +    x*    y*(1-z)*D->Ex[i+1][j+1][k]
                 +(1-x)*(1-y)*    z*D->Ex[i][j][k+1]
                 +    x*(1-y)*    z*D->Ex[i+1][j][k+1]
                 +(1-x)*    y*    z*D->Ex[i][j+1][k+1]
                 +    x*    y*    z*D->Ex[i+1][j+1][k+1];
               Pr=(1-x)*(1-y1)*(1-z)*D->Pr[i][j1-1][k]
                 +    x*(1-y1)*(1-z)*D->Pr[i+1][j1-1][k]
                 +(1-x)*    y1*(1-z)*D->Pr[i][j1][k]
                 +    x*    y1*(1-z)*D->Pr[i+1][j1][k]
                 +(1-x)*(1-y1)*    z*D->Pr[i][j1-1][k+1]
                 +    x*(1-y1)*    z*D->Pr[i+1][j1-1][k+1]
                 +(1-x)*    y1*    z*D->Pr[i][j1][k+1]
                 +    x*    y1*    z*D->Pr[i+1][j1][k+1];
               Pl=(1-x)*(1-y1)*(1-z)*D->Pl[i][j1-1][k]
                 +    x*(1-y1)*(1-z)*D->Pl[i+1][j1-1][k]
                 +(1-x)*    y1*(1-z)*D->Pl[i][j1][k]
                 +    x*    y1*(1-z)*D->Pl[i+1][j1][k]
                 +(1-x)*(1-y1)*    z*D->Pl[i][j1-1][k+1]
                 +    x*(1-y1)*    z*D->Pl[i+1][j1-1][k+1]
                 +(1-x)*    y1*    z*D->Pl[i][j1][k+1]
                 +    x*    y1*    z*D->Pl[i+1][j1][k+1];
               Sr=(1-x)*(1-y)*(1-z1)*D->Sr[i][j][k1-1]
                 +    x*(1-y)*(1-z1)*D->Sr[i+1][j][k1-1]
                 +(1-x)*    y*(1-z1)*D->Sr[i][j+1][k1-1]
                 +    x*    y*(1-z1)*D->Sr[i+1][j+1][k1-1]
                 +(1-x)*(1-y)*    z1*D->Sr[i][j][k1]
                 +    x*(1-y)*    z1*D->Sr[i+1][j][k1]
                 +(1-x)*    y*    z1*D->Sr[i][j+1][k1]
                 +    x*    y*    z1*D->Sr[i+1][j+1][k1];
               Sl=(1-x)*(1-y)*(1-z1)*D->Sl[i][j][k1-1]
                 +    x*(1-y)*(1-z1)*D->Sl[i+1][j][k1-1]
                 +(1-x)*    y*(1-z1)*D->Sl[i][j+1][k1-1]
                 +    x*    y*(1-z1)*D->Sl[i+1][j+1][k1-1]
                 +(1-x)*(1-y)*    z1*D->Sl[i][j][k1]
                 +    x*(1-y)*    z1*D->Sl[i+1][j][k1]
                 +(1-x)*    y*    z1*D->Sl[i][j+1][k1]
                 +    x*    y*    z1*D->Sl[i+1][j+1][k1];

               p->Ex=Ex+extEx; p->Ey=Pr+Pl+extEy; p->Ez=Sr+Sl+extEz;
               p->Bx=Bx+extBx; p->By=Sl-Sr+extBy; p->Bz=Pr-Pl+extBz;

               p=p->next;
               cnt++;
             }
           }		//for(s)        
         }		   //for(i,j)
   }           //End of fieldType=1

}


