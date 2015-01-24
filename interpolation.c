#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void interpolation3D_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,k,ii,jj,kk,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s;
   float x,y,z,Pr,Pl,Sr,Sl,Ex,Bx;
   float totalPr,totalPl,totalSr,totalSl,totalEx,totalBx;
   float Wx[3],Wy[3],Wz[3],WxC[3],WyC[3],WzC[3];
   float extEx,extEy,extEz,extBx,extBy,extBz;
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
             p=particle[i][j][k].head[s]->pt;
             while(p)
             {
               x=p->x;  y=p->y;  z=p->z;
               WxC[0]=0.5*(1-x)*(1-x);
               WxC[1]=0.75-(0.5-x)*(0.5-x);
               WxC[2]=0.5*x*x;
               WyC[0]=0.5*(1-y)*(1-y);
               WyC[1]=0.75-(0.5-y)*(0.5-y);
               WyC[2]=0.5*y*y;
               WzC[0]=0.5*(1-z)*(1-z);
               WzC[1]=0.75-(0.5-z)*(0.5-z);
               WzC[2]=0.5*z*z;

               
               i1=(int)(i+x+0.5);
               j1=(int)(j+y+0.5);
               k1=(int)(k+z+0.5);
               x=i+x-i1;
               y=j+y-j1;
               z=k+z-k1;
               Wx[0]=0.5*(0.5-x)*(0.5-x);
               Wx[1]=0.75-x*x;
               Wx[2]=0.5*(x+0.5)*(x+0.5);
               Wy[0]=0.5*(0.5-y)*(0.5-y);
               Wy[1]=0.75-y*y;
               Wy[2]=0.5*(y+0.5)*(y+0.5);
               Wz[0]=0.5*(0.5-z)*(0.5-z);
               Wz[1]=0.75-z*z;
               Wz[2]=0.5*(z+0.5)*(z+0.5);

               Pr=Pl=Sr=Sl=Ex=Bx=0.0;
               for(kk=0; kk<3; kk++)
                 for(jj=0; jj<3; jj++)
                   for(ii=0; ii<3; ii++)
                   {
                     Pr+=D->Pr[i-1+ii][j-1+jj][k1-1+kk]*WxC[ii]*WyC[jj]*Wz[kk];
                     Pl+=D->Pl[i-1+ii][j-1+jj][k1-1+kk]*WxC[ii]*WyC[jj]*Wz[kk];
                     Sr+=D->Sr[i-1+ii][j1-1+jj][k-1+kk]*WxC[ii]*Wy[jj]*WzC[kk];
                     Sl+=D->Sl[i-1+ii][j1-1+jj][k-1+kk]*WxC[ii]*Wy[jj]*WzC[kk];
                     Ex+=D->Ex[i-1+ii][j1-1+jj][k1-1+kk]*WxC[ii]*Wy[jj]*Wz[kk];
                     Bx+=D->Bx[i-1+ii][j-1+jj][k-1+kk]*WxC[ii]*WyC[jj]*WzC[kk];
                   }


               p->Ex=Ex+extEx; p->Ey=Pr+Pl+extEy; p->Ez=Sr+Sl+extEz;
               p->Bx=Bx+extBx; p->By=Sl-Sr+extBy; p->Bz=Pr-Pl+extBz;
         
               p=p->next;
             }
           }		//for(s)        
         }		   //for(i,j)
   }           //End of fieldType=1

}


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

               Bx=(1-x1)*(1-y1)*(1-z1)*D->Bx[i1-1][j1-1][k1-1]
                 +    x1*(1-y1)*(1-z1)*D->Bx[i1][j1-1][k1-1]
                 +(1-x1)*    y1*(1-z1)*D->Bx[i1-1][j1][k1-1]
                 +    x1*    y1*(1-z1)*D->Bx[i1][j1][k1-1]
                 +(1-x1)*(1-y1)*    z1*D->Bx[i1-1][j1-1][k1]
                 +    x1*(1-y1)*    z1*D->Bx[i1][j1-1][k1]
                 +(1-x1)*    y1*    z1*D->Bx[i1-1][j1][k1]
                 +    x1*    y1*    z1*D->Bx[i1][j1][k1];
               Ex=(1-x1)*(1-y)*(1-z)*D->Ex[i1-1][j][k]
                 +    x1*(1-y)*(1-z)*D->Ex[i1][j][k]
                 +(1-x1)*    y*(1-z)*D->Ex[i1-1][j+1][k]
                 +    x1*    y*(1-z)*D->Ex[i1][j+1][k]
                 +(1-x1)*(1-y)*    z*D->Ex[i1-1][j][k+1]
                 +    x1*(1-y)*    z*D->Ex[i1][j][k+1]
                 +(1-x1)*    y*    z*D->Ex[i1-1][j+1][k+1]
                 +    x1*    y*    z*D->Ex[i1][j+1][k+1];
               Pr=(1-x1)*(1-y1)*(1-z)*D->Pr[i1-1][j1-1][k]
                 +    x1*(1-y1)*(1-z)*D->Pr[i1][j1-1][k]
                 +(1-x1)*y1    *(1-z)*D->Pr[i1-1][j1][k]
                 +    x1*    y1*(1-z)*D->Pr[i1][j1][k]
                 +(1-x1)*(1-y1)*    z*D->Pr[i1-1][j1-1][k+1]
                 +    x1*(1-y1)*    z*D->Pr[i1][j1-1][k+1]
                 +(1-x1)*    y1*    z*D->Pr[i1-1][j1][k+1]
                 +    x1*    y1*    z*D->Pr[i1][j1][k+1];
               Pl=(1-x1)*(1-y1)*(1-z)*D->Pl[i1-1][j1-1][k]
                 +    x1*(1-y1)*(1-z)*D->Pl[i1][j1-1][k]
                 +(1-x1)*y1    *(1-z)*D->Pl[i1-1][j1][k]
                 +    x1*    y1*(1-z)*D->Pl[i1][j1][k]
                 +(1-x1)*(1-y1)*    z*D->Pl[i1-1][j1-1][k+1]
                 +    x1*(1-y1)*    z*D->Pl[i1][j1-1][k+1]
                 +(1-x1)*    y1*    z*D->Pl[i1-1][j1][k+1]
                 +    x1*    y1*    z*D->Pl[i1][j1][k+1];
               Sr=(1-x1)*(1-y)*(1-z1)*D->Sr[i1-1][j][k1-1]
                 +    x1*(1-y)*(1-z1)*D->Sr[i1][j][k1-1]
                 +(1-x1)*    y*(1-z1)*D->Sr[i1-1][j+1][k1-1]
                 +    x1*    y*(1-z1)*D->Sr[i1][j+1][k1-1]
                 +(1-x1)*(1-y)*    z1*D->Sr[i1-1][j][k1]
                 +    x1*(1-y)*    z1*D->Sr[i1][j][k1]
                 +(1-x1)*    y*    z1*D->Sr[i1-1][j+1][k1]
                 +    x1*    y*    z1*D->Sr[i1][j+1][k1];
               Sl=(1-x1)*(1-y)*(1-z1)*D->Sl[i1-1][j][k1-1]
                 +    x1*(1-y)*(1-z1)*D->Sl[i1][j][k1-1]
                 +(1-x1)*    y*(1-z1)*D->Sl[i1-1][j+1][k1-1]
                 +    x1*    y*(1-z1)*D->Sl[i1][j+1][k1-1]
                 +(1-x1)*(1-y)*    z1*D->Sl[i1-1][j][k1]
                 +    x1*(1-y)*    z1*D->Sl[i1][j][k1]
                 +(1-x1)*    y*    z1*D->Sl[i1-1][j+1][k1]
                 +    x1*    y*    z1*D->Sl[i1][j+1][k1];

               p->Ex=Ex+extEx; p->Ey=Pr+Pl+extEy; p->Ez=Sr+Sl+extEz;
               p->Bx=Bx+extBx; p->By=Sl-Sr+extBy; p->Bz=Pr-Pl+extBz;

               p=p->next;
               cnt++;
             }
           }		//for(s)        
         }		   //for(i,j)
   }           //End of fieldType=1

}


