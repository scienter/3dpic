#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include <math.h>
#include "mpi.h"

void updateCurrent3D_DSX_3rd(Domain *D)
{
    int i,j,k,s,n,i1,i2,j1,j2,k1,k2,ii,jj,kk;
//intXc,intYc;
    int istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub;
    double inverDt,dx,dy,dz,xr,yr,zr,x,y,z,dt;
    double Fx[3],Fy[3],Fz[3],Wx[4],Wy[4],Wz[4];
    double xold,xnew,yold,ynew,zold,znew;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    ptclList *p;
    LoadList *LL;   
    Particle ***particle;
    particle=D->particle;

    double maximum(),minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];
  
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
   
    dt=D->dt;
    dx=D->dx; 
    dy=D->dy;
    dz=D->dz;
    inverDt=1.0/D->dt;
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
      for(j=0; j<nySub+5; j++)
        for(k=0; k<nzSub+5; k++)
      {
        D->JxOld[i][j][k]=D->Jx[i][j][k];
        D->JyOld[i][j][k]=D->Jy[i][j][k];
        D->JzOld[i][j][k]=D->Jz[i][j][k];
        D->Jx[i][j][k]=0.0;
        D->Jy[i][j][k]=0.0;
        D->Jz[i][j][k]=0.0;
      }

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;     
         
            while(p) 
            {
              xold=p->oldX;
              xnew=p->x+i;
              yold=p->oldY;
              ynew=p->y+j;
              zold=p->oldZ;
              znew=p->z+k;

              i1=(int)(xold);
              j1=(int)(yold);
              k1=(int)(zold);
              i2=(int)(xnew);
              j2=(int)(ynew);
              k2=(int)(znew);

              if(i1==i2) 
                xr=0.5*(xold+xnew);
              else 
                xr=maximum(i1*1.0,i2*1.0);
              if(j1==j2) 
                yr=0.5*(yold+ynew);
              else 
                yr=maximum(j1*1.0,j2*1.0);
              if(k1==k2) 
                zr=0.5*(zold+znew);
              else 
                zr=maximum(k1*1.0,k2*1.0);

              //calculate old point
              i1=(int)(0.5*(xold+xr));
              j1=(int)(0.5*(yold+yr));
              k1=(int)(0.5*(zold+zr));
              x=(xold+xr)*0.5-i1;
              y=(yold+yr)*0.5-j1;
              z=(zold+zr)*0.5-k1;

              Fx[0]=(xr-xold)*0.5*(1-x)*(1-x);
              Fx[1]=(xr-xold)*(0.75-(0.5-x)*(0.5-x));
              Fx[2]=(xr-xold)*0.5*x*x;
              Fy[0]=(yr-yold)*0.5*(1-y)*(1-y)*dy*inverDt;
              Fy[1]=(yr-yold)*(0.75-(0.5-y)*(0.5-y))*dy*inverDt;
              Fy[2]=(yr-yold)*0.5*y*y*dy*inverDt;
              Fz[0]=(zr-zold)*0.5*(1-z)*(1-z)*dz*inverDt;
              Fz[1]=(zr-zold)*(0.75-(0.5-z)*(0.5-z))*dz*inverDt;
              Fz[2]=(zr-zold)*0.5*z*z*dz*inverDt;
              x1=1+x;
              x2=x;
              x3=1-x;
              x4=2-x;
              Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
              Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
              Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
              Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
              y1=1+y;
              y2=y;
              y3=1-y;
              y4=2-y;
              Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
              Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
              Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
              Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
              z1=1+z;
              z2=z;
              z3=1-z;
              z4=2-z;
              Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
              Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
              Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
              Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;

              for(ii=0; ii<3; ii++)
                for(jj=0; jj<4; jj++)
                  for(kk=0; kk<4; kk++)
                    D->Jx[i1-1+ii][j1-1+jj][k1-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s];

              for(ii=0; ii<4; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<4; kk++)
                    D->Jy[i1-1+ii][j1-1+jj][k1-1+kk]+=Fy[jj]*Wx[ii]*Wz[kk]*coeff[s];

              for(ii=0; ii<4; ii++)
                for(jj=0; jj<4; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jz[i1-1+ii][j1-1+jj][k1-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s];

              //calculate new point
              i1=(int)(0.5*(xnew+xr));
              j1=(int)(0.5*(ynew+yr));
              k1=(int)(0.5*(znew+zr));
              x=(xnew+xr)*0.5-i1;
              y=(ynew+yr)*0.5-j1;
              z=(znew+zr)*0.5-k1;
              Fx[0]=(xnew-xr)*0.5*(1-x)*(1-x);
              Fx[1]=(xnew-xr)*(0.75-(0.5-x)*(0.5-x));
              Fx[2]=(xnew-xr)*0.5*x*x;
              Fy[0]=(ynew-yr)*0.5*(1-y)*(1-y)*dy*inverDt;
              Fy[1]=(ynew-yr)*(0.75-(0.5-y)*(0.5-y))*dy*inverDt;
              Fy[2]=(ynew-yr)*0.5*y*y*dy*inverDt;
              Fz[0]=(znew-zr)*0.5*(1-z)*(1-z)*dz*inverDt;
              Fz[1]=(znew-zr)*(0.75-(0.5-z)*(0.5-z))*dz*inverDt;
              Fz[2]=(znew-zr)*0.5*z*z*dz*inverDt;
              x1=1+x;
              x2=x;
              x3=1-x;
              x4=2-x;
              Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
              Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
              Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
              Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
              y1=1+y;
              y2=y;
              y3=1-y;
              y4=2-y;
              Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
              Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
              Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
              Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
              z1=1+z;
              z2=z;
              z3=1-z;
              z4=2-z;
              Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
              Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
              Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
              Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;

              for(ii=0; ii<3; ii++)
                for(jj=0; jj<4; jj++)
                  for(kk=0; kk<4; kk++)
                    D->Jx[i1-1+ii][j1-1+jj][k1-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s];

              for(ii=0; ii<4; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<4; kk++)
                    D->Jy[i1-1+ii][j1-1+jj][k1-1+kk]+=Fy[jj]*Wx[ii]*Wz[kk]*coeff[s];

              for(ii=0; ii<4; ii++)
                for(jj=0; jj<4; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jz[i1-1+ii][j1-1+jj][k1-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s];

               p=p->next;
             }	//End of while(p)
          }	//End of for(s)     
       }		//End of for(i,j)
}


void updateCurrent3D_DSX_2nd(Domain *D)
{
    int i,j,k,s,n,i1,i2,j1,j2,k1,k2,ii,jj,kk;
    int istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub;
    double inverDt,dx,dt,dy,dz,xr,yr,zr;
    double Fx[2],Fy[2],Fz[2],Wx[3],Wy[3],Wz[3];
    double xold,xnew,yold,ynew,zold,znew,x,y,z;
    ptclList *p;
    LoadList *LL;   
    Particle ***particle;
    particle=D->particle;

    double maximum(),minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
   
    dt=D->dt;
    dx=D->dx; 
    dy=D->dy;
    dz=D->dz;
    inverDt=1.0/D->dt;
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
      for(j=0; j<nySub+5; j++)
        for(k=0; k<nzSub+5; k++)
        {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;     
         
            while(p) 
            {
              xold=p->oldX;
              xnew=p->x+i;
              yold=p->oldY;
              ynew=p->y+j;
              zold=p->oldZ;
              znew=p->z+k;

              i1=(int)(xold+0.5);
              j1=(int)(yold+0.5);
              k1=(int)(zold+0.5);
              i2=(int)(xnew+0.5);
              j2=(int)(ynew+0.5);
              k2=(int)(znew+0.5);

              if(i1==i2) 
                xr=0.5*(xold+xnew);
              else             
                xr=(i1+i2)*0.5;
              if(j1==j2) 
                yr=0.5*(yold+ynew);
              else             
                yr=(j1+j2)*0.5;
              if(k1==k2) 
                zr=0.5*(zold+znew);
              else             
                zr=(k1+k2)*0.5;

              x=(xold+xr)*0.5-i1;
              y=(yold+yr)*0.5-j1;
              z=(zold+zr)*0.5-k1;
              Fx[0]=(xr-xold)*(0.5-x);
              Fx[1]=(xr-xold)*(0.5+x);
              Fy[0]=(yr-yold)*dy*inverDt*(0.5-y);
              Fy[1]=(yr-yold)*dy*inverDt*(0.5+y);
              Fz[0]=(zr-zold)*dz*inverDt*(0.5-z);
              Fz[1]=(zr-zold)*dz*inverDt*(0.5+z);                                 
              Wx[0]=0.5*(0.5-x)*(0.5-x);
              Wx[1]=0.75-x*x;
              Wx[2]=0.5*(0.5+x)*(0.5+x);
              Wy[0]=0.5*(0.5-y)*(0.5-y);
              Wy[1]=0.75-y*y;
              Wy[2]=0.5*(0.5+y)*(0.5+y);
              Wz[0]=0.5*(0.5-z)*(0.5-z);
              Wz[1]=0.75-z*z;
              Wz[2]=0.5*(0.5+z)*(0.5+z);

              for(ii=0; ii<2; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jx[i1-1+ii][j1-1+jj][k1-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s];
              for(ii=0; ii<3; ii++)
                for(jj=0; jj<2; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jy[i1-1+ii][j1-1+jj][k1-1+kk]+=Fy[jj]*Wx[ii]*Wz[kk]*coeff[s];
              for(ii=0; ii<3; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<2; kk++)
                    D->Jz[i1-1+ii][j1-1+jj][k1-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s];


              i1=(int)(xr+0.5);
              j1=(int)(yr+0.5);
              k1=(int)(zr+0.5);
              x=(xnew+xr)*0.5-i1;
              y=(ynew+yr)*0.5-j1;
              z=(znew+zr)*0.5-k1;
              Fx[0]=(xnew-xr)*(0.5-x);
              Fx[1]=(xnew-xr)*(0.5+x);
              Fy[0]=(ynew-yr)*dy*inverDt*(0.5-y);
              Fy[1]=(ynew-yr)*dy*inverDt*(0.5+y);
              Fz[0]=(znew-zr)*dz*inverDt*(0.5-z);
              Fz[1]=(znew-zr)*dz*inverDt*(0.5+z);
                                  
              Wx[0]=0.5*(0.5-x)*(0.5-x);
              Wx[1]=0.75-x*x;
              Wx[2]=0.5*(0.5+x)*(0.5+x);
              Wy[0]=0.5*(0.5-y)*(0.5-y);
              Wy[1]=0.75-y*y;
              Wy[2]=0.5*(0.5+y)*(0.5+y);
              Wz[0]=0.5*(0.5-z)*(0.5-z);
              Wz[1]=0.75-z*z;
              Wz[2]=0.5*(0.5+z)*(0.5+z);

              for(ii=0; ii<2; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jx[i1-1+ii][j1-1+jj][k1-1+kk]+=Fx[ii]*Wy[jj]*Wz[kk]*coeff[s];
              for(ii=0; ii<3; ii++)
                for(jj=0; jj<2; jj++)
                  for(kk=0; kk<3; kk++)
                    D->Jy[i1-1+ii][j1-1+jj][k1-1+kk]+=Fy[jj]*Wx[ii]*Wz[kk]*coeff[s];
              for(ii=0; ii<3; ii++)
                for(jj=0; jj<3; jj++)
                  for(kk=0; kk<2; kk++)
                    D->Jz[i1-1+ii][j1-1+jj][k1-1+kk]+=Fz[kk]*Wx[ii]*Wy[jj]*coeff[s];

              p=p->next;
            }	//End of while(p)
          }		//End of for(s)     
        }		//End of for(i,j)

}


void updateCurrent3D_DSX_1st(Domain *D)
{
    int i,j,k,s,i1,i2,j1,j2,k1,k2;
    int istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub;
    double inverDt,x1,x2,y1,y2,xr,yr,zr,z1,z2;
    double Fx1,Fx2,Wx1,Wx2,Fy1,Fy2,Wy1,Wy2,Fz1,Fz2,Wz1,Wz2,dx,dy,dz,dt;
    ptclList *p;
    LoadList *LL;   
    Particle ***particle;
    particle=D->particle;

    double maximum();
    double minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
   
    dt=D->dt;
    dx=D->dx; 
    dy=D->dy;
    dz=D->dz;
    inverDt=1.0/D->dt;
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
      for(j=0; j<nySub+5; j++)
        for(k=0; k<nzSub+5; k++)
        {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;     
         
            while(p) 
            {
              x2=p->x+i;       y2=p->y+j;       z2=p->z+k;
              x1=p->oldX;      y1=p->oldY;      z1=p->oldZ;
              i1=(int)x1;      j1=(int)y1;      k1=(int)z1;
              i2=(int)x2;      j2=(int)y2;      k2=(int)z2;
              if(i1==i2) 
                xr=0.5*(x1+x2);
              else 
                xr=maximum(i1*1.0,i2*1.0);
              if(j1==j2) 
                yr=0.5*(y1+y2);
              else 
                yr=maximum(j1*1.0,j2*1.0);
              if(z1==z2) 
                zr=0.5*(z1+z2);
              else 
                zr=maximum(k1*1.0,k2*1.0);

              Fx1=(xr-x1);
              Fx2=(x2-xr);
              Fy1=(yr-y1)*dy*inverDt;
              Fy2=(y2-yr)*dy*inverDt;
              Fz1=(zr-z1)*dz*inverDt;
              Fz2=(z2-zr)*dz*inverDt;
              Wx1=0.5*(x1+xr)-i1;
              Wx2=0.5*(xr+x2)-i2;
              Wy1=0.5*(y1+yr)-j1;
              Wy2=0.5*(yr+y2)-j2;
              Wz1=0.5*(z1+zr)-k1;
              Wz2=0.5*(zr+z2)-k2;
 
              D->Jx[i1][j1][k1]    +=Fx1*(1-Wy1)*(1-Wz1)*coeff[s];
              D->Jx[i1][j1+1][k1]  +=Fx1*    Wy1*(1-Wz1)*coeff[s];
              D->Jx[i1][j1][k1+1]  +=Fx1*(1-Wy1)*    Wz1*coeff[s];
              D->Jx[i1][j1+1][k1+1]+=Fx1*    Wy1*    Wz1*coeff[s];
              D->Jy[i1][j1][k1]    +=Fy1*(1-Wx1)*(1-Wz1)*coeff[s];
              D->Jy[i1+1][j1][k1]  +=Fy1*    Wx1*(1-Wz1)*coeff[s];
              D->Jy[i1][j1][k1+1]  +=Fy1*(1-Wx1)*    Wz1*coeff[s];
              D->Jy[i1+1][j1][k1+1]+=Fy1*    Wx1*    Wz1*coeff[s];
              D->Jz[i1][j1][k1]    +=Fz1*(1-Wx1)*(1-Wy1)*coeff[s];
              D->Jz[i1+1][j1][k1]  +=Fz1*    Wx1*(1-Wy1)*coeff[s];
              D->Jz[i1][j1+1][k1]  +=Fz1*(1-Wx1)*    Wy1*coeff[s];
              D->Jz[i1+1][j1+1][k1]+=Fz1*    Wx1*    Wy1*coeff[s];

              D->Jx[i2][j2][k2]    +=Fx2*(1-Wy2)*(1-Wz2)*coeff[s];
              D->Jx[i2][j2+1][k2]  +=Fx2*    Wy2*(1-Wz2)*coeff[s];
              D->Jx[i2][j2][k2+1]  +=Fx2*(1-Wy2)*    Wz2*coeff[s];
              D->Jx[i2][j2+1][k2+1]+=Fx2*    Wy2*    Wz2*coeff[s];
              D->Jy[i2][j2][k2]    +=Fy2*(1-Wx2)*(1-Wz2)*coeff[s];
              D->Jy[i2+1][j2][k2]  +=Fy2*    Wx2*(1-Wz2)*coeff[s];
              D->Jy[i2][j2][k2+1]  +=Fy2*(1-Wx2)*    Wz2*coeff[s];
              D->Jy[i2+1][j2][k2+1]+=Fy2*    Wx2*    Wz2*coeff[s];
              D->Jz[i2][j2][k2]    +=Fz2*(1-Wx2)*(1-Wy2)*coeff[s];
              D->Jz[i2+1][j2][k2]  +=Fz2*    Wx2*(1-Wy2)*coeff[s];
              D->Jz[i2][j2+1][k2]  +=Fz2*(1-Wx2)*    Wy2*coeff[s];
              D->Jz[i2+1][j2+1][k2]+=Fz2*    Wx2*    Wy2*coeff[s];

              p=p->next;
            }	//End of while(p)
          }	   //End of for(s)     
        }      //End of for(i,j)
}


double findR(double x1, double x2,double x3, double x4)
{
  double minimum();
  double maximum();
  double result,result1,result2,result3;

  result1=minimum(x1-0.5,x2-0.5);
  result2=maximum(x1-1.5,x2-1.5);
  result3=maximum(result2,(x3+x4)*0.5);
  result=minimum(result1,result3);

  return result;
}

double maximum(double x1,double x2)
{
   double result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

double minimum(double x1,double x2)
{
   double result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}

int intmaximum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

int intminimum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}
