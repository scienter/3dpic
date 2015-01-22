#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include <math.h>
#include "mpi.h"
/*
void updateCurrent2D_DSX_3rd(Domain *D)
{
    int i,j,s,n,i1,i2,j1,j2,intXc,intYc;
    int istart,iend,jstart,jend;
    int nxSub,nySub;
    double inverDt,dx,dt,dy,xr,yr,x,y;
    double Fx1,Fx2,Fx3,Wx1,Wx2,Wx3,Wx4;
    double Fy1,Fy2,Fy3,Wy1,Wy2,Wy3,Wy4;
    double vz,gamma,xc,yc,wx,wy,xcc,ycc,xold,xnew,yold,ynew;
    double x1,x2,x3,x4,y1,y2,y3,y4;
    ptclList *p;
    FieldDSX **field; 
    LoadList *LL;   
    field=D->fieldDSX;
    Particle **particle;
    particle=D->particle;

    double maximum(),minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];
  
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    nxSub=D->nxSub;
    nySub=D->nySub;
   
    dt=D->dt;
    dx=D->dx; 
    dy=D->dy;
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
      {
        field[i][j].J1Old=field[i][j].J1;
        field[i][j].J2Old=field[i][j].J2;
        field[i][j].J3Old=field[i][j].J3;
        field[i][j].J1=0.0;
        field[i][j].J2=0.0;
        field[i][j].J3=0.0;
      }

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;     
         
          while(p) 
          {
            gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
            vz=p->p3/gamma;
            xold=p->oldX;
            xnew=p->x+i;
            yold=p->oldY;
            ynew=p->y+j;

            xc=0.5*(xold+xnew);
            yc=0.5*(yold+ynew);
            intXc=(int)(xc);
            intYc=(int)(yc);
            xcc=xc-intXc;
            ycc=yc-intYc;

            i1=(int)(xold);
            j1=(int)(yold);
            i2=(int)(xnew);
            j2=(int)(ynew);

            if(i1==i2) 
              xr=0.5*(xold+xnew);
            else 
              xr=maximum(i1*1.0,i2*1.0);
            if(j1==j2) 
              yr=0.5*(yold+ynew);
            else 
              yr=maximum(j1*1.0,j2*1.0);

            //calculate old point
            i1=(int)(0.5*(xold+xr));
            j1=(int)(0.5*(yold+yr));
            x=(xold+xr)*0.5-i1;
            y=(yold+yr)*0.5-j1;
            Fx1=(xr-xold)*0.5*(1-x)*(1-x);
            Fx2=(xr-xold)*(0.75-(0.5-x)*(0.5-x));
            Fx3=(xr-xold)*0.5*x*x;
            y1=1+y;
            y2=y;
            y3=1-y;
            y4=2-y;
            Wy1=(2-y1)*(2-y1)*(2-y1)/6.0;
            Wy2=(4-6*y2*y2+3*y2*y2*y2)/6.0;
            Wy3=(4-6*y3*y3+3*y3*y3*y3)/6.0;
            Wy4=(2-y4)*(2-y4)*(2-y4)/6.0;

            field[i1-1][j1-1].J1+=Fx1*Wy1*coeff[s];
            field[i1-1][j1+0].J1+=Fx1*Wy2*coeff[s];
            field[i1-1][j1+1].J1+=Fx1*Wy3*coeff[s];
            field[i1-1][j1+2].J1+=Fx1*Wy4*coeff[s];
            field[i1+0][j1-1].J1+=Fx2*Wy1*coeff[s];
            field[i1+0][j1+0].J1+=Fx2*Wy2*coeff[s];
            field[i1+0][j1+1].J1+=Fx2*Wy3*coeff[s];
            field[i1+0][j1+2].J1+=Fx2*Wy4*coeff[s];
            field[i1+1][j1-1].J1+=Fx3*Wy1*coeff[s];
            field[i1+1][j1+0].J1+=Fx3*Wy2*coeff[s];
            field[i1+1][j1+1].J1+=Fx3*Wy3*coeff[s];
            field[i1+1][j1+2].J1+=Fx3*Wy4*coeff[s];

            Fy1=(yr-yold)*0.5*(1-y)*(1-y)*dy*inverDt;
            Fy2=(yr-yold)*(0.75-(0.5-y)*(0.5-y))*dy*inverDt;
            Fy3=(yr-yold)*0.5*y*y*dy*inverDt;
            x1=1+x;
            x2=x;
            x3=1-x;
            x4=2-x;
            Wx1=(2-x1)*(2-x1)*(2-x1)/6.0;
            Wx2=(4-6*x2*x2+3*x2*x2*x2)/6.0;
            Wx3=(4-6*x3*x3+3*x3*x3*x3)/6.0;
            Wx4=(2-x4)*(2-x4)*(2-x4)/6.0;

            field[i1-1][j1-1].J2+=Fy1*Wx1*coeff[s];
            field[i1+0][j1-1].J2+=Fy1*Wx2*coeff[s];
            field[i1+1][j1-1].J2+=Fy1*Wx3*coeff[s];
            field[i1+2][j1-1].J2+=Fy1*Wx4*coeff[s];
            field[i1-1][j1+0].J2+=Fy2*Wx1*coeff[s];
            field[i1+0][j1+0].J2+=Fy2*Wx2*coeff[s];
            field[i1+1][j1+0].J2+=Fy2*Wx3*coeff[s];
            field[i1+2][j1+0].J2+=Fy2*Wx4*coeff[s];
            field[i1-1][j1+1].J2+=Fy3*Wx1*coeff[s];
            field[i1+0][j1+1].J2+=Fy3*Wx2*coeff[s];
            field[i1+1][j1+1].J2+=Fy3*Wx3*coeff[s];
            field[i1+2][j1+1].J2+=Fy3*Wx4*coeff[s];

            //calculate new point
            i2=(int)(0.5*(xnew+xr));
            j2=(int)(0.5*(ynew+yr));
            x=(xnew+xr)*0.5-i2;
            y=(ynew+yr)*0.5-j2;
            Fx1=(xnew-xr)*0.5*(1-x)*(1-x);
            Fx2=(xnew-xr)*(0.75-(0.5-x)*(0.5-x));
            Fx3=(xnew-xr)*0.5*x*x;
            y1=1+y;
            y2=y;
            y3=1-y;
            y4=2-y;
            Wy1=(2-y1)*(2-y1)*(2-y1)/6.0;
            Wy2=(4-6*y2*y2+3*y2*y2*y2)/6.0;
            Wy3=(4-6*y3*y3+3*y3*y3*y3)/6.0;
            Wy4=(2-y4)*(2-y4)*(2-y4)/6.0;

            field[i2-1][j2-1].J1+=Fx1*Wy1*coeff[s];
            field[i2-1][j2+0].J1+=Fx1*Wy2*coeff[s];
            field[i2-1][j2+1].J1+=Fx1*Wy3*coeff[s];
            field[i2-1][j2+2].J1+=Fx1*Wy4*coeff[s];
            field[i2+0][j2-1].J1+=Fx2*Wy1*coeff[s];
            field[i2+0][j2+0].J1+=Fx2*Wy2*coeff[s];
            field[i2+0][j2+1].J1+=Fx2*Wy3*coeff[s];
            field[i2+0][j2+2].J1+=Fx2*Wy4*coeff[s];
            field[i2+1][j2-1].J1+=Fx3*Wy1*coeff[s];
            field[i2+1][j2+0].J1+=Fx3*Wy2*coeff[s];
            field[i2+1][j2+1].J1+=Fx3*Wy3*coeff[s];
            field[i2+1][j2+2].J1+=Fx3*Wy4*coeff[s];

            Fy1=(ynew-yr)*0.5*(1-y)*(1-y)*dy*inverDt;
            Fy2=(ynew-yr)*(0.75-(0.5-y)*(0.5-y))*dy*inverDt;
            Fy3=(ynew-yr)*0.5*y*y*dy*inverDt;
            x1=1+x;
            x2=x;
            x3=1-x;
            x4=2-x;
            Wx1=(2-x1)*(2-x1)*(2-x1)/6.0;
            Wx2=(4-6*x2*x2+3*x2*x2*x2)/6.0;
            Wx3=(4-6*x3*x3+3*x3*x3*x3)/6.0;
            Wx4=(2-x4)*(2-x4)*(2-x4)/6.0;

            field[i2-1][j2-1].J2+=Fy1*Wx1*coeff[s];
            field[i2+0][j2-1].J2+=Fy1*Wx2*coeff[s];
            field[i2+1][j2-1].J2+=Fy1*Wx3*coeff[s];
            field[i2+2][j2-1].J2+=Fy1*Wx4*coeff[s];
            field[i2-1][j2+0].J2+=Fy2*Wx1*coeff[s];
            field[i2+0][j2+0].J2+=Fy2*Wx2*coeff[s];
            field[i2+1][j2+0].J2+=Fy2*Wx3*coeff[s];
            field[i2+2][j2+0].J2+=Fy2*Wx4*coeff[s];
            field[i2-1][j2+1].J2+=Fy3*Wx1*coeff[s];
            field[i2+0][j2+1].J2+=Fy3*Wx2*coeff[s];
            field[i2+1][j2+1].J2+=Fy3*Wx3*coeff[s];
            field[i2+2][j2+1].J2+=Fy3*Wx4*coeff[s];

            //calculate z component
            x=xcc;   y=ycc;
            x1=1+x;
            x2=x;
            x3=1-x;
            x4=2-x;
            y1=1+y;
            y2=y;
            y3=1-y;
            y4=2-y;
            Wx1=(2-x1)*(2-x1)*(2-x1)/6.0;
            Wx2=(4-6*x2*x2+3*x2*x2*x2)/6.0;
            Wx3=(4-6*x3*x3+3*x3*x3*x3)/6.0;
            Wx4=(2-x4)*(2-x4)*(2-x4)/6.0;
            Wy1=(2-y1)*(2-y1)*(2-y1)/6.0;
            Wy2=(4-6*y2*y2+3*y2*y2*y2)/6.0;
            Wy3=(4-6*y3*y3+3*y3*y3*y3)/6.0;
            Wy4=(2-y4)*(2-y4)*(2-y4)/6.0;

            field[intXc-1][intYc-1].J3+=Wx1*Wy1*vz*coeff[s];
            field[intXc+0][intYc-1].J3+=Wx2*Wy1*vz*coeff[s];
            field[intXc+1][intYc-1].J3+=Wx3*Wy1*vz*coeff[s];
            field[intXc+2][intYc-1].J3+=Wx4*Wy1*vz*coeff[s];
            field[intXc-1][intYc+0].J3+=Wx1*Wy2*vz*coeff[s];
            field[intXc+0][intYc+0].J3+=Wx2*Wy2*vz*coeff[s];
            field[intXc+1][intYc+0].J3+=Wx3*Wy2*vz*coeff[s];
            field[intXc+2][intYc+0].J3+=Wx4*Wy2*vz*coeff[s];
            field[intXc-1][intYc+1].J3+=Wx1*Wy3*vz*coeff[s];
            field[intXc+0][intYc+1].J3+=Wx2*Wy3*vz*coeff[s];
            field[intXc+1][intYc+1].J3+=Wx3*Wy3*vz*coeff[s];
            field[intXc+2][intYc+1].J3+=Wx4*Wy3*vz*coeff[s];
            field[intXc-1][intYc+2].J3+=Wx1*Wy4*vz*coeff[s];
            field[intXc+0][intYc+2].J3+=Wx2*Wy4*vz*coeff[s];
            field[intXc+1][intYc+2].J3+=Wx3*Wy4*vz*coeff[s];
            field[intXc+2][intYc+2].J3+=Wx4*Wy4*vz*coeff[s];

             p=p->next;
           }	//End of while(p)
        }	//End of for(s)     
     }		//End of for(i,j)

}
*/
/*
void updateCurrent2D_DSX_2nd(Domain *D)
{
    int i,j,s,n,i1,i2,j1,j2,intXc,intYc;
    int istart,iend,jstart,jend;
    int nxSub,nySub;
    double inverDt,dx,dt,dy,xr,yr;
    double Fx1,Fx2,Wx1,Wx2,Wx3,Wi1,Wj1;
    double Fy1,Fy2,Wy1,Wy2,Wy3;
    double vz,gamma,xc,yc,wx,wy,xcc,ycc,xold,xnew,yold,ynew;
    ptclList *p;
    FieldDSX **field; 
    LoadList *LL;   
    field=D->fieldDSX;
    Particle **particle;
    particle=D->particle;

    double maximum(),minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    nxSub=D->nxSub;
    nySub=D->nySub;
   
    dt=D->dt;
    dx=D->dx; 
    dy=D->dy;
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
      {
        field[i][j].J1Old=field[i][j].J1;
        field[i][j].J2Old=field[i][j].J2;
        field[i][j].J3Old=field[i][j].J3;
        field[i][j].J1=0.0;
        field[i][j].J2=0.0;
        field[i][j].J3=0.0;
      }

    for(i=istart; i<=iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;     
         
          while(p) 
          {
            gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
            vz=p->p3/gamma;
            xold=p->oldX;
            xnew=p->x+i;
            yold=p->oldY;
            ynew=p->y+j;

            xc=0.5*(xold+xnew);
            yc=0.5*(yold+ynew);
            intXc=(int)(xc+0.5);
            intYc=(int)(yc+0.5);
            xcc=xc-intXc;
            ycc=yc-intYc;

            i1=(int)(xold+0.5);
            j1=(int)(yold+0.5);
            i2=(int)(xnew+0.5);
            j2=(int)(ynew+0.5);

            if(i1==i2) 
              xr=0.5*(xold+xnew);
            else             
              xr=(i1+i2)*0.5;
            if(j1==j2) 
              yr=0.5*(yold+ynew);
            else             
              yr=(j1+j2)*0.5;

            Wi1=(xold+xr)*0.5-i1;
            Wj1=(yold+yr)*0.5-j1;
            Fx1=(xr-xold)*(0.5-Wi1);
            Fx2=(xr-xold)*(0.5+Wi1);
            Fy1=(yr-yold)*dy*inverDt*(0.5-Wj1);
            Fy2=(yr-yold)*dy*inverDt*(0.5+Wj1);
                                 
            Wx1=0.5*(0.5-Wi1)*(0.5-Wi1);
            Wx2=0.75-Wi1*Wi1;
            Wx3=0.5*(0.5+Wi1)*(0.5+Wi1);
            Wy1=0.5*(0.5-Wj1)*(0.5-Wj1);
            Wy2=0.75-Wj1*Wj1;
            Wy3=0.5*(0.5+Wj1)*(0.5+Wj1);

            field[i1-1][j1-1].J1+=Fx1*Wy1*coeff[s];
            field[i1-1][j1].J1+=Fx1*Wy2*coeff[s];
            field[i1-1][j1+1].J1+=Fx1*Wy3*coeff[s];
            field[i1][j1-1].J1+=Fx2*Wy1*coeff[s];
            field[i1][j1].J1+=Fx2*Wy2*coeff[s];
            field[i1][j1+1].J1+=Fx2*Wy3*coeff[s];

            field[i1-1][j1-1].J2+=Fy1*Wx1*coeff[s];
            field[i1][j1-1].J2+=Fy1*Wx2*coeff[s];
            field[i1+1][j1-1].J2+=Fy1*Wx3*coeff[s];
            field[i1-1][j1].J2+=Fy2*Wx1*coeff[s];
            field[i1][j1].J2+=Fy2*Wx2*coeff[s];
            field[i1+1][j1].J2+=Fy2*Wx3*coeff[s];

            i1=(int)(xr+0.5);
            j1=(int)(yr+0.5);

            Wi1=(xnew+xr)*0.5-i1;
            Wj1=(ynew+yr)*0.5-j1;
            Fx1=(xnew-xr)*(0.5-Wi1);
            Fx2=(xnew-xr)*(0.5+Wi1);
            Fy1=(ynew-yr)*dy*inverDt*(0.5-Wj1);
            Fy2=(ynew-yr)*dy*inverDt*(0.5+Wj1);
                                 
            Wx1=0.5*(0.5-Wi1)*(0.5-Wi1);
            Wx2=0.75-Wi1*Wi1;
            Wx3=0.5*(0.5+Wi1)*(0.5+Wi1);
            Wy1=0.5*(0.5-Wj1)*(0.5-Wj1);
            Wy2=0.75-Wj1*Wj1;
            Wy3=0.5*(0.5+Wj1)*(0.5+Wj1);

            field[i1-1][j1-1].J1+=Fx1*Wy1*coeff[s];
            field[i1-1][j1].J1+=Fx1*Wy2*coeff[s];
            field[i1-1][j1+1].J1+=Fx1*Wy3*coeff[s];
            field[i1][j1-1].J1+=Fx2*Wy1*coeff[s];
            field[i1][j1].J1+=Fx2*Wy2*coeff[s];
            field[i1][j1+1].J1+=Fx2*Wy3*coeff[s];

            field[i1-1][j1-1].J2+=Fy1*Wx1*coeff[s];
            field[i1][j1-1].J2+=Fy1*Wx2*coeff[s];
            field[i1+1][j1-1].J2+=Fy1*Wx3*coeff[s];
            field[i1-1][j1].J2+=Fy2*Wx1*coeff[s];
            field[i1][j1].J2+=Fy2*Wx2*coeff[s];
            field[i1+1][j1].J2+=Fy2*Wx3*coeff[s];

            Wx1=0.5*(xcc-0.5)*(xcc-0.5);
            Wx2=0.75-xcc*xcc;
            Wx3=0.5*(0.5+xcc)*(0.5+xcc);
            Wy1=0.5*(ycc-0.5)*(ycc-0.5);
            Wy2=0.75-ycc*ycc;
            Wy3=0.5*(0.5+ycc)*(0.5+ycc);

            field[intXc-1][intYc-1].J3+=Wx1*Wy1*vz*coeff[s];
            field[intXc-1][intYc+0].J3+=Wx1*Wy2*vz*coeff[s];
            field[intXc-1][intYc+1].J3+=Wx1*Wy3*vz*coeff[s];
            field[intXc+0][intYc-1].J3+=Wx2*Wy1*vz*coeff[s];
            field[intXc+0][intYc+0].J3+=Wx2*Wy2*vz*coeff[s];
            field[intXc+0][intYc+1].J3+=Wx2*Wy3*vz*coeff[s];
            field[intXc+1][intYc-1].J3+=Wx3*Wy1*vz*coeff[s];
            field[intXc+1][intYc+0].J3+=Wx3*Wy2*vz*coeff[s];
            field[intXc+1][intYc+1].J3+=Wx3*Wy3*vz*coeff[s];

            p=p->next;
          }	//End of while(p)
        }	//End of for(s)     
     }		//End of for(i,j)

}
*/

void updateCurrent3D_DSX_1st(Domain *D)
{
    int i,j,k,s,n,i1,i2,j1,j2,k1,k2,intXc,intYc,intZc;
    int istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub;
    double inverDt,x1,x2,y1,y2,xr,yr,zr,z1,z2;
    double Fx1,Fx2,Wx1,Wx2,Fy1,Fy2,Wy1,Wy2,Fz1,Fz2,Wz1,Wz2,dx,dy,dz,dt;
    double vz,gamma,xc,yc,zc,wx,wy,wz,xcc,ycc,zcc;
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
              gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
              vz=p->p3/gamma;
              x2=p->x+i;       y2=p->y+j;       z2=p->z+k;
              x1=p->oldX;      y1=p->oldY;      z1=p->oldZ;
              xc=0.5*(x1+x2);  yc=0.5*(y1+y2);  zc=0.5*(z1+z2);
              intXc=(int)xc;   intYc=(int)yc;   intZc=(int)zc;
              xcc=xc-intXc;    ycc=yc-intYc;    zcc=zc-intZc;
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

              D->Jx[i2][j2][k2]    +=Fx1*(1-Wy1)*(1-Wz1)*coeff[s];
              D->Jx[i2][j2+1][k2]  +=Fx1*    Wy1*(1-Wz1)*coeff[s];
              D->Jx[i2][j2][k2+1]  +=Fx1*(1-Wy1)*    Wz1*coeff[s];
              D->Jx[i2][j2+1][k2+1]+=Fx1*    Wy1*    Wz1*coeff[s];
              D->Jy[i2][j2][k2]    +=Fy1*(1-Wx1)*(1-Wz1)*coeff[s];
              D->Jy[i2+1][j2][k2]  +=Fy1*    Wx1*(1-Wz1)*coeff[s];
              D->Jy[i2][j2][k2+1]  +=Fy1*(1-Wx1)*    Wz1*coeff[s];
              D->Jy[i2+1][j2][k2+1]+=Fy1*    Wx1*    Wz1*coeff[s];
              D->Jz[i2][j2][k2]    +=Fz1*(1-Wx1)*(1-Wy1)*coeff[s];
              D->Jz[i2+1][j2][k2]  +=Fz1*    Wx1*(1-Wy1)*coeff[s];
              D->Jz[i2][j2+1][k2]  +=Fz1*(1-Wx1)*    Wy1*coeff[s];
              D->Jz[i2+1][j2+1][k2]+=Fz1*    Wx1*    Wy1*coeff[s];

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
