#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void loadLaser3D(Domain *D,LaserList *L,double t)
{
   float rU,rD,longitudinal,t0,flat;
   float zR,w0,w,phi,omega,kx,pphi,amp;
   float x,y,z,r2,w2;
   int istart,iend,jstart,jend,kstart,kend;
   int positionX,rank,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=(int)(D->ny*0.5);	//center position
   kC=(int)(D->nz*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*pi*L->omega/D->omega;
   kx=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;


   if(D->fieldType==1 && laserOK==1)
   {
     if(L->polarity==2)
     {
//       D->Pr[positionX][jC-D->minYSub][kC-D->minZSub]=longitudinal*sin(omega*t);  
//       D->Pl[positionX][jC-D->minYSub][kC-D->minZSub]=longitudinal*sin(omega*t);  

       w2=w*w;
       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=kstart; k<kend; k++)
         {
           z=(k-kstart+D->minZSub-kC)*D->dz;
           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           D->Pr[positionX][j][k]=amp;            
           D->Pl[positionX][j][k]=amp;           
         } 
       }
     }  
/*         
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
*/
   }     //End of fieldType=1
}


/*
void boostLoadLaser2D(Domain *D,LaserList *L)
{
   float x,x0,rU,rD,beta,gamma,labWaist,f,yy;
   float zR,w0,z,w,phi,omega,k,y,pphi,amp,unitCompen;
   int i,j,jC;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;

   jC=(int)(D->ny*0.5);	//center position
   gamma=D->gamma;
   beta=D->beta;

   labWaist=L->beamWaist*D->lambda;
   zR=L->rayleighLength/gamma;
//   zR=pi/(L->lambda/D->gamma/(1.0+D->beta))*labWaist*labWaist/gamma/D->lambda;
   w0=labWaist/D->lambda;  
   f=-L->focus/gamma;   
   omega=k=2*pi;
   x0=-2*rU;
   unitCompen=gamma*gamma*(1.0-beta*beta);

   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2)
     {
       for(i=2; i<D->nxSub+2; i++)
         for(j=2; j<D->nySub+2; j++)
         {
           x=(i+D->minXSub-2)*D->dx;
           y=(j-2+D->minYSub-jC)*D->dy;
           z=(f+x+x0);          
           w=w0*sqrt(1.0+z*z/zR/zR);
           if(x>=2*x0)
           {
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
//             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
             field[i][j].Pr=amp;            
//             field[i][j].Pl=amp;            
           }
           else      {
             field[i][j].Pr=0.0;            
//             field[i][j].Pl=0.0;            
           }
         }
     }
     else if(L->polarity==3)
     {
       for(i=1; i<=D->nxSub; i++)
         for(j=1; j<=D->nySub; j++)
         {
           x=(i+D->minXSub-1)*D->dx;
           z=f+x+x0;          
           if(x>=2*x0)
           {
             y=(j-1+D->minYSub-jC)*D->dy;
             w=w0*sqrt(1.0+z*z/zR/zR);
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
             field[i][j].Sr=amp;            
             field[i][j].Sl=amp;            
           }
           else      {
             field[i][j].Sr=0.0;            
             field[i][j].Sl=0.0;            
           }
         }
     }
   }     //End of fieldType=1
}
*/
/*
void loadLaser2D(Domain *D,LaserList *L,double t)
{
   float rU,rD,longitudinal,t0,flat;
   float zR,w0,z,w,phi,omega,k,y,pphi,amp;
   int istart,iend,jstart,jend;
   int positionX,positionY,rank,j,jC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=(int)(D->ny*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   z=-L->focus;   
   w=w0*sqrt(1.0+z*z/zR/zR);
   phi=atan(z/zR);
   omega=2*pi*L->omega/D->omega;
   k=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   positionY=L->loadPointY+jstart-D->minYSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
      laserOK=1;


   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2 && laserOK==1)
     {
//         field[positionX][positionY].Pr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Pl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Pr=amp;            
         field[positionX][j].Pl=amp;
       }

     }           
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
   }     //End of fieldType=1

}
*/
/*
void loadLaserOpp2D(Domain *D,LaserList *L,double t)
{
   float rU,rD,longitudinal,t0,flat,loadPosition;
   float zR,w0,z,w,phi,omega,k,y,pphi,amp;
   int istart,iend,jstart,jend;
   int positionX,positionY,rank,j,jC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=(int)(D->ny*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   loadPosition=L->loadPointX*D->dx;
   z=loadPosition-L->focus;   
   w=w0*sqrt(1.0+z*z/zR/zR);
   phi=atan(z/zR);
   omega=2*pi*L->omega/D->omega;
   k=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub-2;
   positionY=L->loadPointY+jstart-D->minYSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
      laserOK=1;


   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2 && laserOK==1)
     {
//         field[positionX][positionY].Pr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Pl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=-z/zR*y*y/w/w-0.5*phi-k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Pr=amp;            
         field[positionX][j].Pl=amp;
       }

     }           
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=-z/zR*y*y/w/w-0.5*phi-k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
   }     //End of fieldType=1

}
*/
