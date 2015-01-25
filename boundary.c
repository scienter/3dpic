#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void boundary(Domain *D,External *Ext)
{
     int i,j,k,s,rank,rankX,rankY,rankZ;
     int remainX,remainY,remainZ,subX,subY,subZ,tmp;
     int nxSub,nySub,numdataUp,numdataBt,numberData;
     int startj,startk,nxSub1D,nySub2D,nzSub3D;
     int minX,maxX,minY,maxY,minZ,maxZ;
     int myrank, nTasks;
     float ***memoryAsign();
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     D->nxSub=D->nx/D->L;
     subX=D->nxSub;
     remainX=D->nx%D->L;
     minX=maxX=0;
     D->nySub=D->ny/D->M;
     subY=D->nySub;
     remainY=D->ny%D->M;
     minY=maxY=0;
     D->nzSub=D->nz/D->N;
     subZ=D->nzSub;
     remainZ=D->nz%D->N;
     minZ=maxZ=0;

     D->minXSub=0;
     D->maxXSub=D->nxSub;
     for(rankZ=0; rankZ<D->N; rankZ++)
     {
       minY=maxY=0;
       for(rankY=0; rankY<D->M; rankY++)
       {
          rank=rankY+(rankZ*D->M);
          if(rankY<remainY)   tmp=subY+1;
          else                tmp=subY;
          minY=maxY;
          maxY=minY+tmp;
          if(myrank==rank)
          {
             D->minYSub=minY;
             D->maxYSub=maxY;
             D->nySub=tmp;
          }
       }
     }

     for(rankY=0; rankY<D->M; rankY++)
     {
       minZ=maxZ=0;
       for(rankZ=0; rankZ<D->N; rankZ++)
       {
          rank=rankY+(rankZ*D->M);
          if(rankZ<remainZ)   tmp=subZ+1;
          else                tmp=subZ;
          minZ=maxZ;
          maxZ=minZ+tmp;
          if(myrank==rank)
          {
             D->minZSub=minZ;
             D->maxZSub=maxZ;
             D->nzSub=tmp;
          }
       }
     }
     

printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d,minZSub=%d,maxZSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub,D->minZSub,D->maxZSub);

     D->istart=2;
     D->iend=D->nxSub+2;
     D->jstart=0;
     D->jend=1;
     D->kstart=0;
     D->kend=1;
     if(D->dimension>1)  {
        D->jstart=2;
        D->jend=D->nySub+2;
     }
     if(D->dimension>2)  {
        D->kstart=2;
        D->kend=D->nzSub+2;
     }

     // Field setting
     nxSub1D=D->nxSub+5;
     nySub2D=1;
     nzSub3D=1;
     if(D->dimension>1)  
       nySub2D=D->nySub+5;
     if(D->dimension>2)  
       nzSub3D=D->nzSub+5;

     if(D->fieldType==1)
     {
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ExC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->JxOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->JyOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->JzOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
     }

     // Particle setting
     nxSub1D=D->nxSub+3;
     nySub2D=1;
     nzSub3D=1;
     startj=0;
     startk=0;
     if(D->dimension>1) {
       nySub2D=D->nySub+3;
       startj=D->jstart-1;
     }
     if(D->dimension>2) { 
       nzSub3D=D->nzSub+3;
       startk=D->kstart-1;
     }

     D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
     for(i=0; i<nxSub1D; i++) {
       D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
       for(j=startj; j<nySub2D; j++) 
         D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
     }

     // setting up particle's pointer
     for(i=0; i<nxSub1D; i++)	//i starts at 0 because of boost frame
       for(j=startj; j<nySub2D; j++)
         for(k=startk; k<nzSub3D; k++)
         {
           D->particle[i][j][k].rho=0.0;  
           D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
           for(s=0; s<D->nSpecies; s++)
           {
             D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
             D->particle[i][j][k].head[s]->pt = NULL;
           }
         }
   
    // current J trasffering boundary
    numberData=3*3*(D->nx+5)*(D->nzSub+5);
    D->YplusJ=(double *)malloc(numberData*sizeof(double ));
    numberData=2*3*(D->nx+5)*(D->nzSub+5);
    D->YminusJ=(double *)malloc(numberData*sizeof(double ));
    numberData=3*3*(D->nx+5)*(D->nySub+5);
    D->ZplusJ=(double *)malloc(numberData*sizeof(double ));
    numberData=2*3*(D->nx+5)*(D->nySub+5);
    D->ZminusJ=(double *)malloc(numberData*sizeof(double ));
/*
     D->probe = (Probe **)malloc(D->probeNum*sizeof(Probe *));
     for(i=0; i<D->probeNum; i++)
       D->probe[i] = (Probe *)malloc((D->maxStep+1)*sizeof(Probe ));
     for(i=0; i<D->probeNum; i++)
       for(j=0; j<=D->maxStep; j++)
       {
         D->probe[i][j].Ex=0.0;
         D->probe[i][j].Ey=0.0;
         D->probe[i][j].Ez=0.0;
         D->probe[i][j].Bx=0.0;
         D->probe[i][j].By=0.0;
         D->probe[i][j].Bz=0.0;
       }
*/
    //Share Field
    switch(D->dimension) {
    case 1 :
    break;
    case 3 :
      numberData=3*(D->nx+5)*(D->nzSub+5)*1;
      D->plusYC = (float *)malloc(numberData*sizeof(float ));
      numberData=3*(D->nx+5)*(D->nzSub+5)*1;
      D->minusYC = (float *)malloc(numberData*sizeof(float ));
      numberData=3*(D->nx+5)*(D->nySub+5)*1;
      D->plusZC = (float *)malloc(numberData*sizeof(float ));
      numberData=3*(D->nx+5)*(D->nySub+5)*1;
      D->minusZC = (float *)malloc(numberData*sizeof(float ));

      numberData=6*(D->nx+5)*(D->nzSub+5)*2;
      D->plusY = (float *)malloc(numberData*sizeof(float ));
      numberData=6*(D->nx+5)*(D->nzSub+5)*3;
      D->minusY = (float *)malloc(numberData*sizeof(float ));
      numberData=6*(D->nx+5)*(D->nySub+5)*2;
      D->plusZ = (float *)malloc(numberData*sizeof(float ));
      numberData=6*(D->nx+5)*(D->nySub+5)*3;
      D->minusZ = (float *)malloc(numberData*sizeof(float ));
    break;
    }    

}

float ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   float ***field;

   field = (float ***)malloc((nx)*sizeof(float **));
   for(i=0; i<nx; i++)
   {
     field[i] = (float **)malloc((ny)*sizeof(float *));
     for(j=0; j<ny; j++)
       field[i][j] = (float *)malloc((nz)*sizeof(float ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0.0;
       }

   return field;
}
