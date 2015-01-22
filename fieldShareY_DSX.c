#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
/*//
void MPI_TransferF_XplusFilter(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    behindF=(float *)malloc(numberData*sizeof(float )); 
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}
*/


void MPI_TransferJ_DSX_Yplus(Domain *D)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    //first 3 is 3rd interpolation, 2nd 3 is Jx,Jy,Jz;
    numberData=3*3*(nx+5)*(nzSub+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank%D->M;

    //Transferring even ~ odd cores 
    if(rank%2==1)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=0; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=D->M-1)
    {
      start=0; 
      for(j=0; j<3; j++)
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jx[i][jend+j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jy[i][jend+j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jz[i][jend+j][k];
          start+=ibottom;
        }
        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=0; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jstart+j][k]+=D->upJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1 && rank!=D->M-1)
    {
      start=0; 
      for(j=0; j<3; j++)
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jx[i][jend+j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jy[i][jend+j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->upJ[start+i]=D->Jz[i][jend+j][k];
          start+=ibottom;
        }        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }     
    MPI_Barrier(MPI_COMM_WORLD);        
}


void MPI_TransferJ_DSX_Yminus(Domain *D)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank; 

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    //first 2 is 3rd interpolation, 2nd 3 is Jx,Jy,Jz.
    numberData=2*3*(nx+5)*(nzSub+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank%D->M;

    //Transferring even ~ odd cores 
    if(rank%2==1 && rank!=D->M-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==0 && rank!=0)
    {
      start=0; 
      for(j=1; j<3; j++)
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jx[i][jstart-j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jy[i][jstart-j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jz[i][jstart-j][k];
          start+=ibottom;
        }
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(rank%2==0 && rank!=D->M-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<3; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             D->Jx[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jy[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
           for(i=ibegin; i<ibottom; i++)
             D->Jz[i][jend-j][k]+=D->btJ[i+start];
           start+=ibottom;
         }
    }
    else if(rank%2==1)
    {
      start=0; 
      for(j=1; j<3; j++)
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jx[i][jstart-j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jy[i][jstart-j][k];
          start+=ibottom;
          for(i=ibegin; i<ibottom; i++)
            D->btJ[start+i]=D->Jz[i][jstart-j][k];
          start+=ibottom;
        }
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferF_DSX_Yminus(Domain *D,float ***f1,float ***f2,float ***f3)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         

    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    ibegin=0;
    ibottom=nx+5;   

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3*(nx+5)*(nzSub+5)*3; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=0; j<3; j++)	// 3 layers in y direction
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f1[i][jstart+j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f2[i][jstart+j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f3[i][jstart+j][k];
        start+=nx+5;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<3; j++)	// 3 layers in y direction
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusY,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=0; j<3; j++)	// 3 layers in y direction
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f1[i][jstart+j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f2[i][jstart+j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusY[start+i]=f3[i][jstart+j][k];
        start+=nx+5;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
      start=0; 
      for(j=0; j<3; j++)	// 3 layers in y direction
        for(k=0; k<nzSub+5; k++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f2[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f3[i][jend+j][k]=D->minusY[start+i];
          start+=nx+5;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusY,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferF_DSX_Yplus(Domain *D,float ***f1,float ***f2,float ***f3)
{
    int i,j,k,numberData,start,end,nx,nySub,nzSub,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         
   
    nx=D->nx;
    nySub=D->nySub;
    nzSub=D->nzSub;
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3*(nx+5)*(nzSub+5)*2; 
    rank=myrank%D->M;   

    //Transferring even ~ odd cores 
    start=0; 
    for(j=1; j<=2; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f1[i][jend-j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f2[i][jend-j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f3[i][jend-j][k];
        start+=nx+5;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusY,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<=2; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
         }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(D->plusY,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(j=1; j<=2; j++)
      for(k=0; k<nzSub+5; k++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f1[i][jend-j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f2[i][jend-j][k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusY[start+i]=f3[i][jend-j][k];
        start+=nx+5;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusY,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(j=1; j<=2; j++)
         for(k=0; k<nzSub+5; k++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f2[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f3[i][jstart-j][k]=D->plusY[start+i];
           start+=nx+5;
         }
    }
    else if(rank%2==1 && rank!=D->M-1)
       MPI_Send(D->plusY,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}
