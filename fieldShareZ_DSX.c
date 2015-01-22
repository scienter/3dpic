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

/*
void MPI_TransferJ_DSX_Yplus(Domain *D)
{
    int i,n,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    FieldDSX **field;
    field=D->fieldDSX;

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    numberData=3*3*(nx+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //Transferring even ~ odd cores 
    if(myrank%2==1)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<3; n++)
       {
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J1+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J2+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J3+=D->upJ[i+start];
         start+=ibottom;
       }
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
    {
      start=0; 
      for(n=0; n<3; n++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J1;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J2;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J3;
        start+=ibottom;
      }
        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<3; n++)
       {
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J1+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J2+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J3+=D->upJ[i+start];
         start+=ibottom;
       }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
    {
      start=0; 
      for(n=0; n<3; n++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J1;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J2;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J3;
        start+=ibottom;
      }        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }     
    MPI_Barrier(MPI_COMM_WORLD);        
}
*/
/*
void MPI_TransferJ_DSX_Yminus(Domain *D)
{
    int i,n,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    FieldDSX **field;
    field=D->fieldDSX;

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    numberData=8*(nx+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //Transferring even ~ odd cores 
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3Old=D->btJ[i+start];
    }
    else if(myrank%2==0 && myrank!=0)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3Old;
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3Old=D->btJ[i+start];
    }
    else if(myrank%2==1)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3Old;
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);
}
*/

void MPI_TransferF_DSX_Zminus(Domain *D,float ***f1,float ***f2,float ***f3)
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

    numberData=3*(nx+5)*(nySub+5)*3; 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=0; k<3; k++)	// 3 layers in z direction
      for(j=0; j<nySub+5; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f1[i][j][kstart+k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f2[i][j][kstart+k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f3[i][j][kstart+k];
        start+=nx+5;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZ,numberData, MPI_FLOAT, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<3; k++)	// 3 layers in z direction
        for(j=0; j<nySub+5; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusZ,numberData, MPI_FLOAT, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=0; k<3; k++)	// 3 layers in z direction
      for(j=0; j<nySub+5; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f1[i][j][kstart+k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f2[i][j][kstart+k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->minusZ[start+i]=f3[i][j][kstart+k];
        start+=nx+5;
      }
        
    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(D->minusZ,numberData, MPI_FLOAT, myrank+D->M, myrank+D->M, MPI_COMM_WORLD,&status);  
      start=0; 
      for(k=0; k<3; k++)	// 3 layers in z direction
        for(j=0; j<nySub+5; j++)
        {
          for(i=ibegin; i<ibottom; i++)
            f1[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f2[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
          for(i=ibegin; i<ibottom; i++)
            f3[i][j][kend+k]=D->minusZ[start+i];
          start+=nx+5;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusZ,numberData, MPI_FLOAT, myrank-D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferF_DSX_Zplus(Domain *D,float ***f1,float ***f2,float ***f3)
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

    numberData=3*(nx+5)*(nySub+5)*2; 
    rank=(int)(myrank/D->M);   

    //Transferring even ~ odd cores 
    start=0; 
    for(k=1; k<=2; k++)
      for(j=0; j<nySub+5; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f1[i][j][kend-k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f2[i][j][kend-k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f3[i][j][kend-k];
        start+=nx+5;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusZ,numberData, MPI_FLOAT, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<=2; k++)
         for(j=0; j<nySub+5; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
         }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(D->plusZ,numberData, MPI_FLOAT, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(k=1; k<=2; k++)
      for(j=0; j<nySub+5; j++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f1[i][j][kend-k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f2[i][j][kend-k];
        start+=nx+5;
        for(i=ibegin; i<ibottom; i++)
          D->plusZ[start+i]=f3[i][j][kend-k];
        start+=nx+5;
      }
        
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->plusZ,numberData, MPI_FLOAT, myrank-D->M, myrank-D->M, MPI_COMM_WORLD,&status);  
       start=0;
       for(k=1; k<=2; k++)
         for(j=0; j<nySub+5; j++)
         {
           for(i=ibegin; i<ibottom; i++)
             f1[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f2[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
           for(i=ibegin; i<ibottom; i++)
             f3[i][j][kstart-k]=D->plusZ[start+i];
           start+=nx+5;
         }
    }
    else if(rank%2==1 && rank!=D->N-1)
       MPI_Send(D->plusZ,numberData, MPI_FLOAT, myrank+D->M, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}