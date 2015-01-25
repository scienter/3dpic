#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

// Temporary routine to dump and restore data
void saveDump3D(Domain D,int iteration)
{
   FILE *out;
   char name[100];
   int i,j,k,istart,iend,jstart,jend,kstart,kend,n,s,cnt;
   Particle ***particle;
   particle=D.particle;
   Probe **probe;
   probe=D.probe;
   ptclList *p;
   int myrank, nprocs;    

   istart=D.istart;
   iend=D.iend;
   jstart=D.jstart;
   jend=D.jend;
   kstart=D.kstart;
   kend=D.kend;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   out = fopen(name, "w");   

   // Save simulation Domain information
   fwrite(&(D.minXSub),sizeof(int),1,out);
   fwrite(&(D.maxXSub),sizeof(int),1,out);
   fwrite(&(D.minYSub),sizeof(int),1,out);
   fwrite(&(D.maxYSub),sizeof(int),1,out);
   fwrite(&(D.minZSub),sizeof(int),1,out);
   fwrite(&(D.maxZSub),sizeof(int),1,out);
   fwrite(&(D.maxStep),sizeof(int),1,out);

   // Save informations of particles inside the domain
   for(s=0; s<D.nSpecies; s++)
     for (i=istart; i<iend; i++) 
       for (j=jstart; j<jend; j++)
         for (k=kstart; k<kend; k++)
         {
           p = particle[i][j][k].head[s]->pt;
           cnt = 0;
           while(p)  { 
             p=p->next; 
             cnt++;    
           }
           fwrite(&cnt,sizeof(int),1,out);

           p = particle[i][j][k].head[s]->pt;
           while(p)  
           { 
             fwrite(&(p->x),sizeof(float),1,out);   
             fwrite(&(p->y),sizeof(float),1,out);   
             fwrite(&(p->z),sizeof(float),1,out);   
             fwrite(&(p->p1),sizeof(float),1,out);   
             fwrite(&(p->p2),sizeof(float),1,out);   
             fwrite(&(p->p3),sizeof(float),1,out);   
             fwrite(&(p->index),sizeof(float),1,out);   

             p = p->next; 
           }            
         }		//End of for(i,j)

   if(D.fieldType==1)
   {
     for(i=0; i<D.nxSub+5; i++) 
       for(j=0; j<D.nySub+5; j++)
         for(k=0; k<D.nzSub+5; k++)
         {
           fwrite(&(D.Ex[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Pr[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Pl[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Bx[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Sr[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Sl[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.ExC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.PrC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.PlC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.BxC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.SrC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.SlC[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Jx[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Jy[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.Jz[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.JxOld[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.JyOld[i][j][k]),sizeof(float),1,out);   
           fwrite(&(D.JzOld[i][j][k]),sizeof(float),1,out);   
         }
   }
/*  
   //save Probe data
   if(D.probeNum>0)
   {
     for(n=0; n<D.probeNum; n++)
       for(i=0; i<=D.maxStep; i++)
       { 
         fwrite(&(probe[n][i].Pr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Pl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].E1),sizeof(float),1,out);   
         fwrite(&(probe[n][i].B1),sizeof(float),1,out);   
       }
   }
*/
   fclose(out);
}


// Temporary routine to dump and restore data
void restoreData3D(Domain *D, int iteration)
{
   FILE *in;
   char name[100];
   int i,j,k,istart,iend,jstart,jend,kstart,kend,n,s,cnt,maxStep;
   float tmp;
   double tmp1;
   ptclList *p;
   Particle ***particle;
   particle=D->particle;
   Probe **probe;
   probe=D->probe;
   int myrank, nTasks;    

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   in = fopen(name, "r");   

   // restore simulation Domain information
   fread(&(D->minXSub),sizeof(int),1,in);
   fread(&(D->maxXSub),sizeof(int),1,in);
   fread(&(D->minYSub),sizeof(int),1,in);
   fread(&(D->maxYSub),sizeof(int),1,in);
   fread(&(D->minZSub),sizeof(int),1,in);
   fread(&(D->maxZSub),sizeof(int),1,in);
   fread(&(maxStep),sizeof(int),1,in);

   // restore informations of particles inside the domain
   for(s=0; s<D->nSpecies; s++)
     for (i=istart; i<iend; i++) 
       for (j=jstart; j<jend; j++)
         for (k=kstart; k<kend; k++)
         {
           fread(&cnt,sizeof(int),1,in);

           for(n=0; n<cnt; n++)  
           { 
             p = (ptclList *)malloc(sizeof(ptclList)); 
             p->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = p;

             fread(&(p->x),sizeof(float),1,in);   
             fread(&(p->y),sizeof(float),1,in);   
             fread(&(p->z),sizeof(float),1,in);   
             fread(&(p->p1),sizeof(float),1,in);   
             fread(&(p->p2),sizeof(float),1,in);   
             fread(&(p->p3),sizeof(float),1,in);   
             fread(&(p->index),sizeof(float),1,in); 
           }
         }

   if(D->fieldType==1)
   {
     for(i=0; i<D->nxSub+5; i++) 
       for(j=0; j<D->nySub+5; j++)
         for(k=0; k<D->nzSub+5; k++)
         {
           fread(&(D->Ex[i][j][k]),sizeof(float),1,in);
           fread(&(D->Pr[i][j][k]),sizeof(float),1,in);
           fread(&(D->Pl[i][j][k]),sizeof(float),1,in);
           fread(&(D->Bx[i][j][k]),sizeof(float),1,in);
           fread(&(D->Sr[i][j][k]),sizeof(float),1,in);
           fread(&(D->Sl[i][j][k]),sizeof(float),1,in);
           fread(&(D->ExC[i][j][k]),sizeof(float),1,in);
           fread(&(D->PrC[i][j][k]),sizeof(float),1,in);
           fread(&(D->PlC[i][j][k]),sizeof(float),1,in);
           fread(&(D->BxC[i][j][k]),sizeof(float),1,in);
           fread(&(D->SrC[i][j][k]),sizeof(float),1,in);
           fread(&(D->SlC[i][j][k]),sizeof(float),1,in);
           fread(&(D->Jx[i][j][k]),sizeof(float),1,in);
           fread(&(D->Jy[i][j][k]),sizeof(float),1,in);
           fread(&(D->Jz[i][j][k]),sizeof(float),1,in);
           fread(&(D->JxOld[i][j][k]),sizeof(float),1,in);
           fread(&(D->JyOld[i][j][k]),sizeof(float),1,in);
           fread(&(D->JzOld[i][j][k]),sizeof(float),1,in);
         }
   }       //End of fieldType=1
/*
   //restore Probe data
   if(D->probeNum>0)
   {
     for(n=0; n<D->probeNum; n++)
       for(i=0; i<=maxStep; i++)
       { 
         fread(&(probe[n][i].Pr),sizeof(float),1,in);   
         fread(&(probe[n][i].Pl),sizeof(float),1,in);   
         fread(&(probe[n][i].Sr),sizeof(float),1,in);   
         fread(&(probe[n][i].Sl),sizeof(float),1,in);   
         fread(&(probe[n][i].E1),sizeof(float),1,in);   
         fread(&(probe[n][i].B1),sizeof(float),1,in);   
       }
   }
*/   
   fclose(in);
}

