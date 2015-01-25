#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void clean3D(Domain *D)
{
    Particle ***particle;
    particle=D->particle;
    LoadList *LL,*tmpLL;
    LaserList *L, *tmpL;
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s;
    ptclList *p,*tmp;
 
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    //remove particles
    for(i=0; i<=iend; i++)
      for(j=jstart-1; j<=jend; j++)
        for(k=kstart-1; k<=kend; k++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {	
              tmp=p->next;
              particle[i][j][k].head[s]->pt=tmp; 
              p->next=NULL;
              free(p);
              p=particle[i][j][k].head[s]->pt;
            }
            free(particle[i][j][k].head[s]);
          }
          free(particle[i][j][k].head);
        }		
      
    LL=D->loadList;
    while(LL)
    {	
      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;
      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);

    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);

    //remove trans
    free(D->YplusJ);
    free(D->YminusJ);
    free(D->ZplusJ);
    free(D->ZminusJ);
    free(D->plusYC);
    free(D->minusYC);
    free(D->plusY);
    free(D->minusY);
    free(D->plusZC);
    free(D->minusZC);
    free(D->plusZ);
    free(D->minusZ);


    //remove field
    if(D->fieldType==1)
    {
      for(i=0; i<D->nxSub+5; i++)
      {
        for(j=0; j<D->nySub+5; j++)
        {
          free(D->Ex[i][j]);
          free(D->Pr[i][j]);
          free(D->Pl[i][j]);
          free(D->Bx[i][j]);
          free(D->Sr[i][j]);
          free(D->Sl[i][j]);
          free(D->ExC[i][j]);
          free(D->PrC[i][j]);
          free(D->PlC[i][j]);
          free(D->BxC[i][j]);
          free(D->SrC[i][j]);
          free(D->SlC[i][j]);
          free(D->Jx[i][j]);
          free(D->Jy[i][j]);
          free(D->Jz[i][j]);
          free(D->JxOld[i][j]);
          free(D->JyOld[i][j]);
          free(D->JzOld[i][j]);
        }
        free(D->Ex[i]);
        free(D->Pr[i]);
        free(D->Pl[i]);
        free(D->Bx[i]);
        free(D->Sr[i]);
        free(D->Sl[i]);
        free(D->ExC[i]);
        free(D->PrC[i]);
        free(D->PlC[i]);
        free(D->BxC[i]);
        free(D->SrC[i]);
        free(D->SlC[i]);
        free(D->Jx[i]);
        free(D->Jy[i]);
        free(D->Jz[i]);
        free(D->JxOld[i]);
        free(D->JyOld[i]);
        free(D->JzOld[i]);
      }
      free(D->Ex);
      free(D->Pr);
      free(D->Pl);
      free(D->Bx);
      free(D->Sr);
      free(D->Sl);
      free(D->ExC);
      free(D->PrC);
      free(D->PlC);
      free(D->BxC);
      free(D->SrC);
      free(D->SlC);
      free(D->Jx);
      free(D->Jy);
      free(D->Jz);
      free(D->JxOld);
      free(D->JyOld);
      free(D->JzOld);
    }		//End if(fieldType=1)
/*
    //remove probe
    if(D->probeNum>0)
    {
      for(i=0; i<D->probeNum; i++)
       free(D->probe[i]);
      free(D->probe);
      free(D->probeX);
      free(D->probeY);
    }
*/
/*
    //remove boost field
    Boost **boost;
    boost=D->boost;

    for(i=0; i<D->nxSub+5; i++)
      free(boost[i]);
    free(boost);
*/
}

