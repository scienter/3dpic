#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
/*
void removeEdge2DBoost(Domain *D)
{
    int i,j,istart,iend,jstart,jend,s;
    Particle **particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //remove left boundary
    for(i=0; i<istart; i++)
    {
      for(s=0; s<D->nSpecies; s++)
        for(j=jstart-1; j<=jend; j++)
        { 
          p=particle[i][j].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
        }
    }
    

    //remove right boundary
    i=iend;
    for(s=0; s<D->nSpecies; s++)
      for(j=jstart-1; j<=jend; j++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

    //remove up boundary
    j=jend;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

    //remove bottom boundary
    j=jstart-1;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

}
*/

void removeEdge3D(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s;
    Particle ***particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    //remove left boundary
    for(i=0; i<istart-1; i++)
    {
      for(s=0; s<D->nSpecies; s++)
        for(j=jstart-1; j<=jend; j++)
          for(k=kstart-1; k<=kend; k++)
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
          }
    }

    //remove right boundary
    i=iend;
    for(s=0; s<D->nSpecies; s++)
      for(j=jstart-1; j<=jend; j++)
        for(k=kstart-1; k<=kend; k++)
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
        }

    //remove plusZ boundary
    k=kend;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
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
        }

    //remove minusZ boundary
    k=kstart-1;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
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
        }

    //remove plusY boundary
    j=jend;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
        for(k=kstart-1; k<=kend; k++)
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
        }

    //remove minusY boundary
    j=jstart-1;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
        for(k=kstart-1; k<=kend; k++)
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
        }

}
