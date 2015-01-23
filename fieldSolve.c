#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"

void solveField3DC_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    float dx,dy,dz,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(k=kstart; k<kend; k++)
      for(j=jstart; j<jend; j++)
      {
        nowPrC=D->PrC[istart-1][j][k];
        nowSrC=D->SrC[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]+=dt/dy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])+dt/dz*(D->Sr[i][j][k]-D->Sr[i][j][k-1]-D->Sl[i][j][k]+D->Sl[i][j][k-1])-pi*dt*(D->JxOld[i][j][k]+D->Jx[i][j][k]);
          D->BxC[i][j][k]+=dt/dz*(D->Pr[i][j][k+1]-D->Pr[i][j][k]+D->Pl[i][j][k+1]-D->Pl[i][j][k])-dt/dy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
          prevPrC=nowPrC;
          nowPrC=D->PrC[i][j][k];
          D->PrC[i][j][k]=prevPrC+0.25*dt/dz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])+0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          D->PlC[i-1][j][k]=D->PlC[i][j][k]+0.25*dt/dz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])-0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          prevSrC=nowSrC;
          nowSrC=D->SrC[i][j][k];
          D->SrC[i][j][k]=prevSrC-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])+0.25*dt/dz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
          D->SlC[i-1][j][k]=D->SlC[i][j][k]-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-0.25*dt/dz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
        }	//End of i
      }		//End of j,k
}

void solveField3D_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    float dx,dy,dz,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    // Pr,Pl,E1,Sr,Sl,B1
    for(k=kstart; k<kend; k++)
      for(j=jstart; j<jend; j++)
      {
        nowPr=D->Pr[istart-1][j][k];
        nowSr=D->Sr[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]+=dt/dy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])+dt/dz*(D->SrC[i][j][k]-D->SrC[i][j][k-1]-D->SlC[i][j][k]+D->SlC[i][j][k-1])-2*pi*dt*D->Jx[i][j][k];
          D->Bx[i][j][k]+=dt/dz*(D->PrC[i][j][k+1]-D->PrC[i][j][k]+D->PlC[i][j][k+1]-D->PlC[i][j][k])-dt/dy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]);
          prevPr=nowPr;
          nowPr=D->Pr[i][j][k];
          D->Pr[i][j][k]=prevPr+0.25*dt/dz*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j][k-1]-D->BxC[i-1][j][k-1])+0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jy[i][j][k];
          D->Pl[i-1][j][k]=D->Pl[i][j][k]+0.25*dt/dz*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j][k-1]-D->BxC[i-1][j][k-1])-0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jy[i][j][k];
          prevSr=nowSr;
          nowSr=D->Sr[i][j][k];
          D->Sr[i][j][k]=prevSr-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])+0.25*dt/dz*(D->ExC[i][j][k+1]+D->ExC[i-1][j][k+1]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jz[i][j][k];
          D->Sl[i-1][j][k]=D->Sl[i][j][k]-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-0.25*dt/dz*(D->ExC[i][j][k+1]+D->ExC[i-1][j][k+1]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jz[i][j][k];
        }	//End of i
      }		//End of j,k
}

/*
void solveField2D_DSX(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    FieldDSX **field;
    field=D->fieldDSX;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    // Pr,Pl,E1,Sr,Sl,B1
    for(j=jstart; j<jend; j++)
    {
      nowPr=field[istart-1][j].Pr;
      nowSr=field[istart-1][j].Sr;
      for(i=istart; i<iend; i++)
      {
         field[i][j].E1+=0.5*dt/dy*(field[i-1][j].PrC-field[i-1][j].PlC+field[i][j].PrC-field[i][j].PlC-field[i-1][j-1].PrC+field[i-1][j-1].PlC-field[i][j-1].PrC+field[i][j-1].PlC)-pi*dt*(field[i-1][j].J1+field[i][j].J1);
         field[i][j].B1+=-0.5*dt/dy*(field[i-1][j].SrC+field[i-1][j].SlC+field[i][j].SrC+field[i][j].SlC-field[i-1][j-1].SrC-field[i-1][j-1].SlC-field[i][j-1].SrC-field[i][j-1].SlC);

         prevSr=nowSr;
         nowSr=field[i][j].Sr;
         field[i][j].Sr=prevSr-0.5*dt/dy*(field[i][j+1].B1C-field[i][j].B1C)-0.5*pi*dt*(field[i][j].J3+field[i][j+1].J3);
         field[i-1][j].Sl=field[i][j].Sl-0.5*dt/dy*(field[i][j+1].B1C-field[i][j].B1C)-0.5*pi*dt*(field[i][j].J3+field[i][j+1].J3);     
         prevPr=nowPr;
         nowPr=field[i][j].Pr;
         field[i][j].Pr=prevPr+0.5*dt/dy*(field[i][j+1].E1C-field[i][j].E1C)-pi*dt*field[i][j].J2;
         field[i-1][j].Pl=field[i][j].Pl-0.5*dt/dy*(field[i][j+1].E1C-field[i][j].E1C)-pi*dt*field[i][j].J2;
      }
    }
}
*/
