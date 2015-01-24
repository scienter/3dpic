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
          D->ExC[i][j][k]+=0.5*dt/dy*(D->Pr[i][j][k]+D->Pr[i-1][j][k]-D->Pr[i][j-1][k]-D->Pr[i-1][j-1][k]-D->Pl[i][j][k]-D->Pl[i-1][j][k]+D->Pl[i][j-1][k]+D->Pl[i-1][j-1][k])
                          +0.5*dt/dz*(D->Sr[i][j][k]+D->Sr[i-1][j][k]-D->Sr[i][j][k-1]-D->Sr[i-1][j][k-1]-D->Sl[i][j][k]-D->Sl[i-1][j][k]+D->Sl[i][j][k-1]+D->Sl[i-1][j][k-1])
                          -0.5*pi*dt*(D->JxOld[i][j][k]+D->JxOld[i-1][j][k]+D->Jx[i][j][k]+D->Jx[i-1][j][k]);
          D->BxC[i][j][k]+=0.5*dt/dz*(D->Pr[i][j][k+1]+D->Pr[i-1][j][k+1]-D->Pr[i][j][k]-D->Pr[i-1][j][k]+D->Pl[i][j][k+1]+D->Pl[i-1][j][k+1]-D->Pl[i][j][k]-D->Pl[i-1][j][k])
                          -0.5*dt/dy*(D->Sr[i][j+1][k]+D->Sr[i-1][j+1][k]-D->Sr[i][j][k]-D->Sr[i-1][j][k]+D->Sl[i][j+1][k]+D->Sl[i-1][j+1][k]-D->Sl[i][j][k]-D->Sl[i-1][j][k]);
          prevPrC=nowPrC;
          nowPrC=D->PrC[i][j][k];
          D->PrC[i][j][k]=prevPrC
                         +0.5*dt/dz*(D->Bx[i][j][k]-D->Bx[i][j][k-1])
                         +0.5*dt/dy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])
                         -0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          D->PlC[i-1][j][k]=D->PlC[i][j][k]
                         +0.5*dt/dz*(D->Bx[i][j][k]-D->Bx[i][j][k-1])
                         -0.5*dt/dy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])
                         -0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          prevSrC=nowSrC;
          nowSrC=D->SrC[i][j][k];
          D->SrC[i][j][k]=prevSrC
                         -0.5*dt/dy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])
                         +0.5*dt/dz*(D->Ex[i][j][k+1]-D->Ex[i][j][k])
                         -0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
          D->SlC[i-1][j][k]=D->SlC[i][j][k]
                         -0.5*dt/dy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])
                         -0.5*dt/dz*(D->Ex[i][j][k+1]-D->Ex[i][j][k])
                         -0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
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

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(k=kstart; k<kend; k++)
      for(j=jstart; j<jend; j++)
      {
        nowPr=D->Pr[istart-1][j][k];
        nowSr=D->Sr[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]+=0.5*dt/dy*(D->PrC[i][j][k]+D->PrC[i-1][j][k]-D->PrC[i][j-1][k]-D->PrC[i-1][j-1][k]-D->PlC[i][j][k]-D->PlC[i-1][j][k]+D->PlC[i][j-1][k]+D->PlC[i-1][j-1][k])
                         +0.5*dt/dz*(D->SrC[i][j][k]+D->SrC[i-1][j][k]-D->SrC[i][j][k-1]-D->SrC[i-1][j][k-1]-D->SlC[i][j][k]-D->SlC[i-1][j][k]+D->SlC[i][j][k-1]+D->SlC[i-1][j][k-1])
                         -pi*dt*(D->Jx[i][j][k]+D->Jx[i-1][j][k]);
          D->Bx[i][j][k]+=0.5*dt/dz*(D->PrC[i][j][k+1]+D->PrC[i-1][j][k+1]-D->PrC[i][j][k]-D->PrC[i-1][j][k]+D->PlC[i][j][k+1]+D->PlC[i-1][j][k+1]-D->PlC[i][j][k]-D->PlC[i-1][j][k])
                         -0.5*dt/dy*(D->SrC[i][j+1][k]+D->SrC[i-1][j+1][k]-D->SrC[i][j][k]-D->SrC[i-1][j][k]+D->SlC[i][j+1][k]+D->SlC[i-1][j+1][k]-D->SlC[i][j][k]-D->SlC[i-1][j][k]);
          prevPr=nowPr;
          nowPr=D->Pr[i][j][k];
          D->Pr[i][j][k]=prevPr
                        +0.5*dt/dz*(D->BxC[i][j][k]-D->BxC[i][j][k-1])
                        +0.5*dt/dy*(D->ExC[i][j+1][k]-D->ExC[i][j][k])
                        -pi*dt*D->Jy[i][j][k];
          D->Pl[i-1][j][k]=D->Pl[i][j][k]
                        +0.5*dt/dz*(D->BxC[i][j][k]-D->BxC[i][j][k-1])
                        -0.5*dt/dy*(D->ExC[i][j+1][k]-D->ExC[i][j][k])
                        -pi*dt*D->Jy[i][j][k];
          prevSr=nowSr;
          nowSr=D->Sr[i][j][k];
          D->Sr[i][j][k]=prevSr
                        -0.5*dt/dy*(D->BxC[i][j][k]-D->BxC[i][j-1][k])
                        +0.5*dt/dz*(D->ExC[i][j][k+1]-D->ExC[i][j][k])
                        -pi*dt*D->Jz[i][j][k];
          D->Sl[i-1][j][k]=D->Sl[i][j][k]
                        -0.5*dt/dy*(D->BxC[i][j][k]-D->BxC[i][j-1][k])
                        -0.5*dt/dz*(D->ExC[i][j][k+1]-D->ExC[i][j][k])
                        -pi*dt*D->Jz[i][j][k];
        }	//End of i
      }		//End of j,k
}

