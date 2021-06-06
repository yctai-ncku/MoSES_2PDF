#ifndef __SWE_KERNELS_CUH
#define __SWE_KERNELS_CUH

#include <math.h>
#include <stdio.h>

#define MD      3
#define MN      7
#define MU      0.21
// #define GRAVITY 9.8
#define H0      1.0e-6
#define H00     0.0
#define U0      1.0e-6
#define Hmann   0.01
#define BLOCK_SIZE 16
#define CFL     0.1
#define PI     3.14159

// __device__ float CFL = 0.1f;

__global__
void makeTopo1Kernel(double *result,
                     double* bfkt, 
                     double  xmin, double ymin,
                     double dx, double dy,
                     int nxd, int nyd);

__global__
void makeTopo2Kernel(double* result, 
                double* topo,
                double* bfkt,
                double dx, double dy,
                int nxd, int nyd, int nx, int ny);

// __global__
// void makeTopo21Kernel(double* result, 
//                 double* topo,
//                 double* bfkt,
//                 double dx, double dy,
//                 int nxd, int nyd, int nx, int ny);

__global__
void makeTopo3Kernel(double* result, 
    double* bfkt,
    double* posx, double* posy,
    int nxd, int nyd);
    

__global__
void makeTopo4Kernel(double* result, 
    double* bfkt,
    double* posx, double* posy,
    double* dxdxi11, double* dxdxi12,
    double* dxdxi21, double* dxdxi22,
    double* dbdx, double* dbdy,
    int nxd, int nyd);
    

__global__
void makeTopo5Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    int nxd, int nyd);

__global__
void makeTopo6Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    double* svec, 
    int nxd, int nyd);
    
__global__
void makeTopo7Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    double* svec, 
    double* Jacb31, double* Jacb32,
    double* dxdxi11, double* dxdxi12,
    double* dxdxi21, double* dxdxi22,  
    double* dettmp, 
    int nxd, int nyd);

__global__
void makeTopo8Kernel(double* result,
    double* cvalue, double* svec, 
    double* Jacb31, double* Jacb32,
    double* dxdxi11, double* dxdxi12,
    double* dxdxi21, double* dxdxi22,  
    double* dettmp, 
    double* Detmin,
    double* i_ddxi11, double* i_ddxi12, 
    double* i_ddxi21, double* i_ddxi22,  
    int nxd, int nyd);

__global__
void makeTopo9Kernel(double* result,
    double* cvalue, double* svec,     
    double* i_ddxi11, double* i_ddxi12, 
    double* i_ddxi21, double* i_ddxi22,  
    double* invJ11, double* invJ12, double* invJ13,
    double* invJ21, double* invJ22, double* invJ23,
    double* invJ31, double* invJ32, double* invJ33, 
    int nxd, int nyd);

    
__global__
void makeTopo10Kernel(double* result,
    double* cvalue, double* depth, 
    double* u,
    double  phiS0, 
    int nxd, int nyd);


__global__
void makeTopo11Kernel(double* result,
    double* u, 
    double* tande, double delta0,
    int nxd, int nyd);
    

__global__
void Boundary1Kernel(double* result,
    double* dxdxi11, double* dxdxi12, 
    double* dxdxi21, double* dxdxi22,
    double* cvalue , double* Detmin,
    double* svec, 
    double* Jacb31, double* Jacb32,
    double* invJ11, double* invJ12, double* invJ13,
    double* invJ21, double* invJ22, double* invJ23,  
    double* invJ31, double* invJ32, double* invJ33,  
    int nxd, int nyd, int nx, int ny);


__global__
void Boundary2Kernel(double* result,
    double* dxdxi11, double* dxdxi12, 
    double* dxdxi21, double* dxdxi22,
    double* cvalue , double* Detmin,
    double* svec, 
    double* Jacb31, double* Jacb32,
    double* invJ11, double* invJ12, double* invJ13,
    double* invJ21, double* invJ22, double* invJ23,  
    double* invJ31, double* invJ32, double* invJ33,  
    int nxd, int nyd, int nx, int ny);


__global__
void JacobKernel(double* result,
    double* svec, double* cvalue, 
    double* posx, double* posy, 
    double* J13dxi, double* J23dxi, double* J33dxi,  
    double* J13det, double* J23det, double* J33det,  
    int nxd, int nyd);

__global__
void Boundary3Kernel(double* result,
    double* J13dxi, double* J23dxi, double* J33dxi,  
    double* J13det, double* J23det, double* J33det,  
    int nxd, int nyd, int nx, int ny);

__global__
void Boundary4Kernel(double* result,
    double* J13dxi, double* J23dxi, double* J33dxi,  
    double* J13det, double* J23det, double* J33det,  
    int nxd, int nyd, int nx, int ny);


__global__
void UzeroKernel(double* result,
    double* u, double* uzero,
    int nxd, int nyd);

__global__
void Boundary5Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny);

__global__
void Boundary6Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny);

__global__
void Boundary7Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny);

__global__
void Boundary8Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny);

__global__
void Boundary9Kernel(double* result,
    double* Hpx, double* Hpy,
    double* Ppx, double* Ppy, 
    double* PDx, double* PDy,
    double*  ux, double*  uy,    
    double* apEW, double* apSN,    
    double* apFEW, double* apFSN,    
    int nxd, int nyd);

__global__
void TVD1Kernel(double* result,
    double* u,    
    double* dux, double* duy,    
    int nxd, int nyd);


__global__
void TVD2Kernel(double* result,
    double* dux, double* duy,   
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    int nxd, int nyd);


__global__
void TVD3Kernel(double* result,
    double* dux, double* duy,   
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    double* t1x,   double* t2x, 
    double* t1y,   double* t2y, 
    int nxd, int nyd);
    

__global__
void TVD4Kernel(double* result,
    double* t1x, double* t2x,  
    double* t1y, double* t2y, 
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    int nxd, int nyd);

__global__
void TVD5Kernel(double* result,
    double* t1x, double* t2x, 
    double* t1y, double* t2y, 
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy,
    double*  ux, double* uy, 
    int nxd, int nyd);

__global__
void MeanKernel(double* result,
    double* dxdxi11, double* dxdxi21,
    double* dxdxi12, double* dxdxi22,

    double* J13dxi, double* J23dxi, double* J33dxi,
    double* J13det, double* J23det, double* J33det,

    double* invJ11, double* invJ12, double* invJ13, 
    double* invJ21, double* invJ22, double* invJ23, 
    double* invJ31, double* invJ32, double* invJ33, 

    double* Detmin, double* cvalue, double* svec, 

    double* dxdxi11_avgEW, double* dxdxi21_avgEW,
    double* dxdxi12_avgSN, double* dxdxi22_avgSN,

    double* J13dxi_avgEW, double* J23dxi_avgEW, double* J33dxi_avgEW,
    double* J13det_avgEW, double* J23det_avgEW, double* J33det_avgEW,

    double* J13dxi_avgSN, double* J23dxi_avgSN, double* J33dxi_avgSN,
    double* J13det_avgSN, double* J23det_avgSN, double* J33det_avgSN,

    double* invJ11_avgEW, double* invJ12_avgEW, double* invJ13_avgEW, 
    double* invJ21_avgEW, double* invJ22_avgEW, double* invJ23_avgEW, 
    double* invJ31_avgEW, double* invJ32_avgEW, double* invJ33_avgEW, 

    double* invJ11_avgSN, double* invJ12_avgSN, double* invJ13_avgSN, 
    double* invJ21_avgSN, double* invJ22_avgSN, double* invJ23_avgSN, 
    double* invJ31_avgSN, double* invJ32_avgSN, double* invJ33_avgSN, 

    double* Detmin_avgEW, double* Detmin_avgSN, 
    double* cval_avgEW, double* cval_avgSN, 
    double* svec_avgEW, double* svec_avgSN, 
     
    int nxd, int nyd);


__global__
void InterfacesKernel(double* result,
    double* u, 
    double* ux, double* uy,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd);


__global__
void Interfaces2Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd, int nx, int ny);


// __global__
// void Interfaces3Kernel(double* result,
//     double* uE, double* uW,  
//     double* uN, double* uS,       
//     int nxd, int nyd, int nx, int ny);

__global__
void KeepPositivi1Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd);

__global__
void KeepPositivi2Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd);

__global__
void Flux1Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS, 
    double* vexE, double* veyE, 
    double* vexW, double* veyW,
    double* vexFE, double* veyFE, 
    double* vexFW, double* veyFW,

    double* vexN, double* veyN, 
    double* vexS, double* veyS,
    double* vexFN, double* veyFN, 
    double* vexFS, double* veyFS,

    int nxd, int nyd);

__global__
void Flux2Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS, 
    double* vexE, double* veyE, 
    double* vexW, double* veyW,
    double* vexFE, double* veyFE, 
    double* vexFW, double* veyFW,

    double* vexN, double* veyN, 
    double* vexS, double* veyS,
    double* vexFN, double* veyFN, 
    double* vexFS, double* veyFS,

    double* w_wertE, double* w_wertW,
    double* w_wertFE, double* w_wertFW, 

    double* w_wertN, double* w_wertS,
    double* w_wertFN, double* w_wertFS, 

    double* svec, double* cvalue,

    int nxd, int nyd);

__global__
void Flux3Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS, 
    double* vexE, double* veyE, 
    double* vexW, double* veyW,
    double* vexFE, double* veyFE, 
    double* vexFW, double* veyFW,

    double* vexN, double* veyN, 
    double* vexS, double* veyS,
    double* vexFN, double* veyFN, 
    double* vexFS, double* veyFS,

    double* w_wertE, double* w_wertW,
    double* w_wertFE, double* w_wertFW, 

    double* w_wertN, double* w_wertS,
    double* w_wertFN, double* w_wertFS, 

    double* q_xiE, double* q_etE,
    double* q_xiW, double* q_etW,
    double* q_xiFE, double* q_etFE,
    double* q_xiFW, double* q_etFW,

    double* NpressFE, double* NpressFW, double* M11EW,
    
    double* invJ11_avgEW, double* invJ12_avgEW, double* invJ13_avgEW, 
    double* invJ21_avgEW, double* invJ22_avgEW, double* invJ23_avgEW, 
    double* cval_avgEW,

    double* q_xiN, double* q_etN,
    double* q_xiS, double* q_etS,
    double* q_xiFN, double* q_etFN,
    double* q_xiFS, double* q_etFS,

    double* NpressFN, double* NpressFS, double* M22SN,

    double* invJ11_avgSN, double* invJ12_avgSN, double* invJ13_avgSN, 
    double* invJ21_avgSN, double* invJ22_avgSN, double* invJ23_avgSN, 
    double* cval_avgSN,

    int nxd, int nyd);

__global__
void Flux4Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS, 

    double* q_xiE, double* q_xiW,
    double* q_xiFE, double* q_xiFW,
    double* NpressFE, double* NpressFW,  double* invJ11_avgEW,  

    double* apE, double* apW,
    double* apFE, double* apFW, 

    double* q_etN, double* q_etS,
    double* q_etFN, double* q_etFS,

    double* NpressFN, double* NpressFS, double* invJ22_avgSN, 

    double* apN, double* apS,
    double* apFN, double* apFS,

    int nxd, int nyd);

__global__
void Flux5Kernel(double* result,

    double* apE, double* apW,
    double* apFE, double* apFW, 

    double* apEW, double* apFEW,

    double* apN, double* apS,
    double* apFN, double* apFS,

    double* apSN, double* apFSN,

    int nxd, int nyd);

__global__
void Flux6Kernel(double* result,
    double* apEW, double* apFEW,
    double* apSN, double* apFSN,

    double* em_x, double* em_y, 
    double* em_Fx, double* em_Fy, 

    double* czw1x , double* czw2x, 
    double* czwF1x, double* czwF2x,
    double* czw1y , double* czw2y, 
    double* czwF1y, double* czwF2y,

    double* uE, double* uW, 
    double* uN, double* uS,

    double* cval_avgEW, double* cval_avgSN,
    double* Detmin_avgEW, double* Detmin_avgSN,
    double* M11EW, double* M22SN,
    
    int nxd, int nyd);

__global__
void Flux7Kernel(double* result,
    double* FpE, double* FpW,
    double* GpN, double* GpS,

    double* czw1x , double* czw2x, 
    double* czwF1x, double* czwF2x,
    double* czw1y , double* czw2y, 
    double* czwF1y, double* czwF2y,

    double* uE, double* uW, 
    double* uN, double* uS,

    double* Detmin_avgEW, double* Detmin_avgSN,

    double* q_xiE, double* q_xiFE,
    double* q_xiW, double* q_xiFW, 

    double* q_etN, double* q_etFN,
    double* q_etS, double* q_etFS,

    double* dxdxi11_avgEW, double* dxdxi21_avgEW,
    double* dxdxi12_avgSN, double* dxdxi22_avgSN,

    
    int nxd, int nyd);


__global__
void CFL1Kernel(double* result,
    double* em_x, double* em_y,
    double* em_Fx, double* em_Fy,
    double* em_valS, double* em_valF,
    double dx, double dy,
    int nxd, int nyd);

__global__
void CFL2Kernel(double* result,
    double* em_valS, double* em_valF,
    double* Val,
    int nxd, int nyd);

__global__
void CFL3Kernel(double* result,
    double* dtval, double* d_max1, double* totalTime);

__global__
void CFL4Kernel(double* result,
    double* dt, double* dtval,double* totalTime, int io);

__global__
void reduceKernel(double *waveSpeed, double* d_max, int size);


__global__
void Flux8Kernel(double* result,
    double* Hpx, double* Hpy,
    double* Ppx, double* Ppy,
    double* FpE, double* FpW,
    double* GpN, double* GpS,

    double* apEW, double* apFEW, 
    double* apSN, double* apFSN,
    double* uE, double* uW, 
    double* uN, double* uS,
    double* u,
    double* ux, double* uy,

    double* Detmin_avgEW, double* Detmin_avgSN,
    double* cval_avgEW, double* cval_avgSN, 
 
    double* invJ11_avgEW, double* invJ12_avgEW, double* invJ13_avgEW,
    double* invJ21_avgEW, double* invJ22_avgEW, double* invJ23_avgEW,
    double* invJ31_avgEW, double* invJ32_avgEW, double* invJ33_avgEW,

    double* invJ11_avgSN, double* invJ12_avgSN, double* invJ13_avgSN,
    double* invJ21_avgSN, double* invJ22_avgSN, double* invJ23_avgSN,
    double* invJ31_avgSN, double* invJ32_avgSN, double* invJ33_avgSN,

    double* dudxE, double* dvdxE,
    double* dudyE, double* dvdyE,

    double* dudxN, double* dvdxN,
    double* dudyN, double* dvdyN,

    double dx, double dy,

    int nxd, int nyd);

__global__
void Flux9Kernel(double* result,
    double* duxidxix, double* dvetdxix,
    double* duxidetx, double* dvetdetx,
    double* duxidxiy, double* dvetdxiy,
    double* duxidety, double* dvetdety, 
 
    double* invJ11_avgEW, double* invJ12_avgEW, double* invJ13_avgEW,
    double* invJ21_avgEW, double* invJ22_avgEW, double* invJ23_avgEW,
    double* invJ31_avgEW, double* invJ32_avgEW, double* invJ33_avgEW,

    double* invJ11_avgSN, double* invJ12_avgSN, double* invJ13_avgSN,
    double* invJ21_avgSN, double* invJ22_avgSN, double* invJ23_avgSN,
    double* invJ31_avgSN, double* invJ32_avgSN, double* invJ33_avgSN,

    double* dudxE, double* dvdxE,
    double* dudyE, double* dvdyE,

    double* dudxN, double* dvdxN,
    double* dudyN, double* dvdyN,

    int nxd, int nyd);

__global__
void Flux10Kernel(double* result,
    double* PDx, double* PDy,

    double* uE, double* uW, 
    double* uN, double* uS,

    double* Detmin_avgEW, double* Detmin_avgSN,

    double* duxidxix, double* dvetdxix,
    double* duxidetx, double* dvetdetx,
    double* duxidxiy, double* dvetdxiy,
    double* duxidety, double* dvetdety,
 
    double* invJ11_avgEW, double* invJ12_avgEW, 
    double* invJ21_avgEW, double* invJ22_avgEW, 

    double* invJ11_avgSN, double* invJ12_avgSN, 
    double* invJ21_avgSN, double* invJ22_avgSN, 

    int nxd, int nyd);


__global__
void Flux11Kernel(double* result,
    double* u, 
    double* vex, double* vey, 
    double* vexF,double* veyF,
    int nxd, int nyd);

__global__
void Flux12Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* vexF, double* veyF,
    double* svec, double* cvalue,
    int nxd, int nyd);

__global__
void Flux13Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* usw, double* vel,
    int nxd, int nyd);

__global__
void Flux14Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* vexF, double* veyF,
    double* vexw, double* veyw, 
    double* usw, double* vel,
    double* q_xi, double* q_et,
    double* q_xiF,double* q_etF,
    double* invJ11, double* invJ12, double* invJ13, 
    double* invJ21, double* invJ22, double* invJ23, 
    int nxd, int nyd);

__global__
void Flux15Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* vexF, double* veyF, 
    double* q_xi, double* q_et,
    double* q_xiF,double* q_etF,
    double* J13dxi, double* J23dxi, double* J33dxi, 
    double* J13det, double* J23det, double* J33det,
    double* Ac, double* AcF, 
    int nxd, int nyd);

__global__
void Flux16Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* vexF, double* veyF, 
    double* q_xi, double* q_et,
    double* q_xiF,double* q_etF,
    double* J13dxi, double* J23dxi, double* J33dxi, 
    double* J13det, double* J23det, double* J33det,
    double* Ac, double* AcF, 
    int nxd, int nyd);

__global__
void Flux16Kernel(double* result,
    double* Npress1, double* Npress2,  double* NpressF,
    double* Ac, double* AcF, 
    double* Detmin, double* cvalue,
    double* u,
    int nxd, int nyd);

__global__
void Flux17Kernel(double* result,
    double* Npress1, double* Npress2,  double* NpressF,
    double* Ac, double* AcF, 
    double* Detmin, double* svec,
    double* vex, double* vey,
    double* vexw, double* veyw,
    double* vexF, double* veyF,
    double* tande,
    double* u, double* s,
    double Cd, double N_R, double varTheta,
    int nxd, int nyd);

__global__
void Flux18Kernel(double* result,
    double* v, double* Detmin, 
    double* Hpx, double* Hpy, 
    double* Ppx, double* Ppy, 
    double* PDx, double* PDy, 
    double* s, 
    double* u, double* uzero,
    double* dt,
    double dx, double dy,
    double N_R,
    int nxd, int nyd);

__global__
void Flux19Kernel(double* result,
    double* v, double* Detmin, 
    double* u, 
    int nxd, int nyd);  

__global__
void Flux201Kernel(double* result,
    double* u, double* utmp,
    int nxd, int nyd);

__global__
void Flux20Kernel(double* result,
    double* u, double* utmp,
    int nxd, int nyd); 

__global__
void Flux21Kernel(double* result,
    double* u, 
    int nxd, int nyd); 

__global__
void Flux22Kernel(double* result,
    double* u, double* uone,
    int nxd, int nyd); 

__global__
void Flux23Kernel(double* result,
    double* u, double* uzero,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd);

__global__
void Flux24Kernel(double* result,
    double* u, double* uone,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd);

__global__
void Flux25Kernel(double* result,
    double* u, double* utwo,
    int nxd, int nyd);

__global__
void Flux26Kernel(double* result,
    double* u, double* uone,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd);

__global__
void Flux27Kernel(double* result,
    double* u, double* utwo,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd);

__global__
void Flux28Kernel(double* result,
    double* u, double* utwo,
    int nxd, int nyd);

__global__
void Flux29Kernel(double* result,
    double* u, double* utwo, double* uone,
    int nxd, int nyd);

#endif