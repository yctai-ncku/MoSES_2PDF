#include "adNOC_2Pb_kernels.cuh"
#include "adNOC_2Pb_runKernels.cuh"


__global__
void makeTopo1Kernel(double *result,
                     double* bfkt, 
                     double  xmin, double ymin,
                     double dx, double dy,
                     int nxd, int nyd){
    
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
        
        bfkt[0 * nxd * nyd + (row   ) * nxd + (col+MD)] = xmin + col*dx;

        bfkt[1 * nxd * nyd + (row+MD) * nxd + (col   )] = ymin + row*dy; 
    
    }

}

__global__
void makeTopo2Kernel(double* result, 
                double* topo,
                double* bfkt,
                double dx, double dy,
                int nxd, int nyd, int nx, int ny) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
        bfkt[0 * nxd * nyd + row * nxd + 0] = bfkt[0 * nxd * nyd + row * nxd + MD]-(MD-0)*dx;
        bfkt[0 * nxd * nyd + row * nxd + 1] = bfkt[0 * nxd * nyd + row * nxd + MD]-(MD-1)*dx;
        bfkt[0 * nxd * nyd + row * nxd + 2] = bfkt[0 * nxd * nyd + row * nxd + MD]-(MD-2)*dx;
        
        bfkt[1 * nxd * nyd + 0 * nxd + col] = bfkt[1 * nxd * nyd + MD * nxd + col]-(MD-0)*dy;
        bfkt[1 * nxd * nyd + 1 * nxd + col] = bfkt[1 * nxd * nyd + MD * nxd + col]-(MD-1)*dy;
        bfkt[1 * nxd * nyd + 2 * nxd + col] = bfkt[1 * nxd * nyd + MD * nxd + col]-(MD-2)*dy;

        bfkt[0 * nxd * nyd + row * nxd + (nx+MD+0)] = bfkt[0 * nxd * nyd + row * nxd + (nx+MD-1)]+(1+0)*dx;
        bfkt[0 * nxd * nyd + row * nxd + (nx+MD+1)] = bfkt[0 * nxd * nyd + row * nxd + (nx+MD-1)]+(1+1)*dx;
        bfkt[0 * nxd * nyd + row * nxd + (nx+MD+2)] = bfkt[0 * nxd * nyd + row * nxd + (nx+MD-1)]+(1+2)*dx;

        bfkt[1 * nxd * nyd + (ny+MD+0) * nxd + col] = bfkt[1 * nxd * nyd + (ny+MD-1) * nxd + col]+(1+0)*dy;
        bfkt[1 * nxd * nyd + (ny+MD+1) * nxd + col] = bfkt[1 * nxd * nyd + (ny+MD-1) * nxd + col]+(1+1)*dy;
        bfkt[1 * nxd * nyd + (ny+MD+2) * nxd + col] = bfkt[1 * nxd * nyd + (ny+MD-1) * nxd + col]+(1+2)*dy;

        bfkt[2 * nxd * nyd + (row   ) * nxd + (col  )] = topo[row * nxd + col];
    }

    
}

__global__
void makeTopo3Kernel(double* result, 
    double* bfkt,
    double* posx, double* posy,
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
        posx[row * nxd + col] = bfkt[0 * nxd * nyd + row * nxd + col];
        posy[row * nxd + col] = bfkt[1 * nxd * nyd + row * nxd + col];
    }

}

__global__
void makeTopo4Kernel(double* result, 
    double* bfkt,
    double* posx, double* posy,
    double* dxdxi11, double* dxdxi12,
    double* dxdxi21, double* dxdxi22,
    double* dbdx, double* dbdy,
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
        dxdxi11[row * nxd + col] = bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)]*( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(  posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)]))
                                  +bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )]*((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])-( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))
                                  -bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)]*( posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)])*(  posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]));
    
        dxdxi12[row * nxd + col] = bfkt[0 * nxd * nyd + (row-1) * nxd + (col  )]*( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(  posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                  +bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )]*((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])-( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                  -bfkt[0 * nxd * nyd + (row+1) * nxd + (col  )]*( posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])*(  posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]));
 
        dxdxi21[row * nxd + col] = bfkt[1 * nxd * nyd + (row  ) * nxd + (col-1)]*( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(  posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)]))
                                  +bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )]*((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])-( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))
                                  -bfkt[1 * nxd * nyd + (row  ) * nxd + (col+1)]*( posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)])*(  posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]));
    
        dxdxi22[row * nxd + col] = bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )]*( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(  posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                  +bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )]*((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])-( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                  -bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )]*( posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])*(  posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]));
 
           dbdx[row * nxd + col] = bfkt[2 * nxd * nyd + (row  ) * nxd + (col-1)]*( bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)])/((bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )])* ( bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)]))
                                  +bfkt[2 * nxd * nyd + (row  ) * nxd + (col  )]*((bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )])-( bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)]))/((bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )])*(bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)]))
                                  -bfkt[2 * nxd * nyd + (row  ) * nxd + (col+1)]*( bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )])/((bfkt[0 * nxd * nyd + (row  ) * nxd + (col-1)] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)])* ( bfkt[0 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[0 * nxd * nyd + (row  ) * nxd + (col+1)]));
            
           dbdy[row * nxd + col] = bfkt[2 * nxd * nyd + (row-1) * nxd + (col  )]*( bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )])/((bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )])* ( bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )]))
                                  +bfkt[2 * nxd * nyd + (row  ) * nxd + (col  )]*((bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )])-( bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )]))/((bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )])*(bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )]))
                                  -bfkt[2 * nxd * nyd + (row+1) * nxd + (col  )]*( bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )])/((bfkt[1 * nxd * nyd + (row-1) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )])* ( bfkt[1 * nxd * nyd + (row  ) * nxd + (col  )] - bfkt[1 * nxd * nyd + (row+1) * nxd + (col  )]));
   
    }

}

__global__
void makeTopo5Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
        cvalue[row * nxd + col] = 1.0/sqrt(1.0+
                                dbdx[row * nxd + col]*dbdx[row * nxd + col]
                               +dbdy[row * nxd + col]*dbdy[row * nxd + col]);
    }

}

__global__
void makeTopo6Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    double* svec, 
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
        svec[0 * nxd * nyd + row * nxd + col] = cvalue[row * nxd + col]*dbdx[row * nxd + col];
        svec[1 * nxd * nyd + row * nxd + col] = cvalue[row * nxd + col]*dbdy[row * nxd + col];
    }

}

__global__
void makeTopo7Kernel(double* result, 
    double* dbdx, double* dbdy,
    double* cvalue,
    double* svec, 
    double* Jacb31, double* Jacb32,
    double* dxdxi11, double* dxdxi12,
    double* dxdxi21, double* dxdxi22,  
    double* dettmp, 
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
        Jacb31[row * nxd + col] = (svec[0 * nxd * nyd + row * nxd + col]*dxdxi11[row * nxd + col]
                                  +svec[1 * nxd * nyd + row * nxd + col]*dxdxi21[row * nxd + col])
                                  /cvalue[row * nxd + col];

        Jacb32[row * nxd + col] = (svec[0 * nxd * nyd + row * nxd + col]*dxdxi12[row * nxd + col]
                                  +svec[1 * nxd * nyd + row * nxd + col]*dxdxi22[row * nxd + col])
                                  /cvalue[row * nxd + col];

        dettmp[row * nxd + col] = dxdxi11[row * nxd + col]*dxdxi22[row * nxd + col]
                                 -dxdxi21[row * nxd + col]*dxdxi12[row * nxd + col];


    }

}

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
    int nxd, int nyd) 
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
        Detmin[row * nxd + col] = dxdxi11[row * nxd + col]*(dxdxi22[row * nxd + col]*cvalue[row * nxd + col]-(-svec[1 * nxd * nyd + row * nxd + col])*Jacb32[row * nxd + col] )
                                 -dxdxi21[row * nxd + col]*(dxdxi12[row * nxd + col]*cvalue[row * nxd + col]-Jacb32[row * nxd + col] *(-svec[0 * nxd * nyd + row * nxd + col]))
                                 + Jacb31[row * nxd + col]*(dxdxi12[row * nxd + col]*(-svec[1 * nxd * nyd + row * nxd + col])-dxdxi22[row * nxd + col] *(-svec[0 * nxd * nyd + row * nxd + col]));
       

        i_ddxi11[row * nxd + col] =  dxdxi22[row * nxd + col]/dettmp[row * nxd + col];
        i_ddxi12[row * nxd + col] = -dxdxi12[row * nxd + col]/dettmp[row * nxd + col];
        i_ddxi21[row * nxd + col] = -dxdxi21[row * nxd + col]/dettmp[row * nxd + col];
        i_ddxi22[row * nxd + col] =  dxdxi11[row * nxd + col]/dettmp[row * nxd + col];
    }


}

__global__
void makeTopo9Kernel(double* result,
    double* cvalue, double* svec,     
    double* i_ddxi11, double* i_ddxi12, 
    double* i_ddxi21, double* i_ddxi22,  
    double* invJ11, double* invJ12, double* invJ13,
    double* invJ21, double* invJ22, double* invJ23,
    double* invJ31, double* invJ32, double* invJ33, 
    int nxd, int nyd)  
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-2) && col < (nxd-2) && row >1 && col >1 ) {
   
        invJ11[row * nxd + col] =  i_ddxi11[row * nxd + col]-i_ddxi11[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]
                                                            -i_ddxi12[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col];
   
        invJ12[row * nxd + col] =  i_ddxi12[row * nxd + col]-i_ddxi11[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]
                                                            -i_ddxi12[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col];

        invJ21[row * nxd + col] =  i_ddxi21[row * nxd + col]-i_ddxi21[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]
                                                            -i_ddxi22[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col];

        invJ22[row * nxd + col] =  i_ddxi22[row * nxd + col]-i_ddxi21[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]
                                                            -i_ddxi22[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col];

        
        invJ13[row * nxd + col] = cvalue[row * nxd + col]*( svec[0 * nxd * nyd + row * nxd + col]*i_ddxi11[row * nxd + col]
                                                           +svec[1 * nxd * nyd + row * nxd + col]*i_ddxi12[row * nxd + col]);

        invJ23[row * nxd + col] = cvalue[row * nxd + col]*( svec[0 * nxd * nyd + row * nxd + col]*i_ddxi21[row * nxd + col]
                                                           +svec[1 * nxd * nyd + row * nxd + col]*i_ddxi22[row * nxd + col]);

        invJ31[row * nxd + col] = -svec[0 * nxd * nyd + row * nxd + col];

        invJ32[row * nxd + col] = -svec[1 * nxd * nyd + row * nxd + col];

        invJ33[row * nxd + col] = cvalue[row * nxd + col];

    }


}

__global__
void makeTopo10Kernel(double* result,
    double* cvalue, double* depth, 
    double* u,
    double  phis0, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
   
        u[0 * nxd * nyd + row * nxd + col] = depth[row * nxd + col]*phis0*cvalue[row * nxd + col];

        u[1 * nxd * nyd + row * nxd + col] = 0.0;
        u[2 * nxd * nyd + row * nxd + col] = 0.0;

        u[3 * nxd * nyd + row * nxd + col] = depth[row * nxd + col]*(1.0-phis0)*cvalue[row * nxd + col];
        u[4 * nxd * nyd + row * nxd + col] = 0.0;
        u[5 * nxd * nyd + row * nxd + col] = 0.0;
        u[6 * nxd * nyd + row * nxd + col] = phis0;
        

    }


}

__global__
void makeTopo11Kernel(double* result,
    double* u, 
    double* tande, double delta0,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
 
        tande[row * nxd + col] = tan(delta0*PI/180.);

    }

}

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
    int nxd, int nyd, int nx, int ny)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    // int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && row >=0 ) {
 
        dxdxi11[row * nxd + 0] = dxdxi11[row * nxd + MD];
        dxdxi11[row * nxd + 1] = dxdxi11[row * nxd + MD];
        dxdxi11[row * nxd + 2] = dxdxi11[row * nxd + MD];

        dxdxi12[row * nxd + 0] = dxdxi12[row * nxd + MD];
        dxdxi12[row * nxd + 1] = dxdxi12[row * nxd + MD];
        dxdxi12[row * nxd + 2] = dxdxi12[row * nxd + MD];

        dxdxi21[row * nxd + 0] = dxdxi21[row * nxd + MD];
        dxdxi21[row * nxd + 1] = dxdxi21[row * nxd + MD];
        dxdxi21[row * nxd + 2] = dxdxi21[row * nxd + MD];

        dxdxi22[row * nxd + 0] = dxdxi22[row * nxd + MD];
        dxdxi22[row * nxd + 1] = dxdxi22[row * nxd + MD];
        dxdxi22[row * nxd + 2] = dxdxi22[row * nxd + MD];

        cvalue[row * nxd + 0] = cvalue[row * nxd + MD];
        cvalue[row * nxd + 1] = cvalue[row * nxd + MD];
        cvalue[row * nxd + 2] = cvalue[row * nxd + MD];

        Detmin[row * nxd + 0] = Detmin[row * nxd + MD];
        Detmin[row * nxd + 1] = Detmin[row * nxd + MD];
        Detmin[row * nxd + 2] = Detmin[row * nxd + MD];

        svec[0 * nxd * nyd + row * nxd + 0] = svec[0 * nxd * nyd + row * nxd + MD];
        svec[0 * nxd * nyd + row * nxd + 1] = svec[0 * nxd * nyd + row * nxd + MD];
        svec[0 * nxd * nyd + row * nxd + 2] = svec[0 * nxd * nyd + row * nxd + MD];

        svec[1 * nxd * nyd + row * nxd + 0] = svec[1 * nxd * nyd + row * nxd + MD];
        svec[1 * nxd * nyd + row * nxd + 1] = svec[1 * nxd * nyd + row * nxd + MD];
        svec[1 * nxd * nyd + row * nxd + 2] = svec[1 * nxd * nyd + row * nxd + MD];

        Jacb31[row * nxd + 0] = Jacb31[row * nxd + MD];
        Jacb31[row * nxd + 1] = Jacb31[row * nxd + MD];
        Jacb31[row * nxd + 2] = Jacb31[row * nxd + MD];

        Jacb32[row * nxd + 0] = Jacb32[row * nxd + MD];
        Jacb32[row * nxd + 1] = Jacb32[row * nxd + MD];
        Jacb32[row * nxd + 2] = Jacb32[row * nxd + MD];

        invJ11[row * nxd + 0] = invJ11[row * nxd + MD];
        invJ11[row * nxd + 1] = invJ11[row * nxd + MD];
        invJ11[row * nxd + 2] = invJ11[row * nxd + MD];

        invJ12[row * nxd + 0] = invJ12[row * nxd + MD];
        invJ12[row * nxd + 1] = invJ12[row * nxd + MD];
        invJ12[row * nxd + 2] = invJ12[row * nxd + MD];

        invJ13[row * nxd + 0] = invJ13[row * nxd + MD];
        invJ13[row * nxd + 1] = invJ13[row * nxd + MD];
        invJ13[row * nxd + 2] = invJ13[row * nxd + MD];

        invJ21[row * nxd + 0] = invJ21[row * nxd + MD];
        invJ21[row * nxd + 1] = invJ21[row * nxd + MD];
        invJ21[row * nxd + 2] = invJ21[row * nxd + MD];

        invJ22[row * nxd + 0] = invJ22[row * nxd + MD];
        invJ22[row * nxd + 1] = invJ22[row * nxd + MD];
        invJ22[row * nxd + 2] = invJ22[row * nxd + MD];

        invJ23[row * nxd + 0] = invJ23[row * nxd + MD];
        invJ23[row * nxd + 1] = invJ23[row * nxd + MD];
        invJ23[row * nxd + 2] = invJ23[row * nxd + MD];

        invJ31[row * nxd + 0] = invJ31[row * nxd + MD];
        invJ31[row * nxd + 1] = invJ31[row * nxd + MD];
        invJ31[row * nxd + 2] = invJ31[row * nxd + MD];

        invJ32[row * nxd + 0] = invJ32[row * nxd + MD];
        invJ32[row * nxd + 1] = invJ32[row * nxd + MD];
        invJ32[row * nxd + 2] = invJ32[row * nxd + MD];

        invJ33[row * nxd + 0] = invJ33[row * nxd + MD];
        invJ33[row * nxd + 1] = invJ33[row * nxd + MD];
        invJ33[row * nxd + 2] = invJ33[row * nxd + MD];



        dxdxi11[row * nxd + (nx+MD+0)] = dxdxi11[row * nxd + (nx+MD-1)];
        dxdxi11[row * nxd + (nx+MD+1)] = dxdxi11[row * nxd + (nx+MD-1)];
        dxdxi11[row * nxd + (nx+MD+2)] = dxdxi11[row * nxd + (nx+MD-1)];

        dxdxi12[row * nxd + (nx+MD+0)] = dxdxi12[row * nxd + (nx+MD-1)];
        dxdxi12[row * nxd + (nx+MD+1)] = dxdxi12[row * nxd + (nx+MD-1)];
        dxdxi12[row * nxd + (nx+MD+2)] = dxdxi12[row * nxd + (nx+MD-1)];

        dxdxi21[row * nxd + (nx+MD+0)] = dxdxi21[row * nxd + (nx+MD-1)];
        dxdxi21[row * nxd + (nx+MD+1)] = dxdxi21[row * nxd + (nx+MD-1)];
        dxdxi21[row * nxd + (nx+MD+1)] = dxdxi21[row * nxd + (nx+MD-1)];

        dxdxi22[row * nxd + (nx+MD+0)] = dxdxi22[row * nxd + (nx+MD-1)];
        dxdxi22[row * nxd + (nx+MD+1)] = dxdxi22[row * nxd + (nx+MD-1)];
        dxdxi22[row * nxd + (nx+MD+2)] = dxdxi22[row * nxd + (nx+MD-1)];

        cvalue[row * nxd + (nx+MD+0)] = cvalue[row * nxd + (nx+MD-1)];
        cvalue[row * nxd + (nx+MD+1)] = cvalue[row * nxd + (nx+MD-1)];
        cvalue[row * nxd + (nx+MD+2)] = cvalue[row * nxd + (nx+MD-1)];

        Detmin[row * nxd + (nx+MD+0)] = Detmin[row * nxd + (nx+MD-1)];
        Detmin[row * nxd + (nx+MD+1)] = Detmin[row * nxd + (nx+MD-1)];
        Detmin[row * nxd + (nx+MD+2)] = Detmin[row * nxd + (nx+MD-1)];

        svec[0 * nxd * nyd + row * nxd + (nx+MD+0)] = svec[0 * nxd * nyd + row * nxd + (nx+MD-1)];
        svec[0 * nxd * nyd + row * nxd + (nx+MD+1)] = svec[0 * nxd * nyd + row * nxd + (nx+MD-1)];
        svec[0 * nxd * nyd + row * nxd + (nx+MD+2)] = svec[0 * nxd * nyd + row * nxd + (nx+MD-1)];

        svec[1 * nxd * nyd + row * nxd + (nx+MD+0)] = svec[1 * nxd * nyd + row * nxd + (nx+MD-1)];
        svec[1 * nxd * nyd + row * nxd + (nx+MD+1)] = svec[1 * nxd * nyd + row * nxd + (nx+MD-1)];
        svec[1 * nxd * nyd + row * nxd + (nx+MD+2)] = svec[1 * nxd * nyd + row * nxd + (nx+MD-1)];

        Jacb31[row * nxd + (nx+MD+0)] = Jacb31[row * nxd + (nx+MD-1)];
        Jacb31[row * nxd + (nx+MD+1)] = Jacb31[row * nxd + (nx+MD-1)];
        Jacb31[row * nxd + (nx+MD+2)] = Jacb31[row * nxd + (nx+MD-1)];

        Jacb32[row * nxd + (nx+MD+0)] = Jacb32[row * nxd + (nx+MD-1)];
        Jacb32[row * nxd + (nx+MD+1)] = Jacb32[row * nxd + (nx+MD-1)];
        Jacb32[row * nxd + (nx+MD+2)] = Jacb32[row * nxd + (nx+MD-1)];

        invJ11[row * nxd + (nx+MD+0)] = invJ11[row * nxd + (nx+MD-1)];
        invJ11[row * nxd + (nx+MD+1)] = invJ11[row * nxd + (nx+MD-1)];
        invJ11[row * nxd + (nx+MD+2)] = invJ11[row * nxd + (nx+MD-1)];

        invJ12[row * nxd + (nx+MD+0)] = invJ12[row * nxd + (nx+MD-1)];
        invJ12[row * nxd + (nx+MD+1)] = invJ12[row * nxd + (nx+MD-1)];
        invJ12[row * nxd + (nx+MD+2)] = invJ12[row * nxd + (nx+MD-1)];

        invJ13[row * nxd + (nx+MD+0)] = invJ13[row * nxd + (nx+MD-1)];
        invJ13[row * nxd + (nx+MD+1)] = invJ13[row * nxd + (nx+MD-1)];
        invJ13[row * nxd + (nx+MD+2)] = invJ13[row * nxd + (nx+MD-1)];

        invJ21[row * nxd + (nx+MD+0)] = invJ21[row * nxd + (nx+MD-1)];
        invJ21[row * nxd + (nx+MD+1)] = invJ21[row * nxd + (nx+MD-1)];
        invJ21[row * nxd + (nx+MD+2)] = invJ21[row * nxd + (nx+MD-1)];

        invJ22[row * nxd + (nx+MD+0)] = invJ22[row * nxd + (nx+MD-1)];
        invJ22[row * nxd + (nx+MD+1)] = invJ22[row * nxd + (nx+MD-1)];
        invJ22[row * nxd + (nx+MD+2)] = invJ22[row * nxd + (nx+MD-1)];

        invJ23[row * nxd + (nx+MD+0)] = invJ23[row * nxd + (nx+MD-1)];
        invJ23[row * nxd + (nx+MD+1)] = invJ23[row * nxd + (nx+MD-1)];
        invJ23[row * nxd + (nx+MD+2)] = invJ23[row * nxd + (nx+MD-1)];

        invJ31[row * nxd + (nx+MD+0)] = invJ31[row * nxd + (nx+MD-1)];
        invJ31[row * nxd + (nx+MD+1)] = invJ31[row * nxd + (nx+MD-1)];
        invJ31[row * nxd + (nx+MD+2)] = invJ31[row * nxd + (nx+MD-1)];

        invJ32[row * nxd + (nx+MD+0)] = invJ32[row * nxd + (nx+MD-1)];
        invJ32[row * nxd + (nx+MD+1)] = invJ32[row * nxd + (nx+MD-1)];
        invJ32[row * nxd + (nx+MD+2)] = invJ32[row * nxd + (nx+MD-1)];

        invJ33[row * nxd + (nx+MD+0)] = invJ33[row * nxd + (nx+MD-1)];
        invJ33[row * nxd + (nx+MD+1)] = invJ33[row * nxd + (nx+MD-1)];
        invJ33[row * nxd + (nx+MD+2)] = invJ33[row * nxd + (nx+MD-1)];


    }

}


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
    int nxd, int nyd, int nx, int ny)
{
    // int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (col < (nxd) && col >=0 ) {
 
        dxdxi11[0 * nxd + col] = dxdxi11[MD * nxd + col];
        dxdxi11[1 * nxd + col] = dxdxi11[MD * nxd + col];
        dxdxi11[2 * nxd + col] = dxdxi11[MD * nxd + col];

        dxdxi12[0 * nxd + col] = dxdxi12[MD * nxd + col];
        dxdxi12[1 * nxd + col] = dxdxi12[MD * nxd + col];
        dxdxi12[2 * nxd + col] = dxdxi12[MD * nxd + col];

        dxdxi21[0 * nxd + col] = dxdxi21[MD * nxd + col];
        dxdxi21[1 * nxd + col] = dxdxi21[MD * nxd + col];
        dxdxi21[2 * nxd + col] = dxdxi21[MD * nxd + col];

        dxdxi22[0 * nxd + col] = dxdxi22[MD * nxd + col];
        dxdxi22[1 * nxd + col] = dxdxi22[MD * nxd + col];
        dxdxi22[2 * nxd + col] = dxdxi22[MD * nxd + col];

        cvalue[0 * nxd + col] = cvalue[MD * nxd + col];
        cvalue[1 * nxd + col] = cvalue[MD * nxd + col];
        cvalue[2 * nxd + col] = cvalue[MD * nxd + col];

        Detmin[0 * nxd + col] = Detmin[MD * nxd + col];
        Detmin[1 * nxd + col] = Detmin[MD * nxd + col];
        Detmin[2 * nxd + col] = Detmin[MD * nxd + col];

        svec[0 * nxd * nyd + 0 * nxd + col] = svec[0 * nxd * nyd + MD * nxd + col];
        svec[0 * nxd * nyd + 1 * nxd + col] = svec[0 * nxd * nyd + MD * nxd + col];
        svec[0 * nxd * nyd + 2 * nxd + col] = svec[0 * nxd * nyd + MD * nxd + col];

        svec[1 * nxd * nyd + 0 * nxd + col] = svec[1 * nxd * nyd + MD * nxd + col];
        svec[1 * nxd * nyd + 1 * nxd + col] = svec[1 * nxd * nyd + MD * nxd + col];
        svec[1 * nxd * nyd + 2 * nxd + col] = svec[1 * nxd * nyd + MD * nxd + col];

        Jacb31[0 * nxd + col] = Jacb31[MD * nxd + col];
        Jacb31[1 * nxd + col] = Jacb31[MD * nxd + col];
        Jacb31[2 * nxd + col] = Jacb31[MD * nxd + col];

        Jacb32[0 * nxd + col] = Jacb32[MD * nxd + col];
        Jacb32[1 * nxd + col] = Jacb32[MD * nxd + col];
        Jacb32[2 * nxd + col] = Jacb32[MD * nxd + col];

        invJ11[0 * nxd + col] = invJ11[MD * nxd + col];
        invJ11[1 * nxd + col] = invJ11[MD * nxd + col];
        invJ11[2 * nxd + col] = invJ11[MD * nxd + col];

        invJ12[0 * nxd + col] = invJ12[MD * nxd + col];
        invJ12[1 * nxd + col] = invJ12[MD * nxd + col];
        invJ12[2 * nxd + col] = invJ12[MD * nxd + col];

        invJ13[0 * nxd + col] = invJ13[MD * nxd + col];
        invJ13[1 * nxd + col] = invJ13[MD * nxd + col];
        invJ13[2 * nxd + col] = invJ13[MD * nxd + col];

        invJ21[0 * nxd + col] = invJ21[MD * nxd + col];
        invJ21[1 * nxd + col] = invJ21[MD * nxd + col];
        invJ21[2 * nxd + col] = invJ21[MD * nxd + col];

        invJ22[0 * nxd + col] = invJ22[MD * nxd + col];
        invJ22[1 * nxd + col] = invJ22[MD * nxd + col];
        invJ22[2 * nxd + col] = invJ22[MD * nxd + col];

        invJ23[0 * nxd + col] = invJ23[MD * nxd + col];
        invJ23[1 * nxd + col] = invJ23[MD * nxd + col];
        invJ23[2 * nxd + col] = invJ23[MD * nxd + col];

        invJ31[0 * nxd + col] = invJ31[MD * nxd + col];
        invJ31[1 * nxd + col] = invJ31[MD * nxd + col];
        invJ31[2 * nxd + col] = invJ31[MD * nxd + col];

        invJ32[0 * nxd + col] = invJ32[MD * nxd + col];
        invJ32[1 * nxd + col] = invJ32[MD * nxd + col];
        invJ32[2 * nxd + col] = invJ32[MD * nxd + col];

        invJ33[0 * nxd + col] = invJ33[MD * nxd + col];
        invJ33[1 * nxd + col] = invJ33[MD * nxd + col];
        invJ33[2 * nxd + col] = invJ33[MD * nxd + col];



        dxdxi11[(ny+MD+0) * nxd + col] = dxdxi11[(ny+MD-1) * nxd + col];
        dxdxi11[(ny+MD+1) * nxd + col] = dxdxi11[(ny+MD-1) * nxd + col];
        dxdxi11[(ny+MD+2) * nxd + col] = dxdxi11[(ny+MD-1) * nxd + col];

        dxdxi12[(ny+MD+0) * nxd + col] = dxdxi12[(ny+MD-1) * nxd + col];
        dxdxi12[(ny+MD+1) * nxd + col] = dxdxi12[(ny+MD-1) * nxd + col];
        dxdxi12[(ny+MD+2) * nxd + col] = dxdxi12[(ny+MD-1) * nxd + col];

        dxdxi21[(ny+MD+0) * nxd + col] = dxdxi21[(ny+MD-1) * nxd + col];
        dxdxi21[(ny+MD+1) * nxd + col] = dxdxi21[(ny+MD-1) * nxd + col];
        dxdxi21[(ny+MD+2) * nxd + col] = dxdxi21[(ny+MD-1) * nxd + col];

        dxdxi22[(ny+MD+0) * nxd + col] = dxdxi22[(ny+MD-1) * nxd + col];
        dxdxi22[(ny+MD+1) * nxd + col] = dxdxi22[(ny+MD-1) * nxd + col];
        dxdxi22[(ny+MD+2) * nxd + col] = dxdxi22[(ny+MD-1) * nxd + col];

        cvalue[(ny+MD+0) * nxd + col] = cvalue[(ny+MD-1) * nxd + col];
        cvalue[(ny+MD+1) * nxd + col] = cvalue[(ny+MD-1) * nxd + col];
        cvalue[(ny+MD+2) * nxd + col] = cvalue[(ny+MD-1) * nxd + col];

        Detmin[(ny+MD+0) * nxd + col] = Detmin[(ny+MD-1) * nxd + col];
        Detmin[(ny+MD+1) * nxd + col] = Detmin[(ny+MD-1) * nxd + col];
        Detmin[(ny+MD+2) * nxd + col] = Detmin[(ny+MD-1) * nxd + col];

        svec[0 * nxd * nyd + (ny+MD+0) * nxd + col] = svec[0 * nxd * nyd + (ny+MD-1) * nxd + col];
        svec[0 * nxd * nyd + (ny+MD+1) * nxd + col] = svec[0 * nxd * nyd + (ny+MD-1) * nxd + col];
        svec[0 * nxd * nyd + (ny+MD+2) * nxd + col] = svec[0 * nxd * nyd + (ny+MD-1) * nxd + col];

        svec[1 * nxd * nyd + (ny+MD+0) * nxd + col] = svec[1 * nxd * nyd + (ny+MD-1) * nxd + col];
        svec[1 * nxd * nyd + (ny+MD+1) * nxd + col] = svec[1 * nxd * nyd + (ny+MD-1) * nxd + col];
        svec[1 * nxd * nyd + (ny+MD+2) * nxd + col] = svec[1 * nxd * nyd + (ny+MD-1) * nxd + col];

        Jacb31[(ny+MD+0) * nxd + col] = Jacb31[(ny+MD-1) * nxd + col];
        Jacb31[(ny+MD+1) * nxd + col] = Jacb31[(ny+MD-1) * nxd + col];
        Jacb31[(ny+MD+2) * nxd + col] = Jacb31[(ny+MD-1) * nxd + col];

        Jacb32[(ny+MD+0) * nxd + col] = Jacb32[(ny+MD-1) * nxd + col];
        Jacb32[(ny+MD+1) * nxd + col] = Jacb32[(ny+MD-1) * nxd + col];
        Jacb32[(ny+MD+2) * nxd + col] = Jacb32[(ny+MD-1) * nxd + col];

        invJ11[(ny+MD+0) * nxd + col] = invJ11[(ny+MD-1) * nxd + col];
        invJ11[(ny+MD+1) * nxd + col] = invJ11[(ny+MD-1) * nxd + col];
        invJ11[(ny+MD+2) * nxd + col] = invJ11[(ny+MD-1) * nxd + col];

        invJ12[(ny+MD+0) * nxd + col] = invJ12[(ny+MD-1) * nxd + col];
        invJ12[(ny+MD+1) * nxd + col] = invJ12[(ny+MD-1) * nxd + col];
        invJ12[(ny+MD+2) * nxd + col] = invJ12[(ny+MD-1) * nxd + col];

        invJ13[(ny+MD+0) * nxd + col] = invJ13[(ny+MD-1) * nxd + col];
        invJ13[(ny+MD+1) * nxd + col] = invJ13[(ny+MD-1) * nxd + col];
        invJ13[(ny+MD+2) * nxd + col] = invJ13[(ny+MD-1) * nxd + col];

        invJ21[(ny+MD+0) * nxd + col] = invJ21[(ny+MD-1) * nxd + col];
        invJ21[(ny+MD+1) * nxd + col] = invJ21[(ny+MD-1) * nxd + col];
        invJ21[(ny+MD+2) * nxd + col] = invJ21[(ny+MD-1) * nxd + col];

        invJ22[(ny+MD+0) * nxd + col] = invJ22[(ny+MD-1) * nxd + col];
        invJ22[(ny+MD+1) * nxd + col] = invJ22[(ny+MD-1) * nxd + col];
        invJ22[(ny+MD+2) * nxd + col] = invJ22[(ny+MD-1) * nxd + col];

        invJ23[(ny+MD+0) * nxd + col] = invJ23[(ny+MD-1) * nxd + col];
        invJ23[(ny+MD+1) * nxd + col] = invJ23[(ny+MD-1) * nxd + col];
        invJ23[(ny+MD+2) * nxd + col] = invJ23[(ny+MD-1) * nxd + col];

        invJ31[(ny+MD+0) * nxd + col] = invJ31[(ny+MD-1) * nxd + col];
        invJ31[(ny+MD+1) * nxd + col] = invJ31[(ny+MD-1) * nxd + col];
        invJ31[(ny+MD+2) * nxd + col] = invJ31[(ny+MD-1) * nxd + col];

        invJ32[(ny+MD+0) * nxd + col] = invJ32[(ny+MD-1) * nxd + col];
        invJ32[(ny+MD+1) * nxd + col] = invJ32[(ny+MD-1) * nxd + col];
        invJ32[(ny+MD+2) * nxd + col] = invJ32[(ny+MD-1) * nxd + col];

        invJ33[(ny+MD+0) * nxd + col] = invJ33[(ny+MD-1) * nxd + col];
        invJ33[(ny+MD+1) * nxd + col] = invJ33[(ny+MD-1) * nxd + col];
        invJ33[(ny+MD+2) * nxd + col] = invJ33[(ny+MD-1) * nxd + col];


    }
    // result[row * nxd + col] = u[6 * nxd * nyd + row * nxd + col];

}

__global__
void JacobKernel(double* result,
    double* svec, double* cvalue, 
    double* posx, double* posy, 
    double* J13dxi, double* J23dxi, double* J33dxi,  
    double* J13det, double* J23det, double* J33det,  
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-MD) && col < (nxd-MD) && row >2 && col >2 ) {
   
        J13dxi[row * nxd + col] = (-svec[0 * nxd * nyd + (row  ) * nxd + (col-1)])*( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(  posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)]))
                                 +(-svec[0 * nxd * nyd + (row  ) * nxd + (col  )])*((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])-( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))
                                 -(-svec[0 * nxd * nyd + (row  ) * nxd + (col+1)])*( posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)])*(  posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]));
        
        J23dxi[row * nxd + col] = (-svec[1 * nxd * nyd + (row  ) * nxd + (col-1)])*( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(  posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)]))
                                 +(-svec[1 * nxd * nyd + (row  ) * nxd + (col  )])*((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])-( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))
                                 -(-svec[1 * nxd * nyd + (row  ) * nxd + (col+1)])*( posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)])*(  posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]));
        
        J33dxi[row * nxd + col] = (cvalue[(row  ) * nxd + (col-1)])*( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(  posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)]))
                                 +(cvalue[(row  ) * nxd + (col  )])*((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])-( posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])*(posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]))
                                 -(cvalue[(row  ) * nxd + (col+1)])*( posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col  )])/((posx[(row  ) * nxd + (col-1)] - posx[(row  ) * nxd + (col+1)])*(  posx[(row  ) * nxd + (col  )] - posx[(row  ) * nxd + (col+1)]));
        
        J13det[row * nxd + col] = (-svec[0 * nxd * nyd + (row-1) * nxd + (col  )])*( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(  posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 +(-svec[0 * nxd * nyd + (row  ) * nxd + (col  )])*((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])-( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 -(-svec[0 * nxd * nyd + (row+1) * nxd + (col  )])*( posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])*(  posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]));
         
        J23det[row * nxd + col] = (-svec[1 * nxd * nyd + (row-1) * nxd + (col  )])*( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(  posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 +(-svec[1 * nxd * nyd + (row  ) * nxd + (col  )])*((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])-( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 -(-svec[1 * nxd * nyd + (row+1) * nxd + (col  )])*( posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])*(  posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]));
        
        J33det[row * nxd + col] = (cvalue[(row-1) * nxd + (col  )])*( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(  posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 +(cvalue[(row  ) * nxd + (col  )])*((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])-( posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))/((posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])*(posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]))
                                 -(cvalue[(row+1) * nxd + (col  )])*( posy[(row-1) * nxd + (col  )] - posy[(row  ) * nxd + (col  )])/((posy[(row-1) * nxd + (col  )] - posy[(row+1) * nxd + (col  )])*(  posy[(row  ) * nxd + (col  )] - posy[(row+1) * nxd + (col  )]));
        
    }

}

__global__
void Boundary3Kernel(double* result,
    double* J13dxi, double* J23dxi, double* J33dxi,  
    double* J13det, double* J23det, double* J33det,  
    int nxd, int nyd, int nx, int ny)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && row >=0 && col < (nxd) && col >=0 ) {
   
        J13dxi[row * nxd + 0] = J13dxi[row * nxd + MD];
        J13dxi[row * nxd + 1] = J13dxi[row * nxd + MD];
        J13dxi[row * nxd + 2] = J13dxi[row * nxd + MD]; 

        J23dxi[row * nxd + 0] = J23dxi[row * nxd + MD];
        J23dxi[row * nxd + 1] = J23dxi[row * nxd + MD];
        J23dxi[row * nxd + 2] = J23dxi[row * nxd + MD];

        J33dxi[row * nxd + 0] = J33dxi[row * nxd + MD];
        J33dxi[row * nxd + 1] = J33dxi[row * nxd + MD];
        J33dxi[row * nxd + 2] = J33dxi[row * nxd + MD];

        J13det[row * nxd + 0] = J13det[row * nxd + MD];
        J13det[row * nxd + 1] = J13det[row * nxd + MD];
        J13det[row * nxd + 2] = J13det[row * nxd + MD]; 

        J23det[row * nxd + 0] = J23det[row * nxd + MD];
        J23det[row * nxd + 1] = J23det[row * nxd + MD];
        J23det[row * nxd + 2] = J23det[row * nxd + MD];

        J33det[row * nxd + 0] = J33det[row * nxd + MD];
        J33det[row * nxd + 1] = J33det[row * nxd + MD];
        J33det[row * nxd + 2] = J33det[row * nxd + MD];


        J13dxi[row * nxd + (nx+MD+0)] = J13dxi[row * nxd + (nx+MD-1)];
        J13dxi[row * nxd + (nx+MD+1)] = J13dxi[row * nxd + (nx+MD-1)];
        J13dxi[row * nxd + (nx+MD+2)] = J13dxi[row * nxd + (nx+MD-1)]; 

        J23dxi[row * nxd + (nx+MD+0)] = J23dxi[row * nxd + (nx+MD-1)];
        J23dxi[row * nxd + (nx+MD+1)] = J23dxi[row * nxd + (nx+MD-1)];
        J23dxi[row * nxd + (nx+MD+2)] = J23dxi[row * nxd + (nx+MD-1)]; 

        J33dxi[row * nxd + (nx+MD+0)] = J33dxi[row * nxd + (nx+MD-1)];
        J33dxi[row * nxd + (nx+MD+1)] = J33dxi[row * nxd + (nx+MD-1)];
        J33dxi[row * nxd + (nx+MD+2)] = J33dxi[row * nxd + (nx+MD-1)]; 

        J13det[row * nxd + (nx+MD+0)] = J13det[row * nxd + (nx+MD-1)];
        J13det[row * nxd + (nx+MD+1)] = J13det[row * nxd + (nx+MD-1)];
        J13det[row * nxd + (nx+MD+2)] = J13det[row * nxd + (nx+MD-1)]; 

        J23det[row * nxd + (nx+MD+0)] = J23det[row * nxd + (nx+MD-1)];
        J23det[row * nxd + (nx+MD+1)] = J23det[row * nxd + (nx+MD-1)];
        J23det[row * nxd + (nx+MD+2)] = J23det[row * nxd + (nx+MD-1)]; 

        J33det[row * nxd + (nx+MD+0)] = J33det[row * nxd + (nx+MD-1)];
        J33det[row * nxd + (nx+MD+1)] = J33det[row * nxd + (nx+MD-1)];
        J33det[row * nxd + (nx+MD+2)] = J33det[row * nxd + (nx+MD-1)]; 

        J13dxi[0 * nxd + col] = J13dxi[MD * nxd + col];
        J13dxi[1 * nxd + col] = J13dxi[MD * nxd + col];
        J13dxi[2 * nxd + col] = J13dxi[MD * nxd + col]; 

        J23dxi[0 * nxd + col] = J23dxi[MD * nxd + col];
        J23dxi[1 * nxd + col] = J23dxi[MD * nxd + col];
        J23dxi[2 * nxd + col] = J23dxi[MD * nxd + col];

        J33dxi[0 * nxd + col] = J33dxi[MD * nxd + col];
        J33dxi[1 * nxd + col] = J33dxi[MD * nxd + col];
        J33dxi[2 * nxd + col] = J33dxi[MD * nxd + col];

        J13det[0 * nxd + col] = J13det[MD * nxd + col];
        J13det[1 * nxd + col] = J13det[MD * nxd + col];
        J13det[2 * nxd + col] = J13det[MD * nxd + col]; 

        J23det[0 * nxd + col] = J23det[MD * nxd + col];
        J23det[1 * nxd + col] = J23det[MD * nxd + col];
        J23det[2 * nxd + col] = J23det[MD * nxd + col];

        J33det[0 * nxd + col] = J33det[MD * nxd + col];
        J33det[1 * nxd + col] = J33det[MD * nxd + col];
        J33det[2 * nxd + col] = J33det[MD * nxd + col];


        J13dxi[(ny+MD+0) * nxd + col] = J13dxi[(ny+MD-1) * nxd + col];
        J13dxi[(ny+MD+1) * nxd + col] = J13dxi[(ny+MD-1) * nxd + col];
        J13dxi[(ny+MD+2) * nxd + col] = J13dxi[(ny+MD-1) * nxd + col]; 

        J23dxi[(ny+MD+0) * nxd + col] = J23dxi[(ny+MD-1) * nxd + col];
        J23dxi[(ny+MD+1) * nxd + col] = J23dxi[(ny+MD-1) * nxd + col];
        J23dxi[(ny+MD+2) * nxd + col] = J23dxi[(ny+MD-1) * nxd + col]; 

        J33dxi[(ny+MD+0) * nxd + col] = J33dxi[(ny+MD-1) * nxd + col];
        J33dxi[(ny+MD+1) * nxd + col] = J33dxi[(ny+MD-1) * nxd + col];
        J33dxi[(ny+MD+2) * nxd + col] = J33dxi[(ny+MD-1) * nxd + col]; 

        J13det[(ny+MD+0) * nxd + col] = J13det[(ny+MD-1) * nxd + col];
        J13det[(ny+MD+1) * nxd + col] = J13det[(ny+MD-1) * nxd + col];
        J13det[(ny+MD+2) * nxd + col] = J13det[(ny+MD-1) * nxd + col]; 

        J23det[(ny+MD+0) * nxd + col] = J23det[(ny+MD-1) * nxd + col];
        J23det[(ny+MD+1) * nxd + col] = J23det[(ny+MD-1) * nxd + col];
        J23det[(ny+MD+2) * nxd + col] = J23det[(ny+MD-1) * nxd + col]; 

        J33det[(ny+MD+0) * nxd + col] = J33det[(ny+MD-1) * nxd + col];
        J33det[(ny+MD+1) * nxd + col] = J33det[(ny+MD-1) * nxd + col];
        J33det[(ny+MD+2) * nxd + col] = J33det[(ny+MD-1) * nxd + col]; 

    }

}


__global__
void Inflow1Kernel(double* result,
    double* inflow, double* loc, 
    double* u,double* cvalue, 
    double phiS0,
    int locflowlen, int Iniflowlen, int inflowCount,
    int* dire,
    int nxd, int nyd)
{

    int xloc, yloc;
    int direX, direY;

    for(int i=0;i<locflowlen;i++){

        direX = dire[i * 3 +1];
        direY = dire[i * 3 +2];
        xloc  =  loc[i * 3 +1];
        yloc  =  loc[i * 3 +2];

        __syncthreads();

        for(int mm=0;mm<(direX+MD+3);mm++){
            for(int nn=0;nn<(direY+MD+3);nn++){
                u[0 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[0 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(1-phiS0)/phiS0;
                u[3 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[0 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(1.0-phiS0)*cvalue[(yloc+nn) * nxd + (xloc+mm)];
            
            }

        }
    }


}

__global__
void Inflow2Kernel(double* result,
    double* inflow, double* loc, 
    double* u, 
    int locflowlen, int Iniflowlen, int inflowCount,
    int* dire,
    int nxd, int nyd)
{

    int xloc, yloc;
    int direX, direY;

    for(int i=0;i<locflowlen;i++){

        direX = dire[i * 3 +1];
        direY = dire[i * 3 +2];
        xloc  =  loc[i * 3 +1];
        yloc  =  loc[i * 3 +2];
        __syncthreads();
        for(int mm=0;mm<(direX+MD+3);mm++){
            for(int nn=0;nn<(direY+MD+3);nn++){
                u[1 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[1 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(u[0 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)]+H0);
                u[2 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[2 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(u[0 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)]+H0);
                
                u[4 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[1 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(u[3 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)]+H0);
                u[5 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = inflow[2 * Iniflowlen * locflowlen + inflowCount*locflowlen + i ]*(u[3 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)]+H0);
                 

                u[6 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] = u[0 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)]/(u[0 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] +u[3 * nxd * nyd + (yloc+nn) * nxd + (xloc+mm)] +H0);
            }
        }

    }

}


__global__
void UzeroKernel(double* result,
    double* u, double* uzero,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-MD) && col < (nxd-MD) && row >2 && col >2 ) {
   
        for(int i=0; i<MN; i++){

            uzero[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + row * nxd + col];
        
        }

    }

}

__global__
void Boundary5Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    // int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && row >=0 ) {
        
        for(int i=0; i<MN; i++){

            u[i * nxd * nyd + row * nxd + (MD-0-1)] = 0.5*(u[i * nxd * nyd + row * nxd + (MD-0)] + u[i * nxd * nyd + row * nxd + (MD-0+1)]);
            u[i * nxd * nyd + row * nxd + (MD-1-1)] = 0.5*(u[i * nxd * nyd + row * nxd + (MD-1)] + u[i * nxd * nyd + row * nxd + (MD-1+1)]);
            u[i * nxd * nyd + row * nxd + (MD-2-1)] = 0.5*(u[i * nxd * nyd + row * nxd + (MD-2)] + u[i * nxd * nyd + row * nxd + (MD-2+1)]);
        
            u[i * nxd * nyd + row * nxd + (nx+MD+0)] = 0.5*(u[i * nxd * nyd + row * nxd + (nx+MD+0-1)] + u[i * nxd * nyd + row * nxd + (nx+MD+0-2)]);
            u[i * nxd * nyd + row * nxd + (nx+MD+1)] = 0.5*(u[i * nxd * nyd + row * nxd + (nx+MD+1-1)] + u[i * nxd * nyd + row * nxd + (nx+MD+1-2)]);
            u[i * nxd * nyd + row * nxd + (nx+MD+2)] = 0.5*(u[i * nxd * nyd + row * nxd + (nx+MD+2-1)] + u[i * nxd * nyd + row * nxd + (nx+MD+2-2)]);
        
        }
    }

}

__global__
void Boundary6Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny)
{
    // int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (col < (nxd) && col >=0 ) {
        
        for(int i=0; i<MN; i++){

            u[i * nxd * nyd + (MD-0-1) * nxd + col] = 0.5*(u[i * nxd * nyd + (MD-0) * nxd + col] + u[i * nxd * nyd + (MD-0+1) * nxd + col]);
            u[i * nxd * nyd + (MD-1-1) * nxd + col] = 0.5*(u[i * nxd * nyd + (MD-1) * nxd + col] + u[i * nxd * nyd + (MD-1+1) * nxd + col]);
            u[i * nxd * nyd + (MD-2-1) * nxd + col] = 0.5*(u[i * nxd * nyd + (MD-2) * nxd + col] + u[i * nxd * nyd + (MD-2+1) * nxd + col]);
        
            u[i * nxd * nyd + (ny+MD+0) * nxd + col] = 0.5*(u[i * nxd * nyd + (ny+MD+0-1) * nxd + col] +u[i * nxd * nyd + (ny+MD+0-2) * nxd + col]);
            u[i * nxd * nyd + (ny+MD+1) * nxd + col] = 0.5*(u[i * nxd * nyd + (ny+MD+1-1) * nxd + col] +u[i * nxd * nyd + (ny+MD+1-2) * nxd + col]);
            u[i * nxd * nyd + (ny+MD+2) * nxd + col] = 0.5*(u[i * nxd * nyd + (ny+MD+2-1) * nxd + col] +u[i * nxd * nyd + (ny+MD+2-2) * nxd + col]);
        
        
        }
    }

}

__global__
void Boundary7Kernel(double* result,
    double* u,  
    int nxd, int nyd, int nx, int ny)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && row >=0 && col < (nxd) && col >=0) {
        for(int i=0; i<MD; i++){

            u[1 * nxd * nyd + row * nxd + (MD-i-1)] = 0.0;
            u[2 * nxd * nyd + row * nxd + (MD-i-1)] = 0.0;
            u[4 * nxd * nyd + row * nxd + (MD-i-1)] = 0.0;
            u[5 * nxd * nyd + row * nxd + (MD-i-1)] = 0.0;
            
            u[1 * nxd * nyd + row * nxd + (nx+MD+i)] = 0.0;
            u[2 * nxd * nyd + row * nxd + (nx+MD+i)] = 0.0;
            u[4 * nxd * nyd + row * nxd + (nx+MD+i)] = 0.0;
            u[5 * nxd * nyd + row * nxd + (nx+MD+i)] = 0.0;

            u[1 * nxd * nyd + (MD-i-1) * nxd + col] = 0.0;
            u[2 * nxd * nyd + (MD-i-1) * nxd + col] = 0.0;
            u[4 * nxd * nyd + (MD-i-1) * nxd + col] = 0.0;
            u[5 * nxd * nyd + (MD-i-1) * nxd + col] = 0.0;
            
            u[1 * nxd * nyd + (ny+MD+i) * nxd + col] = 0.0;
            u[2 * nxd * nyd + (ny+MD+i) * nxd + col] = 0.0;
            u[4 * nxd * nyd + (ny+MD+i) * nxd + col] = 0.0;
            u[5 * nxd * nyd + (ny+MD+i) * nxd + col] = 0.0;
    
        }
        
    }

}


__global__
void Boundary9Kernel(double* result,
    double* Hpx, double* Hpy,
    double* Ppx, double* Ppy, 
    double* PDx, double* PDy,
    double*  ux, double*  uy,    
    double* apEW, double* apSN,    
    double* apFEW, double* apFSN,    
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd) && col < (nxd) && row >=0 && col >=0 ) {
   
        for(int i=0; i<MN; i++){
            Hpx[i * nxd * nyd + row * nxd + col] = 0.0;
            Hpy[i * nxd * nyd + row * nxd + col] = 0.0;

            Ppx[i * nxd * nyd + row * nxd + col] = 0.0;
            Ppy[i * nxd * nyd + row * nxd + col] = 0.0;

            PDx[i * nxd * nyd + row * nxd + col] = 0.0;
            PDy[i * nxd * nyd + row * nxd + col] = 0.0;

            ux[i * nxd * nyd + row * nxd + col] = 0.0;
            uy[i * nxd * nyd + row * nxd + col] = 0.0;
            
        }

         apEW[row * nxd + col] = 0.0;
         apSN[row * nxd + col] = 0.0;
        apFEW[row * nxd + col] = 0.0;
        apFSN[row * nxd + col] = 0.0;

    }

}


__global__
void TVD1Kernel(double* result,
    double* u,    
    double* dux, double* duy,    
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-1) && col < (nxd-1) && row >=0 && col >=0 ) {
        for(int i=0; i<MN; i++){
            dux[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row  ) * nxd + (col+1)]-u[i * nxd * nyd + (row  ) * nxd + (col  )];
            duy[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row+1) * nxd + (col  )]-u[i * nxd * nyd + (row  ) * nxd + (col  )];
        }
    }

}

__global__
void TVD2Kernel(double* result,
    double* dux, double* duy,   
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-1) && col < (nxd-1) && row >0 && col >0 ) {
        for(int i=0; i<MN; i++){
            sgnAx[i * nxd * nyd + row * nxd + col] = ((dux[i * nxd * nyd + (row  ) * nxd + (col-1)]) > 0.0) ? 1.0 : -1.0;
            sgnBx[i * nxd * nyd + row * nxd + col] = ((dux[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;

            sgnAy[i * nxd * nyd + row * nxd + col] = ((duy[i * nxd * nyd + (row-1) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;
            sgnBy[i * nxd * nyd + row * nxd + col] = ((duy[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;

        }
    }

}

__global__
void TVD3Kernel(double* result,
    double* dux, double* duy,   
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    double* t1x,   double* t2x, 
    double* t1y,   double* t2y, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-1) && col < (nxd-1) && row >0 && col >0 ) {
        for(int i=0; i<MN; i++){
            t1x[i * nxd * nyd + row * nxd + col] = 0.5*((sgnAx[i * nxd * nyd + row * nxd + col]) + (sgnBx[i * nxd * nyd + row * nxd + col]))
                                                   *min((dux[i * nxd * nyd + (row  ) * nxd + (col-1)]*sgnAx[i * nxd * nyd + row * nxd + col]), (dux[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnBx[i * nxd * nyd + row * nxd + col]));
            
            t2x[i * nxd * nyd + row * nxd + col] = 0.5*(dux[i * nxd * nyd + (row  ) * nxd + (col-1)] + dux[i * nxd * nyd + (row  ) * nxd + (col  )]);

            t1y[i * nxd * nyd + row * nxd + col] = 0.5*((sgnAy[i * nxd * nyd + row * nxd + col]) + (sgnBy[i * nxd * nyd + row * nxd + col]))
                                                   *min((duy[i * nxd * nyd + (row-1) * nxd + (col  )]*sgnAy[i * nxd * nyd + row * nxd + col]), (duy[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnBy[i * nxd * nyd + row * nxd + col]));
            
            t2y[i * nxd * nyd + row * nxd + col] = 0.5*(duy[i * nxd * nyd + (row-1) * nxd + (col  )] + duy[i * nxd * nyd + (row  ) * nxd + (col  )]);

        }
    }

}

__global__
void TVD4Kernel(double* result,
    double* t1x, double* t2x,  
    double* t1y, double* t2y, 
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-1) && col < (nxd-1) && row >0 && col >0 ) {
        for(int i=0; i<MN; i++){
            sgnAx[i * nxd * nyd + row * nxd + col] = ((t1x[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;
            sgnBx[i * nxd * nyd + row * nxd + col] = ((t2x[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;

            sgnAy[i * nxd * nyd + row * nxd + col] = ((t1y[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;
            sgnBy[i * nxd * nyd + row * nxd + col] = ((t2y[i * nxd * nyd + (row  ) * nxd + (col  )]) > 0.0) ? 1.0 : -1.0;

        }
    }

}

__global__
void TVD5Kernel(double* result,
    double* t1x, double* t2x, 
    double* t1y, double* t2y, 
    double* sgnAx, double* sgnBx, 
    double* sgnAy, double* sgnBy,
    double*  ux, double* uy, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd-1) && col < (nxd-1) && row >0 && col >0 ) {
        for(int i=0; i<MN; i++){
            ux[i * nxd * nyd + row * nxd + col] = 0.5*((sgnAx[i * nxd * nyd + row * nxd + col]) + (sgnBx[i * nxd * nyd + row * nxd + col]))
                                                   *min((t1x[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnAx[i * nxd * nyd + row * nxd + col]), (t2x[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnBx[i * nxd * nyd + row * nxd + col]));
            
            uy[i * nxd * nyd + row * nxd + col] = 0.5*((sgnAy[i * nxd * nyd + row * nxd + col]) + (sgnBy[i * nxd * nyd + row * nxd + col]))
                                                   *min((t1y[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnAy[i * nxd * nyd + row * nxd + col]), (t2y[i * nxd * nyd + (row  ) * nxd + (col  )]*sgnBy[i * nxd * nyd + row * nxd + col]));
            
          
        }

    }

}

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
     
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-2) && col < (nxd-2) && row >0 && col >0 ) {
        
        dxdxi11_avgEW[row * nxd + col] = 0.5*(dxdxi11[row * nxd + col] + dxdxi11[(row  ) * nxd + (col+1)]);
        dxdxi21_avgEW[row * nxd + col] = 0.5*(dxdxi21[row * nxd + col] + dxdxi21[(row  ) * nxd + (col+1)]);
        
        dxdxi12_avgSN[row * nxd + col] = 0.5*(dxdxi12[row * nxd + col] + dxdxi12[(row+1) * nxd + (col  )]);
        dxdxi22_avgSN[row * nxd + col] = 0.5*(dxdxi22[row * nxd + col] + dxdxi22[(row+1) * nxd + (col  )]);

        J13dxi_avgEW[row * nxd + col] = 0.5*(J13dxi[row * nxd + col] + J13dxi[(row  ) * nxd + (col+1)]);
        J23dxi_avgEW[row * nxd + col] = 0.5*(J23dxi[row * nxd + col] + J23dxi[(row  ) * nxd + (col+1)]);
        J33dxi_avgEW[row * nxd + col] = 0.5*(J33dxi[row * nxd + col] + J33dxi[(row  ) * nxd + (col+1)]);
        
        J13det_avgEW[row * nxd + col] = 0.5*(J13det[row * nxd + col] + J13det[(row  ) * nxd + (col+1)]);
        J23det_avgEW[row * nxd + col] = 0.5*(J23det[row * nxd + col] + J23det[(row  ) * nxd + (col+1)]);
        J33det_avgEW[row * nxd + col] = 0.5*(J33det[row * nxd + col] + J33det[(row  ) * nxd + (col+1)]);
        
        J13dxi_avgSN[row * nxd + col] = 0.5*(J13dxi[row * nxd + col] + J13dxi[(row+1) * nxd + (col  )]);
        J23dxi_avgSN[row * nxd + col] = 0.5*(J23dxi[row * nxd + col] + J23dxi[(row+1) * nxd + (col  )]);
        J33dxi_avgSN[row * nxd + col] = 0.5*(J33dxi[row * nxd + col] + J33dxi[(row+1) * nxd + (col  )]);
        
        J13det_avgSN[row * nxd + col] = 0.5*(J13det[row * nxd + col] + J13det[(row+1) * nxd + (col  )]);
        J23det_avgSN[row * nxd + col] = 0.5*(J23det[row * nxd + col] + J23det[(row+1) * nxd + (col  )]);
        J33det_avgSN[row * nxd + col] = 0.5*(J33det[row * nxd + col] + J33det[(row+1) * nxd + (col  )]);

        invJ11_avgEW[row * nxd + col] = 0.5*(invJ11[row * nxd + col] + invJ11[(row  ) * nxd + (col+1)]);
        invJ12_avgEW[row * nxd + col] = 0.5*(invJ12[row * nxd + col] + invJ12[(row  ) * nxd + (col+1)]);
        invJ13_avgEW[row * nxd + col] = 0.5*(invJ13[row * nxd + col] + invJ13[(row  ) * nxd + (col+1)]);

        invJ21_avgEW[row * nxd + col] = 0.5*(invJ21[row * nxd + col] + invJ21[(row  ) * nxd + (col+1)]);
        invJ22_avgEW[row * nxd + col] = 0.5*(invJ22[row * nxd + col] + invJ22[(row  ) * nxd + (col+1)]);
        invJ23_avgEW[row * nxd + col] = 0.5*(invJ23[row * nxd + col] + invJ23[(row  ) * nxd + (col+1)]);

        invJ31_avgEW[row * nxd + col] = 0.5*(invJ31[row * nxd + col] + invJ31[(row  ) * nxd + (col+1)]);
        invJ32_avgEW[row * nxd + col] = 0.5*(invJ32[row * nxd + col] + invJ32[(row  ) * nxd + (col+1)]);
        invJ33_avgEW[row * nxd + col] = 0.5*(invJ33[row * nxd + col] + invJ33[(row  ) * nxd + (col+1)]);

        invJ11_avgSN[row * nxd + col] = 0.5*(invJ11[row * nxd + col] + invJ11[(row+1) * nxd + (col  )]);
        invJ12_avgSN[row * nxd + col] = 0.5*(invJ12[row * nxd + col] + invJ12[(row+1) * nxd + (col  )]);
        invJ13_avgSN[row * nxd + col] = 0.5*(invJ13[row * nxd + col] + invJ13[(row+1) * nxd + (col  )]);

        invJ21_avgSN[row * nxd + col] = 0.5*(invJ21[row * nxd + col] + invJ21[(row+1) * nxd + (col  )]);
        invJ22_avgSN[row * nxd + col] = 0.5*(invJ22[row * nxd + col] + invJ22[(row+1) * nxd + (col  )]);
        invJ23_avgSN[row * nxd + col] = 0.5*(invJ23[row * nxd + col] + invJ23[(row+1) * nxd + (col  )]);

        invJ31_avgSN[row * nxd + col] = 0.5*(invJ31[row * nxd + col] + invJ31[(row+1) * nxd + (col  )]);
        invJ32_avgSN[row * nxd + col] = 0.5*(invJ32[row * nxd + col] + invJ32[(row+1) * nxd + (col  )]);
        invJ33_avgSN[row * nxd + col] = 0.5*(invJ33[row * nxd + col] + invJ33[(row+1) * nxd + (col  )]);

        Detmin_avgEW[row * nxd + col] = 0.5*(Detmin[row * nxd + col] + Detmin[(row  ) * nxd + (col+1)]);
        Detmin_avgSN[row * nxd + col] = 0.5*(Detmin[row * nxd + col] + Detmin[(row+1) * nxd + (col  )]);

        cval_avgEW[row * nxd + col] = 0.5*(cvalue[row * nxd + col] + cvalue[(row  ) * nxd + (col+1)]);
        cval_avgSN[row * nxd + col] = 0.5*(cvalue[row * nxd + col] + cvalue[(row+1) * nxd + (col  )]);

        svec_avgEW[0 * nxd * nyd + row * nxd + col] = 0.5*(svec[0 * nxd * nyd + row * nxd + col] + svec[0 * nxd * nyd + (row  ) * nxd + (col+1)]);
        svec_avgEW[1 * nxd * nyd + row * nxd + col] = 0.5*(svec[1 * nxd * nyd + row * nxd + col] + svec[1 * nxd * nyd + (row  ) * nxd + (col+1)]);
        
        svec_avgSN[0 * nxd * nyd + row * nxd + col] = 0.5*(svec[0 * nxd * nyd + row * nxd + col] + svec[0 * nxd * nyd + (row+1) * nxd + (col  )]);
        svec_avgSN[1 * nxd * nyd + row * nxd + col] = 0.5*(svec[1 * nxd * nyd + row * nxd + col] + svec[1 * nxd * nyd + (row+1) * nxd + (col  )]);
        
    }

}

__global__
void InterfacesKernel(double* result,
    double* u, 
    double* ux, double* uy,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row < (nyd-1) && col < (nxd-1) && row >=0 && col >=0 ) {

        for(int i=0; i<MN; i++){
            uE[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row  ) * nxd + (col+1)] - 0.5*ux[i * nxd * nyd + (row  ) * nxd + (col+1)];
            uW[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row  ) * nxd + (col  )] + 0.5*ux[i * nxd * nyd + (row  ) * nxd + (col  )];
        
            uN[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row+1) * nxd + (col  )] - 0.5*uy[i * nxd * nyd + (row+1) * nxd + (col  )];
            uS[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + (row  ) * nxd + (col  )] + 0.5*uy[i * nxd * nyd + (row  ) * nxd + (col  )];
        
        }
         

    }

}

__global__
void Interfaces2Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd, int nx, int ny)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row < (nyd) && row >=0  && col < (nxd) && col >=0) {

        for(int i=0; i<MN; i++){
            uE[i * nxd * nyd + row * nxd + 0] = 0.5*(uE[i * nxd * nyd + row * nxd + MD] + uE[i * nxd * nyd + row * nxd + (MD+1)]);
            uE[i * nxd * nyd + row * nxd + 1] = 0.5*(uE[i * nxd * nyd + row * nxd + MD] + uE[i * nxd * nyd + row * nxd + (MD+1)]);
            uE[i * nxd * nyd + row * nxd + 2] = 0.5*(uE[i * nxd * nyd + row * nxd + MD] + uE[i * nxd * nyd + row * nxd + (MD+1)]);
            
            uW[i * nxd * nyd + row * nxd + 0] = 0.5*(uW[i * nxd * nyd + row * nxd + MD] + uW[i * nxd * nyd + row * nxd + (MD+1)]);
            uW[i * nxd * nyd + row * nxd + 1] = 0.5*(uW[i * nxd * nyd + row * nxd + MD] + uW[i * nxd * nyd + row * nxd + (MD+1)]);
            uW[i * nxd * nyd + row * nxd + 2] = 0.5*(uW[i * nxd * nyd + row * nxd + MD] + uW[i * nxd * nyd + row * nxd + (MD+1)]);
            
            uN[i * nxd * nyd + row * nxd + 0] = 0.5*(uN[i * nxd * nyd + row * nxd + MD] + uN[i * nxd * nyd + row * nxd + (MD+1)]);
            uN[i * nxd * nyd + row * nxd + 1] = 0.5*(uN[i * nxd * nyd + row * nxd + MD] + uN[i * nxd * nyd + row * nxd + (MD+1)]);
            uN[i * nxd * nyd + row * nxd + 2] = 0.5*(uN[i * nxd * nyd + row * nxd + MD] + uN[i * nxd * nyd + row * nxd + (MD+1)]);
            
            uS[i * nxd * nyd + row * nxd + 0] = 0.5*(uS[i * nxd * nyd + row * nxd + MD] + uS[i * nxd * nyd + row * nxd + (MD+1)]);
            uS[i * nxd * nyd + row * nxd + 1] = 0.5*(uS[i * nxd * nyd + row * nxd + MD] + uS[i * nxd * nyd + row * nxd + (MD+1)]);
            uS[i * nxd * nyd + row * nxd + 2] = 0.5*(uS[i * nxd * nyd + row * nxd + MD] + uS[i * nxd * nyd + row * nxd + (MD+1)]);

            uE[i * nxd * nyd + row * nxd + (nx+MD+0)] = 0.5*(uE[i * nxd * nyd + row * nxd + (nx+MD-1)] + uE[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uE[i * nxd * nyd + row * nxd + (nx+MD+1)] = 0.5*(uE[i * nxd * nyd + row * nxd + (nx+MD-1)] + uE[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uE[i * nxd * nyd + row * nxd + (nx+MD+2)] = 0.5*(uE[i * nxd * nyd + row * nxd + (nx+MD-1)] + uE[i * nxd * nyd + row * nxd + (nx+MD-2)]);

            uW[i * nxd * nyd + row * nxd + (nx+MD+0)] = 0.5*(uW[i * nxd * nyd + row * nxd + (nx+MD-1)] + uW[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uW[i * nxd * nyd + row * nxd + (nx+MD+1)] = 0.5*(uW[i * nxd * nyd + row * nxd + (nx+MD-1)] + uW[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uW[i * nxd * nyd + row * nxd + (nx+MD+2)] = 0.5*(uW[i * nxd * nyd + row * nxd + (nx+MD-1)] + uW[i * nxd * nyd + row * nxd + (nx+MD-2)]);

            uN[i * nxd * nyd + row * nxd + (nx+MD+0)] = 0.5*(uN[i * nxd * nyd + row * nxd + (nx+MD-1)] + uN[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uN[i * nxd * nyd + row * nxd + (nx+MD+1)] = 0.5*(uN[i * nxd * nyd + row * nxd + (nx+MD-1)] + uN[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uN[i * nxd * nyd + row * nxd + (nx+MD+2)] = 0.5*(uN[i * nxd * nyd + row * nxd + (nx+MD-1)] + uN[i * nxd * nyd + row * nxd + (nx+MD-2)]);

            uS[i * nxd * nyd + row * nxd + (nx+MD+0)] = 0.5*(uS[i * nxd * nyd + row * nxd + (nx+MD-1)] + uS[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uS[i * nxd * nyd + row * nxd + (nx+MD+1)] = 0.5*(uS[i * nxd * nyd + row * nxd + (nx+MD-1)] + uS[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            uS[i * nxd * nyd + row * nxd + (nx+MD+2)] = 0.5*(uS[i * nxd * nyd + row * nxd + (nx+MD-1)] + uS[i * nxd * nyd + row * nxd + (nx+MD-2)]);
            

            uE[i * nxd * nyd + 0 * nxd + col] = 0.5*(uE[i * nxd * nyd + MD * nxd + col] + uE[i * nxd * nyd + (MD+1) * nxd + col]);
            uE[i * nxd * nyd + 1 * nxd + col] = 0.5*(uE[i * nxd * nyd + MD * nxd + col] + uE[i * nxd * nyd + (MD+1) * nxd + col]);
            uE[i * nxd * nyd + 2 * nxd + col] = 0.5*(uE[i * nxd * nyd + MD * nxd + col] + uE[i * nxd * nyd + (MD+1) * nxd + col]);

            uW[i * nxd * nyd + 0 * nxd + col] = 0.5*(uW[i * nxd * nyd + MD * nxd + col] + uW[i * nxd * nyd + (MD+1) * nxd + col]);
            uW[i * nxd * nyd + 1 * nxd + col] = 0.5*(uW[i * nxd * nyd + MD * nxd + col] + uW[i * nxd * nyd + (MD+1) * nxd + col]);
            uW[i * nxd * nyd + 2 * nxd + col] = 0.5*(uW[i * nxd * nyd + MD * nxd + col] + uW[i * nxd * nyd + (MD+1) * nxd + col]);

            uN[i * nxd * nyd + 0 * nxd + col] = 0.5*(uN[i * nxd * nyd + MD * nxd + col] + uN[i * nxd * nyd + (MD+1) * nxd + col]);
            uN[i * nxd * nyd + 1 * nxd + col] = 0.5*(uN[i * nxd * nyd + MD * nxd + col] + uN[i * nxd * nyd + (MD+1) * nxd + col]);
            uN[i * nxd * nyd + 2 * nxd + col] = 0.5*(uN[i * nxd * nyd + MD * nxd + col] + uN[i * nxd * nyd + (MD+1) * nxd + col]);

            uS[i * nxd * nyd + 0 * nxd + col] = 0.5*(uS[i * nxd * nyd + MD * nxd + col] + uS[i * nxd * nyd + (MD+1) * nxd + col]);
            uS[i * nxd * nyd + 1 * nxd + col] = 0.5*(uS[i * nxd * nyd + MD * nxd + col] + uS[i * nxd * nyd + (MD+1) * nxd + col]);
            uS[i * nxd * nyd + 2 * nxd + col] = 0.5*(uS[i * nxd * nyd + MD * nxd + col] + uS[i * nxd * nyd + (MD+1) * nxd + col]);
          

            uE[i * nxd * nyd + (ny+MD+0) * nxd + col] = 0.5*(uE[i * nxd * nyd + (ny+MD-1) * nxd + col] + uE[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uE[i * nxd * nyd + (ny+MD+1) * nxd + col] = 0.5*(uE[i * nxd * nyd + (ny+MD-1) * nxd + col] + uE[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uE[i * nxd * nyd + (ny+MD+2) * nxd + col] = 0.5*(uE[i * nxd * nyd + (ny+MD-1) * nxd + col] + uE[i * nxd * nyd + (ny+MD-2) * nxd + col]);

            uW[i * nxd * nyd + (ny+MD+0) * nxd + col] = 0.5*(uW[i * nxd * nyd + (ny+MD-1) * nxd + col] + uW[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uW[i * nxd * nyd + (ny+MD+1) * nxd + col] = 0.5*(uW[i * nxd * nyd + (ny+MD-1) * nxd + col] + uW[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uW[i * nxd * nyd + (ny+MD+2) * nxd + col] = 0.5*(uW[i * nxd * nyd + (ny+MD-1) * nxd + col] + uW[i * nxd * nyd + (ny+MD-2) * nxd + col]);

            uN[i * nxd * nyd + (ny+MD+0) * nxd + col] = 0.5*(uN[i * nxd * nyd + (ny+MD-1) * nxd + col] + uN[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uN[i * nxd * nyd + (ny+MD+1) * nxd + col] = 0.5*(uN[i * nxd * nyd + (ny+MD-1) * nxd + col] + uN[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uN[i * nxd * nyd + (ny+MD+2) * nxd + col] = 0.5*(uN[i * nxd * nyd + (ny+MD-1) * nxd + col] + uN[i * nxd * nyd + (ny+MD-2) * nxd + col]);

            uS[i * nxd * nyd + (ny+MD+0) * nxd + col] = 0.5*(uS[i * nxd * nyd + (ny+MD-1) * nxd + col] + uS[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uS[i * nxd * nyd + (ny+MD+1) * nxd + col] = 0.5*(uS[i * nxd * nyd + (ny+MD-1) * nxd + col] + uS[i * nxd * nyd + (ny+MD-2) * nxd + col]);
            uS[i * nxd * nyd + (ny+MD+2) * nxd + col] = 0.5*(uS[i * nxd * nyd + (ny+MD-1) * nxd + col] + uS[i * nxd * nyd + (ny+MD-2) * nxd + col]);

        }
         

    }

}



__global__
void KeepPositivi1Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if ( row < (nyd) && col < (nxd) && col >=0 && row >=0 ) {

        if(uE[0 * nxd * nyd + row * nxd + col] <= H00){
            uE[0 * nxd * nyd + row * nxd + col] = 0.0;
            uE[1 * nxd * nyd + row * nxd + col] = 0.0;
            uE[2 * nxd * nyd + row * nxd + col] = 0.0;
        }
        
        if(uE[3 * nxd * nyd + row * nxd + col] <= H00){
            uE[3 * nxd * nyd + row * nxd + col] = 0.0;
            uE[4 * nxd * nyd + row * nxd + col] = 0.0;
            uE[5 * nxd * nyd + row * nxd + col] = 0.0;
        }

        if(uW[0 * nxd * nyd + row * nxd + col] <= H00){
            uW[0 * nxd * nyd + row * nxd + col] = 0.0;
            uW[1 * nxd * nyd + row * nxd + col] = 0.0;
            uW[2 * nxd * nyd + row * nxd + col] = 0.0;
        }
        
        if(uW[3 * nxd * nyd + row * nxd + col] <= H00){
            uW[3 * nxd * nyd + row * nxd + col] = 0.0;
            uW[4 * nxd * nyd + row * nxd + col] = 0.0;
            uW[5 * nxd * nyd + row * nxd + col] = 0.0;
        }

        if(uN[0 * nxd * nyd + row * nxd + col] <= H00){
            uN[0 * nxd * nyd + row * nxd + col] = 0.0;
            uN[1 * nxd * nyd + row * nxd + col] = 0.0;
            uN[2 * nxd * nyd + row * nxd + col] = 0.0;
        }
        
        if(uN[3 * nxd * nyd + row * nxd + col] <= H00){
            uN[3 * nxd * nyd + row * nxd + col] = 0.0;
            uN[4 * nxd * nyd + row * nxd + col] = 0.0;
            uN[5 * nxd * nyd + row * nxd + col] = 0.0;
        }

        if(uS[0 * nxd * nyd + row * nxd + col] <= H00){
            uS[0 * nxd * nyd + row * nxd + col] = 0.0;
            uS[1 * nxd * nyd + row * nxd + col] = 0.0;
            uS[2 * nxd * nyd + row * nxd + col] = 0.0;
        }
        
        if(uS[3 * nxd * nyd + row * nxd + col] <= H00){
            uS[3 * nxd * nyd + row * nxd + col] = 0.0;
            uS[4 * nxd * nyd + row * nxd + col] = 0.0;
            uS[5 * nxd * nyd + row * nxd + col] = 0.0;
        }



    }

}

__global__
void KeepPositivi2Kernel(double* result,
    double* uE, double* uW,  
    double* uN, double* uS,       
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if ( row < nyd && col < (nxd) &&  col >=0 && row >=0 ) {

        uE[6 * nxd * nyd + row * nxd + col] = uE[0 * nxd * nyd + row * nxd + col]/(uE[0 * nxd * nyd + row * nxd + col] + uE[3 * nxd * nyd + row * nxd + col]+H0);
        uW[6 * nxd * nyd + row * nxd + col] = uW[0 * nxd * nyd + row * nxd + col]/(uW[0 * nxd * nyd + row * nxd + col] + uW[3 * nxd * nyd + row * nxd + col]+H0);
        uS[6 * nxd * nyd + row * nxd + col] = uS[0 * nxd * nyd + row * nxd + col]/(uS[0 * nxd * nyd + row * nxd + col] + uS[3 * nxd * nyd + row * nxd + col]+H0);
        uN[6 * nxd * nyd + row * nxd + col] = uN[0 * nxd * nyd + row * nxd + col]/(uN[0 * nxd * nyd + row * nxd + col] + uN[3 * nxd * nyd + row * nxd + col]+H0);

    }

}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-1) && col < (nxd-1) &&  col >=0 && row >=0 ) {

        vexE[row * nxd + col] =  uE[1 * nxd * nyd + row * nxd + col]/(uE[0 * nxd * nyd + row * nxd + col]+H0);
        veyE[row * nxd + col] =  uE[2 * nxd * nyd + row * nxd + col]/(uE[0 * nxd * nyd + row * nxd + col]+H0);
        
        vexW[row * nxd + col] =  uW[1 * nxd * nyd + row * nxd + col]/(uW[0 * nxd * nyd + row * nxd + col]+H0);
        veyW[row * nxd + col] =  uW[2 * nxd * nyd + row * nxd + col]/(uW[0 * nxd * nyd + row * nxd + col]+H0);

        vexFE[row * nxd + col] =  uE[4 * nxd * nyd + row * nxd + col]/(uE[3 * nxd * nyd + row * nxd + col]+H0);
        veyFE[row * nxd + col] =  uE[5 * nxd * nyd + row * nxd + col]/(uE[3 * nxd * nyd + row * nxd + col]+H0);
        
        vexFW[row * nxd + col] =  uW[4 * nxd * nyd + row * nxd + col]/(uW[3 * nxd * nyd + row * nxd + col]+H0);
        veyFW[row * nxd + col] =  uW[5 * nxd * nyd + row * nxd + col]/(uW[3 * nxd * nyd + row * nxd + col]+H0);
     
        vexN[row * nxd + col] =  uN[1 * nxd * nyd + row * nxd + col]/(uN[0 * nxd * nyd + row * nxd + col]+H0);
        veyN[row * nxd + col] =  uN[2 * nxd * nyd + row * nxd + col]/(uN[0 * nxd * nyd + row * nxd + col]+H0);
        
        vexS[row * nxd + col] =  uS[1 * nxd * nyd + row * nxd + col]/(uS[0 * nxd * nyd + row * nxd + col]+H0);
        veyS[row * nxd + col] =  uS[2 * nxd * nyd + row * nxd + col]/(uS[0 * nxd * nyd + row * nxd + col]+H0);

        vexFN[row * nxd + col] =  uN[4 * nxd * nyd + row * nxd + col]/(uN[3 * nxd * nyd + row * nxd + col]+H0);
        veyFN[row * nxd + col] =  uN[5 * nxd * nyd + row * nxd + col]/(uN[3 * nxd * nyd + row * nxd + col]+H0);
        
        vexFS[row * nxd + col] =  uS[4 * nxd * nyd + row * nxd + col]/(uS[3 * nxd * nyd + row * nxd + col]+H0);
        veyFS[row * nxd + col] =  uS[5 * nxd * nyd + row * nxd + col]/(uS[3 * nxd * nyd + row * nxd + col]+H0);
     
    }

}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {

        w_wertE[row * nxd + col] =  (svec[0 * nxd * nyd + (row  ) * nxd + (col+1)]*vexE[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col+1)]*veyE[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col+1)];

        w_wertW[row * nxd + col] =  (svec[0 * nxd * nyd + (row  ) * nxd + (col  )]*vexW[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col  )]*veyW[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col  )];

        w_wertFE[row * nxd + col] = (svec[0 * nxd * nyd + (row  ) * nxd + (col+1)]*vexFE[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col+1)]*veyFE[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col+1)];

        w_wertFW[row * nxd + col] = (svec[0 * nxd * nyd + (row  ) * nxd + (col  )]*vexFW[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col  )]*veyFW[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col  )];

        w_wertN[row * nxd + col] =  (svec[0 * nxd * nyd + (row+1) * nxd + (col  )]*vexN[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row+1) * nxd + (col  )]*veyN[row * nxd + col])
                                  /cvalue[(row+1) * nxd + (col  )];

        w_wertS[row * nxd + col] =  (svec[0 * nxd * nyd + (row  ) * nxd + (col  )]*vexS[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col  )]*veyS[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col  )];

        w_wertFN[row * nxd + col] = (svec[0 * nxd * nyd + (row+1) * nxd + (col  )]*vexFN[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row+1) * nxd + (col  )]*veyFN[row * nxd + col])
                                  /cvalue[(row+1) * nxd + (col  )];

        w_wertFS[row * nxd + col] = (svec[0 * nxd * nyd + (row  ) * nxd + (col  )]*vexFS[row * nxd + col] 
                                   + svec[1 * nxd * nyd + (row  ) * nxd + (col  )]*veyFS[row * nxd + col])
                                  /cvalue[(row  ) * nxd + (col  )];

       
    }

}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {

        q_xiE[row * nxd + col] = invJ11_avgEW[row * nxd + col]*vexE[row * nxd + col] 
                               + invJ12_avgEW[row * nxd + col]*veyE[row * nxd + col] 
                               + invJ13_avgEW[row * nxd + col]*w_wertE[row * nxd + col];

        q_etE[row * nxd + col] = invJ21_avgEW[row * nxd + col]*vexE[row * nxd + col] 
                               + invJ22_avgEW[row * nxd + col]*veyE[row * nxd + col] 
                               + invJ23_avgEW[row * nxd + col]*w_wertE[row * nxd + col];

        q_xiW[row * nxd + col] = invJ11_avgEW[row * nxd + col]*vexW[row * nxd + col] 
                               + invJ12_avgEW[row * nxd + col]*veyW[row * nxd + col] 
                               + invJ13_avgEW[row * nxd + col]*w_wertW[row * nxd + col];

        q_etW[row * nxd + col] = invJ21_avgEW[row * nxd + col]*vexW[row * nxd + col] 
                               + invJ22_avgEW[row * nxd + col]*veyW[row * nxd + col] 
                               + invJ23_avgEW[row * nxd + col]*w_wertW[row * nxd + col];

        q_xiFE[row * nxd + col] = invJ11_avgEW[row * nxd + col]*vexFE[row * nxd + col] 
                                + invJ12_avgEW[row * nxd + col]*veyFE[row * nxd + col] 
                                + invJ13_avgEW[row * nxd + col]*w_wertFE[row * nxd + col];

        q_etFE[row * nxd + col] = invJ21_avgEW[row * nxd + col]*vexFE[row * nxd + col] 
                                + invJ22_avgEW[row * nxd + col]*veyFE[row * nxd + col] 
                                + invJ23_avgEW[row * nxd + col]*w_wertFE[row * nxd + col];

        q_xiFW[row * nxd + col] = invJ11_avgEW[row * nxd + col]*vexFW[row * nxd + col] 
                                + invJ12_avgEW[row * nxd + col]*veyFW[row * nxd + col] 
                                + invJ13_avgEW[row * nxd + col]*w_wertFW[row * nxd + col];

        q_etFW[row * nxd + col] = invJ21_avgEW[row * nxd + col]*vexFW[row * nxd + col] 
                                + invJ22_avgEW[row * nxd + col]*veyFW[row * nxd + col] 
                                + invJ23_avgEW[row * nxd + col]*w_wertFW[row * nxd + col];

        NpressFE[row * nxd + col] = (uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col]) 
                                   * cval_avgEW[row * nxd + col];

        NpressFW[row * nxd + col] = (uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col]) 
                                   * cval_avgEW[row * nxd + col];

        M11EW[row * nxd + col] = invJ11_avgEW[row * nxd + col]*invJ11_avgEW[row * nxd + col] 
                               + invJ12_avgEW[row * nxd + col]*invJ12_avgEW[row * nxd + col] 
                               + invJ13_avgEW[row * nxd + col]*invJ13_avgEW[row * nxd + col];


        q_xiN[row * nxd + col] = invJ11_avgSN[row * nxd + col]*vexN[row * nxd + col] 
                               + invJ12_avgSN[row * nxd + col]*veyN[row * nxd + col] 
                               + invJ13_avgSN[row * nxd + col]*w_wertN[row * nxd + col];

        q_etN[row * nxd + col] = invJ21_avgSN[row * nxd + col]*vexN[row * nxd + col] 
                               + invJ22_avgSN[row * nxd + col]*veyN[row * nxd + col] 
                               + invJ23_avgSN[row * nxd + col]*w_wertN[row * nxd + col];

        q_xiS[row * nxd + col] = invJ11_avgSN[row * nxd + col]*vexS[row * nxd + col] 
                               + invJ12_avgSN[row * nxd + col]*veyS[row * nxd + col] 
                               + invJ13_avgSN[row * nxd + col]*w_wertS[row * nxd + col];

        q_etS[row * nxd + col] = invJ21_avgSN[row * nxd + col]*vexS[row * nxd + col] 
                               + invJ22_avgSN[row * nxd + col]*veyS[row * nxd + col] 
                               + invJ23_avgSN[row * nxd + col]*w_wertS[row * nxd + col];

        q_xiFN[row * nxd + col] = invJ11_avgSN[row * nxd + col]*vexFN[row * nxd + col] 
                                + invJ12_avgSN[row * nxd + col]*veyFN[row * nxd + col] 
                                + invJ13_avgSN[row * nxd + col]*w_wertFN[row * nxd + col];

        q_etFN[row * nxd + col] = invJ21_avgSN[row * nxd + col]*vexFN[row * nxd + col] 
                                + invJ22_avgSN[row * nxd + col]*veyFN[row * nxd + col] 
                                + invJ23_avgSN[row * nxd + col]*w_wertFN[row * nxd + col];

        q_xiFS[row * nxd + col] = invJ11_avgSN[row * nxd + col]*vexFS[row * nxd + col] 
                                + invJ12_avgSN[row * nxd + col]*veyFS[row * nxd + col] 
                                + invJ13_avgSN[row * nxd + col]*w_wertFS[row * nxd + col];

        q_etFS[row * nxd + col] = invJ21_avgSN[row * nxd + col]*vexFS[row * nxd + col] 
                                + invJ22_avgSN[row * nxd + col]*veyFS[row * nxd + col] 
                                + invJ23_avgSN[row * nxd + col]*w_wertFS[row * nxd + col];

        NpressFN[row * nxd + col] = (uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col]) 
                                   * cval_avgSN[row * nxd + col];

        NpressFS[row * nxd + col] = (uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col]) 
                                   * cval_avgSN[row * nxd + col];

        M22SN[row * nxd + col] = invJ21_avgSN[row * nxd + col]*invJ21_avgSN[row * nxd + col] 
                               + invJ22_avgSN[row * nxd + col]*invJ22_avgSN[row * nxd + col] 
                               + invJ23_avgSN[row * nxd + col]*invJ23_avgSN[row * nxd + col];

    }

}



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


    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
        
    double Drtio1 = 1.0 - (1.420/2.60);
   

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        apE[row * nxd + col] =  fabs(q_xiE[row * nxd + col]) 
                            + invJ11_avgEW[row * nxd + col]*sqrt(NpressFE[row * nxd + col]*Drtio1);
        
        apW[row * nxd + col] =  fabs(q_xiW[row * nxd + col]) 
                            + invJ11_avgEW[row * nxd + col]*sqrt(NpressFW[row * nxd + col]*Drtio1);

        apFE[row * nxd + col] =  fabs(q_xiFE[row * nxd + col]) 
                            + invJ11_avgEW[row * nxd + col]*sqrt(NpressFE[row * nxd + col]/(1.0 - uE[6 * nxd * nyd + row * nxd + col] + H0));
        
        apFW[row * nxd + col] =  fabs(q_xiFW[row * nxd + col]) 
                            + invJ11_avgEW[row * nxd + col]*sqrt(NpressFW[row * nxd + col]/(1.0 - uW[6 * nxd * nyd + row * nxd + col] + H0));
        
        apN[row * nxd + col] =  fabs(q_etN[row * nxd + col]) 
                            + invJ22_avgSN[row * nxd + col]*sqrt(NpressFN[row * nxd + col]*Drtio1);
        
        apS[row * nxd + col] =  fabs(q_etS[row * nxd + col]) 
                            + invJ22_avgSN[row * nxd + col]*sqrt(NpressFS[row * nxd + col]*Drtio1);
        
        apFN[row * nxd + col] =  fabs(q_etFN[row * nxd + col]) 
                            + invJ22_avgSN[row * nxd + col]*sqrt(NpressFN[row * nxd + col]/(1.0 - uN[6 * nxd * nyd + row * nxd + col] + H0));
       
        apFS[row * nxd + col] =  fabs(q_etFS[row * nxd + col]) 
                            + invJ22_avgSN[row * nxd + col]*sqrt(NpressFS[row * nxd + col]/(1.0 - uS[6 * nxd * nyd + row * nxd + col] + H0));
       
    }

}


__global__
void Flux5Kernel(double* result,
    double* apE, double* apW,
    double* apFE, double* apFW, 

    double* apEW, double* apFEW,

    double* apN, double* apS,
    double* apFN, double* apFS,

    double* apSN, double* apFSN,

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
         apEW[row * nxd + col] = max( apE[row * nxd + col],  apW[row * nxd + col]);
        apFEW[row * nxd + col] = max(apFE[row * nxd + col], apFW[row * nxd + col]);
        
         apSN[row * nxd + col] = max( apN[row * nxd + col],  apS[row * nxd + col]);
        apFSN[row * nxd + col] = max(apFN[row * nxd + col], apFS[row * nxd + col]);

    }

}

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
    
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    double Drtio1 = 1.0 - (1.420/2.60);
    double epsilon0 = 1.0;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        em_x[row * nxd + col] = max(1.0e-10, apEW[row * nxd + col]);
        em_y[row * nxd + col] = max(1.0e-10, apSN[row * nxd + col]);
        
        em_Fx[row * nxd + col] = max(1.0e-10, apFEW[row * nxd + col]);
        em_Fy[row * nxd + col] = max(1.0e-10, apFSN[row * nxd + col]);

        czw1x[row * nxd + col] = uE[0 * nxd * nyd + row * nxd + col] 
                              * (uE[0 * nxd * nyd + row * nxd + col] + uE[3 * nxd * nyd + row * nxd + col] ) 
                              * ( 0.5*epsilon0*cval_avgEW[row * nxd + col] *Detmin_avgEW[row * nxd + col]*M11EW[row * nxd + col] ) 
                              * Drtio1;

        czw2x[row * nxd + col] = uW[0 * nxd * nyd + row * nxd + col] 
                              * (uW[0 * nxd * nyd + row * nxd + col] + uW[3 * nxd * nyd + row * nxd + col] ) 
                              * ( 0.5*epsilon0*cval_avgEW[row * nxd + col] *Detmin_avgEW[row * nxd + col]*M11EW[row * nxd + col] ) 
                              * Drtio1;

        czwF1x[row * nxd + col] =  (uE[0 * nxd * nyd + row * nxd + col] + uE[3 * nxd * nyd + row * nxd + col] ) 
                                 * (uE[0 * nxd * nyd + row * nxd + col] + uE[3 * nxd * nyd + row * nxd + col] ) 
                                 * ( 0.5*epsilon0*cval_avgEW[row * nxd + col] *Detmin_avgEW[row * nxd + col]*M11EW[row * nxd + col] );
        
        czwF2x[row * nxd + col] =  (uW[0 * nxd * nyd + row * nxd + col] + uW[3 * nxd * nyd + row * nxd + col] ) 
                                 * (uW[0 * nxd * nyd + row * nxd + col] + uW[3 * nxd * nyd + row * nxd + col] ) 
                                 * ( 0.5*epsilon0*cval_avgEW[row * nxd + col] *Detmin_avgEW[row * nxd + col]*M11EW[row * nxd + col] );
        
        czw1y[row * nxd + col] = uN[0 * nxd * nyd + row * nxd + col] 
                              * (uN[0 * nxd * nyd + row * nxd + col] + uN[3 * nxd * nyd + row * nxd + col] ) 
                              * ( 0.5*epsilon0*cval_avgSN[row * nxd + col] *Detmin_avgSN[row * nxd + col]*M22SN[row * nxd + col] ) 
                              * Drtio1;

        czw2y[row * nxd + col] = uS[0 * nxd * nyd + row * nxd + col] 
                              * (uS[0 * nxd * nyd + row * nxd + col] + uS[3 * nxd * nyd + row * nxd + col] ) 
                              * ( 0.5*epsilon0*cval_avgSN[row * nxd + col] *Detmin_avgSN[row * nxd + col]*M22SN[row * nxd + col] ) 
                              * Drtio1;

        czwF1y[row * nxd + col] =  (uN[0 * nxd * nyd + row * nxd + col] + uN[3 * nxd * nyd + row * nxd + col] ) 
                                 * (uN[0 * nxd * nyd + row * nxd + col] + uN[3 * nxd * nyd + row * nxd + col] ) 
                                 * ( 0.5*epsilon0*cval_avgSN[row * nxd + col] *Detmin_avgSN[row * nxd + col]*M22SN[row * nxd + col] );
        
        czwF2y[row * nxd + col] =  (uS[0 * nxd * nyd + row * nxd + col] + uS[3 * nxd * nyd + row * nxd + col] ) 
                                 * (uS[0 * nxd * nyd + row * nxd + col] + uS[3 * nxd * nyd + row * nxd + col] ) 
                                 * ( 0.5*epsilon0*cval_avgSN[row * nxd + col] *Detmin_avgSN[row * nxd + col]*M22SN[row * nxd + col] );
        

    }

}

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

    
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        FpE[0 * nxd * nyd + row * nxd + col] = uE[0 * nxd * nyd + row * nxd + col]*q_xiE[row * nxd + col]*Detmin_avgEW[row * nxd + col];
        
        FpE[1 * nxd * nyd + row * nxd + col] = uE[1 * nxd * nyd + row * nxd + col]*q_xiE[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czw1x[row * nxd + col]*dxdxi11_avgEW[row * nxd + col];
        
        FpE[2 * nxd * nyd + row * nxd + col] = uE[2 * nxd * nyd + row * nxd + col]*q_xiE[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czw1x[row * nxd + col]*dxdxi21_avgEW[row * nxd + col];

        FpE[3 * nxd * nyd + row * nxd + col] = uE[3 * nxd * nyd + row * nxd + col]*q_xiFE[row * nxd + col]*Detmin_avgEW[row * nxd + col];
        
        FpE[4 * nxd * nyd + row * nxd + col] = uE[4 * nxd * nyd + row * nxd + col]*q_xiFE[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czwF1x[row * nxd + col]*dxdxi11_avgEW[row * nxd + col];
        
        FpE[5 * nxd * nyd + row * nxd + col] = uE[5 * nxd * nyd + row * nxd + col]*q_xiFE[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czwF1x[row * nxd + col]*dxdxi21_avgEW[row * nxd + col];
     

        FpW[0 * nxd * nyd + row * nxd + col] = uW[0 * nxd * nyd + row * nxd + col]*q_xiW[row * nxd + col]*Detmin_avgEW[row * nxd + col];
        
        FpW[1 * nxd * nyd + row * nxd + col] = uW[1 * nxd * nyd + row * nxd + col]*q_xiW[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czw2x[row * nxd + col]*dxdxi11_avgEW[row * nxd + col];
        
        FpW[2 * nxd * nyd + row * nxd + col] = uW[2 * nxd * nyd + row * nxd + col]*q_xiW[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czw2x[row * nxd + col]*dxdxi21_avgEW[row * nxd + col];

        FpW[3 * nxd * nyd + row * nxd + col] = uW[3 * nxd * nyd + row * nxd + col]*q_xiFW[row * nxd + col]*Detmin_avgEW[row * nxd + col];
        
        FpW[4 * nxd * nyd + row * nxd + col] = uW[4 * nxd * nyd + row * nxd + col]*q_xiFW[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czwF2x[row * nxd + col]*dxdxi11_avgEW[row * nxd + col];
        
        FpW[5 * nxd * nyd + row * nxd + col] = uW[5 * nxd * nyd + row * nxd + col]*q_xiFW[row * nxd + col]*Detmin_avgEW[row * nxd + col] + czwF2x[row * nxd + col]*dxdxi21_avgEW[row * nxd + col];
        

        GpN[0 * nxd * nyd + row * nxd + col] = uN[0 * nxd * nyd + row * nxd + col]*q_etN[row * nxd + col]*Detmin_avgSN[row * nxd + col];
        
        GpN[1 * nxd * nyd + row * nxd + col] = uN[1 * nxd * nyd + row * nxd + col]*q_etN[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czw1y[row * nxd + col]*dxdxi12_avgSN[row * nxd + col];
        
        GpN[2 * nxd * nyd + row * nxd + col] = uN[2 * nxd * nyd + row * nxd + col]*q_etN[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czw1y[row * nxd + col]*dxdxi22_avgSN[row * nxd + col];
        
        GpN[3 * nxd * nyd + row * nxd + col] = uN[3 * nxd * nyd + row * nxd + col]*q_etFN[row * nxd + col]*Detmin_avgSN[row * nxd + col];
        
        GpN[4 * nxd * nyd + row * nxd + col] = uN[4 * nxd * nyd + row * nxd + col]*q_etFN[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czwF1y[row * nxd + col]*dxdxi12_avgSN[row * nxd + col];
        
        GpN[5 * nxd * nyd + row * nxd + col] = uN[5 * nxd * nyd + row * nxd + col]*q_etFN[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czwF1y[row * nxd + col]*dxdxi22_avgSN[row * nxd + col];
   

        GpS[0 * nxd * nyd + row * nxd + col] = uS[0 * nxd * nyd + row * nxd + col]*q_etS[row * nxd + col]*Detmin_avgSN[row * nxd + col];
        
        GpS[1 * nxd * nyd + row * nxd + col] = uS[1 * nxd * nyd + row * nxd + col]*q_etS[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czw2y[row * nxd + col]*dxdxi12_avgSN[row * nxd + col];
        
        GpS[2 * nxd * nyd + row * nxd + col] = uS[2 * nxd * nyd + row * nxd + col]*q_etS[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czw2y[row * nxd + col]*dxdxi22_avgSN[row * nxd + col];

        GpS[3 * nxd * nyd + row * nxd + col] = uS[3 * nxd * nyd + row * nxd + col]*q_etFS[row * nxd + col]*Detmin_avgSN[row * nxd + col];
        
        GpS[4 * nxd * nyd + row * nxd + col] = uS[4 * nxd * nyd + row * nxd + col]*q_etFS[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czwF2y[row * nxd + col]*dxdxi12_avgSN[row * nxd + col];
        
        GpS[5 * nxd * nyd + row * nxd + col] = uS[5 * nxd * nyd + row * nxd + col]*q_etFS[row * nxd + col]*Detmin_avgSN[row * nxd + col] + czwF2y[row * nxd + col]*dxdxi22_avgSN[row * nxd + col];
   
    }

}

__global__
void CFL1Kernel(double* result,
    double* em_x, double* em_y,
    double* em_Fx, double* em_Fy,
    double* em_valS, double* em_valF,
    double dx, double dy,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        em_valS[row * nxd + col] =  max( em_x[row * nxd + col]/dx ,  em_y[row * nxd + col]/dy);
        em_valF[row * nxd + col] =  max(em_Fx[row * nxd + col]/dx , em_Fy[row * nxd + col]/dy);
    
    }

}

__global__
void CFL2Kernel(double* result,
    double* em_valS, double* em_valF,
    double* Val,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        Val[row * nxd + col] =  max(em_valS[row * nxd + col] ,  em_valF[row * nxd + col]);
        
    }

}


__global__
void reduceKernel(double *waveSpeed, double* d_max, int size)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    double maxTmp = 0.0;

    __shared__ double cache[1024];

    //reduce multiple elements per thread
    for (int i = index; i < size; i += stride) {
        maxTmp = fmaxf(maxTmp, waveSpeed[i]);
    }

    cache[threadIdx.x] = maxTmp;
    __syncthreads();

    for (int i = blockDim.x / 2; i > 0; i /= 2) {
        if (threadIdx.x < i) {
            cache[threadIdx.x] = fmaxf(cache[threadIdx.x], cache[threadIdx.x + i]);
        }   
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        d_max[blockIdx.x] = cache[0];
    }
}

__global__
void CFL3Kernel(double* result,
    double* dtval, double* d_max, double* totalTime)
{

    (*dtval) = CFL / (*d_max); 
    __syncthreads();

    *dtval = min(2.0,*dtval);
    __syncthreads();

    *dtval = max(0.002,*dtval);
    __syncthreads();

}

__global__
void CFL4Kernel(double* result,
    double* dt, double* dtval, double* totalTime, int io)
{          

    *dt   = *dtval/(2.0-io);
    __syncthreads();

    if(io==1){
        (*totalTime) = (*dt) + (*totalTime);
        __syncthreads();
    }
        __syncthreads();
    
}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    // double Drtio1 = 1.0 - (1.420/2.60);
    double epsilon0 = 1.0;
    double Dratio   = 1.420/2.60;

    if ( row<(nyd-1) && col < (nxd-1) && col >0 && row >0 ) {

        Hpx[0 * nxd * nyd + row * nxd + col] = 0.5*(FpE[0 * nxd * nyd + row * nxd + col] + FpW[0 * nxd * nyd + row * nxd + col]
                                                 - apEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[0 * nxd * nyd + row * nxd + col] -  uW[0 * nxd * nyd + row * nxd + col]));

        Hpx[1 * nxd * nyd + row * nxd + col] = 0.5*(FpE[1 * nxd * nyd + row * nxd + col] + FpW[1 * nxd * nyd + row * nxd + col]
                                                 - apEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[1 * nxd * nyd + row * nxd + col] -  uW[1 * nxd * nyd + row * nxd + col]));

        Hpx[2 * nxd * nyd + row * nxd + col] = 0.5*(FpE[2 * nxd * nyd + row * nxd + col] + FpW[2 * nxd * nyd + row * nxd + col]
                                                 - apEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[2 * nxd * nyd + row * nxd + col] -  uW[2 * nxd * nyd + row * nxd + col]));
        
        Hpx[3 * nxd * nyd + row * nxd + col] = 0.5*(FpE[3 * nxd * nyd + row * nxd + col] + FpW[3 * nxd * nyd + row * nxd + col]
                                                - apFEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[3 * nxd * nyd + row * nxd + col] -  uW[3 * nxd * nyd + row * nxd + col]));               
                                                   
        Hpx[4 * nxd * nyd + row * nxd + col] = 0.5*(FpE[4 * nxd * nyd + row * nxd + col] + FpW[4 * nxd * nyd + row * nxd + col]
                                                - apFEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[4 * nxd * nyd + row * nxd + col] -  uW[4 * nxd * nyd + row * nxd + col]));               
        
        Hpx[5 * nxd * nyd + row * nxd + col] = 0.5*(FpE[5 * nxd * nyd + row * nxd + col] + FpW[5 * nxd * nyd + row * nxd + col]
                                                - apFEW[row * nxd + col]*Detmin_avgEW[row * nxd + col]
                                                   *(uE[5 * nxd * nyd + row * nxd + col] -  uW[5 * nxd * nyd + row * nxd + col]));               
        
                                                   
        Hpy[0 * nxd * nyd + row * nxd + col] =  0.5*(GpN[0 * nxd * nyd + row * nxd + col] + GpS[0 * nxd * nyd + row * nxd + col]
                                                  - apSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[0 * nxd * nyd + row * nxd + col] -  uS[0 * nxd * nyd + row * nxd + col]));

        Hpy[1 * nxd * nyd + row * nxd + col] =  0.5*(GpN[1 * nxd * nyd + row * nxd + col] + GpS[1 * nxd * nyd + row * nxd + col]
                                                  - apSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[1 * nxd * nyd + row * nxd + col] -  uS[1 * nxd * nyd + row * nxd + col]));

        Hpy[2 * nxd * nyd + row * nxd + col] =  0.5*(GpN[2 * nxd * nyd + row * nxd + col] + GpS[2 * nxd * nyd + row * nxd + col]
                                                  - apSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[2 * nxd * nyd + row * nxd + col] -  uS[2 * nxd * nyd + row * nxd + col]));
    
        Hpy[3 * nxd * nyd + row * nxd + col] =  0.5*(GpN[3 * nxd * nyd + row * nxd + col] + GpS[3 * nxd * nyd + row * nxd + col]
                                                 - apFSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[3 * nxd * nyd + row * nxd + col] -  uS[3 * nxd * nyd + row * nxd + col]));
 
        Hpy[4 * nxd * nyd + row * nxd + col] =  0.5*(GpN[4 * nxd * nyd + row * nxd + col] + GpS[4 * nxd * nyd + row * nxd + col]
                                                 - apFSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[4 * nxd * nyd + row * nxd + col] -  uS[4 * nxd * nyd + row * nxd + col]));

        Hpy[5 * nxd * nyd + row * nxd + col] =  0.5*(GpN[5 * nxd * nyd + row * nxd + col] + GpS[5 * nxd * nyd + row * nxd + col]
                                                 - apFSN[row * nxd + col]*Detmin_avgSN[row * nxd + col]
                                                    *(uN[5 * nxd * nyd + row * nxd + col] -  uS[5 * nxd * nyd + row * nxd + col]));

 
        Ppx[0 * nxd * nyd + row * nxd + col] =  0.0;

        Ppx[1 * nxd * nyd + row * nxd + col] = -Dratio*(0.5*epsilon0*cval_avgEW[row * nxd + col]*Detmin_avgEW[row * nxd + col])
                                             *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                             *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                             *invJ11_avgEW[row * nxd + col]; 

        Ppx[2 * nxd * nyd + row * nxd + col] = -Dratio*(0.5*epsilon0*cval_avgEW[row * nxd + col]*Detmin_avgEW[row * nxd + col])
                                             *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                             *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                             *invJ12_avgEW[row * nxd + col]; 

        Ppx[3 * nxd * nyd + row * nxd + col] = 0.0;
        
        Ppx[4 * nxd * nyd + row * nxd + col] = (0.5*epsilon0*cval_avgEW[row * nxd + col]*Detmin_avgEW[row * nxd + col])
                                              *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                              *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                              *invJ11_avgEW[row * nxd + col];
  
        Ppx[5 * nxd * nyd + row * nxd + col] = (0.5*epsilon0*cval_avgEW[row * nxd + col]*Detmin_avgEW[row * nxd + col])
                                              *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                              *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))
                                              *invJ12_avgEW[row * nxd + col];

    
        Ppy[0 * nxd * nyd + row * nxd + col] = 0.0;

        Ppy[1 * nxd * nyd + row * nxd + col] = -Dratio*(0.5*epsilon0*cval_avgSN[row * nxd + col]*Detmin_avgSN[row * nxd + col])
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *invJ21_avgSN[row * nxd + col];


        Ppy[2 * nxd * nyd + row * nxd + col] = -Dratio*(0.5*epsilon0*cval_avgSN[row * nxd + col]*Detmin_avgSN[row * nxd + col])
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *invJ22_avgSN[row * nxd + col];

        Ppy[3 * nxd * nyd + row * nxd + col] = 0.0;

        Ppy[4 * nxd * nyd + row * nxd + col] = (0.5*epsilon0*cval_avgSN[row * nxd + col]*Detmin_avgSN[row * nxd + col])
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *invJ21_avgSN[row * nxd + col];

        Ppy[5 * nxd * nyd + row * nxd + col] = (0.5*epsilon0*cval_avgSN[row * nxd + col]*Detmin_avgSN[row * nxd + col])
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col])))
                                                *invJ22_avgSN[row * nxd + col];


        dudxE[row * nxd + col] = ( u[4 * nxd * nyd + (row  ) * nxd + (col+1)]/(u[3 * nxd * nyd + (row  ) * nxd + (col+1)]+H0) 
                                -  u[4 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dx;

        dvdxE[row * nxd + col] = ( u[5 * nxd * nyd + (row  ) * nxd + (col+1)]/(u[3 * nxd * nyd + (row  ) * nxd + (col+1)]+H0) 
                                -  u[5 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dx;

        dudyE[row * nxd + col] = 0.5*(uy[4 * nxd * nyd + (row  ) * nxd + (col+1)]/(u[3 * nxd * nyd + (row  ) * nxd + (col+1)]+H0)
                                    + uy[4 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dy;
          
        dvdyE[row * nxd + col] = 0.5*(uy[5 * nxd * nyd + (row  ) * nxd + (col+1)]/(u[3 * nxd * nyd + (row  ) * nxd + (col+1)]+H0)
                                    + uy[5 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dy;

        
        dudxN[row * nxd + col] = 0.5*(ux[4 * nxd * nyd + (row+1) * nxd + (col  )]/(u[3 * nxd * nyd + (row+1) * nxd + (col  )]+H0)
                                    + ux[4 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dx;
        
        dvdxN[row * nxd + col] = 0.5*(ux[5 * nxd * nyd + (row+1) * nxd + (col  )]/(u[3 * nxd * nyd + (row+1) * nxd + (col  )]+H0)
                                    + ux[5 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dx;

        dudyN[row * nxd + col] = ( u[4 * nxd * nyd + (row+1) * nxd + (col  )]/(u[3 * nxd * nyd + (row+1) * nxd + (col  )]+H0)
                                -  u[4 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dy;
        
        dvdyN[row * nxd + col] = ( u[5 * nxd * nyd + (row+1) * nxd + (col  )]/(u[3 * nxd * nyd + (row+1) * nxd + (col  )]+H0)
                                -  u[5 * nxd * nyd + (row  ) * nxd + (col  )]/(u[3 * nxd * nyd + (row  ) * nxd + (col  )]+H0))/dy;



    }

}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-1) && col < (nxd-1) && col >0 && row >0 ) {
        
        duxidxix[row * nxd + col] = invJ11_avgEW[row * nxd + col]*dudxE[row * nxd + col] 
                                  + invJ12_avgEW[row * nxd + col]*dvdxE[row * nxd + col] 
                                  - invJ13_avgEW[row * nxd + col]*(invJ31_avgEW[row * nxd + col]*dudxE[row * nxd + col]+invJ32_avgEW[row * nxd + col]*dvdxE[row * nxd + col])/(invJ33_avgEW[row * nxd + col]);

        dvetdxix[row * nxd + col] = invJ21_avgEW[row * nxd + col]*dudxE[row * nxd + col] 
                                  + invJ22_avgEW[row * nxd + col]*dvdxE[row * nxd + col] 
                                  - invJ23_avgEW[row * nxd + col]*(invJ31_avgEW[row * nxd + col]*dudxE[row * nxd + col]+invJ32_avgEW[row * nxd + col]*dvdxE[row * nxd + col])/invJ33_avgEW[row * nxd + col];

        duxidetx[row * nxd + col] = invJ11_avgEW[row * nxd + col]*dudyE[row * nxd + col] 
                                  + invJ12_avgEW[row * nxd + col]*dvdyE[row * nxd + col] 
                                  - invJ13_avgEW[row * nxd + col]*(invJ31_avgEW[row * nxd + col]*dudyE[row * nxd + col]+invJ32_avgEW[row * nxd + col]*dvdyE[row * nxd + col])/invJ33_avgEW[row * nxd + col];
    
        dvetdetx[row * nxd + col] = invJ21_avgEW[row * nxd + col]*dudyE[row * nxd + col] 
                                  + invJ22_avgEW[row * nxd + col]*dvdyE[row * nxd + col] 
                                  - invJ23_avgEW[row * nxd + col]*(invJ31_avgEW[row * nxd + col]*dudyE[row * nxd + col]+invJ32_avgEW[row * nxd + col]*dvdyE[row * nxd + col])/invJ33_avgEW[row * nxd + col];
         
                                  
        duxidxiy[row * nxd + col] = invJ11_avgSN[row * nxd + col]*dudxN[row * nxd + col] 
                                  + invJ12_avgSN[row * nxd + col]*dvdxN[row * nxd + col] 
                                  - invJ13_avgSN[row * nxd + col]*(invJ31_avgSN[row * nxd + col]*dudxN[row * nxd + col]+invJ32_avgSN[row * nxd + col]*dvdxN[row * nxd + col])/invJ33_avgSN[row * nxd + col];
         
        dvetdxiy[row * nxd + col] = invJ21_avgSN[row * nxd + col]*dudxN[row * nxd + col] 
                                  + invJ22_avgSN[row * nxd + col]*dvdxN[row * nxd + col] 
                                  - invJ23_avgSN[row * nxd + col]*(invJ31_avgSN[row * nxd + col]*dudxN[row * nxd + col]+invJ32_avgSN[row * nxd + col]*dvdxN[row * nxd + col])/invJ33_avgSN[row * nxd + col];
        
        duxidety[row * nxd + col] = invJ11_avgSN[row * nxd + col]*dudyN[row * nxd + col] 
                                  + invJ12_avgSN[row * nxd + col]*dvdyN[row * nxd + col] 
                                  - invJ13_avgSN[row * nxd + col]*(invJ31_avgSN[row * nxd + col]*dudyN[row * nxd + col]+invJ32_avgSN[row * nxd + col]*dvdyN[row * nxd + col])/invJ33_avgSN[row * nxd + col];
                                 
        dvetdety[row * nxd + col] = invJ21_avgSN[row * nxd + col]*dudyN[row * nxd + col] 
                                  + invJ22_avgSN[row * nxd + col]*dvdyN[row * nxd + col] 
                                  - invJ23_avgSN[row * nxd + col]*(invJ31_avgSN[row * nxd + col]*dudyN[row * nxd + col]+invJ32_avgSN[row * nxd + col]*dvdyN[row * nxd + col])/invJ33_avgSN[row * nxd + col];
                                 

    }

}

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

    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-1) && col < (nxd-1) && col >0 && row >0 ) {
        
        PDx[0 * nxd * nyd + row * nxd + col] = 0.0;

        PDx[1 * nxd * nyd + row * nxd + col] = 0.0;

        PDx[2 * nxd * nyd + row * nxd + col] = 0.0;

        PDx[3 * nxd * nyd + row * nxd + col] = 0.0;

        PDx[4 * nxd * nyd + row * nxd + col] = 2.0*Detmin_avgEW[row * nxd + col]
                                             *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col]))) 
                                             *(invJ11_avgEW[row * nxd + col]*duxidxix[row * nxd + col]
                                             + invJ21_avgEW[row * nxd + col]*duxidetx[row * nxd + col]);
        
        PDx[5 * nxd * nyd + row * nxd + col] = Detmin_avgEW[row * nxd + col]
                                            *(0.5*((uE[0 * nxd * nyd + row * nxd + col]+uE[3 * nxd * nyd + row * nxd + col])+(uW[0 * nxd * nyd + row * nxd + col]+uW[3 * nxd * nyd + row * nxd + col])))  
                                            *(invJ12_avgEW[row * nxd + col]*duxidxix[row * nxd + col] + invJ22_avgEW[row * nxd + col]*duxidetx[row * nxd + col]
			                                + invJ11_avgEW[row * nxd + col]*dvetdxix[row * nxd + col] + invJ21_avgEW[row * nxd + col]*dvetdetx[row * nxd + col]);
	                 
        
        PDy[0 * nxd * nyd + row * nxd + col] = 0.0;

        PDy[1 * nxd * nyd + row * nxd + col] = 0.0;

        PDy[2 * nxd * nyd + row * nxd + col] = 0.0;

        PDy[3 * nxd * nyd + row * nxd + col] = 0.0;

        PDy[4 * nxd * nyd + row * nxd + col] = Detmin_avgSN[row * nxd + col]
                                            *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col]))) 
                                            *(invJ12_avgSN[row * nxd + col]*duxidxiy[row * nxd + col] + invJ22_avgSN[row * nxd + col]*duxidety[row * nxd + col]
                                            + invJ11_avgSN[row * nxd + col]*dvetdxiy[row * nxd + col] + invJ21_avgSN[row * nxd + col]*dvetdety[row * nxd + col]);

        PDy[5 * nxd * nyd + row * nxd + col] = 2.0* Detmin_avgSN[row * nxd + col]
                                            *(0.5*((uN[0 * nxd * nyd + row * nxd + col]+uN[3 * nxd * nyd + row * nxd + col])+(uS[0 * nxd * nyd + row * nxd + col]+uS[3 * nxd * nyd + row * nxd + col]))) 
                                            *(invJ12_avgSN[row * nxd + col]*dvetdxiy[row * nxd + col] 
                                            + invJ22_avgSN[row * nxd + col]*dvetdety[row * nxd + col]);


    }

}

__global__
void Flux11Kernel(double* result,
    double* u, 
    double* vex, double* vey, 
    double* vexF,double* veyF, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        if( u[0 * nxd * nyd + row * nxd + col] <= H00 )
        {  
            vex[row * nxd + col] = 0.0; 
            vey[row * nxd + col] = 0.0;
        }
        else
        { 
            vex[row * nxd + col] = u[1 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0); 
            vey[row * nxd + col] = u[2 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0);
        
        }

        if( u[3 * nxd * nyd + row * nxd + col] <= H00 )
		{  
            vexF[row * nxd + col] = 0.0; 
            veyF[row * nxd + col] = 0.0;
        }
	    else
        {  
            vexF[row * nxd + col] = u[4 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
            veyF[row * nxd + col] = u[5 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
        }

    }

}

__global__
void Flux12Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* vexF, double* veyF, 
    double* svec, double* cvalue,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ){
        
        
        w_wert[row * nxd + col] = (svec[0 * nxd * nyd + row * nxd + col]*vex[row * nxd + col] 
                                  +svec[1 * nxd * nyd + row * nxd + col]*vey[row * nxd + col] )/cvalue[row * nxd + col];

        w_wertF[row * nxd + col] = (svec[0 * nxd * nyd + row * nxd + col]*vexF[row * nxd + col] 
                                   +svec[1 * nxd * nyd + row * nxd + col]*veyF[row * nxd + col] )/cvalue[row * nxd + col];

    }

}


__global__
void Flux13Kernel(double* result,
    double* w_wert, double* w_wertF, 
    double* vex, double* vey, 
    double* usw, double* vel,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        
        usw[row * nxd + col] =  vex[row * nxd + col]*vex[row * nxd + col]
                               +vey[row * nxd + col]*vey[row * nxd + col]
                            +w_wert[row * nxd + col]*w_wert[row * nxd + col];

        vel[row * nxd + col] =  sqrt(vex[row * nxd + col]*vex[row * nxd + col]
                                    +vey[row * nxd + col]*vey[row * nxd + col]
                                 +w_wert[row * nxd + col]*w_wert[row * nxd + col]);

    }

}

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
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;


    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        if(usw[row * nxd + col] <= U0*U0)
		{ 
		  
          vexw[row * nxd + col] = vex[row * nxd + col]/U0; 
          veyw[row * nxd + col] = vey[row * nxd + col]/U0;
		 
		  vel[row * nxd + col] = U0;
		}
	      else
		{
          vexw[row * nxd + col] = vex[row * nxd + col]/vel[row * nxd + col]; 
          veyw[row * nxd + col] = vey[row * nxd + col]/vel[row * nxd + col];
		  
        }
        
        q_xi[row * nxd + col] = invJ11[row * nxd + col]*vex[row * nxd + col]
                               +invJ12[row * nxd + col]*vey[row * nxd + col]
                               +invJ13[row * nxd + col]*w_wert[row * nxd + col];

        q_et[row * nxd + col] = invJ21[row * nxd + col]*vex[row * nxd + col]
                               +invJ22[row * nxd + col]*vey[row * nxd + col]
                               +invJ23[row * nxd + col]*w_wert[row * nxd + col];

        q_xiF[row * nxd + col] = invJ11[row * nxd + col]*vexF[row * nxd + col]
                                +invJ12[row * nxd + col]*veyF[row * nxd + col]
                                +invJ13[row * nxd + col]*w_wertF[row * nxd + col];

        q_etF[row * nxd + col] = invJ21[row * nxd + col]*vexF[row * nxd + col]
                                +invJ22[row * nxd + col]*veyF[row * nxd + col]
                                +invJ23[row * nxd + col]*w_wertF[row * nxd + col];
        

    }


}

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
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
                
        Ac[row * nxd + col] = vex[row * nxd + col]*J13dxi[row * nxd + col]*q_xi[row * nxd + col]  
                            + vey[row * nxd + col]*J23dxi[row * nxd + col]*q_xi[row * nxd + col] 
                         + w_wert[row * nxd + col]*J33dxi[row * nxd + col]*q_xi[row * nxd + col]
                            + vex[row * nxd + col]*J13det[row * nxd + col]*q_et[row * nxd + col] 
                            + vey[row * nxd + col]*J23det[row * nxd + col]*q_et[row * nxd + col]  
                         + w_wert[row * nxd + col]*J33det[row * nxd + col]*q_et[row * nxd + col];

        AcF[row * nxd + col] = vexF[row * nxd + col]*J13dxi[row * nxd + col]*q_xiF[row * nxd + col]
                             + veyF[row * nxd + col]*J23dxi[row * nxd + col]*q_xiF[row * nxd + col] 
                          + w_wertF[row * nxd + col]*J33dxi[row * nxd + col]*q_xiF[row * nxd + col]
                             + vexF[row * nxd + col]*J13det[row * nxd + col]*q_etF[row * nxd + col] 
                             + veyF[row * nxd + col]*J23det[row * nxd + col]*q_etF[row * nxd + col] 
                          + w_wertF[row * nxd + col]*J33det[row * nxd + col]*q_etF[row * nxd + col];

    }

}

__global__
void Flux16Kernel(double* result,
    double* Npress1, double* Npress2, double* NpressF,
    double* Ac, double* AcF, 
    double* Detmin, double* cvalue,
    double* u,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    double Drtio1 = 1.0 - (1.420/2.60);
    double Dratio   = 1.420/2.60;

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        
        Npress1[row * nxd + col] = max(0.0, 
                                       (Detmin[row * nxd + col]*u[0 * nxd * nyd + row * nxd + col]*(Drtio1*cvalue[row * nxd + col] - Ac[row * nxd + col] + AcF[row * nxd + col]*Dratio)));
        
        Npress2[row * nxd + col] = max(0.0, 
                                       (Detmin[row * nxd + col]*u[0 * nxd * nyd + row * nxd + col]*(cvalue[row * nxd + col] - Ac[row * nxd + col])));
        
        NpressF[row * nxd + col] = max(0.0, 
                                       (Detmin[row * nxd + col]*u[3 * nxd * nyd + row * nxd + col]*(cvalue[row * nxd + col] - AcF[row * nxd + col])));


    }

}

__global__
void Flux17Kernel(double* result,
    double* Npress1, double* Npress2, double* NpressF,
    double* Ac, double* AcF, 
    double* Detmin, double* svec,
    double* vex, double* vey,
    double* vexw, double* veyw,
    double* vexF, double* veyF,
    double* tande,
    double* u, double* s,
    double Cd, double N_R, double varTheta,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    double Dratio   = 1.420/2.60;
    double epsilon0 = 1.0;

    if (row<(nyd-1) && col < (nxd-1) && col >=0 && row >=0 ) {
        
        
        s[0 * nxd * nyd + row * nxd + col] = 0.0;

        s[1 * nxd * nyd + row * nxd + col] = - Npress2[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col]
                                             - Npress1[row * nxd + col]*tande[row * nxd + col]*vexw[row * nxd + col] 
                                             + Dratio*Detmin[row * nxd + col]*Cd*u[6 * nxd * nyd + row * nxd + col]*(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *(vexF[row * nxd + col] - vex[row * nxd + col]);

        s[2 * nxd * nyd + row * nxd + col] = - Npress2[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col]
                                             - Npress1[row * nxd + col]*tande[row * nxd + col]*veyw[row * nxd + col]
                                             + Dratio*Detmin[row * nxd + col]*Cd*u[6 * nxd * nyd + row * nxd + col]*(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *(veyF[row * nxd + col] - vey[row * nxd + col]);


        s[3 * nxd * nyd + row * nxd + col] = 0.0;

        s[4 * nxd * nyd + row * nxd + col]=  -NpressF[row * nxd + col]*svec[0 * nxd * nyd + row * nxd + col] 
                                             - Detmin[row * nxd + col]*Cd*u[6 * nxd * nyd + row * nxd + col]*(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *(vexF[row * nxd + col] - vex[row * nxd + col]) 
                                             + Detmin[row * nxd + col]
                                             *(-(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *varTheta *vexF[row * nxd + col]/epsilon0)/N_R;
        
        s[5 * nxd * nyd + row * nxd + col]=  -NpressF[row * nxd + col]*svec[1 * nxd * nyd + row * nxd + col] 
                                             - Detmin[row * nxd + col]*Cd*u[6 * nxd * nyd + row * nxd + col]*(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *(veyF[row * nxd + col] - vey[row * nxd + col]) 
                                             + Detmin[row * nxd + col]
                                             *(-(1.0 - u[6 * nxd * nyd + row * nxd + col])
                                             *(u[0 * nxd * nyd + row * nxd + col] + u[3 * nxd * nyd + row * nxd + col])
                                             *varTheta *veyF[row * nxd + col]/epsilon0)/N_R; 
        // result[row * nxd + col] = Npress2[row * nxd + col];


    }

}

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
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    double epsilon0 = 1.0;

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        for (int i=0;i<(MN-1);i++){
            v[i * nxd * nyd + row * nxd + col] = uzero[i * nxd * nyd + row * nxd + col]*Detmin[row * nxd + col]
                                                - (*dt/dx)*( Hpx[i * nxd * nyd + row * nxd + col] -  Hpx[i * nxd * nyd + (row  ) * nxd + (col-1)] )
                                                - (*dt/dy)*( Hpy[i * nxd * nyd + row * nxd + col] -  Hpy[i * nxd * nyd + (row-1) * nxd + (col  )] )
                                                + (*dt/dx)*( Ppx[i * nxd * nyd + row * nxd + col] -  Ppx[i * nxd * nyd + (row  ) * nxd + (col-1)] )*(epsilon0*u[6 * nxd * nyd + row * nxd + col])
                                                + (*dt/dy)*( Ppy[i * nxd * nyd + row * nxd + col] -  Ppy[i * nxd * nyd + (row-1) * nxd + (col  )] )*(epsilon0*u[6 * nxd * nyd + row * nxd + col])
                                                + (*dt/dx)*( PDx[i * nxd * nyd + row * nxd + col] -  PDx[i * nxd * nyd + (row  ) * nxd + (col-1)] )*(epsilon0*(1.0-u[6 * nxd * nyd + row * nxd + col])/N_R)
                                                + (*dt/dy)*( PDy[i * nxd * nyd + row * nxd + col] -  PDy[i * nxd * nyd + (row-1) * nxd + (col  )] )*(epsilon0*(1.0-u[6 * nxd * nyd + row * nxd + col])/N_R)
                                                +    (*dt)*s[i * nxd * nyd + row * nxd + col];

        }
       
        // result[row * nxd + col] = s[5 * nxd * nyd + row * nxd + col];


    }

}

__global__
void Flux19Kernel(double* result,
    double* v, double* Detmin, 
    double* u, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (col < (nxd-MD) && row<(nyd-MD) && col >2 && row >2 ) {

        for (int i=0;i<(MN-1);i++){
            u[i * nxd * nyd + row * nxd + col] = v[i * nxd * nyd + row * nxd + col]/Detmin[row * nxd + col];
        // result[row * nxd + col] =  v[0 * nxd * nyd + row * nxd + col];

        }
    }

}



__global__
void Flux20Kernel(double* result,
    double* u, double* utmp,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        if (u[0 * nxd * nyd + row * nxd + col] <= H00)
        {
            u[0 * nxd * nyd + row * nxd + col] = 0.0;
            u[1 * nxd * nyd + row * nxd + col] = 0.0;
            u[2 * nxd * nyd + row * nxd + col] = 0.0;
        }

        if (u[3 * nxd * nyd + row * nxd + col] <= H00)
		{		
            u[3 * nxd * nyd + row * nxd + col] = 0.0;
            u[4 * nxd * nyd + row * nxd + col] = 0.0;
            u[5 * nxd * nyd + row * nxd + col] = 0.0;
		}


    }

}

__global__
void Flux21Kernel(double* result,
    double* u, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        u[6 * nxd * nyd + row * nxd + col] = u[0 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+u[3 * nxd * nyd + row * nxd + col]+H0);

    }

}



__global__
void Flux22Kernel(double* result,
    double* u, double* uone,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        for (int i=0;i<(MN);i++){

            uone[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + row * nxd + col];
        
        }

    }

}

__global__
void Flux23Kernel(double* result,
    double* u, double* uzero,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {
      
        usxnew[row * nxd + col] = u[1 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0);
        
        ufxnew[row * nxd + col] = u[4 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
        
        usxold[row * nxd + col] = uzero[1 * nxd * nyd + row * nxd + col]/(uzero[0 * nxd * nyd + row * nxd + col]+H0);

		ufxold[row * nxd + col] = uzero[4 * nxd * nyd + row * nxd + col]/(uzero[3 * nxd * nyd + row * nxd + col]+H0);
        
        usynew[row * nxd + col] = u[2 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0);
        
        ufynew[row * nxd + col] = u[5 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
        
        usyold[row * nxd + col] = uzero[2 * nxd * nyd + row * nxd + col]/(uzero[0 * nxd * nyd + row * nxd + col]+H0);

        ufyold[row * nxd + col] = uzero[5 * nxd * nyd + row * nxd + col]/(uzero[3 * nxd * nyd + row * nxd + col]+H0);

    }

}

__global__
void Flux24Kernel(double* result,
    double* u, double* uone,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (col < (nxd-MD) && row<(nyd-MD) && col >2 && row >2 ) {

        if ((usxnew[row * nxd + col] - ufxnew[row * nxd + col])*(usxold[row * nxd + col] - ufxold[row * nxd + col]) >= 0.0){
            uone[4 * nxd * nyd + row * nxd + col] = u[4 * nxd * nyd + row * nxd + col];
        }
        else{
            uone[4 * nxd * nyd + row * nxd + col] = usxnew[row * nxd + col] * u[3 * nxd * nyd + row * nxd + col];
        }

        if ((usynew[row * nxd + col] - ufynew[row * nxd + col])*(usyold[row * nxd + col] - ufyold[row * nxd + col]) >= 0.0){
            uone[5 * nxd * nyd + row * nxd + col] = u[5 * nxd * nyd + row * nxd + col];
        }
        else{
            uone[5 * nxd * nyd + row * nxd + col] = usynew[row * nxd + col] * u[3 * nxd * nyd + row * nxd + col];
        }
        
        uone[6 * nxd * nyd + row * nxd + col] = u[6 * nxd * nyd + row * nxd + col];


    }

}



__global__
void Flux25Kernel(double* result,
    double* u, double* utwo,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    
    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        for (int i=0;i<(MN);i++){
            utwo[i * nxd * nyd + row * nxd + col] = u[i * nxd * nyd + row * nxd + col];
        
        }


    }

}

__global__
void Flux26Kernel(double* result,
    double* u, double* uone,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        usxnew[row * nxd + col] = u[1 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0);
        
        ufxnew[row * nxd + col] = u[4 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
        
        usxold[row * nxd + col] = uone[1 * nxd * nyd + row * nxd + col]/(uone[0 * nxd * nyd + row * nxd + col]+H0);

		ufxold[row * nxd + col] = uone[4 * nxd * nyd + row * nxd + col]/(uone[3 * nxd * nyd + row * nxd + col]+H0);
        
        usynew[row * nxd + col] = u[2 * nxd * nyd + row * nxd + col]/(u[0 * nxd * nyd + row * nxd + col]+H0);
        
        ufynew[row * nxd + col] = u[5 * nxd * nyd + row * nxd + col]/(u[3 * nxd * nyd + row * nxd + col]+H0);
        
        usyold[row * nxd + col] = uone[2 * nxd * nyd + row * nxd + col]/(uone[0 * nxd * nyd + row * nxd + col]+H0);

        ufyold[row * nxd + col] = uone[5 * nxd * nyd + row * nxd + col]/(uone[3 * nxd * nyd + row * nxd + col]+H0);


    }

}

__global__
void Flux27Kernel(double* result,
    double* u, double* utwo,
    double* usxnew, double* ufxnew, 
    double* usxold, double* ufxold, 
    double* usynew, double* ufynew, 
    double* usyold, double* ufyold, 
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        if ((usxnew[row * nxd + col] - ufxnew[row * nxd + col])*(usxold[row * nxd + col] - ufxold[row * nxd + col]) >= 0.0){
            utwo[4 * nxd * nyd + row * nxd + col] = u[4 * nxd * nyd + row * nxd + col];
        }
        else{
            utwo[4 * nxd * nyd + row * nxd + col] = usxnew[row * nxd + col] * u[3 * nxd * nyd + row * nxd + col];
        }

        if ((usynew[row * nxd + col] - ufynew[row * nxd + col])*(usyold[row * nxd + col] - ufyold[row * nxd + col]) >= 0.0){
            utwo[5 * nxd * nyd + row * nxd + col] = u[5 * nxd * nyd + row * nxd + col];
        }
        else{
            utwo[5 * nxd * nyd + row * nxd + col] = usynew[row * nxd + col] * u[3 * nxd * nyd + row * nxd + col];
        }
        
        utwo[6 * nxd * nyd + row * nxd + col] = u[6 * nxd * nyd + row * nxd + col];


    }

}

__global__
void Flux28Kernel(double* result,
    double* u, double* utwo,
    int nxd, int nyd)
{
    int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
    int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    

    if (row<(nyd-MD) && col < (nxd-MD) && col >2 && row >2 ) {

        for (int i=0;i<(MN);i++){

            u[i * nxd * nyd + row * nxd + col] = utwo[i * nxd * nyd + row * nxd + col];
        
        }

    }

}
