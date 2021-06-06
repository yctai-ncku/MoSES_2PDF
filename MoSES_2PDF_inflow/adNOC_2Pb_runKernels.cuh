#ifndef __SWE_RUNKERNELS_CUH
#define __SWE_RUNKERNELS_CUH

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream> /* for ofstream */ 
#include <time.h>
#include <vector>
#include <cstring>
#include <string.h>

#include "matalloc.h"
#include "adNOC_2Pb_kernels.cuh"

#include "cuda_runtime.h"

#define PI 3.14159
#define MD 3
#define MN 7


using namespace std;

class RunKernels {
public:
    double *TotalStep_h, *dt_h;

    RunKernels(char *fname);
    clock_t run();

private:
    int Ttype, NX, NY, nx, ny, nxd, nyd, nxyd, arraySize, inflowSize;//, Total_step;
	char InitData[200], DEMData[200], locData[200], IniHData[200], IniUData[200], IniVData[200];
	double dx, dy;
	double MINX, MAXX, MINY, MAXY, delta0, Cd, N_R, varTheta, phiS0;
	double TotalSim, EachOut;
	double xllcorner, yllcorner;

	vector<string> Topodata, Initdata, inputTopoTmp, inputInitTmp;
	vector<string> IniUdata, IniVdata, IniHdata, locdata, inputIniUTmp, inputIniVTmp, inputIniHTmp, inputlocTmp; 

	string TopoS, InitS;
	string IniUS, IniVS, IniHS, locS;
	int StartTopo, StartInit;

	int locflowlen, Iniflowlen;

	double afa_rho = 0.9;

	double **topo, **depth, ***inputFlow, **inputLoc, *inflowTime;

	double *depth_h, *topo_h, *speed_h;
	double *resultHs_h, *resultHf_h, *resultUs_h, *resultVs_h,*resultUf_h,*resultVf_h, *resultphi_h;

	double *bfkt_h, *cvalue_h, *svec_h;

	double *result_h;

	double *inflow_h, *loc_h;

	int **dire;
	int *dire_h; 

	int* dev_dire =  nullptr;

	double* dev_inflow =  nullptr;
	double* dev_loc    =  nullptr;

	double* dev_topo   =  nullptr;
	double* dev_depth  = nullptr;
	double* dev_result = nullptr;

	double* dev_bfkt   = nullptr;
	double* dev_posx   = nullptr;
	double* dev_posy   = nullptr;

	double* dev_dxdxi11 = nullptr;
	double* dev_dxdxi12 = nullptr;
	double* dev_dxdxi21 = nullptr;
	double* dev_dxdxi22 = nullptr;

	double* dev_dbdx    = nullptr;
	double* dev_dbdy    = nullptr;
	double* dev_cvalue  = nullptr;

	double* dev_svec    = nullptr;
	double* dev_Jacb31  = nullptr;
	double* dev_Jacb32  = nullptr;

	double* dev_dettmp  = nullptr;
	double* dev_Detmin  = nullptr;

	double* dev_i_ddxi11 = nullptr;
	double* dev_i_ddxi12 = nullptr;
	double* dev_i_ddxi21 = nullptr;
	double* dev_i_ddxi22 = nullptr;

	double* dev_invJ11 = nullptr;
	double* dev_invJ12 = nullptr;
	double* dev_invJ13 = nullptr;
	double* dev_invJ21 = nullptr;
	double* dev_invJ22 = nullptr;
	double* dev_invJ23 = nullptr;
	double* dev_invJ31 = nullptr;
	double* dev_invJ32 = nullptr;
	double* dev_invJ33 = nullptr;

	double* dev_u      = nullptr;
	double* dev_uzero  = nullptr;

	double* dev_tande  = nullptr;

	double* dev_J13dxi  = nullptr;
	double* dev_J23dxi  = nullptr;
	double* dev_J33dxi  = nullptr;
	double* dev_J13det  = nullptr;
	double* dev_J23det  = nullptr;
	double* dev_J33det  = nullptr;

	double* dev_ux      = nullptr;
	double* dev_uy      = nullptr;

	double* dev_apEW    = nullptr;
	double* dev_apSN    = nullptr;
	double* dev_apFEW   = nullptr;
	double* dev_apFSN   = nullptr;

	double* dev_dux     = nullptr;
	double* dev_duy     = nullptr;

	double* dev_t1x     = nullptr;
	double* dev_t2x     = nullptr;
	double* dev_t1y     = nullptr;
	double* dev_t2y     = nullptr;

	double* dev_sgnAx   = nullptr;
	double* dev_sgnBx   = nullptr;
	double* dev_sgnAy   = nullptr;
	double* dev_sgnBy   = nullptr;

	double* dev_sgnAx2  = nullptr;
	double* dev_sgnBx2  = nullptr;
	double* dev_sgnAy2  = nullptr;
	double* dev_sgnBy2  = nullptr;

	double* dev_dxdxi11_avgEW   = nullptr;
	double* dev_dxdxi21_avgEW   = nullptr;

	double* dev_dxdxi12_avgSN   = nullptr;
	double* dev_dxdxi22_avgSN   = nullptr;

	double* dev_J13dxi_avgEW   = nullptr;
	double* dev_J23dxi_avgEW   = nullptr;
	double* dev_J33dxi_avgEW   = nullptr;

	double* dev_J13det_avgEW   = nullptr;
	double* dev_J23det_avgEW   = nullptr;
	double* dev_J33det_avgEW   = nullptr;

	double* dev_J13dxi_avgSN   = nullptr;
	double* dev_J23dxi_avgSN   = nullptr;
	double* dev_J33dxi_avgSN   = nullptr;

	double* dev_J13det_avgSN   = nullptr;
	double* dev_J23det_avgSN   = nullptr;
	double* dev_J33det_avgSN   = nullptr;

	double* dev_invJ11_avgEW   = nullptr;
	double* dev_invJ12_avgEW   = nullptr;
	double* dev_invJ13_avgEW   = nullptr;

	double* dev_invJ21_avgEW   = nullptr;
	double* dev_invJ22_avgEW   = nullptr;
	double* dev_invJ23_avgEW   = nullptr;

	double* dev_invJ31_avgEW   = nullptr;
	double* dev_invJ32_avgEW   = nullptr;
	double* dev_invJ33_avgEW   = nullptr;

	double* dev_invJ11_avgSN   = nullptr;
	double* dev_invJ12_avgSN   = nullptr;
	double* dev_invJ13_avgSN   = nullptr;

	double* dev_invJ21_avgSN   = nullptr;
	double* dev_invJ22_avgSN   = nullptr;
	double* dev_invJ23_avgSN   = nullptr;

	double* dev_invJ31_avgSN   = nullptr;
	double* dev_invJ32_avgSN   = nullptr;
	double* dev_invJ33_avgSN   = nullptr;

	double* dev_Detmin_avgEW   = nullptr;
	double* dev_Detmin_avgSN   = nullptr;
	
	double* dev_cval_avgEW   = nullptr;
	double* dev_cval_avgSN   = nullptr;
	
	double* dev_svec_avgEW   = nullptr;
	double* dev_svec_avgSN   = nullptr;

	double* dev_uE = nullptr;
	double* dev_uW = nullptr;
	double* dev_uN = nullptr;
	double* dev_uS = nullptr;

	double* dev_vexE = nullptr;
	double* dev_vexW = nullptr;
	double* dev_veyE = nullptr;
	double* dev_veyW = nullptr;

	double* dev_w_wertE = nullptr;
	double* dev_w_wertW = nullptr;

	double* dev_vexFE = nullptr;
	double* dev_vexFW = nullptr;
	double* dev_veyFE = nullptr;
	double* dev_veyFW = nullptr;

	double* dev_w_wertFE = nullptr;
	double* dev_w_wertFW = nullptr;

	double* dev_vexN = nullptr;
	double* dev_vexS = nullptr;
	double* dev_veyN = nullptr;
	double* dev_veyS = nullptr;

	double* dev_w_wertN = nullptr;
	double* dev_w_wertS = nullptr;

	double* dev_vexFN = nullptr;
	double* dev_vexFS = nullptr;
	double* dev_veyFN = nullptr;
	double* dev_veyFS = nullptr;

	double* dev_w_wertFN = nullptr;
	double* dev_w_wertFS = nullptr;

	double* dev_q_xiE = nullptr;
	double* dev_q_etE = nullptr;
	double* dev_q_xiW = nullptr;
	double* dev_q_etW = nullptr;
	
	double* dev_q_xiFE = nullptr;
	double* dev_q_etFE = nullptr;
	double* dev_q_xiFW = nullptr;
	double* dev_q_etFW = nullptr;

	double* dev_NpressFE = nullptr;
	double* dev_NpressFW = nullptr;
	
	double* dev_M11EW    = nullptr;

	double* dev_q_xiN = nullptr;
	double* dev_q_etN = nullptr;
	double* dev_q_xiS = nullptr;
	double* dev_q_etS = nullptr;
	
	double* dev_q_xiFN = nullptr;
	double* dev_q_etFN = nullptr;
	double* dev_q_xiFS = nullptr;
	double* dev_q_etFS = nullptr;

	double* dev_NpressFN = nullptr;
	double* dev_NpressFS = nullptr;
	
	double* dev_M22SN    = nullptr;

	double* dev_apE = nullptr;
	double* dev_apW = nullptr;
	double* dev_apFE = nullptr;
	double* dev_apFW = nullptr;

	double* dev_apN = nullptr;
	double* dev_apS = nullptr;
	double* dev_apFN = nullptr;
	double* dev_apFS = nullptr;

	double* dev_em_x  = nullptr;
	double* dev_em_y  = nullptr;
	double* dev_em_Fx = nullptr;
	double* dev_em_Fy = nullptr;

	double* dev_FpE   = nullptr;
	double* dev_FpW   = nullptr;
	double* dev_GpN   = nullptr;
	double* dev_GpS   = nullptr;

	double* dev_czw1x  = nullptr;
	double* dev_czw2x  = nullptr;
	double* dev_czwF1x  = nullptr;
	double* dev_czwF2x  = nullptr;
	
	double* dev_czw1y  = nullptr;
	double* dev_czw2y  = nullptr;
	double* dev_czwF1y  = nullptr;
	double* dev_czwF2y  = nullptr;

	double* dev_em_valS = nullptr;
	double* dev_em_valF = nullptr;
	double* dev_Val     = nullptr;

	double* dev_Hpx     =  nullptr;
	double* dev_Hpy     =  nullptr;
	double* dev_Ppx     =  nullptr;
	double* dev_Ppy     =  nullptr;
	double* dev_PDx     =  nullptr;
	double* dev_PDy     =  nullptr;

	double* dev_dudxE = nullptr;
	double* dev_dvdxE = nullptr;
	double* dev_dudyE = nullptr;
	double* dev_dvdyE = nullptr;

	double* dev_dudxN = nullptr;
	double* dev_dvdxN = nullptr;
	double* dev_dudyN = nullptr;
	double* dev_dvdyN = nullptr;

	double* dev_duxidxix = nullptr;
	double* dev_dvetdxix = nullptr;
	double* dev_duxidetx = nullptr;
	double* dev_dvetdetx = nullptr;

	double* dev_duxidxiy = nullptr;
	double* dev_dvetdxiy = nullptr;
	double* dev_duxidety = nullptr;
	double* dev_dvetdety = nullptr;

	double* dev_vex      = nullptr;
	double* dev_vey      = nullptr;
	double* dev_vexF     = nullptr;
	double* dev_veyF     = nullptr;
	

	double* dev_w_wert   = nullptr;
	double* dev_w_wertF  = nullptr;
	double* dev_usw      = nullptr;
	double* dev_vel      = nullptr;

	double* dev_vexw     = nullptr;
	double* dev_veyw     = nullptr;
	double* dev_q_xi     = nullptr;
	double* dev_q_et     = nullptr;
	double* dev_q_xiF    = nullptr;
	double* dev_q_etF    = nullptr;

	double* dev_Ac       = nullptr;
	double* dev_AcF      = nullptr;

	double* dev_Npress1  = nullptr;
	double* dev_Npress2  = nullptr;
	double* dev_NpressF  = nullptr;

	double* dev_s = nullptr;
	double* dev_v = nullptr;

	double* dev_uone   = nullptr;
	double* dev_utwo   = nullptr;
	
	double* dev_usxnew = nullptr;
	double* dev_ufxnew = nullptr;
	double* dev_usxold = nullptr;
	double* dev_ufxold = nullptr;

	double* dev_usynew = nullptr;
	double* dev_ufynew = nullptr;
	double* dev_usyold = nullptr;
	double* dev_ufyold = nullptr;

	double* dev_waveSpeed = nullptr;
	double* dev_TotalTime = nullptr;
	double* dev_max  = nullptr;
	double* dev_maxW = nullptr;
	double* dev_dt   = nullptr;

	double* dev_dtval = nullptr;
	double* dev_utmp  = nullptr;


    void readinParamFile(char *fname);
    void readParamLine(FILE *fid, char *line, int len);
    void dataInitialization();

    void memoryMalloc();
    void kernelStep();
    void freeMemory();
	void outputFile();
	void split(const string& s, vector<string>& sv, const char* delim);
};

#endif