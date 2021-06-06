#include "adNOC_2Pb_runKernels.cuh"
#include <assert.h>

RunKernels::RunKernels(char *fname)
{
    readinParamFile(fname);
    
    dataInitialization();
}

void RunKernels::readinParamFile(char *fname)
{
	FILE* FID =  fopen(fname, "r");
	if(!FID) {
		printf("Parameter file %s is missing. Abort...\n", fname);
		exit(1);
	}

	char line[1024];

	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &TotalSim);
  
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &EachOut);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &delta0);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &Cd);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &N_R);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &varTheta);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%lf", &phiS0);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%s", DEMData);
	
	readParamLine(FID, line, 1024);
	sscanf(line, "%s", InitData);
	
}

void RunKernels::readParamLine(FILE *fid, char *line, int len)
{
	char *c;
	fgets(line, len, fid);

	// Remove the things after '#', or end of line ('\n')
	for(c = line; *c; c++) {
		int br = 0;
		switch(*c) {
		case '#':
		case '\n':
		*c = 0;
		br = 1;
		break;
		default:
		break;
		}
		if(br) break;
	}

	return;
}

void RunKernels::split(const string& s, vector<string>& sv, const char* delim)
{
	sv.clear();
	char* buffer = new char[s.size() + 1];
	buffer[s.size()] = '\0';
	copy(s.begin(), s.end(), buffer);

	char* p = std::strtok(buffer, delim);
	do {
		sv.push_back(p);
	} while ((p = strtok(NULL, delim)));

	delete[] buffer;

	return;
}

void RunKernels::dataInitialization()
{

	fprintf(stdout,"\n\t***********************************************************************************\n");
	fprintf(stdout,"\t\t* * * * *     MoSES_2PDF:  landslide (debris flow) mode     * * * * * \n");
	fprintf(stdout,"\t-----------------------------------------------------------------------------------\n");
  	fprintf(stdout,"\t\t2D Central Scheme  (adNOC) Mixture Code in CUDA: landslide mode\n");
	fprintf(stdout,"\t\tAuthor : Chi-Jyun Ko, Po-Chih Chen, Hock-Kiet Wong and Yih-Chin Tai\n");
	fprintf(stdout,"\t\tLab for Computer Simulation and Visualization (CSV Lab), NCKU, Taiwan\n");  
  	fprintf(stdout,"\t***********************************************************************************\n\n");
	

	ifstream inputFileTopo;
	ifstream inputFileInit;
  
	inputFileTopo.open(DEMData,  ios::in);
	inputFileInit.open(InitData, ios::in);

	fprintf(stdout,"\t\tinput Topo File    : %s\n", DEMData);
	fprintf(stdout,"\t\tinput Initial File : %s\n\n", InitData);

	// read Topo data
	if (inputFileTopo.fail() ){
	  cout << "\nError can't open Topo file.\n"<< endl;
	  assert(0);
	}

	// read init data
	if (inputFileInit.fail() ){
		cout << "\nError can't open Initial file.\n"<< endl;
		assert(0);
	}

	
	while(getline(inputFileTopo,TopoS)){
		inputTopoTmp.push_back(TopoS);
	}

	while(getline(inputFileInit,InitS)){
		inputInitTmp.push_back(InitS);
	}

	// read Topo file
	for(int i=0;i<6;i++){
		split(inputTopoTmp[i], Topodata, " ");
		if(i==0){
			NX = stoi(Topodata[1]);
		}
		else if(i==1){
			NY = stoi(Topodata[1]);
		}
		else if(i==2){
			xllcorner = stof(Topodata[1]);
		}else if(i==3){
			yllcorner = stof(Topodata[1]);
		}
		else if(i==4){
			dx = stof(Topodata[1]);
			dy = stof(Topodata[1]);
		}
		else if(i==5){
			if(Topodata[0] == "NODATA_value"){
				StartTopo = 0;
			}
			else{
				StartTopo = 1;
			}
		}
	}

	// read Init File
	for(int i=0;i<6;i++){
		split(inputInitTmp[i], Initdata, " ");
		if(i==5){
			if(Initdata[0] == "NODATA_value"){
				StartInit = 0;
			}
			else{
				StartInit = 1;
			}
		}
	}


	nx = NX;
	ny = NY;

	if(nx*ny>100000){
		fprintf(stdout,"\n\t***********************************\n");
  		fprintf(stdout,"\t     Number of grids exceeded\n");
  		fprintf(stdout,"\t***********************************\n\n");

		assert(0);
	}

	dx = dx*10;
	dy = dy*10;

	// dx = (MAXX - MINX)/(nx-1);
	// dy = (MAXY - MINY)/(ny-1);
	MINX = 0.0;
	MINY = 0.0;

	MAXX = dx*(nx-1);
	MAXY = dy*(ny-1);

	nxd = nx + 2*MD;
	nyd = ny + 2*MD;

	nxyd = max(nxd,nyd);

	arraySize = nxd * nyd;

	fprintf(stdout,"\t\tData points          : %d,%d\n", nx, ny);
  	fprintf(stdout,"\t\tDomain [dm]          : (%5.3f,%5.3f)(%5.3f,%5.3f)\n",
		MINX, MAXX, MINY, MAXY);
  	fprintf(stdout,"\t\tGrid size            : %6.2f,%6.2f (%d,%d)\n",dx,dy,nx,ny);
  	/*  fprintf(stdout,"\tMaximum time         : %5.3f\n", tf);*/
	fprintf(stdout,"\t\tCFL number           : %5.3f\n", CFL);
	fprintf(stdout,"\t\tdelta0               : %5.3f\n", delta0);
	fprintf(stdout,"\t\tCd                   : %5.3f\n", Cd);
	fprintf(stdout,"\t\tN_R                  : %5.3f\n", N_R);
	fprintf(stdout,"\t\tvarTheta             : %5.3f\n", varTheta);
	fprintf(stdout,"\t\tinitial value of solid volume fraction : %5.3f\n\n", phiS0);

	fprintf(stdout,"\t\tTotal simulation time (sec) : %5.3f\n", TotalSim);
	fprintf(stdout,"\t\tEach output time (sec)      : %5.3f\n\n", EachOut);
	

	NEW_MATRIX(topo, double, nxd, nyd);
	NEW_MATRIX(depth, double, nxd, nyd);

	// input Topo to matrix
	if(StartTopo==1){
		for(int j=0; j<NY; j++) {
			split(inputTopoTmp[j+5], Topodata, " ");
			
			for(int i=0; i<NX; i++) {
				topo[i+MD][j+MD] = stod(Topodata[i]);
				topo[i+MD][j+MD] = 10.0*topo[i+MD][j+MD];

				if(topo[i+MD][j+MD]<0){
					topo[i+MD][j+MD] = 0;
				}
			}
	
		}
	}
	else{
		for(int j=0; j<NY; j++) {
			split(inputTopoTmp[j+6], Topodata, " ");
			
			for(int i=0; i<NX; i++) {
				topo[i+MD][j+MD] = stod(Topodata[i]);
				topo[i+MD][j+MD] = 10.0*topo[i+MD][j+MD];

				if(topo[i+MD][j+MD]<0){
					topo[i+MD][j+MD] = 0;
				}
			}
	
		}
	}
	
	// input Init to matrix
	if(StartInit==1){
		for(int j=0; j<NY; j++) {
			split(inputInitTmp[j+5], Initdata, " ");
			
			for(int i=0; i<NX; i++) {
				depth[i+MD][j+MD] = stod(Initdata[i]);
				depth[i+MD][j+MD] = 10.0*depth[i+MD][j+MD];

				if(depth[i+MD][j+MD]<0){
					depth[i+MD][j+MD] = 0;
				}
			}
	
		}
	
	}
	else{
		for(int j=0; j<NY; j++) {
			split(inputInitTmp[j+6], Initdata, " ");
			
			for(int i=0; i<NX; i++) {
				depth[i+MD][j+MD] = stod(Initdata[i]);
				depth[i+MD][j+MD] = 10.0*depth[i+MD][j+MD];

				if(depth[i+MD][j+MD]<0){
					depth[i+MD][j+MD] = 0;
				}
			}
	
		}
	}


	// B.C. initioal condition
	for(int i=0; i<MD; i++) {
		for(int j=MD; j<(ny+MD); j++) {
			topo[i      ][j] = topo[MD     ][j];
			topo[nx+MD+i][j] = topo[nx+MD-1][j];

			depth[i      ][j] = depth[MD     ][j];
			depth[nx+MD+i][j] = depth[nx+MD-1][j];
		
		}
	}

	for(int j=0; j<MD; j++) {
		for(int i=0; i<(nx+2*MD); i++) {
			topo[i][j      ] = topo[i][MD     ];
			topo[i][ny+MD+j] = topo[i][ny+MD-1];

			depth[i][j      ] = depth[i][MD     ];
			depth[i][ny+MD+j] = depth[i][ny+MD-1];
		}
	}
	

}

clock_t RunKernels::run()
{
	cudaMallocHost((void **)&TotalStep_h, sizeof(double));
	cudaMallocHost((void **)&dt_h, sizeof(double));

	cudaMallocHost((void **)&depth_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&topo_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&speed_h, sizeof(double) * arraySize);

	cudaMallocHost((void **)&resultHs_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&resultHf_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&resultUs_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&resultVs_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&resultUf_h, sizeof(double) * arraySize);
	cudaMallocHost((void **)&resultVf_h, sizeof(double) * arraySize);

	cudaMallocHost((void **)&resultphi_h, sizeof(double) * arraySize);
	
	cudaMallocHost((void **)&bfkt_h,   sizeof(double) * arraySize * 3);
	cudaMallocHost((void **)&svec_h,   sizeof(double) * arraySize * 2);
	cudaMallocHost((void **)&cvalue_h, sizeof(double) * arraySize);
		
	for(int i = 0; i < nxd; i++){
		for(int j = 0; j < nyd; j++){
		
			 topo_h[j*nxd+i] =  topo[i][j];
			depth_h[j*nxd+i] = depth[i][j];

		}
	}
    
	clock_t start, end;

	start = clock(); //cuda start
	memoryMalloc();
	
	kernelStep();
	
	freeMemory();
	end = clock(); //cuda stop

	// outputFile();
	
	return end - start;
}

void RunKernels::kernelStep()
{

	cudaMemset(dev_TotalTime, 0, sizeof(double));
	cudaMemset(dev_dt, 0.0, sizeof(double));
	cudaMemset(dev_dtval, 0.0, sizeof(double));
	cudaMemset(dt_h, 0.0, sizeof(double));
	
	cudaMemcpy(dev_topo,   topo_h, sizeof(double) * arraySize, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_depth, depth_h, sizeof(double) * arraySize, cudaMemcpyHostToDevice);

	cudaMallocHost((void **)&result_h, sizeof(double) * arraySize);

	int bx = (nxd + BLOCK_SIZE - 1) / BLOCK_SIZE;
	int by = (nyd  + BLOCK_SIZE - 1) / BLOCK_SIZE;
	dim3 blocksPerGrid(bx, by);
	dim3 threadsPerBlock(BLOCK_SIZE, BLOCK_SIZE);


	double Htmp, hstmp, hftmp, ustmp, uftmp, vstmp, vftmp, phitmp;

	int outputStep = TotalSim/EachOut;
	double tf[outputStep+1] = {0.0};

	for(int ii=1;ii<=outputStep;ii++){
		tf[ii] = ii*EachOut*10;
	}

	double outtime[outputStep+1]={0};
	int iter = 1, nt, io;
	int nstop = 0, schreiben =0;
	// int Totalnt = 0;
	int outsteplen = sizeof(tf)/sizeof(tf)[0];

	makeTopo1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_bfkt,
		MINX, MINY,
		dx, dy,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_topo,
		dev_bfkt,
		dx, dy,
		nxd, nyd, nx, ny);
	cudaDeviceSynchronize();


	makeTopo3Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_bfkt,
		dev_posx, dev_posy,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo4Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_bfkt,
		dev_posx, dev_posy,
		dev_dxdxi11, dev_dxdxi12,
		dev_dxdxi21, dev_dxdxi22,
		dev_dbdx, dev_dbdy,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo5Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dbdx, dev_dbdy,
		dev_cvalue,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo6Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dbdx, dev_dbdy,
		dev_cvalue,
		dev_svec,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo7Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dbdx, dev_dbdy,
		dev_cvalue,
		dev_svec,
		dev_Jacb31, dev_Jacb32,
		dev_dxdxi11, dev_dxdxi12, 
		dev_dxdxi21, dev_dxdxi22,
		dev_dettmp,  
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo8Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_cvalue, dev_svec,
		dev_Jacb31, dev_Jacb32,
		dev_dxdxi11, dev_dxdxi12, 
		dev_dxdxi21, dev_dxdxi22,
		dev_dettmp,  
		dev_Detmin,
		dev_i_ddxi11, dev_i_ddxi12, 
		dev_i_ddxi21, dev_i_ddxi22, 
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo9Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_cvalue, dev_svec,
		dev_i_ddxi11, dev_i_ddxi12, 
		dev_i_ddxi21, dev_i_ddxi22,
		dev_invJ11, dev_invJ12, dev_invJ13,  
		dev_invJ21, dev_invJ22, dev_invJ23,
		dev_invJ31, dev_invJ32, dev_invJ33,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo10Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_cvalue, dev_depth,
		dev_u,
		phiS0,
		nxd, nyd);
	cudaDeviceSynchronize();

	makeTopo11Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_u,
		dev_tande, delta0,
		nxd, nyd);
	cudaDeviceSynchronize();

	Boundary1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dxdxi11, dev_dxdxi12,
		dev_dxdxi21, dev_dxdxi22,
		dev_cvalue, dev_Detmin,
		dev_svec,
		dev_Jacb31, dev_Jacb32,
		dev_invJ11, dev_invJ12, dev_invJ13,
		dev_invJ21, dev_invJ22, dev_invJ23,
		dev_invJ31, dev_invJ32, dev_invJ33,  
		nxd, nyd, nx ,ny);
	cudaDeviceSynchronize();

	Boundary2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dxdxi11, dev_dxdxi12,
		dev_dxdxi21, dev_dxdxi22,
		dev_cvalue, dev_Detmin,
		dev_svec,
		dev_Jacb31, dev_Jacb32,
		dev_invJ11, dev_invJ12, dev_invJ13,
		dev_invJ21, dev_invJ22, dev_invJ23,
		dev_invJ31, dev_invJ32, dev_invJ33,  
		nxd, nyd, nx ,ny);
	cudaDeviceSynchronize();

	JacobKernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_svec, dev_cvalue,
		dev_posx, dev_posy, 
		dev_J13dxi, dev_J23dxi, dev_J33dxi, 
		dev_J13det, dev_J23det, dev_J33det,  
		nxd, nyd);
	cudaDeviceSynchronize();

	Boundary3Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_J13dxi, dev_J23dxi, dev_J33dxi, 
		dev_J13det, dev_J23det, dev_J33det,  
		nxd, nyd, nx ,ny);
	cudaDeviceSynchronize();

	Boundary4Kernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_J13dxi, dev_J23dxi, dev_J33dxi, 
		dev_J13det, dev_J23det, dev_J33det,  
		nxd, nyd, nx ,ny);
	cudaDeviceSynchronize();

	cudaMemcpy(bfkt_h,   &dev_bfkt[0],   sizeof(double)* arraySize * 3, cudaMemcpyDeviceToHost);
	cudaMemcpy(svec_h,   &dev_svec[0],   sizeof(double)* arraySize * 2, cudaMemcpyDeviceToHost);
	cudaMemcpy(cvalue_h, &dev_cvalue[0], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
				

	char check_Topo[100]; sprintf(check_Topo,"./result2Pb/DEM.dat");
	FILE  *fpTopo;

	fpTopo=fopen(check_Topo, "w");
	fprintf(fpTopo, "VARIABLES = \"x\", \"y\", \"z\", \"c\", \"S1\", \"S2\"\n ");
	for (int i=MD;i<nxd-MD;i++) {
	  for (int j=MD;j<nyd-MD;j++) {
		fprintf(fpTopo, "%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",bfkt_h[0 * nxd * nyd + j * nxd + i]*0.1, bfkt_h[1 * nxd * nyd + j * nxd + i]*0.1, bfkt_h[2 * nxd * nyd + j * nxd + i]*0.1, cvalue_h[j * nxd + i],svec_h[0 * nxd * nyd + j * nxd + i],svec_h[1 * nxd * nyd + j * nxd + i]);
	  } 
	}	
	fclose(fpTopo);

	FILE  *fpInit;

	if ((fpInit=fopen("./result2Pb/001.dat", "w")) == NULL)
    {
      fclose(fpInit);
      exit(0);
    }
	fprintf(fpInit, "VARIABLES = \"H\", \"phi\", \"Us\", \"Uf\", \"Vs\", \"Vf\"\n ");
	for (int i=MD;i<nxd-MD;i++) {
		for (int j=MD;j<nyd-MD;j++) {
			fprintf(fpInit, "%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",0.1*depth[i][j],phiS0, 0.0, 0.0, 0.0, 0.0);
		} 
	}
	fclose(fpInit);

	FILE  *fpInfo;
	if ((fpInfo=fopen("./result2Pb/Info.dat", "w")) == NULL)
    {
      fclose(fpInfo);
      exit(0);
	}
	fprintf(fpInfo, "VARIABLES = \"x-point\", \"y-point\", \"dx\", \"dy\", \"xllcorner\", \"yllcorner\", \"TotalStep\"\n ");
	fprintf(fpInfo, "\t%d\t\t %d\t %10.2f\t %10.2f\t %10.4f\t %10.4f\t %d", NX, NY, (dx*0.1), (dy*0.1), xllcorner, yllcorner,(outputStep+1));
	fclose(fpInfo);

	cudaFree(dev_i_ddxi11);
	cudaFree(dev_i_ddxi12);
	cudaFree(dev_i_ddxi21);
	cudaFree(dev_i_ddxi22);

	cudaFree(dev_dettmp);
	cudaFree(dev_bfkt);

	MeanKernel <<<blocksPerGrid, threadsPerBlock>>>(
		dev_result, 
		dev_dxdxi11, dev_dxdxi21,
		dev_dxdxi12, dev_dxdxi22,

		dev_J13dxi, dev_J23dxi, dev_J33dxi, 
		dev_J13det, dev_J23det, dev_J33det, 
		
		dev_invJ11, dev_invJ12, dev_invJ13, 
		dev_invJ21, dev_invJ22, dev_invJ23, 
		dev_invJ31, dev_invJ32, dev_invJ33, 

		dev_Detmin, dev_cvalue, dev_svec, 

		dev_dxdxi11_avgEW, dev_dxdxi21_avgEW,
		dev_dxdxi12_avgSN, dev_dxdxi22_avgSN, 
		
		dev_J13dxi_avgEW, dev_J23dxi_avgEW, dev_J33dxi_avgEW,
		dev_J13det_avgEW, dev_J23det_avgEW, dev_J33det_avgEW,

		dev_J13dxi_avgSN, dev_J23dxi_avgSN, dev_J33dxi_avgSN,
		dev_J13det_avgSN, dev_J23det_avgSN, dev_J33det_avgSN,

		dev_invJ11_avgEW, dev_invJ12_avgEW, dev_invJ13_avgEW,
		dev_invJ21_avgEW, dev_invJ22_avgEW, dev_invJ23_avgEW,
		dev_invJ31_avgEW, dev_invJ32_avgEW, dev_invJ33_avgEW,

		dev_invJ11_avgSN, dev_invJ12_avgSN, dev_invJ13_avgSN,
		dev_invJ21_avgSN, dev_invJ22_avgSN, dev_invJ23_avgSN,
		dev_invJ31_avgSN, dev_invJ32_avgSN, dev_invJ33_avgSN,

		dev_Detmin_avgEW, dev_Detmin_avgSN,
		dev_cval_avgEW,   dev_cval_avgSN,
		dev_svec_avgEW,   dev_svec_avgSN,  

		nxd, nyd);
	cudaDeviceSynchronize();

	

	for (nt = 1; (!nstop) && (nt<100000); nt++){
		for (io=0; io<2; io++){

			if(io == 0){
				UzeroKernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_uzero, 
					nxd, nyd);
				cudaDeviceSynchronize();
			}

			cudaDeviceSynchronize();

			Boundary5Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u,  
				nxd, nyd, nx ,ny);
			cudaDeviceSynchronize();

			Boundary6Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u,  
				nxd, nyd, nx ,ny);
			cudaDeviceSynchronize();

			Boundary7Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u,  
				nxd, nyd, nx ,ny);
			cudaDeviceSynchronize();

			Boundary8Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u,  
				nxd, nyd, nx ,ny);
			cudaDeviceSynchronize();

			Boundary9Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_Hpx, dev_Hpy,
				dev_Ppx, dev_Ppy,
				dev_PDx, dev_PDy,
				dev_ux,  dev_uy,
				dev_apEW, dev_apSN,
				dev_apFEW, dev_apFSN,   
				nxd, nyd);
			cudaDeviceSynchronize();

			TVD1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u, 
				dev_dux, dev_duy,
				nxd, nyd);
			cudaDeviceSynchronize();

			TVD2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_dux, dev_duy,
				dev_sgnAx, dev_sgnBx,
				dev_sgnAy, dev_sgnBy, 
				nxd, nyd);
			cudaDeviceSynchronize();

			TVD3Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_dux, dev_duy,
				dev_sgnAx, dev_sgnBx,
				dev_sgnAy, dev_sgnBy, 
				dev_t1x, dev_t2x,
				dev_t1y, dev_t2y, 
				nxd, nyd);
			cudaDeviceSynchronize();

			TVD4Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_t1x, dev_t2x,
				dev_t1y, dev_t2y,
				dev_sgnAx, dev_sgnBx,
				dev_sgnAy, dev_sgnBy, 
				nxd, nyd);
			cudaDeviceSynchronize();

			TVD5Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_t1x, dev_t2x,
				dev_t1y, dev_t2y,
				dev_sgnAx, dev_sgnBx,
				dev_sgnAy, dev_sgnBy, 
				dev_ux, dev_uy, 
				nxd, nyd);
			cudaDeviceSynchronize();

			InterfacesKernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u, 
				dev_ux, dev_uy,
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				nxd, nyd);
			cudaDeviceSynchronize();

			Interfaces2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				nxd, nyd,nx ,ny);
			cudaDeviceSynchronize();


			KeepPositivi1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				nxd, nyd);
			cudaDeviceSynchronize();

			KeepPositivi2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				dev_vexE, dev_veyE, 
				dev_vexW, dev_veyW,
				dev_vexFE, dev_veyFE, 
				dev_vexFW, dev_veyFW,

				dev_vexN, dev_veyN, 
				dev_vexS, dev_veyS,
				dev_vexFN, dev_veyFN, 
				dev_vexFS, dev_veyFS,

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				dev_vexE, dev_veyE, 
				dev_vexW, dev_veyW,
				dev_vexFE, dev_veyFE, 
				dev_vexFW, dev_veyFW,

				dev_vexN, dev_veyN, 
				dev_vexS, dev_veyS,
				dev_vexFN, dev_veyFN, 
				dev_vexFS, dev_veyFS,

				dev_w_wertE, dev_w_wertW,
				dev_w_wertFE, dev_w_wertFW,

				dev_w_wertN, dev_w_wertS,
				dev_w_wertFN, dev_w_wertFS,
				
				dev_svec, dev_cvalue,

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux3Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,  
				dev_vexE, dev_veyE, 
				dev_vexW, dev_veyW,
				dev_vexFE, dev_veyFE, 
				dev_vexFW, dev_veyFW,

				dev_vexN, dev_veyN, 
				dev_vexS, dev_veyS,
				dev_vexFN, dev_veyFN, 
				dev_vexFS, dev_veyFS,

				dev_w_wertE, dev_w_wertW,
				dev_w_wertFE, dev_w_wertFW,

				dev_w_wertN, dev_w_wertS,
				dev_w_wertFN, dev_w_wertFS,
				
				dev_q_xiE , dev_q_etE,
				dev_q_xiW , dev_q_etW,
				dev_q_xiFE, dev_q_etFE,
				dev_q_xiFW, dev_q_etFW,

				dev_NpressFE, dev_NpressFW, dev_M11EW,
				
				dev_invJ11_avgEW, dev_invJ12_avgEW, dev_invJ13_avgEW,
				dev_invJ21_avgEW, dev_invJ22_avgEW, dev_invJ23_avgEW, 
				dev_cval_avgEW,

				dev_q_xiN , dev_q_etN,
				dev_q_xiS , dev_q_etS,
				dev_q_xiFN, dev_q_etFN,
				dev_q_xiFS, dev_q_etFS,

				dev_NpressFN, dev_NpressFS, dev_M22SN,
				
				dev_invJ11_avgSN, dev_invJ12_avgSN, dev_invJ13_avgSN,
				dev_invJ21_avgSN, dev_invJ22_avgSN, dev_invJ23_avgSN, 
				dev_cval_avgSN,

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux4Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS, 
				
				dev_q_xiE , dev_q_xiW,
				dev_q_xiFE, dev_q_xiFW,

				dev_NpressFE, dev_NpressFW, dev_invJ11_avgEW,
				
				dev_apE, dev_apW,
				dev_apFE, dev_apFW, 

				dev_q_etN , dev_q_etS,
				dev_q_etFN, dev_q_etFS,

				dev_NpressFN, dev_NpressFS, dev_invJ22_avgSN,

				dev_apN, dev_apS,
				dev_apFN, dev_apFS, 

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux5Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_apE, dev_apW,
				dev_apFE, dev_apFW, 

				dev_apEW , dev_apFEW,

				dev_apN, dev_apS,
				dev_apFN, dev_apFS, 

				dev_apSN, dev_apFSN,

				nxd, nyd);
			cudaDeviceSynchronize();


			Flux6Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_apEW, dev_apFEW,
				dev_apSN, dev_apFSN,

				dev_em_x ,  dev_em_y,
				dev_em_Fx, dev_em_Fy,

				dev_czw1x, dev_czw2x,
				dev_czwF1x,dev_czwF2x, 
				dev_czw1y, dev_czw2y,
				dev_czwF1y,dev_czwF2y,
				
				dev_uE, dev_uW, 
				dev_uN, dev_uS,
				
				dev_cval_avgEW, dev_cval_avgSN,
				dev_Detmin_avgEW, dev_Detmin_avgSN,
				dev_M11EW, dev_M22SN,

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux7Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_FpE, dev_FpW,
				dev_GpN, dev_GpS,

				dev_czw1x, dev_czw2x,
				dev_czwF1x,dev_czwF2x, 
				dev_czw1y, dev_czw2y,
				dev_czwF1y,dev_czwF2y,
				
				dev_uE, dev_uW, 
				dev_uN, dev_uS,
				
				dev_Detmin_avgEW, dev_Detmin_avgSN,
				
				dev_q_xiE, dev_q_xiFE,
				dev_q_xiW, dev_q_xiFW,

				dev_q_etN, dev_q_etFN,
				dev_q_etS, dev_q_etFS,

				dev_dxdxi11_avgEW, dev_dxdxi21_avgEW,  
				dev_dxdxi12_avgSN, dev_dxdxi22_avgSN, 

				nxd, nyd);
			cudaDeviceSynchronize();


			CFL1Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_em_x, dev_em_y,
				dev_em_Fx, dev_em_Fy,
				dev_em_valS, dev_em_valF,
				dx, dy,
				nxd, nyd);
			cudaDeviceSynchronize();

			CFL2Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_em_valS, dev_em_valF,
				dev_Val,
				nxd, nyd);
			cudaDeviceSynchronize();
			
			if(io==0){

				int threads = 256;
				int blocks = 256;//min((arraySize + threads - 1) / threads, 256);
				
				reduceKernel <<<blocks, threads>>> (dev_Val, dev_max, arraySize);
				cudaDeviceSynchronize();

				reduceKernel <<<1, blocks>>> (dev_max, dev_maxW, blocks);
				cudaDeviceSynchronize();

				CFL3Kernel <<<1, threadsPerBlock>>>(
					dev_result, 
					dev_dtval, dev_maxW, dev_TotalTime);
				cudaDeviceSynchronize();
				cudaMemcpy(dt_h, dev_dtval, sizeof(double)* 1, cudaMemcpyDeviceToHost);
				cudaMemcpy(TotalStep_h, dev_TotalTime, sizeof(double)* 1, cudaMemcpyDeviceToHost);
				
				if((*TotalStep_h + *dt_h) >= tf[iter]){
					*dt_h =  tf[iter] - *TotalStep_h;
					schreiben = 1;
					cudaMemcpy(dev_dtval, dt_h, sizeof(double)* 1, cudaMemcpyHostToDevice);
					// cout << *TotalStep_h <<endl; 
					iter++;
				}
				cudaDeviceSynchronize();
				// printf("%10.5f", *dt_h);

			}



			CFL4Kernel <<<1, threadsPerBlock>>>(
				dev_result, 
				dev_dt, dev_dtval,dev_TotalTime,io);
			cudaDeviceSynchronize();


			Flux8Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_Hpx, dev_Hpy,
				dev_Ppx, dev_Ppy,
				dev_FpE, dev_FpW,
				dev_GpN, dev_GpS,

				dev_apEW, dev_apFEW,
				dev_apSN, dev_apFSN, 
				dev_uE, dev_uW, 
				dev_uN, dev_uS,
				dev_u,
				dev_ux, dev_uy,
				
				dev_Detmin_avgEW, dev_Detmin_avgSN,
				dev_cval_avgEW, dev_cval_avgSN,
				
				dev_invJ11_avgEW, dev_invJ12_avgEW, dev_invJ13_avgEW,  
				dev_invJ21_avgEW, dev_invJ22_avgEW, dev_invJ23_avgEW,  
				dev_invJ31_avgEW, dev_invJ32_avgEW, dev_invJ33_avgEW,  
				
				dev_invJ11_avgSN, dev_invJ12_avgSN, dev_invJ13_avgSN,  
				dev_invJ21_avgSN, dev_invJ22_avgSN, dev_invJ23_avgSN,  
				dev_invJ31_avgSN, dev_invJ32_avgSN, dev_invJ33_avgSN,  
				
				dev_dudxE, dev_dvdxE, 
				dev_dudyE, dev_dvdyE,

				dev_dudxN, dev_dvdxN, 
				dev_dudyN, dev_dvdyN,
				
				dx, dy,

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux9Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_duxidxix, dev_dvetdxix,
				dev_duxidetx, dev_dvetdetx,
				dev_duxidxiy, dev_dvetdxiy,
				dev_duxidety, dev_dvetdety,
				
				dev_invJ11_avgEW, dev_invJ12_avgEW, dev_invJ13_avgEW,  
				dev_invJ21_avgEW, dev_invJ22_avgEW, dev_invJ23_avgEW,  
				dev_invJ31_avgEW, dev_invJ32_avgEW, dev_invJ33_avgEW,  
				
				dev_invJ11_avgSN, dev_invJ12_avgSN, dev_invJ13_avgSN,  
				dev_invJ21_avgSN, dev_invJ22_avgSN, dev_invJ23_avgSN,  
				dev_invJ31_avgSN, dev_invJ32_avgSN, dev_invJ33_avgSN,  
				
				dev_dudxE, dev_dvdxE, 
				dev_dudyE, dev_dvdyE,

				dev_dudxN, dev_dvdxN, 
				dev_dudyN, dev_dvdyN,

				nxd, nyd);
			cudaDeviceSynchronize();


			Flux10Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_PDx, dev_PDy,

				dev_uE, dev_uW, 
				dev_uN, dev_uS,
				
				dev_Detmin_avgEW, dev_Detmin_avgSN,

				dev_duxidxix, dev_dvetdxix,
				dev_duxidetx, dev_dvetdetx,
				dev_duxidxiy, dev_dvetdxiy,
				dev_duxidety, dev_dvetdety,
				
				dev_invJ11_avgEW, dev_invJ12_avgEW,
				dev_invJ21_avgEW, dev_invJ22_avgEW, 

				dev_invJ11_avgSN, dev_invJ12_avgSN,   
				dev_invJ21_avgSN, dev_invJ22_avgSN, 

				nxd, nyd);
			cudaDeviceSynchronize();

			Flux11Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u, 
				dev_vex, dev_vey,
				dev_vexF, dev_veyF, 
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux12Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_w_wert, dev_w_wertF,
				dev_vex, dev_vey,
				dev_vexF, dev_veyF, 
				dev_svec, dev_cvalue,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux13Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_w_wert, dev_w_wertF,
				dev_vex, dev_vey,
				dev_usw, dev_vel, 
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux14Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_w_wert, dev_w_wertF,
				dev_vex, dev_vey,
				dev_vexF, dev_veyF,
				dev_vexw, dev_veyw, 
				dev_usw, dev_vel, 
				dev_q_xi, dev_q_et,
				dev_q_xiF, dev_q_etF,
				dev_invJ11, dev_invJ12, dev_invJ13,
				dev_invJ21, dev_invJ22, dev_invJ23,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux15Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_w_wert, dev_w_wertF,
				dev_vex, dev_vey,
				dev_vexF, dev_veyF,
				dev_q_xi, dev_q_et,
				dev_q_xiF, dev_q_etF,
				dev_J13dxi, dev_J23dxi, dev_J33dxi,
				dev_J13det, dev_J23det, dev_J33det,
				dev_Ac, dev_AcF,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux16Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_Npress1, dev_Npress2, dev_NpressF,
				dev_Ac, dev_AcF,
				dev_Detmin, dev_cvalue,
				dev_u,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux17Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_Npress1, dev_Npress2, dev_NpressF,
				dev_Ac, dev_AcF,
				dev_Detmin, dev_svec,
				dev_vex, dev_vey,
				dev_vexw, dev_veyw,
				dev_vexF, dev_veyF,
				dev_tande,
				dev_u, dev_s,
				Cd, N_R, varTheta,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux18Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_v, dev_Detmin, 
				dev_Hpx, dev_Hpy,
				dev_Ppx, dev_Ppy,
				dev_PDx, dev_PDy,
				dev_s,
				dev_u, dev_uzero,
				dev_dt,
				dx, dy,
				N_R,
				nxd, nyd);
			cudaDeviceSynchronize();

			Flux19Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_v, dev_Detmin, 
				dev_u, 
				nxd, nyd);
			cudaDeviceSynchronize();


			Flux20Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u, dev_utmp,
				nxd, nyd);
			cudaDeviceSynchronize();


			Flux21Kernel <<<blocksPerGrid, threadsPerBlock>>>(
				dev_result, 
				dev_u, 
				nxd, nyd);
			cudaDeviceSynchronize();
			
			if(io == 0)
			{
			
				Flux22Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_uone,
					nxd, nyd);
				cudaDeviceSynchronize();

				Flux23Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_uzero,
					dev_usxnew, dev_ufxnew, 
					dev_usxold, dev_ufxold, 
					dev_usynew, dev_ufynew, 
					dev_usyold, dev_ufyold, 
					nxd, nyd);
				cudaDeviceSynchronize();

				Flux24Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_uone,
					dev_usxnew, dev_ufxnew, 
					dev_usxold, dev_ufxold, 
					dev_usynew, dev_ufynew, 
					dev_usyold, dev_ufyold, 
					nxd, nyd);
				cudaDeviceSynchronize();
			
			}
			else{
				Flux25Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_utwo, 
					nxd, nyd);
				cudaDeviceSynchronize();

				Flux26Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_uone,
					dev_usxnew, dev_ufxnew, 
					dev_usxold, dev_ufxold, 
					dev_usynew, dev_ufynew, 
					dev_usyold, dev_ufyold, 
					nxd, nyd);
				cudaDeviceSynchronize();

				Flux27Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_utwo,
					dev_usxnew, dev_ufxnew, 
					dev_usxold, dev_ufxold, 
					dev_usynew, dev_ufynew, 
					dev_usyold, dev_ufyold, 
					nxd, nyd);
				cudaDeviceSynchronize();

				Flux28Kernel <<<blocksPerGrid, threadsPerBlock>>>(
					dev_result, 
					dev_u, dev_utwo, 
					nxd, nyd);
				cudaDeviceSynchronize();
			}


				if(io)
				{
					if(schreiben==1){
							
						cudaMemcpy(resultHs_h, &dev_u[0 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultHf_h, &dev_u[3 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultUs_h, &dev_u[1 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultVs_h, &dev_u[2 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultUf_h, &dev_u[4 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultVf_h, &dev_u[5 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
						cudaMemcpy(resultphi_h,&dev_u[6 * nyd * nxd], sizeof(double)* arraySize, cudaMemcpyDeviceToHost);
					
						
						char outfile_Web[100]; sprintf(outfile_Web,"./result2Pb/%03d.dat",iter);

						outtime[iter-1] = *TotalStep_h;
						FILE *fpTmp;
						fpTmp=fopen("./result2Pb/Time.dat", "w");
						for (int nn=0;nn<(iter);nn++){
							fprintf(fpTmp, "%20.4f", outtime[nn]);
						}
						fclose(fpTmp);


						FILE  *fpout;

						fpout=fopen(outfile_Web, "w");
						fprintf(fpout, "VARIABLES = \"H\", \"phi\", \"Us\", \"Uf\", \"Vs\", \"Vf\"\n ");
						for (int i=MD;i<nxd-MD;i++) {
							for (int j=MD;j<nyd-MD;j++) {
								Htmp  = 0.1*(resultHs_h[j * nxd + i] + resultHf_h[j * nxd + i]);
								hstmp = resultHs_h[j * nxd + i];
								hftmp = resultHf_h[j * nxd + i];

								if (hstmp > 0.00001){
									vstmp = resultVs_h[j * nxd + i]/resultHs_h[j * nxd + i];
									ustmp = resultUs_h[j * nxd + i]/resultHs_h[j * nxd + i];
									phitmp = resultphi_h[j * nxd + i];
								}else{
									vstmp = 0.0;
									ustmp = 0.0;
									phitmp = 0.0;
								}

								if (hftmp > 0.00001){
									vftmp = resultVf_h[j * nxd + i]/resultHf_h[j * nxd + i];
									uftmp = resultUf_h[j * nxd + i]/resultHf_h[j * nxd + i];
								}else{
									vftmp = 0.0;
									uftmp = 0.0;
								}

								fprintf(fpout, "%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",Htmp,phitmp,ustmp,uftmp,vstmp,vftmp);
							} 
						}
						fclose(fpout);

						
						if(iter == outsteplen){
							nstop = 1;
							// Totalnt = nt;
						}
						// }
						schreiben = 0;
					}
				}

			}
			
			
			cudaDeviceSynchronize();

	}

	if (cudaPeekAtLastError() != cudaSuccess) 
	{
			cout << cudaGetErrorString(cudaPeekAtLastError()) << endl;
	}

	cudaMemcpy(TotalStep_h, dev_TotalTime, sizeof(double)* 1, cudaMemcpyDeviceToHost);
	cout << "\nTotal time : " << *TotalStep_h/10 << " sec  ";
	// cout << "\nTotal step : " << Totalnt      << "   \n";
	fprintf(stdout, "\nTotal number of steps: %d\n", nt);
}


void RunKernels::memoryMalloc()
{
	cudaMalloc((void **)&dev_topo, sizeof(double) * arraySize); 
	cudaMalloc((void **)&dev_depth, sizeof(double) * arraySize); 
	cudaMalloc((void **)&dev_result, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_bfkt, sizeof(double) * arraySize * 3);

	cudaMalloc((void **)&dev_posx, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_posy, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dxdxi11, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dxdxi12, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dxdxi21, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dxdxi22, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dbdx  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dbdy  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_cvalue, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_svec  , sizeof(double) * arraySize * 2);
	cudaMalloc((void **)&dev_Jacb31, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Jacb32, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dettmp, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Detmin, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_i_ddxi11, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_i_ddxi12, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_i_ddxi21, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_i_ddxi22, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ11, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ12, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ13, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ21, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ22, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ23, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ31, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ32, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ33, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_u     , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_uzero , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_Hpx   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_Hpy   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_Ppx   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_Ppy   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_PDx   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_PDy   , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_ux    , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_uy    , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_dux   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_duy   , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_t1x   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_t2x   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_t1y   , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_t2y   , sizeof(double) * arraySize * 7);
	
	cudaMalloc((void **)&dev_sgnAx , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_sgnBx , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_sgnAy , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_sgnBy , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_uE , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_uW , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_uN , sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_uS , sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_tande , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_J13dxi , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23dxi , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33dxi , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J13det , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23det , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33det , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_apEW  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apSN  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dxdxi11_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dxdxi21_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dxdxi12_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dxdxi22_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_J13dxi_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23dxi_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33dxi_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_J13det_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23det_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33det_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_J13dxi_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23dxi_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33dxi_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_J13det_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J23det_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_J33det_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ11_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ12_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ13_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ21_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ22_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ23_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ31_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ32_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ33_avgEW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ11_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ12_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ13_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ21_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ22_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ23_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_invJ31_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ32_avgSN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_invJ33_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_Detmin_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Detmin_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_cval_avgEW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_cval_avgSN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_svec_avgEW , sizeof(double) * arraySize * 2);
	cudaMalloc((void **)&dev_svec_avgSN , sizeof(double) * arraySize * 2);

	cudaMalloc((void **)&dev_vexE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vexW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_w_wertE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_w_wertW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_vexFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vexFW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyFW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_w_wertFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_w_wertFW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_vexN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vexS , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_w_wertN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_w_wertS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_vexFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vexFS , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyFS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_w_wertFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_w_wertFS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_q_xiE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xiW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_q_xiFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xiFW , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etFW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_NpressFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_NpressFW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_M11EW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_q_xiN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xiS , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_q_xiFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xiFS , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etFS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_NpressFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_NpressFS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_M22SN , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_apE  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apW  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFE , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFW , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_apN  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apS  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFN , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_apFS , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_em_x , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_em_y , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_em_Fx, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_em_Fy, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_FpE, sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_FpW, sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_GpN, sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_GpS, sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_czw1x , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czw2x , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czwF1x, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czwF2x, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_czw1y , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czw2y , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czwF1y, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_czwF2y, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_em_valS, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_em_valF, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Val    , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dudxE, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvdxE, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dudyE, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvdyE, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_dudxN, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvdxN, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dudyN, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvdyN, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_duxidxix, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvetdxix, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_duxidetx, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvetdetx, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_duxidxiy, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvetdxiy, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_duxidety, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_dvetdety, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_vex    , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vey    , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vexF   , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyF   , sizeof(double) * arraySize);
	
	cudaMalloc((void **)&dev_w_wert , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_w_wertF, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_usw    , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_vel    , sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_vexw   , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_veyw   , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xi   , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_et   , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_xiF  , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_q_etF  , sizeof(double) * arraySize);
	
	cudaMalloc((void **)&dev_Ac     , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_AcF    , sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Npress1, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_Npress2, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_NpressF, sizeof(double) * arraySize);
	
	cudaMalloc((void **)&dev_s, sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_v, sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_uone, sizeof(double) * arraySize * 7);
	cudaMalloc((void **)&dev_utwo, sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_usxnew, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_ufxnew, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_usxold, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_ufxold, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_usynew, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_ufynew, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_usyold, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_ufyold, sizeof(double) * arraySize);

	cudaMalloc((void **)&dev_utmp, sizeof(double) * arraySize * 7);

	cudaMalloc((void **)&dev_waveSpeed, sizeof(double) * arraySize);
	cudaMalloc((void **)&dev_max, sizeof(double) * 256);
	cudaMalloc((void **)&dev_maxW, sizeof(double) * 1);
	cudaMalloc((void **)&dev_TotalTime, sizeof(double) * 1);
	cudaMalloc((void **)&dev_dt, sizeof(double) * 1);
	cudaMalloc((void **)&dev_dtval, sizeof(double) * 1);

}

void RunKernels::freeMemory()
{
	cudaFree(dev_topo);
	cudaFree(dev_depth);
	
	cudaFree(dev_result);

	// cudaFree(dev_bfkt);
	cudaFree(dev_posx);
	cudaFree(dev_posy);

	cudaFree(dev_dxdxi11);
	cudaFree(dev_dxdxi12);
	cudaFree(dev_dxdxi21);
	cudaFree(dev_dxdxi22);

	cudaFree(dev_dbdx);
	cudaFree(dev_dbdy);
	cudaFree(dev_cvalue);

	cudaFree(dev_svec);
	cudaFree(dev_Jacb31);
	cudaFree(dev_Jacb32);

	// cudaFree(dev_dettmp);
	cudaFree(dev_Detmin);

	// cudaFree(dev_i_ddxi11);
	// cudaFree(dev_i_ddxi12);
	// cudaFree(dev_i_ddxi21);
	// cudaFree(dev_i_ddxi22);

	cudaFree(dev_invJ11);
	cudaFree(dev_invJ12);
	cudaFree(dev_invJ13);

	cudaFree(dev_invJ21);
	cudaFree(dev_invJ22);
	cudaFree(dev_invJ23);

	cudaFree(dev_invJ31);
	cudaFree(dev_invJ32);
	cudaFree(dev_invJ33);

	cudaFree(dev_u);
	cudaFree(dev_uzero);

	cudaFree(dev_tande);

	cudaFree(dev_J13dxi);
	cudaFree(dev_J23dxi);
	cudaFree(dev_J33dxi);
	cudaFree(dev_J13det);
	cudaFree(dev_J23det);
	cudaFree(dev_J33det);

	cudaFree(dev_Hpx);
	cudaFree(dev_Hpy);
	cudaFree(dev_Ppx);
	cudaFree(dev_Ppy);
	cudaFree(dev_PDx);
	cudaFree(dev_PDy);
	cudaFree(dev_ux);
	cudaFree(dev_uy);

	cudaFree(dev_apEW);
	cudaFree(dev_apSN);
	cudaFree(dev_apFEW);
	cudaFree(dev_apFSN);

	cudaFree(dev_dux);
	cudaFree(dev_duy);

	cudaFree(dev_t1x);
	cudaFree(dev_t2x);
	cudaFree(dev_t1y);
	cudaFree(dev_t2y);

	cudaFree(dev_sgnAx);
	cudaFree(dev_sgnBx);
	cudaFree(dev_sgnAy);
	cudaFree(dev_sgnBy);

	cudaFree(dev_dxdxi11_avgEW);
	cudaFree(dev_dxdxi21_avgEW);

	cudaFree(dev_dxdxi12_avgSN);
	cudaFree(dev_dxdxi22_avgSN);

	cudaFree(dev_J13dxi_avgEW);
	cudaFree(dev_J23dxi_avgEW);
	cudaFree(dev_J33dxi_avgEW);

	cudaFree(dev_J13det_avgEW);
	cudaFree(dev_J23det_avgEW);
	cudaFree(dev_J33det_avgEW);

	cudaFree(dev_J13dxi_avgSN);
	cudaFree(dev_J23dxi_avgSN);
	cudaFree(dev_J33dxi_avgSN);

	cudaFree(dev_J13det_avgSN);
	cudaFree(dev_J23det_avgSN);
	cudaFree(dev_J33det_avgSN);

	cudaFree(dev_invJ11_avgEW);
	cudaFree(dev_invJ12_avgEW);
	cudaFree(dev_invJ13_avgEW);

	cudaFree(dev_invJ21_avgEW);
	cudaFree(dev_invJ22_avgEW);
	cudaFree(dev_invJ23_avgEW);

	cudaFree(dev_invJ31_avgEW);
	cudaFree(dev_invJ32_avgEW);
	cudaFree(dev_invJ33_avgEW);

	cudaFree(dev_invJ11_avgSN);
	cudaFree(dev_invJ12_avgSN);
	cudaFree(dev_invJ13_avgSN);

	cudaFree(dev_invJ21_avgSN);
	cudaFree(dev_invJ22_avgSN);
	cudaFree(dev_invJ23_avgSN);

	cudaFree(dev_invJ31_avgSN);
	cudaFree(dev_invJ32_avgSN);
	cudaFree(dev_invJ33_avgSN);

	cudaFree(dev_Detmin_avgEW);
	cudaFree(dev_Detmin_avgSN);

	cudaFree(dev_cval_avgEW);
	cudaFree(dev_cval_avgSN);

	cudaFree(dev_svec_avgEW);
	cudaFree(dev_svec_avgSN);

	cudaFree(dev_uE);
	cudaFree(dev_uW);
	cudaFree(dev_uN);
	cudaFree(dev_uS);

	cudaFree(dev_vexE);
	cudaFree(dev_vexW);
	cudaFree(dev_veyE);
	cudaFree(dev_veyW);

	cudaFree(dev_w_wertE);
	cudaFree(dev_w_wertW);

	cudaFree(dev_vexFE);
	cudaFree(dev_vexFW);
	cudaFree(dev_veyFE);
	cudaFree(dev_veyFW);

	cudaFree(dev_w_wertFE);
	cudaFree(dev_w_wertFW);

	cudaFree(dev_vexN);
	cudaFree(dev_vexS);
	cudaFree(dev_veyN);
	cudaFree(dev_veyS);

	cudaFree(dev_w_wertFN);
	cudaFree(dev_w_wertFS);

	cudaFree(dev_q_xiE);
	cudaFree(dev_q_etE);
	cudaFree(dev_q_xiW);
	cudaFree(dev_q_etW);

	cudaFree(dev_q_xiFE);
	cudaFree(dev_q_etFE);
	cudaFree(dev_q_xiFW);
	cudaFree(dev_q_etFW);

	cudaFree(dev_NpressFE);
	cudaFree(dev_NpressFW);

	cudaFree(dev_M11EW);

	cudaFree(dev_q_xiN);
	cudaFree(dev_q_etN);
	cudaFree(dev_q_xiS);
	cudaFree(dev_q_etS);

	cudaFree(dev_q_xiFN);
	cudaFree(dev_q_etFN);
	cudaFree(dev_q_xiFS);
	cudaFree(dev_q_etFS);

	cudaFree(dev_NpressFN);
	cudaFree(dev_NpressFS);

	cudaFree(dev_M22SN);

	cudaFree(dev_apE);
	cudaFree(dev_apW);
	cudaFree(dev_apFE);
	cudaFree(dev_apFW);

	cudaFree(dev_apN);
	cudaFree(dev_apS);
	cudaFree(dev_apFN);
	cudaFree(dev_apFS);

	cudaFree(dev_em_x);
	cudaFree(dev_em_y);
	cudaFree(dev_em_Fx);
	cudaFree(dev_em_Fy);

	cudaFree(dev_FpE);
	cudaFree(dev_FpW);
	cudaFree(dev_GpN);
	cudaFree(dev_GpS);


	cudaFree(dev_s);
	cudaFree(dev_v);


	cudaFree(dev_waveSpeed);
	cudaFree(dev_max);
	cudaFree(dev_maxW);
	cudaFree(dev_TotalTime);
	cudaFree(dev_dt);


	
}