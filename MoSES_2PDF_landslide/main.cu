#include <stdio.h>
#include <stdlib.h>

#include "adNOC_2Pb_runKernels.cuh"

using namespace std;

void output_spec(const cudaDeviceProp sDevProp)
{
	printf( "Device name: %s\n", sDevProp.name );
	printf( "Device memory: %lu\n", sDevProp.totalGlobalMem );
	printf( " Memory per-block: %lu\n", sDevProp.sharedMemPerBlock );
	printf( " Register per-block: %d\n", sDevProp.regsPerBlock );
	printf( " Warp size: %d\n", sDevProp.warpSize );
	printf( " Memory pitch: %lu\n", sDevProp.memPitch );
	printf( " Constant Memory: %lu\n", sDevProp.totalConstMem );
	printf( "Max thread per-block: %d\n", sDevProp.maxThreadsPerBlock );
	printf( "Max thread dim: ( %d, %d, %d )\n", sDevProp.maxThreadsDim[0], sDevProp.maxThreadsDim[1], sDevProp.maxThreadsDim[2] );
	printf( "Max grid size: ( %d, %d, %d )\n", sDevProp.maxGridSize[0], sDevProp.maxGridSize[1], sDevProp.maxGridSize[2] );
	printf( "Ver: %d.%d\n", sDevProp.major, sDevProp.minor );
	printf( "Clock: %d\n", sDevProp.clockRate );
	printf( "textureAlignment: %lu\n\n", sDevProp.textureAlignment );
}

void cuda_information() 
{
	int  iDeviceCount = 0;
	cudaGetDeviceCount( &iDeviceCount );
	// printf( "Number of GPU: %d\n", iDeviceCount );

	if( iDeviceCount == 0 )
	{
		printf( "No supported GPU\n" );
		return;
	}

	for( int i = 0; i < iDeviceCount; ++ i )
	{
		// printf( "\n=== Device %i ===\n", i );
		cudaDeviceProp  sDeviceProp;
		cudaGetDeviceProperties( &sDeviceProp, i );
		// output_spec( sDeviceProp );
	}
}

int main(int argc,char *argv[])
{
	cuda_information();

	RunKernels *kernel = new RunKernels((char*) "./par_List");
	clock_t t = kernel->run();

	double sec = (double) t / CLOCKS_PER_SEC;
	printf("Time used on GPU: %.2f sec\n", sec);
	// printf("Total step time : %f sec\n", (*kernel->TotalStep_h));

    return 0;
}