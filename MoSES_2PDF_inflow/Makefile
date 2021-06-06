all:
	@nvcc -o MoSES_2PDF_inflow main.cu adNOC_2Pb_runKernels.cu adNOC_2Pb_kernels.cu
	@./MoSES_2PDF_inflow
clean:
	@rm MoSES_2PDF_inflow
