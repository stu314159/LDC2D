//C/C++ includes
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <cstdio>

//CUDA includes
#include <cuda_runtime.h>

// My Includes
#include "lbm_lib.h"
#include "vtk_lib.h"


// a convenient macro for error checking, if desired...
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


using namespace std;

int main(int argc, char * argv[]){
  
  // the user may forget an input; remind him/her
  if(argc<6){
    cout << "Fewer than 5 input arguments detected!" << endl;
    cout << "Usage: >>ldc2D [Re] [N] [TS] [omega] [dataFreq] where:" << endl;
    cout << "[Re] = flow Reynolds number." << endl;
    cout << "[N] = Number of lattice points along the cavity." << endl;
    cout << "[TS] = Number of time steps to perform." << endl;
    cout << "[omega] = relaxation parameter." << endl;
    cout << "[dataFreq] = # time steps between data outputs" << endl;
    cout << "Exiting the program.  Please try again." << endl;
    exit(1);
  } 
  
  float Re = (float)atof(argv[1]);
  float N = (int)atoi(argv[2]);
  float numTs = (int)atoi(argv[3]);
  float omega = (float)atof(argv[4]);
  int dataFrequency = (int)atoi(argv[5]);

  cout << "Re = " << Re << endl;
  cout << "N = " << N << endl;
  cout << "numTs = " << numTs << endl;
  cout << "omega = " << omega << endl;


  //for this problem, we know that there will be N*N lattice points
  // with 9 lattice directions per point.  
  // make use of this knowledge to simplify the code:
  const int nnodes = N*N;
  const int numSpd = 9;

   //basic host-side data arrays
  float* fDD = new float[nnodes*numSpd];
  int* ndType = new int[nnodes];
  float* ux = new float[nnodes];
  float* uy = new float[nnodes];
  float* pressure = new float[nnodes];
  float* uMag = new float[nnodes];

  float* xCoord = new float[nnodes];//lattice coordinates
  float* yCoord = new float[nnodes];

  // to simplify your life, I also do not allow you to pick the fluid.
  const float rho = 965.3; // density
  const float nu = 0.06/rho; // kinematic viscosity

  // get coordinate information; populate xCoord and yCoord
  LDC2D_getGeometry(xCoord,yCoord,1.0,1.0,N);

  // host side scaling and BC variables
  float u;
  float u_conv;
  float t_conv;
  float p_conv;
  float l_conv;

  // call setup function to initialize fDD and ndType
  // as well as to get scaling data
  LDC2D_setup(Re,N,omega,rho,nu,
		  fDD,ndType,u,
		  u_conv,t_conv,p_conv);
  l_conv = u_conv*t_conv;



  cout << "LBM mach number = " << u << endl;
  // user should check that this is no greater than ~0.1.
  cout << "Do you wish to continue? (y/[n]) " << endl;;
  char cFlag;
  //cin >> cFlag;
 cFlag = 'y'; 
 if((cFlag=='y') | (cFlag=='Y')){
	  cout << "Ok!  Let's go!" << endl;
  }else{
	cout << "Ok!  Terminating Application.  Better luck next time." << endl;
	exit(1);
  }

  //declare device-side data arrays
  float * fEven;
  float * fOdd;
  int * ndType_d;
  float * ux_d;
  float * uy_d;
  float * pressure_d;

  // Uncomment to select a random GPU from the local server.
  // to avoid everybody just using device 0
 int numDevices;
 gpuErrchk(cudaGetDeviceCount(&numDevices));
 srand(time(NULL));
 int myDevice = rand()%numDevices;
 gpuErrchk(cudaSetDevice(myDevice));



  //allocate memory for these objects
  gpuErrchk(cudaMalloc((void**)&fEven,nnodes*numSpd*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&fOdd,nnodes*numSpd*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&ndType_d,nnodes*sizeof(int)));
  gpuErrchk(cudaMalloc((void**)&ux_d,nnodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&uy_d,nnodes*sizeof(float)));
  gpuErrchk(cudaMalloc((void**)&pressure_d,nnodes*sizeof(float)));

  //copy relevant input data to the device
  gpuErrchk(cudaMemcpy(fEven,fDD,nnodes*numSpd*sizeof(float),cudaMemcpyHostToDevice));
  //cudaMemcpy(fOdd,fDD,nnodes*numSpd*sizeof(float),cudaMemcpyHostToDevice);
  gpuErrchk(cudaMemcpy(ndType_d,ndType,nnodes*sizeof(int),cudaMemcpyHostToDevice));


  // define some utility variables for visualization
  string densityFileStub("pressure");
  string velocityFileStub("velocityMagnitude");
  string vectorVelFileStub("velocity");
  string velocityCSVFileStub("velocityMagCSV");
  string fileSuffix(".vtk");
  string fileSuffixCSV(".csv");
  stringstream ts_ind;
  string ts_ind_str;
  int vtk_ts = 0;
  string fileName1;
  string fileName2;
  string fileName3;
  string fileNameCSV;
  string dataName1("pressure");
  string dataName2("velocityMagnitude");
  string dataName3("velocity");
  int dims[3];
  dims[0]=N; dims[1]=N; dims[2]=1;
  float origin[3];
  origin[0]=0.; origin[1]=0.; origin[2]=0.;
  float spacing[3];
  spacing[0]=l_conv; spacing[1]=l_conv; spacing[2]=l_conv;

  // zero w-values for the visualization.
  float * w = new float[nnodes];
  float * z = new float[nnodes];
  for(int nd =0; nd<nnodes; nd++){
	  w[nd]=0.;
	  z[nd]=0.;
  }

  //execute time steps.  Periodically, write output data for further processing.
  for(int ts = 0; ts<numTs;ts++){

	  // say something comforting about your progress
	  if(ts%1000==0){
		cout << "Executing time step " << ts << endl;
	  }

	  if(ts%2==0){
		  LDC2D_timestep(fOdd,fEven,omega,ndType_d,u,N);


	  }else{
		  LDC2D_timestep(fEven,fOdd,omega,ndType_d,u,N);

	  }


	  // for the first data set and for every [dataFrequency] steps, record macroscopic dependent variables
	  // for later visualization.
	  if((ts%dataFrequency)==0){
		  // compute macroscopic dependent variables for output
		  LDC2D_getVelocityAndDensity(ux_d,uy_d,pressure_d,u_conv,p_conv,fEven,N);

		  // transfer data to the host
		  gpuErrchk(cudaMemcpy(ux,ux_d,nnodes*sizeof(float),cudaMemcpyDeviceToHost));
		  gpuErrchk(cudaMemcpy(uy,uy_d,nnodes*sizeof(float),cudaMemcpyDeviceToHost));
		  gpuErrchk(cudaMemcpy(pressure,pressure_d,nnodes*sizeof(float),cudaMemcpyDeviceToHost));

		  // compute velocity magnitude (on the host)
		  LDC2D_getVelocityMagnitude(uMag,ux,uy,N);
		  // set pressure relative to the central lattice point (on the host)
		  LDC2D_getRelativePressure(pressure,N);

		  // set file names
		  ts_ind << vtk_ts; vtk_ts++;
		  fileName1 = densityFileStub+ts_ind.str()+fileSuffix;
		  fileName2 = velocityFileStub+ts_ind.str()+fileSuffix;
		  fileName3 = vectorVelFileStub+ts_ind.str()+fileSuffix;
		  fileNameCSV= velocityCSVFileStub+ts_ind.str()+fileSuffixCSV;
		  ts_ind.str("");

		  // output data file.
		  SaveVTKImageData_ascii(pressure,fileName1,dataName1,origin,spacing,dims);
		  SaveVTKImageData_ascii(uMag,fileName2,dataName2,origin,spacing,dims);
		  SaveVTKStructuredGridVectorAndMagnitude_ascii(ux,uy,w,xCoord,yCoord,z,fileName3,dataName3,dims);
		  writeCSV(uMag,N,N,fileNameCSV);
	  }
  }

  // be a good leader; free your memory
  cudaFree(fEven);
  cudaFree(fOdd);
  cudaFree(ndType_d);
  cudaFree(ux_d);
  cudaFree(uy_d);
  cudaFree(pressure_d);

  delete [] fDD;
  delete [] ndType;
  delete [] ux;
  delete [] uy;
  delete [] pressure;
  delete [] uMag;
  delete [] xCoord;
  delete [] yCoord;
  delete [] w;
  delete [] z;
  return 0;

}
