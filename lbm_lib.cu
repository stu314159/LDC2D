#define TPB 128

#include "lbm_lib.h"


using namespace std;
// define kernels here
__global__ void ldc2D_ts_kernel(float * fOut, const float * fIn, const float omega,
		const int * ndType, const float u_bc, const int N){
	int nID = threadIdx.x+blockIdx.x*blockDim.x; // my prototype code used nID and I don't want to change it.
	const int nnodes = N*N;
	if(nID<nnodes){
		float f0,f1,f2,f3,f4,f5,f6,f7,f8;
		// execute the time step.
		// get data from fIn
		f0 = fIn[nID];
		f1 = fIn[nID+1*nnodes];
		f2 = fIn[nID+2*nnodes];
		f3 = fIn[nID+3*nnodes];
		f4 = fIn[nID+4*nnodes];
		f5 = fIn[nID+5*nnodes];
		f6 = fIn[nID+6*nnodes];
		f7 = fIn[nID+7*nnodes];
		f8 = fIn[nID+8*nnodes];

		// find out what kind of node I am:
		int nT = ndType[nID]; // 0 = regular interior node; 1 = solid node; 2 = velocity node

		// compute macroscopic variables whether or not we need them
		float rho,ux,uy;
		rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
		ux = (f1-f3+f5-f6-f7+f8)/rho;
		uy = (f2-f4+f5+f6-f7-f8)/rho;

		// if I am a velocity node, I need to update the macroscopic variables and "push" the microscopic variables.
		if(nT==2){
			//get difference between current velocity and specified velocity
			float dx = u_bc - ux; float dy = -uy;
			// update density distribution functions
			//f(spd) = w(spd)*rho*(3*(ex(spd)*dx + ey(spd)*dy))
			f1 += (1./9.)*(rho)*(3.*dx);
			f2 += (1./9.)*(rho)*(3.*dy);
			f3 += (1./9.)*(rho)*(3.*-dx);
			f4 += (1./9.)*(rho)*(3.*-dy);
			f5 += (1./36.)*(rho)*(3.*(dx + dy));
			f6 += (1./36.)*(rho)*(3.*(-dx+dy));
			f7 += (1./36.)*(rho)*(3.*(-dx-dy));
			f8 += (1./36.)*(rho)*(3.*(dx-dy));

			// now set macroscopic velocities for fEq calculation
			ux = u_bc; uy = 0.;
		}

		// if I am not a solid node, I need to compute equilibrium and relax towards it;
		if(nT!=1){ // node types 0 and 2 will execute this.
			float fe0,fe1,fe2,fe3,fe4,fe5,fe6,fe7,fe8;
			float cu; // utility variable.

			// compute fEq
        	// e0 = (0,0)
        	// cu = 3.*(ex(spd)*ux + ey(spd)*uy)
			cu = 0.;
			fe0 = (4./9.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));

			//e1 = (1,0)
			cu = 3.*(ux);
			fe1 = (1./9.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));
			//e2 = (0,1)
			cu = 3.*(uy);
			fe2 = (1./9.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));
			//e3 = (-1,0);
			cu = 3.*(-ux);
			fe3 = (1./9.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));
			//e4 = (0,-1);
			cu = 3.*(-uy);
			fe4 = (1./9.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));

			//e5 = (1,1)
			cu = 3.*(ux+uy);
			fe5 = (1./36.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));

			//e6 = (-1,1)
			cu = 3.*(-ux+uy);
			fe6 = (1./36.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));
			//e7 = (-1,-1)
			cu = 3.*(-ux-uy);
			fe7 = (1./36.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));
			//e8 = (1,-1)
			cu = 3.*(ux-uy);
			fe8 = (1./36.)*rho*(1.+cu+0.5*cu*cu - (3./2.)*(ux*ux+uy*uy));

			// relax towards equilibrium
			f0 -= (f0 - fe0)*omega;
			f1 -= (f1 - fe1)*omega;
			f2 -= (f2 - fe2)*omega;
			f3 -= (f3 - fe3)*omega;
			f4 -= (f4 - fe4)*omega;
			f5 -= (f5 - fe5)*omega;
			f6 -= (f6 - fe6)*omega;
			f7 -= (f7 - fe7)*omega;
			f8 -= (f8 - fe8)*omega;

		}else{ // if I am node type 1: bounce back
			// bounce-back
			float tmp;
			// 1 -> 3, 3 -> 1
			tmp = f1; f1 = f3; f3 = tmp;
			// 2 -> 4, 4 -> 2
			tmp = f2; f2 = f4; f4 = tmp;
			// 5 -> 7, 7 -> 5
			tmp = f5; f5 = f7; f7 = tmp;
			// 6 -> 8, 8 -> 6
			tmp = f6; f6 = f8; f8 = tmp;

		}


		 // stream the result...clunky, but it works.
		    //compute the local stream vector...
		    int x;
		    int y;
		    int yn;
		    int ys;
		    int xe;
		    int xw;

		    //int dir;
		    int dof_num; //int f_num;
		    x = nID%N+1;
		    y = ((nID+1)-x+1)/N + 1;

		    yn = y%N+1;
		    xe = x%N+1;

		    if(y==1){
		      ys = N;
		    }else{
		      ys = y-1;
		    }
		    if(x==1){
		      xw=N;
		    }else{
		      xw=x-1;
		    }

		    dof_num = N*(y-1)+x;
		    fOut[dof_num-1]=f0;

		    dof_num=N*(y-1)+xe;
		    fOut[nnodes+dof_num-1]=f1;

		    dof_num=N*(yn-1)+x;
		    fOut[2*nnodes+dof_num-1]=f2;

		    dof_num=N*(y-1)+xw;
		    fOut[3*nnodes+dof_num-1]=f3;

		    dof_num=N*(ys-1)+x;
		    fOut[4*nnodes+dof_num-1]=f4;

		    dof_num=N*(yn-1)+xe;
		    fOut[5*nnodes+dof_num-1]=f5;

		    dof_num=N*(yn-1)+xw;
		    fOut[6*nnodes+dof_num-1]=f6;

		    dof_num=N*(ys-1)+xw;
		    fOut[7*nnodes+dof_num-1]=f7;

		    dof_num=N*(ys-1)+xe;
		    fOut[8*nnodes+dof_num-1]=f8;


	}


}

__global__ void D2Q9_getVelocityAndDensity(float * ux, float * uy, float * p,
		const float ucf, const float pcf, const float * f, const int N){
	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	const int nnodes=N*N;
	float rho = 0;

	float f0,f1,f2,f3,f4,f5,f6,f7,f8;
	if(tid<nnodes){
		f0 = f[tid];
		f1 = f[tid+nnodes];
		f2 = f[tid+2*nnodes];
		f3 = f[tid+3*nnodes];
		f4 = f[tid+4*nnodes];
		f5 = f[tid+5*nnodes];
		f6 = f[tid+6*nnodes];
		f7 = f[tid+7*nnodes];
		f8 = f[tid+8*nnodes];

		rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
		p[tid]=rho*pcf;
		ux[tid]=((f1-f3+f5-f6-f7+f8)/rho)*ucf;
		uy[tid]=((f2-f4+f5+f6-f7-f8)/rho)*ucf;

	}
}


// define host-side GPU interface functions here
void LDC2D_timestep(float * fOut, const float * fIn, const float omega,
		const int * ndType, const float u_bc, const int N){
	dim3 BLOCKS(TPB,1,1);
	dim3 GRIDS((N*N+TPB-1)/TPB,1,1);
	ldc2D_ts_kernel<<<GRIDS,BLOCKS>>>(fOut,fIn,omega,ndType,u_bc,N);


}

void LDC2D_getVelocityAndDensity(float * ux, float * uy, float * p,
		const float ucf, const float pcf,const float * f, const int N){

	dim3 BLOCKS(TPB,1,1);
	dim3 GRIDS((N*N+TPB-1)/TPB,1,1);
	D2Q9_getVelocityAndDensity<<<GRIDS,BLOCKS>>>(ux,uy,p,ucf,pcf,f,N);

}

void LDC2D_getRelativePressure(float * pressure,const int N){
	float p_ref = pressure[N*N/2];
	for(int nd=0;nd<N*N;nd++){
		pressure[nd]-=p_ref;
	}
}

void LDC2D_getVelocityMagnitude(float * uMag, const float * ux, const float * uy,const int N){
	float u,v;
	for(int nd=0;nd<N*N;nd++){
		u = ux[nd]; v = uy[nd];
		uMag[nd]=sqrtf(u*u+v*v);
	}

}

// define host-side LBM setup and utility functions here

void LDC2D_setup(const float Re, const int N,
		 const float omega, 
		 const float rho, const float nu,
		 float* fDD, int* ndType, float &u_lbm,
		 float &u_conv, float &t_conv, float &p_conv){

	// perform non-dimensionalization and scaling
	// set reference length, velocity and time
	float L_o = 1.0;
	float V_o = Re*nu/L_o;
	float T_o = L_o/V_o;

	// convert to dimensionless units
	//float T_d = 1.0;
	//float L_d = 1.0;

	float U_d = (T_o/L_o)*V_o;
	//float nu_d = 1./Re;

	//convert to LBM units
	float dx = 1./(N - 1.);
	float dt = (dx*dx)*Re*((1./3.)*((1./omega)-0.5));
	u_lbm = (dt/dx)*U_d; // this is also the lbm mach number

	//set conversion factors
	u_conv = (dx/dt)*(L_o/T_o); // multiply LBM velocity by u_conv to get physical velocity
	t_conv = (dt*T_o);
	float l_conv = dx*L_o;
	p_conv = ((l_conv/t_conv)*(l_conv/t_conv)*(1./3.));

	initializeF_toZero(N,rho,fDD);
	initializeNdType(N,ndType);



}

void initializeF_toZero(const int N, const float rho,
			float* fDD){
	const int nnodes = N*N;
	const int numSpd = 9;

	const float w[9] = {4./9.,1./9.,1./9.,1./9.,1./9.,
			1./36.,1./36.,1./36.,1./36.};

	for(int spd=0;spd<numSpd;spd++){
		for(int i=0;i<nnodes;i++){
			fDD[i+spd*nnodes]=rho*w[spd];
		}
	}
}

void initializeNdType(const int N, int* ndType){

	// node types: 	0 = regular internal node
	//             	1 = solid node
	//				2 = velocity boundary node

	//since I did not initialize the array, I will have to touch every entry
	const int nnodes = N*N;

	// this disgraceful display of hackery actually works.
	for(int nd = 0; nd<nnodes;nd++){
		if(nd<N){
			ndType[nd]=1; //solid on the bottom
		}else if(nd%N==0){
			ndType[nd]=1; //solid on the left side
		}else if(nd%N==(N-1)){
			ndType[nd]=1; // solid on the right side
		}else if(nd>=(nnodes-(N-1))){
			ndType[nd]=2; // velocity on the interior of the top
		}else{
			ndType[nd]=0;
		}
	}
}

// get the x and y coordinates of each lattice point.
void LDC2D_getGeometry(float * x, float * y, const float Lx, const float Ly,const int N){
	float dX = Lx/(N-1); float dY = Ly/(N-1);
	float xt = 0.; float yt = 0.;
	for(int i=0;i<N;i++){
		xt = 0;
		for(int j=0;j<N;j++){
			x[j+i*N]=xt;
			y[j+i*N]=yt;
			xt+=dX;
		}
		yt+=dY;
	}
}

// a function to create output csv data files.
// please forgive this hackery, but it should work.
void writeCSV(const float* data, const int Nrows, const int Nentries, const string fileName){

	// need to be able to map each entry of data to an integer between
	// 0 and 256. read all of the values to find the minimum and maximum value.
	float minVal=1e6;//so they will be overwritten immediately
	float maxVal =-1e6;
	float val;
	int scaledVal;
	const int scaleMin = 0;
	const int scaleMax = 256;
	for(int rw=0;rw<Nrows;rw++){
		for(int col=0;col<Nentries;col++){
			val = data[col+rw*Nentries];
			if (val > maxVal){
				maxVal = val;
			}
			if (val < minVal){
				minVal = val;
			}
		}
	}

	// potential problem if all of the values are the same -- such as
	// is the case initially.

	// now each scaled value:
	// sv = data[col+rw*Nentries]

	ofstream output_dat(fileName.c_str(),ios::out);
	if (maxVal != minVal){
		for(int rw=0;rw<Nrows;rw++){
			for(int col=0;col<Nentries;col++){
				val = data[col+rw*Nentries];
				scaledVal = (int)((val - minVal)/(maxVal-minVal)*(scaleMax - scaleMin)+scaleMin);
				output_dat << scaledVal;
				if(col < (Nentries-1)){
					output_dat << ", ";
				}
			}
			output_dat << endl;
		}
	} else { // if minVal == maxVal, set all outputs to zero.
		scaledVal = 0;
		for(int rw=0;rw<Nrows;rw++){
			for(int col=0;col<Nentries;col++){
				output_dat << scaledVal;
				if(col < (Nentries-1)){
					output_dat << ", ";
				}
			}
			output_dat << endl;
		}

	}

	output_dat.close();
}
