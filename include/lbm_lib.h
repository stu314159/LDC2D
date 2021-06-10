#ifndef LBM_LIB_H
#define LBM_LIB_H

#include <tgmath.h>
#include <iostream>
#include <string>
#include <fstream>



void LDC2D_setup(const float Re, const int N,
		 const float omega, 
		 const float rho, const float nu,
		 float* fDD, int* ndType, float &u,
		 float &u_conv, float &t_conv, float &p_conv);



void initializeF_toZero(const int N, const float rho,
			float* fDD);

void initializeNdType(const int N, int* ndType);

void LDC2D_timestep(float * fOut, const float * fIn, const float omega,
		const int * ndType, const float u_bc, const int N);

void LDC2D_getVelocityAndDensity(float * ux, float * uy, float * p,
		const float ucf, const float pcf,const float * f, const int N);

void LDC2D_getGeometry(float * x, float * y, const float Lx, const float Ly,const int N);

void LDC2D_getVelocityMagnitude(float * uMag, const float * ux, const float * uy,const int N);

void LDC2D_getRelativePressure(float * pressure,const int N);

void writeCSV(const float* data, const int Nrows, const int Nentries, const std::string fileName);

#endif
