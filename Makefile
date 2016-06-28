CC=CC
CC_FLAGS=
CC_LIBS=-lcudart -lstdc++
CC_CUL=
CC_CUI=

CU_CC=nvcc
CU_FLAGS=-O3 -arch compute_35    
CU_DEFS=

SOURCE=ldc2D.cxx
TARGET=ldc2D

all: $(TARGET)

$(TARGET): $(SOURCE) vtk_lib.o lbm_lib.o
	$(CU_CC) -o $(TARGET) $(SOURCE) $(CC_FLAGS) -L$(CC_CUL) -I$(CC_CUI) $(CC_LIBS) vtk_lib.o lbm_lib.o


lbm_lib.o: lbm_lib.cu
	$(CU_CC) -c lbm_lib.cu $(CU_FLAGS) $(CU_DEFS)

vtk_lib.o: vtk_lib.cxx
	$(CC) -c vtk_lib.cxx $(CC_FLAGS)

clean:
	rm *.o $(TARGET) 
