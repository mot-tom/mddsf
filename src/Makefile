CC		= gcc
CXX		= g++
CFLAGS		= -O4 -Wall -I/usr/local/include
DEST		= ~/bin
INCLUDE_1	= -I/opt/intel/mkl/include/
LDFLAGS_1	= -L/usr/local/lib -L/opt/intel/mkl/lib/intel64/
LIBS_1		= -lmkl_intel_lp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -fopenmp
OBJS_1		= main.o dsf.o dos.o dispersion.o dsf-omp.o dsf-lmp.o dsf-lmp-omp.o sed-lmp-omp.o dos-lmp.o dos-lmp-omp.o dmat.o dmat2.o hessian-lmp.o hist.o inv-lmp-omp.o xyz-vector.o
SOURCE_1	= main.cpp dsf.cpp dos.cpp dispersion.cpp dsf-omp.cpp dsf-lmp.cpp dsf-lmp-omp.cpp sed-lmp-omp.cpp dos-lmp.cpp dos-lmp-omp.cpp dmat.cpp dmat2.cpp hessian-lmp.cpp hist.cpp inv-lmp-omp.cpp xyz-vector.cpp
PROGRAM_1	= mddsf
INCLUDE_2	= -I/usr/include/python2.7/
LDFLAGS_2	= -L/usr/local/lib -L/usr/lib/python2.7/
LIBS_2		= -lpython2.7
OBJS_2		= density.o
SOURCE_2	= density.cpp
PROGRAM_2	= mddsf-g
OBJS_2_2	= density2.o
SOURCE_2_2	= density2.cpp
PROGRAM_2_2	= mddsf-g2
OBJS_3	= vib-vector.o
SOURCE_3	= vib-vector.cpp
PROGRAM_3	= mddsf-vib
INCLUDE_3	= -I/usr/include/python3.4m/
LDFLAGS_3	= -L/usr/local/lib -L/usr/lib/python3.4/
LIBS_3		= -lpython3.4

all:		$(PROGRAM_1) $(PROGRAM_2) $(PROGRAM_2_2)

$(PROGRAM_1):	$(OBJS_1)
		$(CXX) $(OBJS_1) $(INCLUDE_1) $(LDFLAGS_1) $(LIBS_1) -o $(PROGRAM_1)

$(PROGRAM_2):	$(OBJS_2)
		$(CXX) $(OBJS_2) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2) -o $(PROGRAM_2)
$(PROGRAM_2_2):	$(OBJS_2_2)
		$(CXX) $(OBJS_2_2) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2) -o $(PROGRAM_2_2)
$(PROGRAM_3):	$(OBJS_3)
		$(CXX) $(OBJS_3) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2) -o $(PROGRAM_3)

$(OBJS_1):	$(SOURCE_1)
		$(CXX) -c $(SOURCE_1) $(INCLUDE_1) $(LDFLAGS_1) $(LIBS_1)

$(OBJS_2):	$(SOURCE_2)
		$(CXX) -c $(SOURCE_2) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2)
$(OBJS_2_2):	$(SOURCE_2_2)
		$(CXX) -c $(SOURCE_2_2) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2)
$(OBJS_3):	$(SOURCE_3)
		$(CXX) -c $(SOURCE_3) $(INCLUDE_2) $(LDFLAGS_2) $(LIBS_2)

clean:;	rm -f *.o *~ $(PROGRAM_1) $(PROGRAM_2) $(PROGRAM_2_2) $(PROGRAM_3)

install:	$(PROGRAM_1) $(PROGRAM_2) $(PROGRAM_2_2) $(PROGRAM_3)
		install -s $(PROGRAM_1) $(DEST)
		install -s $(PROGRAM_2) $(DEST)
		install -s $(PROGRAM_2_2) $(DEST)
		install -s $(PROGRAM_3) $(DEST)

