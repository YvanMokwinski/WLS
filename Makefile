TEST_DIRECTORY=tests
APP_DIRECTORY=main

source_app=$(wildcard $(APP_DIRECTORY)/*.cpp)

source_test=$(wildcard $(TEST_DIRECTORY)/*.tests.cpp)
mpi_source_test=$(wildcard $(TEST_DIRECTORY)/*.test.mpi.cpp)

default:$(source_app:%.cpp=%.exe) \
	$(source_test:%.cpp=%.debug.exe) \
	$(source_test:%.cpp=%.exe) \
	$(source_test:%.cpp=%.valgrind.debug.exe) \
	$(source_test:%.cpp=%.valgrind.exe) \
	$(mpi_source_test:%.cpp=%.debug.exe) \
	$(mpi_source_test:%.cpp=%.exe) \
	$(mpi_source_test:%.cpp=%.valgrind.debug.exe) \
	$(mpi_source_test:%.cpp=%.valgrind.exe)

header=$(wildcard include/*.h)
LIBBLAS= -L/home/p901195/ifort/lib/intel64  -liomp5 -L/home/p901195/ifort/mkl/include/../lib/em64t  -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_intel_thread -lmkl_core    -lpthread #-lmkl_solver_ilp64
LIBBLAS=-fopenmp -lgomp -L/home/p901195/ifort/mkl/include/../lib/em64t  -lmkl_gf_ilp64 -lmkl_lapack95_ilp64  -lmkl_gnu_thread -lmkl_core    -lpthread \
-L/home/yvan/Simlab/libExtern/sparsekit2 -lskit
LIBBLAS=-lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core -lpthread
help:
	@echo "make test"	

compile: \
	$(source_test:%.cpp=%.debug.exe) \
	$(source_test:%.cpp=%.exe) \
	$(source_test:%.cpp=%.valgrind.debug.exe) \
	$(source_test:%.cpp=%.valgrind.exe) \
	$(mpi_source_test:%.cpp=%.debug.exe) \
	$(mpi_source_test:%.cpp=%.exe) \
	$(mpi_source_test:%.cpp=%.valgrind.debug.exe) \
	$(mpi_source_test:%.cpp=%.valgrind.exe)

HEADERS=-I./ -Iinclude -I../WCOMMON/include -I../
CFLAGS_DEBUG=-g -Wall -Werror 
CFLAGS= -O3 -DNDEBUG -Wall -Werror
CPP=g++ -std=c++11
MPI=mpic++ -DSIMLAB_MPI -std=c++11

# TO COMPILE
%.exe:%.cpp $(header)
	$(CPP) $(CFLAGS) $< -o $@ $(HEADERS) $(LIBBLAS)

%.debug.exe:%.cpp $(header)
	$(CPP) $(CFLAGS_DEBUG) $< -o $@ $(HEADERS)  $(LIBBLAS)

%.valgrind.exe:%.cpp $(header)
	$(CPP) $(CFLAGS) $< -o $@   -DVALGRIND_TESTING $(HEADERS)  $(LIBBLAS)

%.valgrind.debug.exe:%.cpp $(header)
	$(CPP) $(CFLAGS_DEBUG) $< -o $@  -DVALGRIND_TESTING $(HEADERS)  $(LIBBLAS)


%.mpi.exe:%.mpi.cpp $(header)
	$(MPI) $(CFLAGS) $< -o $@ $(HEADERS) $(LIBBLAS)

%.mpi.debug.exe:%.mpi.cpp $(header)
	$(MPI) $(CFLAGS_DEBUG) $< -o $@ $(HEADERS)  $(LIBBLAS)

%.mpi.valgrind.exe:%.mpi.cpp $(header)
	$(MPI) $(CFLAGS) $< -o $@   -DVALGRIND_TESTING $(HEADERS)  $(LIBBLAS)

%.mpi.valgrind.debug.exe:%.mpi.cpp $(header)
	$(MPI) $(CFLAGS_DEBUG) $< -o $@  -DVALGRIND_TESTING $(HEADERS)  $(LIBBLAS)



# TO TEST
VALGRIND=valgrind --leak-check=full --error-exitcode=1

%.valgrind.debug.exe.log:%.valgrind.debug.exe
	@echo "Testing Valgrind "$<
	@$(VALGRIND) $< > $@ 2>&1

%.valgrind.exe.log:%.valgrind.exe
	@echo "Testing Valgrind "$<
	@$(VALGRIND) $< > $@ 2>&1


%.debug.exe.log:%.debug.exe
	@echo "Testing "$<
	@$< > $@ 2>&1

%.exe.log:%.exe
	@echo "Testing "$<
	@$< > $@ 2>&1

test:
	@echo "Testing ..."
	@+make perform_test
	@echo "Tests passed."

perform_test:$(source_test:%.cpp=%.exe.log) $(source_test:%.cpp=%.debug.exe.log) \
	$(source_test:%.cpp=%.valgrind.exe.log) $(source_test:%.cpp=%.valgrind.debug.exe.log)
#	@\rm tests/*.log


clean:
	\rm -f tests/*.exe tests/*.log tests/*~ include/*~ src/*~

