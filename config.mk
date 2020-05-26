PLATFORM=Linux
CC=gcc
CPP=g++

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=0   
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=0 
LIBS=-L/usr/lib/openblas/lib -lopenbla

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=0    -fno-stack-protector
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=0  -fno-stack-protector
LIBS=-L/usr/lib/ -llapack -lblas -lgfortran


CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=1     -I/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=1   -I/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/
LIBS= -L/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64  -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core  -lpthread

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=1     -DWLS_WITH_MKL
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWLS_ILP64=1   -DWLS_WITH_MKL
LIBS=  -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core  -lpthread -lm

CFLAGS+=-std=c99
CPPFLAGS+=-std=c++11


