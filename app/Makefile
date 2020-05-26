VPATH=$(DISTRIB)/app

INCLUDES=-I$(DISTRIB)/include -I$(DISTRIB)/src/include
cpp_source_lib=$(wildcard $(VPATH)/*.cpp)
CPP_OBJ=$(notdir $(cpp_source_lib:%.cpp=%.cpp.o))
CPP_EXE=$(notdir $(cpp_source_lib:%.cpp=%.cpp.exe))

include $(DISTRIB)/config.mk

install: $(CPP_OBJ) $(CPP_EXE)

%.cpp.exe:%.cpp.o
	$(CPP) $(CPPFLAGS) $(FLAGS) $< -o $@  $(LIBS)
	ln -fs $(DISTRIB)/$(PLATFORM)-$(CC)$(DEBUG)/$@ $(DISTRIB)/bin/$*$(DEBUG)

%.cpp.o:%.cpp
	$(CPP) $(CPPFLAGS) $(FLAGS) -c  $< -o $@ $(INCLUDES)

