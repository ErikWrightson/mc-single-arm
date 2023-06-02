########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
EXECFILE    := pbmodel
LIBFILE     := libF1F209.so

########################################################################
# obj file list
LIBOBJS     := F1F209.o F1F209Wrapper.o
OBJS        := $(LIBOBJS) Main.o
OBJDIR      := obj.$(ARCH)

########################################################################
# Compiler
AR          := ar
CC          := gcc
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
CFLAGS      := -Wall -fPIC -Wno-unused-variable -O3 $(MODE)
CXXFLAGS    := -Wall -fPIC -Wno-unused-variable -O3 $(MODE)
FFLAGS      := -Wall -fPIC -Wno-unused-variable -O3 $(MODE)
ifeq ($(MYOS),Darwin) 
#in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
               -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
LDFLAGS     := -O3 $(MODE)

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean test

########################################################################
all: exe lib

########################################################################
exe: $(OBJS)
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJS) $(SYSLIBS)
	@echo "Linking $(EXECFILE) ...... done!"

########################################################################	
lib: $(LIBOBJS)
	@$(LD) -shared $(LDFLAGS) -o $(LIBFILE) $(LIBOBJS) $(SYSLIBS)
	@echo "Linking $(LIBFILE) ...... done!"

########################################################################
%.o: %.c
	@echo Compiling $< ......
	@$(CC) -c $< -o $@  $(CFLAGS)

%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

########################################################################
clean:
	@rm -f *.o
	@rm -f *~ *# */*~ */*#
	@rm -f test *.so

test:
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo \\CFLAGS\:$(CFLAGS)
	@echo \\CXXFLAGS\:$(CXXFLAGS)
	@echo \\FFLAGS\:$(FFLAGS)
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo \\objs\: $(OBJS)

help: test

env: test
