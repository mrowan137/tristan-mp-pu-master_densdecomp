
#
# Tristan makefile
#

SHELL=/bin/sh

cc= gcc
FC= h5pfc
LD= h5pfc 

#this is the user file used to initialize the problem and provide user boundary conditions
#change the name of the user file for your problem; it will look fore $USER.F90 file. 
# 
USER_FILE ?= user/user_weibel
USR_DIR = $(dir $(USER_FILE))

# Executable name
EXE_FILE = tristan-mp2d

#Performance options
PERFORMANCE= -O3
#-xhost -qopt-report -g -ipo

#Precompiler options for enabling different algorithms
CUSTOM0= -DMPI -DHDF5 -DserIO -Ddd3 -unroll
#CUSTOM0= $(CUSTOM)


ifdef 3D
 EXE_FILE = tristan-mp3d
 CUSTOM = $(CUSTOM0)
else
 EXE_FILE = tristan-mp2d
 CUSTOM = $(CUSTOM0) -DtwoD
endif

# -DMPI use MPI 
# -DHDF5 use HDF5
# -DserIO use serial IO, all processors send chunks to rank 0 which saves into file
#         omitting this option will use parallel IO; performance may vary
# -DtwoD 2D code (ignorable z coordinate), omitting this defaults to 3D
# -DDEBUG when enabled produces a lot of debugging information

# undocumented options: do not use
# -Dfilter2  
# -DABSORB 
# -DLOCALSCRATCH 
# -DLOCALRESTART
# -Dzzag use zigzag current deposit
# -Ddd1 use 1st order density decomposition
# -Ddd2 use 2nd order density decomposition
# -Ddd3 use 3rd order density decomposition

# flags normaly used for C and Fortran compilers
CFLAGS= $(PERFORMANCE)
FFLAGS= $(CUSTOM) $(PERFORMANCE)

OBJ_DIR=obj/
SRC_DIR=code/
BACKUP_DIR=backup/
EXE_DIR=exec/

# Executable name
EXECUTABLE = $(EXE_DIR)$(EXE_FILE)

#INC_PATH=-I$(OBJ_DIR)

# Objects
OBJS= $(OBJ_DIR)system.o $(OBJ_DIR)systemf.o $(OBJ_DIR)par.o $(OBJ_DIR)inputparser.o $(OBJ_DIR)inputparserf.o $(OBJ_DIR)fparser.o $(OBJ_DIR)globaldata.o $(OBJ_DIR)aux.o $(OBJ_DIR)communications.o $(OBJ_DIR)fields.o $(OBJ_DIR)fieldboundaries.o $(OBJ_DIR)particles.o $(OBJ_DIR)domain.o $(OBJ_DIR)dynamic_domain.o $(USER_FILE).o $(OBJ_DIR)output.o $(OBJ_DIR)restart.o $(OBJ_DIR)particles_movedeposit.o $(OBJ_DIR)overload.o  $(OBJ_DIR)initialize.o $(OBJ_DIR)tristanmainloop.o $(OBJ_DIR)tristan.o

###########################
print-%  : ; @echo $* = $($*)

all: dirs $(EXECUTABLE)

dirs : $(OBJ_DIR) $(BACKUP_DIR) $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BACKUP_DIR):
	mkdir -p $(BACKUP_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR)%.o:$(SRC_DIR)%.c
	$(cc) $(CFLAGS) $(INC_PATH) -c $< -o $@

$(USR_DIR)%.o:$(USR_DIR)%.F90
	$(FC) $(FFLAGS) $(INC_PATH) -module $(OBJ_DIR) -c $< -o $@

$(OBJ_DIR)%.o:$(SRC_DIR)%.F90
	$(FC) $(FFLAGS) $(INC_PATH) -module $(OBJ_DIR) -c $< -o $@


$(EXECUTABLE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

	tar --exclude=$@ --exclude=$(OBJ_DIR) --exclude=$(OBJ_DIR)* --exclude='*tar' --exclude='*.o' --exclude='*.mod' --exclude=$(BACKUP_DIR) --exclude=$(BACKUP_DIR)* --exclude='.git/*' --exclude='.git' -cf $(BACKUP_DIR)code.tar --exclude=$(EXE_DIR) --exclude=$(EXE_DIR)* .

clean: 
	rm -f $(OBJS)
	rm -f $(OBJ_DIR)*.mod
	find . -name "*.o" -delete
	rm -f $(EXE_DIR)*
#	rm -f $(EXECUTABLE)

