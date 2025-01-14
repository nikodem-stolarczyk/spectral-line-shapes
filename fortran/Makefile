F90 = gfortran  # compile with gfortran

MOD_DIR   = mod
SRC_DIR   = mHT
SRC_FILES = ${SRC_DIR}/cpf.o ${SRC_DIR}/profile.o

# List of usage examples:
ABSORPTION_OBJ = example_absorption.o
DISPERSION_OBJ = example_dispersion.o
CPF_OBJ = example_cpf.o
OPTIONAL_PARAMETERS_OBJ = example_mHT_optional_parameters.o
PROFILES_OBJ = example_profiles.o

ABSORPTION = $(SRC_FILES) $(ABSORPTION_OBJ)
DISPERSION = $(SRC_FILES) $(DISPERSION_OBJ)
CPF = $(SRC_FILES) $(CPF_OBJ)
OPTIONAL_PARAMETERS = $(SRC_FILES) $(OPTIONAL_PARAMETERS_OBJ)
PROFILES = $(SRC_FILES) $(PROFILES_OBJ)

CMPLFLG = -g -c -O3 -I${MOD_DIR}

.SUFFIXES: .f90 .o

.PHONY: all absorption dispersion cpf optional_parameters profiles clean remove_executables

all: absorption dispersion cpf optional_parameters profiles clean

absorption: $(ABSORPTION)
	$(F90) $(ABSORPTION) -o example_absorption.x

dispersion: $(DISPERSION)
	$(F90) $(DISPERSION) -o example_dispersion.x

cpf: $(CPF)
	$(F90) $(CPF) -o example_cpf.x

optional_parameters: $(OPTIONAL_PARAMETERS)
	$(F90) $(OPTIONAL_PARAMETERS) -o example_mHT_optional_parameters.x

profiles: $(PROFILES)
	$(F90) $(PROFILES) -o example_profiles.x

clean:
	rm -rf *.o ${MOD_DIR}/*.mod ${SRC_DIR}/*.o

remove_executables:
	rm -rf *.x

%.o: %.f90
	mkdir -p $(MOD_DIR)
	$(F90) $(CMPLFLG) -J$(MOD_DIR) -c -o $@ $<
