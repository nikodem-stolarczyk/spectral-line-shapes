F90 = gfortran  # compile with gfortran
#
MOD_DIR   = mod
SRC_DIR   = src
SRC_FILES = ${SRC_DIR}/cpf.o ${SRC_DIR}/profile.o
PROGRAM   = usage_example.o
#
TOTAL = $(SRC_FILES) $(PROGRAM)
#
CMPLFLG = -g -c -O3 -fdefault-real-8 -I${MOD_DIR}
#
.SUFFIXES: .f90 .o
#
all : $(TOTAL)
	$(F90) $(TOTAL) -o test.x
	rm -r *.o ${MOD_DIR} ${SRC_DIR}/*.o
#
clean:
	rm -rf *.x *.o ${MOD_DIR} ${SRC_DIR}/*.o
#
%.o: %.f90
	mkdir -p $(MOD_DIR)
	$(F90) $(CMPLFLG) -J$(MOD_DIR) -c -o $@ $<
