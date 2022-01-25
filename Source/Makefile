# -----------------------------------------------------------------

SHELL = /bin/sh

CC       = gcc
F90      = gfortran
LNKR     = gfortran
CPPFLAGS = -cpp
CFLAGS   = -fopenmp -DOPENMP
FFLAGS   = -fopenmp -DOPENMP -ffree-line-length-0
LDFLAGS  = -fopenmp
LIBS     = -lm

CFLAGS  += -O2 -Wall
FFLAGS  += -O2
LDFLAGS += -O2

# Debug options
#CFLAGS  += -g
#FFLAGS  += -g -debug -C -traceback #-check arg_temp_created -check bounds -fp-stack-check -warn alignments -warn interfaces
#LDFLAGS += -g -debug -C -traceback #-check arg_temp_created -check bounds -fp-stack-check -warn alignments -warn interfaces

# PDR code options
FFLAGS += -DDUST
FFLAGS += -DCALC_TDUST
FFLAGS += -DGUESS_TEMP
FFLAGS += -DPHOTO_CROSS_SECTIONS
#FFLAGS += -DUSE_ALI
#CFLAGS += -DLOG_ODES

# Modify the following two lines to point to the /include and /lib directories of your SUNDIALS installation
INCLUDES  = -I/home/jon/Documents/zupcx4/software/sundials/include
LIBRARIES = -L/home/jon/Documents/zupcx4/software/sundials/lib

SOURCE  = definitions healpix_types modules convert_number
SOURCE += main setup_logfiles count_lines count_substring swap_function spline
SOURCE += read_parameters read_particles read_field read_species read_reactions read_coolants read_cross_sections
SOURCE += evaluation_points dust_temperature guess_temperature column_densities line_opacities
SOURCE += reaction_rates h2_formation_rate photoreaction_rates collisional_rates heating_rates
SOURCE += update_abundances update_coolants update_linewidths ortho_para_ratio lte_populations
#SOURCE += update_coolants update_linewidths ortho_para_ratio lte_populations
SOURCE += level_populations transition_matrix escape_probability # lambda_operator
SOURCE += check_chemistry_convergence check_population_convergence fix_bad_cooling_rates
SOURCE += nearest_neighbours update_temperature write_output analyse_chemistry
SOURCE += photo_cross_sections xray_cross_sections xray_ionization_rates
SOURCE += odes jacobian calculate_abundances

OBJECTS = ${SOURCE:=.o}

LIBRARIES += -lsundials_cvode -lsundials_nvecserial -lm

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c .f90

.c.o :
	${CC} ${CFLAGS} ${INCLUDES} -c $<

.f90.o :
	${F90} ${CPPFLAGS} ${FFLAGS} -c $<

# -----------------------------------------------------------------------------------------

UCL-PDR: ${OBJECTS}
	${LNKR} -o UCL-PDR ${OBJECTS} ${LDFLAGS} ${LIBRARIES}
	mv UCL-PDR ../

install:
	mv UCL-PDR ../

clean:
	rm -f *.o
	rm -f *.mod
	rm -f UCL-PDR

cleandebug:
	rm -f *__genmod.f90
	rm -f *__genmod.mod

# -----------------------------------------------------------------------------------------
