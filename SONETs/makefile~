# INCLUDES = -I/usr/nld/netcdf-4.1.3/include
CC = g++
C_OPTIMIZE_SWITCH = -O2 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
LIBS = -L/opt/local/lib -lgsl -lgslcblas -lnetcdf_c++ -lnetcdf

CFLAGS = -Wall ${C_OPTIMIZE_SWITCH} ${INCLUDES}

run_secorder: run_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o outfile.o
	${CC} run_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o outfile.o -o run_secorder ${LIBS}

run_secorder_gephi: run_secorder_gephi.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o
	${CC} run_secorder_gephi.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o outfile.o -o run_secorder_gephi ${LIBS}

gen_mult_secorder: gen_mult_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o
	${CC} gen_mult_secorder.o secorder_rec_1p.o calc_sqrtcov_rec_1p.o calc_rhos.o calc_stats_1p.o outfile.o -o gen_mult_secorder ${LIBS}

run_secorder_2p: run_secorder_2p.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o
	${CC} run_secorder_2p.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o outfile.o -o run_secorder_2p ${LIBS}

run_ei_balanced: run_ei_balanced.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o
	${CC} run_ei_balanced.o secorder_rec_2p.o calc_sqrtcov_rec_2p.o calc_rhos.o calc_stats_2p.o outfile.o -o run_ei_balanced ${LIBS}

run_secorder_gen: run_secorder_gen.o secorder_gen.o calc_sqrtcov_gen.o calc_rhos.o calc_stats_2p.o calc_stats_1p.o
	${CC} run_secorder_gen.o secorder_gen.o calc_sqrtcov_gen.o calc_rhos.o calc_stats_2p.o calc_stats_1p.o outfile.o -o run_secorder_gen ${LIBS}

testsoinput: testsoinput.o calc_rhos.o secorder_input.o
	${CC} testsoinput.o calc_rhos.o secorder_input.o outfile.o -o testsoinput ${LIBS}

testsoinput.o: secorder_input.hpp calc_rhos.o
secorder_input.o: secorder_input.hpp calc_rhos.o

run_secorder.o: secorder_rec_1p.hpp calc_stats_1p.hpp outfile.hpp
run_secorder_gephi.o: secorder_rec_1p.hpp calc_stats_1p.hpp
gen_mult_secorder.o: secorder_rec_1p.hpp calc_stats_1p.hpp
run_secorder_2p.o: secorder_rec_2p.hpp calc_stats_2p.hpp
run_ei_balanced.o: secorder_rec_2p.hpp calc_stats_2p.hpp
run_secorder_gen.o: secorder_gen.hpp calc_stats_2p.hpp
secorder_rec_1p.o: secorder_rec_1p.hpp calc_sqrtcov_rec_1p.hpp calc_rhos.hpp
secorder_rec_2p.o: secorder_rec_2p.hpp calc_sqrtcov_rec_2p.hpp calc_rhos.hpp calc_stats_2p.hpp
secorder_gen.o: secorder_gen.hpp calc_sqrtcov_gen.hpp calc_rhos.hpp
calc_sqrtcov_rec_1p.o: calc_sqrtcov_rec_1p.hpp
calc_sqrtcov_rec_2p.o: calc_sqrtcov_rec_2p.hpp
calc_sqrtcov_gen.o: calc_sqrtcov_gen.hpp
calc_rhos.o: calc_rhos.hpp
calc_stats_1p.o: calc_stats_1p.hpp
calc_stats_2p.o: calc_stats_2p.hpp
outfile.o: outfile.hpp


%.o : %.cpp
	${CC} -c ${CFLAGS} $<

clean:
	\rm -f *.o
