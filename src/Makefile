CC=gcc
CFLAGS= -ffast-math -std=c11 -fopenmp

VLASOV_1D_OBJ = vlasov_1d.o setup_IC.o io.o advection.o mp5.o poisson.o
VLASOV_1D_DEP := $(VLASOV_1D_OBJ)

CFLAGS += -D__PP__    # switch for the positivity limiter

# switch for the Poisson solver (choose one from three options below)
#CFLAGS += -D__POISSON_2ND__
#CFLAGS += -D__POISSON_4TH__
CFLAGS += -D__POISSON_6TH__

ALL = vlasov_1d

all: $(ALL)

vlasov_1d: $(VLASOV_1D_DEP)
	$(CC) $(CFLAGS) -o $@ $(VLASOV_1D_OBJ) -lm

clean:
	-rm -rf *.o
	-rm -rf *~

distclean: clean
	-rm -rf $(ALL)
