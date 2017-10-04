CC= g++
SUFFIXES += .d

#FLAGS = -g -Wall  -MMD -pedantic
FLAGS = -Wall -MMD -std=c++11
ifdef debug
FLAGS += -g
else
FLAGS += -Ofast
endif


ifdef bersim
FLAGS += -DBERSIM=1
else
FLAGS += -DBERSIM=0   
endif


OBJ = main.o OSDDecoder.o math_funcs.o gf2_algebra.o ecc.o scc.o ldpc_class.o timer.o

all: main

# Add hs10.h dependency to all objects
# $(OBJ): 

main.o: main.c OSDDecoder.h math_funcs.h gf2_algebra.h ecc_c_headers.h ecc.h

OSDDecoder.o: OSDDecoder.c OSDDecoder.h math_funcs.h gf2_algebra.h ecc_c_headers.h ecc.h

math_funcs.o: math_funcs.c math_funcs.h ecc_c_headers.h ecc.h

gf2_algebra.o : gf2_algebra.c gf2_algebra.h math_funcs.h

ecc.o: ecc.c ecc.h scc.h

scc.o: scc.c ecc.h scc.h

ldpc_class.o: ldpc_class.cpp ldpc_class.h ldpc_simulator_data_def.h ecc.h

timer.o: timer.h timer.cpp

.c.o:
	$(CC) -c $(FLAGS) $< -o $@

.cpp.o:
	$(CC) -c $(FLAGS) $< -o $@	

main: $(OBJ)
	$(CC) -lpthread -o main $(OBJ) -lm


clean:
	rm -rf main *.o *.d

cleanprofile:
	rm *.da *.dyn
