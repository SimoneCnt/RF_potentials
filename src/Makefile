CXX = g++
CPPFLAGS = -O3 -g
STATIC_LIBS = -lcrypto -lz 


calc_energy : main_simple_energy.o energy.o aminoacid.o chain.o residue.o
	$(CXX) -o $@ $(LDFLAGS) main_simple_energy.o energy.o aminoacid.o chain.o residue.o $(STATIC_LIBS)

clean:
	rm -f *.o calc_energy
