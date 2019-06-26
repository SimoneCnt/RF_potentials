CXX = g++
CPPFLAGS = -O3
LDFLAGS =  -static
STATIC_LIBS = /usr/lib/libcrypto.a /usr/lib/libz.a


calc_energy : main_simple_energy.o energy.o aminoacid.o chain.o residue.o
	$(CXX) -o $@ $(LDFLAGS) main_simple_energy.o energy.o aminoacid.o chain.o residue.o $(STATIC_LIBS)
