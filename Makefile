all: minBiasPythia.cc ZPythiaGenerator.cc
ZPythiaGenerator: ZPythiaGenerator.cc
	g++ ZPythiaGenerator.cc $(PYTHIA8)/lib/libpythia8.a -o ZPythiaGenerator  -I$(PYTHIA8)/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,$(PYTHIA8)/lib -ldl `root-config --libs --cflags`
minBiasPythia: minBiasPythia.cc
	g++ minBiasPythia.cc $(PYTHIA8)/lib/libpythia8.a -o minBiasPythia  -I$(PYTHIA8)/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,$(PYTHIA8)/lib -ldl `root-config --libs --cflags`
