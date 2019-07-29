all: minBiasPythia ZPythiaGenerator ZStop
ZPythiaGenerator: ZPythiaGenerator.cc
	g++ ZPythiaGenerator.cc $(PYTHIA_PATH)/lib/libpythia8.a -o ZPythiaGenerator  -I$(PYTHIA_PATH)/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,$(PYTHIA_PATH)/lib -ldl `root-config --libs --cflags` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -lAtlas
minBiasPythia: minBiasPythia.cc
	g++ minBiasPythia.cc $(PYTHIA_PATH)/lib/libpythia8.a -o minBiasPythia -I$(PYTHIA_PATH)/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,$(PYTHIA_PATH)/lib -ldl `root-config --libs --cflags` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -lAtlas

ZStop: ZStop.cc
	g++ ZStop.cc $(PYTHIA_PATH)/lib/libpythia8.a -o ZStop -I$(PYTHIA_PATH)/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,$(PYTHIA_PATH)/lib -ldl `root-config --libs --cflags` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -lAtlas
