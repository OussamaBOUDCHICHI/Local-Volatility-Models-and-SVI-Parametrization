#
# TODO : Create project structure
#

name := svi
SOURCES := $(shell find  . -maxdepth 1 -name "*.cpp")


all:

	g++ -std=c++17 -g $(SOURCES) -o $(name) -lnlopt -lm -lboost_iostreams -lboost_system -lboost_filesystem

clean:
	rm -R $(name)