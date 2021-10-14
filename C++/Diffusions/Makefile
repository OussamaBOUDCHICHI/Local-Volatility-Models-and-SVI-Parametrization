# 
# TODO : .o files to build
#	 app to ./bin


appname := BlSc

CXX := clang++
CXXFLAGS := -std=c++17

PATH := /usr/local/include
SRCDIR   := src
BUILDIR  := build
BIN      := bin

#INC  	 := -I include

SOURCES := $(shell find src . -maxdepth 1 -name "*.cpp")

#objects := $(patsubst $(SRCDIR)/%,$(BUILDIR)/%,$(SRCDIR:.cpp=.o))

objects  := $(patsubst $(SRCDIR)%.cpp, $(BUILDIR)%.o, $(SOURCES))

all: makedir $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BIN)/$(appname) $(objects)


$(BUILDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

makedir:
	mkdir -p $(BIN)
	mkdir -p $(BUILDIR)

clean:
	rm -f $(objects)

cleanbuild:
	rm -R $(BUILDIR)

move:
	mkdir -p $(PATH)/CRR
	mkdir -p $(PATH)/CRR/src
	cp -R include $(PATH)/CRR
	cp -R src/CRR.cpp $(PATH)/CRR/src
	
