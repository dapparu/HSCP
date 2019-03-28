#Path
SRC_PATH	=./src
INT_PATH	=./interface
BIN_PATH	=./bin
LIB_PATH	=./lib
TEST_PATH	=./test
OBJ_PATH	=./obj


#Excutable 
EXE		=$(BIN_PATH)/main


#ROOT flags
ROOTCFLAGS	=$(shell root-config --cflags)
ROOTGLIBS	=$(shell root-config --glibs)

CXX		=g++

CFLAGS		=$(ROOTCFLAGS) $(ROOTGLIBS) -O2 -Wall -g

#Code lists
INT		= $(INT_PATH)/Estimator.h \
		  $(INT_PATH)/Builder.h \
		  $(INT_PATH)/Track.h \
		  $(INT_PATH)/Cluster.h \
		  $(INT_PATH)/SimHit.h \
		  $(INT_PATH)/ClusterStrip.h

SRC		= $(SRC_PATH)/Estimator.cc \
		  $(SRC_PATH)/Builder.cc \
		  $(SRC_PATH)/Track.cc \
		  $(SRC_PATH)/Cluster.cc \
		  $(SRC_PATH)/SimHit.cc \
		  $(SRC_PATH)/ClusterStrip.cc

OBJ		= $(OBJ_PATH)/Estimator.o \
		  $(OBJ_PATH)/Builder.o \
		  $(OBJ_PATH)/Track.o \
		  $(OBJ_PATH)/Cluster.o \
		  $(OBJ_PATH)/SimHit.o \
		  $(OBJ_PATH)/ClusterStrip.o \

LIB 		= $(LIB_PATH)/libHSCP.so

TEST		= $(TEST_PATH)/main.cpp



all: 		$(EXE) 

#Creation of the library 
$(LIB): 	$(SRC) $(INT)
	$(CXX) $(CFLAGS) $(SRC) -c -fPIC 
	mv *.o $(OBJ_PATH)/.
	$(CXX) $(CFLAGS) -shared $(OBJ) -o $(LIB)
	
#Creation of the executable
$(EXE):		$(LIB) $(TEST)
	$(CXX) $(CFLAGS) -c $(TEST) -o obj/main.o
	$(CXX) $(CFLAGS) obj/main.o -o $(EXE) -L $(LIB_PATH) -lHSCP

clean:
	rm -rf $(OBJ) $(EXE)
