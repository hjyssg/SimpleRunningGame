#Junyang
#makefile for cs1566 final project
#this makefile only works for mac
#too busy to check other platforms

IDIR =. 

#compiler
CC=g++

#compile flga
CFLAGS= -I$(IDIR) -c -Wall

#opengl tag for mac compiling
LIB = -framework OpenGl -framework GLUT 

LDFLAGS=
SOURCES=RunningGameController.cpp RunningGameEngine.cpp TGALoader.cpp tiny_obj_loader.cc HJY_GL.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=JungleRun

all: $(SOURCES) $(EXECUTABLE) clean
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIB)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

##make "all" depend on clean
##rm .o files
.PHONY: clean

clean:
		rm -f *.o