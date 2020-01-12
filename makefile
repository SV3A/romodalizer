DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS = -std=c++11 -g -Wall
else
	CXXFLAGS = -std=c++11 -O3
endif

# Define program name, directories and objects:
TARGET = ./bin
SRCDIR = ./src
OBJDIR = ./objs
OBJ = main.o ModalAnalysis.o Mesh.o Util.o

VPATH = $(SRCDIR):$(OBJDIR)

# Define compiler and suffixes
CC = clang++

# Compile rule
.cc.o:
	@mkdir -p $(OBJDIR)
	$(CC) $(CXXFLAGS) -c -o $(OBJDIR)/$@ $<

# Link rule:
$(TARGET): $(OBJ)
	$(CC) $(CXXFLAGS) -o $(TARGET) $(addprefix $(OBJDIR)/, $(OBJ))
	chmod go+rx $(TARGET)

# Dependencies:
ModalAnalysis.o: Util.o
main.o: Mesh.o ModalAnalysis.o Util.o

.PHONY: clean

# Clean-up rule:
clean:
	rm -rf $(OBJDIR)
