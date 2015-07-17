CC= g++
# GM_DIR = /home/adam/Documents/Green-Marl/apps/output_cpp/gm_graph
LIB = $(GM_DIR)/lib
INC = $(GM_DIR)/inc
BIN = bin
SRC = src
OBJ = obj
CFLAGS = -std=c++11 -O3 -g -I$(INC) -fopenmp -Wall 
OUT = $(BIN)/graphConverter

# where to find source codes
vpath %.cc $(SRC)

all: $(OUT) Makefile


$(BIN)/graphConverter: $(OBJ)/graphConverter.o $(LIB)/libgmgraph.a $(OBJ)/mmio.o
	if [ ! -d bin ]; then mkdir bin; fi
	$(CC) $(CFLAGS) $(OBJ)/graphConverter.o $(OBJ)/mmio.o $(LFLAGS) -L$(LIB) -lgmgraph -o $@	

$(OBJ)/%.o: %.cc 
	if [ ! -d obj ]; then mkdir obj; fi
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)/*.o $(OUT) 
