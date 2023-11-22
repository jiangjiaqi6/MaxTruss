OBJ = main.o fileLoader.o graph.o log.o MyFile.o
TARGET = main
INCL = -Iinclude
CPP = g++
# -Dopt

CFLAGS := -std=c++17 -Iinclude -fopenmp -w -g -mpopcnt -O3 -DDegSort 

$(TARGET): $(OBJ)
	$(CPP) $(OBJ) $(CFLAGS) $(INCL) -o $(TARGET)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $< -o $@

.PHONY: clean  
clean:
	rm -rf $(OBJ) $(TARGET) 
	# rm -rf graphInfo/*

	