OBJ = main.o fileLoader.o graph.o log.o MyFile.o
TARGET = main
INCL = -Iinclude
CPP = g++


CFLAGS := -std=c++17 -Iinclude -w -g -mpopcnt -O3 -DDegSort  -DLazyUpdate -DMaintenance
# CFLAGS := -std=c++17 -Iinclude -w -g -mpopcnt -O3 -DDegSort -DsemiBinary


$(TARGET): $(OBJ)
	$(CPP) $(OBJ) $(CFLAGS) $(INCL) -o $(TARGET)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $< -o $@

.PHONY: clean  
clean:
	rm -rf $(OBJ) $(TARGET) 
	# rm -rf graphInfoCopy/*
	# rm -rf linearList/*
	# rm -rf kCoreInfo/*
	