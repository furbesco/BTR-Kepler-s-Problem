# Compiler and flags
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Target executable
TARGET = orbit

# Source files
#SRC = keplerian.cpp
SRC = quasi-kepler_1PN.cpp
#SRC = quasi-kepler_2PN.cpp

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean rule
clean:
	rm -f $(TARGET)