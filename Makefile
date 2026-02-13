# Compiler and flags
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Target executable
TARGET = orbit

# Source files
SRC = keplerian.cpp

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean rule
clean:
	rm -f $(TARGET)