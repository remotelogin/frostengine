# Compiler
CC = g++

# Compiler flags
CFLAGS = -g

# X11 library flags
LDFLAGS = -lX11

# Source files
SRCS = main.cpp  # Add your source files here

# Output executable name
TARGET = app

# Default rule to build the application
all: $(TARGET)

# Rule to link the object files and create the executable
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

# Rule to clean up the output
clean:
	rm -f $(TARGET) *.o

# Rule to recompile the application
rebuild: clean all

.PHONY: all clean rebuild
