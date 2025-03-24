# Compiler and flags
CC = gcc
CFLAGS = -Wall -g -std=c99 -lm

# Target output
TARGET = Transmitter

# Source files
SRCS = data.c Transmitter.c short_preamble.c long_preamble.c data_payload.c correlation.c pkt_selection.c

# Object files
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Create the executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(CFLAGS)

# Compile source files to object files
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Clean object and executable files
clean:
	rm -f $(OBJS) $(TARGET)

# Rebuild everything
rebuild: clean all
