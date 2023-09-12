CC = mpicc
CFLAGS = -Wall 
LDFLAGS =

SRC = Ex1Prog1.c Ex1Prog2.c Ex1.c
OBJ = $(SRC:.c=.o)
TARGET = Ex1Prog1 Ex1Prog2 Ex1

all: $(TARGET)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create executables
%: %.o
	$(CC) $(LDFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(PROGRAMS)