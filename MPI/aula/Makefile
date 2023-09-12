CC = mpicc
CFLAGS = -Wall 
LDLIBS = -lm

SRC = ex1.c ex2.c Ex2_working.c
OBJ = $(SRC:.c=.o)
TARGET = ex1 ex2 Ex2_working
all: $(TARGET)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create executables
%: %.o
	$(CC) $< -o $@ $(LDLIBS)

clean:
	rm -f $(OBJS) $(PROGRAMS)