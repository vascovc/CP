CC = mpicc
CFLAGS = -Wall 
LDLIBS = -lm

SRC = original.c alinea_a.c alinea_b.c alinea_c.c alinea_d.c
OBJ = $(SRC:.c=.o)
TARGET = original alinea_a alinea_b alinea_c alinea_d
all: $(TARGET)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create executables
%: %.o
	$(CC) $< -o $@ $(LDLIBS)

clean:
	rm -f $(OBJS) $(PROGRAMS)