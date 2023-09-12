CC = mpicc
CFLAGS = -Wall 

SRC = exer1.c
OBJ = $(SRC:.c=.o)
TARGET = exer1

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)