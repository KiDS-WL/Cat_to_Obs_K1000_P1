# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# define any directories containing header files other than /usr/include
# We're expecting the bjutils to be in this directory
whereami=`pwd`
INCLUDES = -I/$(whereami)/bjutils/include

# the build target executable:
TARGET = xi2bandpow

all: $(TARGET)

$(TARGET): $(TARGET).c
      $(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

clean:
      $(RM) $(TARGET)
