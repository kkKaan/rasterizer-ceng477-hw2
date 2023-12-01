CC = g++
CFLAGS = -O3
TARGET = rasterizer
DIR = source/src
SOURCES = $(DIR)/*.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -rf $(TARGET)