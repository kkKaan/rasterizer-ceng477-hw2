CC = g++
CFLAGS = -O3
TARGET = rasterizer
SOURCES = *.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -rf $(TARGET)