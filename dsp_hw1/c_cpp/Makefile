.PHONY: all clean

CFLAGS+=-std=c99 -O2
LDFLAGS+=-lm     # link to math library
TARGET:train test check
all: $(TARGET)
clean:
	rm -rf train test check model_*   # type make clean to remove the compiled file
