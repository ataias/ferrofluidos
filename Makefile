CC=clang++
CFLAGS=-c -Wall -I /Users/ataias/opt/ -I ./include
LDFLAGS=
SOURCES=poisson.cc sparsePD.cc
OBJECTS=$(SOURCES:.cc=.o)
VPATH = src
BUILDDIR = build
EXECUTABLE=poisson

define cc-command
$(CC) $(CFLAGS) $< -o $@
endef

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILDDIR)/%.o: %.cc
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o