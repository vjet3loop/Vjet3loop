# zgqq: build the library and the command-line driver
#   make                 libzgqq.a and the executable zgqq
#   make test            the test against the reference values (tests/reference.txt)
#   make BOOST=/path/to/boost/include CXX=g++
CXX ?= c++
BOOST ?= /opt/homebrew/include

# double-double arithmetic relies on exact rounding: no fast-math, no fused contractions
CXXFLAGS = -std=c++17 -O2 -ffp-contract=off -Iinclude -I$(BOOST) -DNDEBUG -DBOOST_UBLAS_NDEBUG

OBJS = src/diffeqs.o src/zgqq.o

all: zgqq

libzgqq.a: $(OBJS)
	ar rcs $@ $^

zgqq: src/run.o libzgqq.a
	$(CXX) -o $@ $^

tests/test: tests/test.o libzgqq.a
	$(CXX) -o $@ $^

test: tests/test
	cd tests && ./test

src/%.o: src/%.cpp include/*.hpp src/*.hpp
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

tests/%.o: tests/%.cpp include/*.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# the release: library, executable sources, input files, example, C++ test (no Python);
# without macOS metadata (no ._ files, no extended attributes)
VERSION = v1
DIST = zgqq-$(VERSION)
dist:
	rm -rf $(DIST) $(DIST).tar.gz
	mkdir -p $(DIST)/src $(DIST)/tests
	cp -R README.md Makefile COPYING include data examples $(DIST)/
	cp src/*.cpp src/*.hpp $(DIST)/src/
	cp tests/test.cpp tests/reference.txt $(DIST)/tests/
	COPYFILE_DISABLE=1 tar --no-xattrs -czf $(DIST).tar.gz $(DIST)
	rm -rf $(DIST)
	@echo "written $(DIST).tar.gz"

clean:
	rm -f src/*.o tests/*.o libzgqq.a zgqq tests/test

.PHONY: all test dist clean
