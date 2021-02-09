# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2021 Paulo Rodrigues.

# @makefile

export CXX:=g++-mp-10
export CXXFLAGS:=-std=c++20 -Wfatal-errors -I../include

vpath %.hh ../include/gyronimo
vpath %.cc ../src:../apps
lib_srcs = $(shell basename ../src/*.cc)
lib_objs = $(lib_srcs:.cc=.o)
lib_deps = $(lib_srcs:.cc=.d)

all: lib dox
lib: $(lib_objs) $(lib_deps)
	ar cr libgyronimo.a $(lib_objs)
	ranlib libgyronimo.a
dox: config.dox
	doxygen config.dox
%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
sinclude $(lib_deps)
.PHONY: clean
clean:
	rm -r $(lib_objs) $(lib_deps)
	[[ -f libgyronimo.a ]] && rm libgyronimo.a
	[[ -d gyronimo_dox ]] && rm -r gyronimo_dox
