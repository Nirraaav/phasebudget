# Builds the native dynamic program into the package directory, where
# `phasebudget._native` looks for it. -ffp-contract=off keeps the arithmetic
# identical to the Python pass; do not add -ffast-math.

CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -fPIC -ffp-contract=off -Wall -Wextra
TARGET := ../src/phasebudget/_pbdp.so

$(TARGET): dp.cpp
	$(CXX) $(CXXFLAGS) -shared -o $@ $<

clean:
	rm -f $(TARGET)

.PHONY: clean
