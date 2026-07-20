CXX ?= g++
CXXFLAGS ?= -std=c++17 -O0 -pthread -Wall -Wextra
CPPFLAGS ?= -Iinclude
TARGET := vira
SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

.PHONY: all clean release

all: $(TARGET)

release:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="-std=c++17 -O3 -pthread -Wall -Wextra"

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

src/%.o: src/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
