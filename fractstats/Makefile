TARGET = lbound
OBJS += lib/parse_trace.o
OBJS += lbound_fast.o

CXX = g++
CXXFLAGS += -std=c++11 -fopenmp #‑mveclibabi=svml #-stdlib=libc++ #non-linux
CXXFLAGS += -MMD -MP
CXXFLAGS += -I..
CXXFLAGS += -O3 -ffast-math -march=native -funroll-loops -ffast-math
CXXFLAGS += -Wall -Werror # -ftree-vectorize -fopt-info-vec-missed -fopt-info-vec-optimized
LDFLAGS += $(LIBS)

all:		$(TARGET)

debug: CXXFLAGS += -ggdb  -D_GLIBCXX_DEBUG # debug flags
debug: $(TARGET)

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

DEPS = $(OBJS:%.o=%.d)
-include $(DEPS)

clean:
	-rm $(TARGET) $(OBJS) $(DEPS)
