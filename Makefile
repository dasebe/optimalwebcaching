#TARGET = offline
OBJS += offline.o

CXX = g++ #clang++ #OSX
CXXFLAGS += -std=c++11 #-stdlib=libc++ #non-linux
CXXFLAGS += -I.
CXXFLAGS += -Wall -Werror 
CXXFLAGS += -O3 # release flags
LDFLAGS += $(LIBS)

TARGET_CS = capacityscaling
TARGET_NS = networksimplex

$(TARGET_CS): TARGET = $(TARGET_CS)
$(TARGET_CS): CXXFLAGS += -MMD -MP # dependency tracking flags
$(TARGET_CS): $(OBJS)	
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

debug: CXXFLAGS += -ggdb  -D_GLIBCXX_DEBUG # debug flags
debug:	$(TARGET_CS)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

DEPS = $(OBJS:%.o=%.d)
-include $(DEPS)

clean:
	-rm $(TARGET_CS) $(OBJS) $(DEPS) $(TARGET_NS)

$(TARGET_NS): TARGET = $(TARGET_NS)
$(TARGET_NS): CXXFLAGS += -DNETWORKSIMPLEX
$(TARGET_NS): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
