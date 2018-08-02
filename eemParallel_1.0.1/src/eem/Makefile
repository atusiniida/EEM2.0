-include ../local_settings.mk

CXXFLAGS += -I../utility
LDFLAGS += -L../utility
LIBS += -lutility -lm

TARGET = eemParallel

CXXSRCS = \
	AbstractEEMsearch.cpp \
	AbstractGeneSetAnalysis.cpp \
	CoherenceBasedEEM.cpp \
	CoherenceBasedEEMsearch.cpp \
	DistConverter.cpp \
	EEM.cpp \
	EEMsearch.cpp \
	ExpressionModule.cpp \
	ExpressionModuleSet.cpp \
	GeneSetAnalysis.cpp \
	Main.cpp \


INST_HDRS = \


NOINST_HDRS = \
	AbstractEEMsearch.h \
	AbstractGeneSetAnalysis.h \
	CoherenceBasedEEM.h \
	CoherenceBasedEEMsearch.h \
	DistConverter.h \
	EEM.h \
	EEMsearch.h \
	ExpressionModule.h \
	ExpressionModuleSet.h \
	GeneSetAnalysis.h \


SRCS = $(CXXSRCS)

.SUFFIXES: .o .cpp
CXXOBJS = $(CXXSRCS:.cpp=.o)
OBJS = $(CXXOBJS)

all: $(TARGET)

$(TARGET): $(OBJS) ../utility/libutility.a
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

install:
	mkdir -p $(BINDIR)
	cp $(TARGET) $(BINDIR)

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	$(CXX) -MM -MG $(SRCS) $(CXXFLAGS) > Makefile.depend
	cat Makefile.depend

-include Makefile.depend
