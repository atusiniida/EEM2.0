-include ../local_settings.mk

TARGET = libutility.a

CXXSRCS = \
	Date.cpp \
	Dist.cpp \
	HypergeometricDistributionImpl.cpp \
	MyException.cpp \
	MyFunc.cpp \
	MyMat.cpp \
	StopWatch.cpp \

INST_HDRS = \


NOINST_HDRS = \
	Date.h \
	Dist.h \
	HypergeometricDistributionImpl.h \
	MyException.h \
	MyFunc.h \
	MyMat.h \
	StopWatch.h \


SRCS = $(CXXSRCS)

.SUFFIXES: .o .cpp
CXXOBJS = $(CXXSRCS:.cpp=.o)
OBJS = $(CXXOBJS)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#install:
#	mkdir -p $(INCDIR)
#	cp $(INSTALL_HDRS) $(INCDIR)
#	mkdir -p $(LIBDIR)
#	cp $(TARGET) $(LIBDIR)

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	$(CXX) -MM -MG $(SRCS) $(CXXFLAGS) > Makefile.depend
	cat Makefile.depend

-include Makefile.depend
