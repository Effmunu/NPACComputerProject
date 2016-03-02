# Compile this file (which is in this directory)
PROG      = main
FILENAME  = AnaCalib MappingTool LineShape FitterBase FitterBetaAlpha
LANG = en_US

# File suffixes
SRCSUFFIX  = cpp
SRCSUFFIX2 = C
OBJSUFFIX  = o
DEPSUFFIX  = d

# File directories
SRCDIR     = src
OBJDIR     = .o
DEPDIR     = .d
EXEDIR     = bin

# make sure the directories are there
__dummy := $(shell mkdir -p $(OBJDIR) $(DEPDIR) $(EXEDIR))

# ROOT
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)  
ROOTGLIBS  = $(shell root-config --glibs)
ROOTLDFLAGS= $(shell root-config --ldflags)

# C++ compiler
CXX        = g++ 
LD         = g++
OPTFLAGS = -g -O2
CXXFLAGS   = $(OPTFLAGS) $(ROOTCFLAGS)
LDFLAGS    = $(OPTFLAGS) $(ROOTLDFLAGS) -lMinuit
LIBS       = $(ROOTLIBS) 
SHLIBFLAGS = -shared
OUTPUTOPT  = -o

INCLUDES   = -I. -I$(SRCDIR) $(ROOTCFLAGS)
CXXFLAGS  += $(INCLUDES)
LIBS       = $(ROOTLIBS) 
THEOBJ     = $(addprefix ${OBJDIR}/, $(addsuffix .$(OBJSUFFIX), $(FILENAME)))
THEPROGOBJ = $(addprefix ${OBJDIR}/, $(addsuffix .$(OBJSUFFIX), $(PROG)))
THEEXE     = $(addprefix ${EXEDIR}/, $(PROG))
DEPS       = $(addprefix ${DEPDIR}/, $(addsuffix .$(DEPSUFFIX), $(FILENAME) $(PROG)))


all:	$(THEEXE)
	@echo "in all"

$(THEEXE):   $(THEOBJ) $(THEPROGOBJ)
	@echo "Linking $^ to $@"
	$(LD) $(OUTPUTOPT) $@ $^ $(LDFLAGS) $(LIBS)
	@echo "Successfully completed $@"

clean:
	rm -f $(DEPS) $(THEOBJ) $(THEEXE) $(THEPROGOBJ)
	rm -f .o/*


${SRCDIR}/RootTupleDict.cpp: ${SRCDIR}/RootTuple.h ${SRCDIR}/RootTupleLinkDef.h
	rootcint -f ${SRCDIR}/RootTupleDict.cpp -c ${SRCDIR}/RootTuple.h ${SRCDIR}/RootTupleLinkDef.h

${OBJDIR}/%.$(OBJSUFFIX) : ${SRCDIR}/%.$(SRCSUFFIX) ${DEPDIR}/%.$(DEPSUFFIX) 
	@echo "Compiling $< to $@"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

${DEPDIR}/%.$(DEPSUFFIX): ${SRCDIR}/%.$(SRCSUFFIX)
	@echo "Making $@"
	@$(CXX) -MM $(CXXFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|'$(OBJDIR)'/\1.o '$(DEPDIR)'/\1.d : |g' \
	> $@

# ${OBJDIR}/%.$(OBJSUFFIX) : ${DEPDIR}/%.$(DEPSUFFIX) 
# 	@echo "Compiling to create $@"
# 	@$(CXX) $(CXXFLAGS) -c $(wildcard $(addprefix ${SRCDIR}/,$(subst .$(OBJSUFFIX),.$(SRCSUFFIX),$(notdir $@))) $(addprefix ${SRCDIR}/,$(subst .$(OBJSUFFIX),.$(SRCSUFFIX2),$(notdir $@)))) -o $@


# ${DEPDIR}/%.$(DEPSUFFIX):
# 	@echo "Making $@"
# 	@$(CXX) -MM $(CXXFLAGS) $(wildcard $(addprefix ${SRCDIR}/,$(subst .$(DEPSUFFIX),.$(SRCSUFFIX),$(notdir $@))) $(addprefix ${SRCDIR}/,$(subst .$(DEPSUFFIX),.$(SRCSUFFIX2),$(notdir $@))) ) \
# 	| sed 's|\($*\)\.o[ :]*|'$(OBJDIR)'/\1.o '$(DEPDIR)'/\1.d : |g' \
# 	> $@

-include $(DEPS)
