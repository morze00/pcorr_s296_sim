#CXX= g++
CXX= clang++
#LINKER= g++
LINKER= clang++
LDFLAGS= `root-config --libs`
CPPFLAGS= `root-config --cflags`
CXXFLAGS= -Wno-deprecated -g
LNKFLAGS= -g
#XBDIR=xbtools
XBDIR=xbtools
#DIR_INC=-I.
DIR_INC=-I$(XBDIR)
EXEC=p_corr

OBJ=  p_corr.o track_clone.o xbtools.o

SRC=  libs.hh	\
	track_R3Bsim.h	\
	xbtools.h	\
	constants.hh	\
	xb_list_energy.hh

#DEPS = $(patsubst %,$(IDIR)/%,$(SRC))

%.o : %.cc $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} $(DIR_INC) -c $< -o $@
	@echo "	CXX $@"

%.o : %.C $(DEPS)
	@$(MAKEDEPEND)
	@${CXX} ${CPPFLAGS} ${CXXFLAGS} $(DIR_INC) -c $< -o $@
	@echo "	CXX $@"

default: $(OBJ)
	@${LINKER} ${LDFLAGS} ${LNKFLAGS} -o $(EXEC) $(DIR_INC) $(OBJ)
	@echo "	COMP $(EXEC)"

clean:
	rm -f $(OBJ)
	rm -f $(EXEC)
