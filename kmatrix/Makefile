#Misha Mikhasenko
#14.07.2015

lSRC    := ../src
lDIR    := $(lSRC)
CXFLAGS := -std=c++11 -O2 #-pg -g -std=c++0x  
CFLAGS  := -std=gnu99 -O2 #-pg -g

compile_flags       := -Wall # -MMD -pipe
link_flags          := -s -pipe

####### *.o files
OBJS = \
deflib.o \
MIsobar.o MCoupledChannelIsobar.o NoverD.o MmatrixK.o \
mstructures.o mintegrate.o \
MModel.o MRelationHolder.o MParKeeper.o MChannelPhysics.o MProductionPhysics.o
# MStructureHolder.o 

SET :=  $(addprefix obj/, $(OBJS))

####### libconfig
LIBS += /usr/lib64/libconfig++.so

####### root- files

ROOT_DIR  := -I$(shell root-config --incdir)
ROOT_LIBS := $(shell root-config --libs)# -lMathMore -lRooFit -lRooFitCore -lMinuit

ROOT = $(ROOT_DIR) $(ROOT_LIBS)

####### bin- files
EXEC := test compass_2-+_3wave
ALL: $(EXES)

 #--------------------- E x e c u t a b l e  p r o g r a m s ----------------------#
 #---------------------------------------------------------------------------------#

%: %.cc $(SET)
	g++ $(CXFLAGS) -o $@  $(LIBS) \
	$(addprefix -I,${lDIR}) $^ \
	$(ROOT)

test: test.cc $(SET)
	g++ $(CXFLAGS) -o $@  $(LIBS) \
	$(addprefix -I,${lDIR}) $^ \
	$(ROOT)

 #--------------------------- O b j e c t s ---------------------------------------#
 #---------------------------------------------------------------------------------#

#MClasses
obj/%.o: $(lSRC)/%.cc
	g++ $(CXFLAGS) -o $@ -c $< $(compile_flags) \
	$(addprefix -I,${lDIR}) $(ROOT_DIR)


#Clean
clean:
	rm -f $(EXES) $(wildcard obj/*)

include $(wildcard obj/*.d)

