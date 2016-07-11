#Misha Mikhasenko
#14.07.2015

lSRC    := ../src
lPWA    := $(lSRC)/rpwa
lDIR    := $(lSRC) $(lPWA)
CXFLAGS := -std=c++11 -O2 #-pg -g -std=c++0x  
CFLAGS  := -std=gnu99 -O2 #-pg -g

compile_flags       := -Wall # -MMD -pipe
link_flags          := # -llzma # -s -pipe

####### *.o files
OBJS = \
deflib.o \
mstructures.o mintegrate.o \
MChannel.o MTwoBodyChannel.o MIsobarChannel.o \
MIsobar.o MCoupledChannelIsobar.o MmatrixK.o \
MRelationHolder.o MParKeeper.o \
MChannelPhysics.o MProductionPhysics.o \
MDeck.o

SET :=  $(addprefix obj/, $(OBJS))

####### libconfig
CONF_LIBS += -L/usr/lib/x86_64-linux-gnu -lconfig++

####### root- files

ROOT_DIR  := -I$(shell root-config --incdir)
ROOT_LIBS := $(shell root-config --libs) -lMathMore # -lRooFit -lRooFitCore -lMinuit

ROOT = $(ROOT_DIR) $(ROOT_LIBS)

####### bin- files
EXEC := test read_model_settings \
	test_MIsobar \
	test_MmatrixK \
	test_MProduction \
	test_Secondsheet \
	test_Errors \
	test_Fit \
	test_Damp test_Deck test_MDeck \
	test_Unitarisation
# compass_2-+_3wave

ALL: $(EXES)

 #--------------------- E x e c u t a b l e  p r o g r a m s ----------------------#
 #---------------------------------------------------------------------------------#

%: %.cc $(SET)
	g++ $(link_flags) $(compile_flags) $(CXFLAGS) -o $@ \
	$(addprefix -I,${lDIR}) $^ \
	$(CONF_LIBS) $(ROOT)

test: test.cc $(SET)
	g++ $(link_flags) $(compile_flags) $(CXFLAGS) -o $@ \
	$(addprefix -I,${lDIR}) $^ \
	$(CONF_LIBS) $(ROOT)

 #--------------------------- O b j e c t s ---------------------------------------#
 #---------------------------------------------------------------------------------#

#MClasses
obj/%.o: $(lSRC)/%.cc $(lSRC)/%.h
	g++ $(CXFLAGS) -o $@ -c $< $(compile_flags) \
	$(addprefix -I,${lDIR}) $(ROOT_DIR)

#TBigner
obj/%.o: $(lWGN)/%.cxx
	g++ $(CXFLAGS) -o $@ -c $< $(compile_flags) \
	$(addprefix -I,${lDIR}) $(ROOT_DIR)

#rwpa FunctionsW
obj/%.o: $(lWGN)/%.cxx
	g++ $(CXFLAGS) -o $@ -c $< $(compile_flags) \
	$(addprefix -I,${lDIR}) $(ROOT_DIR)


#Clean
clean:
	rm -f $(EXES) $(wildcard obj/*)

include $(wildcard obj/*.d)

