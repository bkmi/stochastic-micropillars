## MAKEFILE FUER DEN SOLVER ##

### Hier kannst du deinen Compiler wählen. Ich mag clang, der hat nämlich sinnvolle Fehlerausgaben bei Errors.
CXX:=clang++
### Optimize Level, während des Testens kann hier ruhig O0 stehen, sonst ist O3 eine gute Wahl für speed!
OPT?=-O3 -g 
### Die Zeile nicht ändern, inkludiert die Header Files vom Solver und den externen Modulen. 
INCLUDE=-I/$$PWD/solver_inc -I/$$PWD/add_inc/ -I/$$PWD/add_inc/armadillo-7.800.2/include -DARMA_DONT_USE_WRAPPER
### C++ Flags, nicht anfassen wenn du nicht weißt was die machen
CPPFLAGS:=-Wall -pedantic $(OPT) $(INCLUDE)
### C++ Standard, same.
CXXFLAGS:=-std=c++11
### Gibt den Pfad zu den externen C++ (.cc) Files an. Nicht ändern, es sei denn du willst noch mehr CC files dazu tun und hast die Module 
#######  in einem weiteren Ordner. Dann auch die INCLUDE Zeile bearbeiten.
LOCALLIBDIR=$$PWD/add_src/
### Hier sind Standard Bibliotheken, wie BOOST und FFTW3 abgelegt, nicht ändern, es sei denn das Dateisystem ändert sich.
EXTERNALLIBDIR=/usr/lib/x86_64-linux-gnu/
### Nix ändern hier.
LDFLAGS= -L$(LOCALLIBDIR) -L$(EXTERNALLIBDIR)
### Linkt die Libaries: FFTW3, Math, Boost::System, Boost::Filesystem. Wenn du noch mehr libaries linken musst, das ist der Platz dafür. =)
LDLIBS= -lfftw3 -lm -lboost_system -lboost_filesystem -lblas -llapack

### Nix zu ändern auf diesen Zeilen, es werden alle .cc Files aus dem Ordner add_src compiliert
### und die Maschinencodes in den obj/ Ordner gepackt. Sie werden durch $(ADD_OBJ) mit deinem Programm gelinkt.

ADD_SRC= $(wildcard add_src/*.cc)					
ADD_OBJ= $(addprefix obj/,$(notdir $(ADD_SRC:.cc=.o)))			


###compile und link
### DAS HIER IST MAGIE! NICHT ANFASSEN =)
compile_and_dependencies=$(CXX) $(CPPFLAGS) $(CXXFLAGS) \
						 -MMD -MF $(@:.o=.d) -o $@ -c $<

link=$(CXX) $(OPT) -o $@ $^ $(LDFLAGS) $(LDLIBS)


### make, hier musst du alle "EXAMPLE" durch "DEINPROGRAMM_NAME" ersetzen. Optional kannst du auch EXAMPLE lassen
### und "DEINPROGRAMM_NAME" hinter die EXAMPLE.o kopieren. Dann auch YOURPROGAMM_NAME mit in das Target "all:" einfügen.

all: 	testben muPillar_FB #muPillar_Pulse
								
testben: testben.o $(ADD_OBJ)						
	$(call link)

testben.o: testben.cc  
	$(call compile_and_dependencies)

### BEGIN YOUR MODELS

### Feedback Model
muPillar_FB: muPillar_FB.o $(ADD_OBJ)						
	$(call link)

muPillar_FB.o: muPillar_FB.cc  
	$(call compile_and_dependencies)

### Pulse Injection Model
#imuPillar_Pulse: muPillar_Pulse.o $(ADD_OBJ)					#	
#	$(call link)
#
#muPillar_Pulse.o: muPillar_Pulse.cc  
#	$(call compile_and_dependencies)
#



### END YOUR MODELS
obj/%.o: add_src/%.cc
	$(call compile_and_dependencies)


.PHONY: clean
clean:
	rm *.o *.d 
	rm obj/*.o obj/*.d 

