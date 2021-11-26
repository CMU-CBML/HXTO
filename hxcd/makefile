PETSC_DIR=/ocean/projects/eng170006p/xliang4/petsc-3.14.2
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test
Compiler=mpicxx
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX}  -std=c++11 -fopenmp -Wall
#-D_GLIBCXX_USE_CXX11_ABI=0

all: transport

transport: BasicDataStructure.o UserSetting.o Transport.o main.o 
	$(Compiler) BasicDataStructure.o UserSetting.o Transport.o main.o -o transport ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
UserSetting.o: UserSetting.cpp
	$(Compiler)	-c UserSetting.cpp $(CFLAGS)
	
Transport.o: Transport.cpp
	$(Compiler)	-c Transport.cpp $(CFLAGS)
	
main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *o transport

