hn ?= $(shell hostname)
username ?= $(shell whoami)

ifeq "$(hn)" "monolith.mechse.illinois.edu"
include make.monolith
endif

ifeq "$(hn)" "wirelessprv-10-194-150-85.near.illinois.edu"
include make.mattiaLaptop
endif

ifeq "$(hn)" "cedric"
include make.cedric
endif


ifeq "$(config)" "production"
CCFLAGS += -DNDEBUG
endif

SRC_DIR = ./
OBJ_DIR = ./

TEST_OBJS = \
	test.o

# Linking stage
test: ${TEST_OBJS}
	${CC} ${CPPSETTINGS} ${TEST_OBJS} -o $@ ${LIBS} ${LDFLAGS}

# Compiling stage
${OBJ_DIR}/%.o: ${SRC_DIR}/%.cpp 
	${CC} ${CPPSETTINGS} ${INCLUDES} ${CCFLAGS} -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.s rm -f test
