# Makefile for TSI
                                     
# Compiler settings (gcc, icc, win32)
COMPILER := gcc

ifeq ($(COMPILER), icc)
CC       := icc
CPP      := icc
CFLAGS   := -O3
LDFLAGS  := -lm -lc -lpthread -lstdc++ -lgcc
OPTS     :=
DEBUG    :=

elif ($(COMPILER), win32)
CC       := cl
CPP      := cl
CFLAGS   :=
LDFLAGS  :=
OPTS     :=
DEBUG    :=

elif ($(COMPILER), mpicc)
CC       := mpicc
CPP      := mpiCC
CFLAGS   := -pipe -march=pentium2
LDFLAGS  :=  -lm -lc -lpthread -lstdc++ -lgcc
OPTS     := -DTSI_MPI -pipe -O3 -ffast-math -m32 -fomit-frame-pointer
DEBUG    := -DTSI_MPI -g -pg -ggdb -DTSI_DEBUG

else  #gcc
CC       := gcc-4.0
CPP      := g++-4.0
CFLAGS   := -pipe -march=pentium2
LDFLAGS  := -lm -lc -lpthread -lstdc++ -lgcc
OPTS     := -O3 -ffast-math -m32 -fomit-frame-pointer -DTSI_DEBUG2
#OPTS += -fthread-jumps -fcrossjumping -foptimize-sibling-calls
#OPTS += -fcse-follow-jumps  -fcse-skip-blocks -fgcse  -fgcse-lm
#OPTS += -fexpensive-optimizations -fstrength-reduce -frerun-cse-after-loop
#OPTS += -frerun-loop-opt -fcaller-saves -fforce-mem -fpeephole2
#OPTS += -fschedule-insns  -fschedule-insns2 -fsched-interblock
#OPTS += -fsched-spec -fregmove -fstrict-aliasing -fdelete-null-pointer-checks
#OPTS += -freorder-blocks -freorder-functions -funit-at-a-time
#OPTS += -falign-functions  -falign-jumps -falign-loops  -falign-labels
#OPTS += -ftree-pre
#OPTS += -finline-functions -funswitch-loops -fgcse-after-reload
DEBUG    := -g -pg -ggdb -DTSI_DEBUG 
#DEBUG    += -Wall -pedantic -ansi -Wextra -Wcast-qual -Wcast-align -Wconversion
#DEBUG    += -Winit-self -Wswitch-default -Wswitch-enum  -Wfloat-equal -Wshadow
#DEBUG    += -Wmissing-field-initializers -Wunreachable-code -Wdisabled-optimization
#DEBUG    += -Wmissing-prototypes -Wmissing-declarations -Wdeclaration-after-statement 
endif


# Include files
INCLUDE_DIR := include
INCLUDE_FILES := $(wildcard $(INCLUDE_DIR)/*.h)

# Source files
SOURCE_DIR := src
SOURCE_FILES := $(wildcard $(SOURCE_DIR)/*.c)
#DSS_DIR := $(SOURCE_DIR)/dss
#DSS_FILES := $(wildcard $(DSS_DIR)/*.c)
#SI_DIR  := $(SOURCE_DIR)/si
#SI_FILES  := $(wildcard $(SI_DIR)/*.c)
#COMP_DIR  := $(SOURCE_DIR)/compare
#COMP_FILES  := $(wildcard $(COMP_DIR)/*.c)

# Build locations
RELEASE_DIR := obj/release
SRC_DEPS := $(SOURCE_FILES:.c=.o)
SRC_OBJS := $(subst $(SOURCE_DIR),$(RELEASE_DIR),$(SRC_DEPS))
#SI_DEPS := $(SI_FILES:.c=.o)
#SI_OBJS := $(subst $(SI_DIR),$(RELEASE_DIR),$(SI_DEPS))
#DSS_DEPS := $(DSS_FILES:.c=.o)
#DSS_OBJS := $(subst $(DSS_DIR),$(RELEASE_DIR),$(DSS_DEPS))
#COMP_DEPS := $(COMP_FILES:.c=.o)
#COMP_OBJS := $(subst $(COMP_DIR),$(RELEASE_DIR),$(COMP_DEPS))

DEBUG_DIR := obj/debug
SRC_DEPS_DEBUG := $(SRC_FILES:.c=.o)
SRC_OBJS_DEBUG := $(subst $(SOURCE_DIR),$(DEBUG_DIR),$(SRC_DEPS_DEBUG))
#SI_DEPS_DEBUG := $(SI_FILES:.c=.o)
#SI_OBJS_DEBUG := $(subst $(SI_DIR),$(DEBUG_DIR),$(SI_DEPS_DEBUG))
#DSS_DEPS_DEBUG := $(DSS_FILES:.c=.o)
#DSS_OBJS_DEBUG := $(subst $(DSS_DIR),$(DEBUG_DIR),$(DSS_DEPS_DEBUG))
#COMP_DEPS_DEBUG := $(COMP_FILES:.c=.o)
#COMP_OBJS_DEBUG := $(subst $(COMP_DIR),$(DEBUG_DIR),$(COMP_DEPS_DEBUG))

# Targets
RELEASE_TARGET := tsi
DEBUG_TARGET   := tsi-debug


# RELEASE RULES
#release: $(RELEASE_DIR)/$(RELEASE_TARGET)
release: $(RELEASE_TARGET)

#$(RELEASE_DIR)/$(RELEASE_TARGET): $(DSS_OBJS) $(SI_OBJS) $(COMP_OBJS)
#	$(CC) $(LDFLAGS) $(OPTS) -o $@ -O $(DSS_OBJS) $(SI_OBJS)

#$(RELEASE_DIR)/$(RELEASE_TARGET): $(SRC_OBJS)
$(RELEASE_TARGET): $(SRC_OBJS)
	$(CC) $(LDFLAGS) $(OPTS) -o $@ -O $(SRC_OBJS)

$(SRC_OBJS): $(RELEASE_DIR)/%.o: $(SOURCE_DIR)/%.c $(INCLUDE_FILES)
	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)

#$(DSS_OBJS): $(RELEASE_DIR)/%.o: $(DSS_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)

#$(SI_OBJS): $(RELEASE_DIR)/%.o: $(SI_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)

#$(COMP_OBJS): $(RELEASE_DIR)/%.o: $(COMP_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)


# DEBUG RULES
debug: $(DEBUG_TARGET)

#$(DEBUG_DIR)/$(DEBUG_TARGET): $(DSS_OBJS_DEBUG) $(SI_OBJS_DEBUG) $(COMP_OBJS_DEBUG)
#	$(CC) $(LDFLAGS) $(DEBUG) -o $@ -O $(DSS_OBJS_DEBUG) $(SI_OBJS_DEBUG)

#$(DEBUG_DIR)/$(DEBUG_TARGET): $(SRC_OBJS_DEBUG)
$(DEBUG_TARGET): $(SRC_OBJS_DEBUG)
	$(CC) $(LDFLAGS) $(DEBUG) -o $@ -O $(SRC_OBJS_DEBUG)

$(SRC_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(SOURCE_DIR)/%.c $(INCLUDE_FILES)
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

#$(DSS_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(DSS_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

#$(SI_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(SI_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

#$(COMP_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(COMP_DIR)/%.c $(INCLUDE_FILES)
#	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

# OTHER
clean: FORCE
	rm -f $(RELEASE_DIR)/*.o $(DEBUG_DIR)/*.o $(DEBUG_TARGET) $(RELEASE_TARGET)

all: release debug

rebuild: clean all

FORCE:
#end of Makefile
