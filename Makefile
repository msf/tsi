# Makefile for TSI

#                                     
# Compiler settings (gcc, icc, win32)

COMPILER := gcc
ifeq ($(CC), mpicc)
COMPILER := mpicc
endif

ifeq ($(CC), icc)
COMPILER := icc
endif


# Targets
RELEASE_TARGET := tsi
DEBUG_TARGET   := tsid


# Default settings for gcc
CC      := gcc#-3.4
CPP     := g++#-3.4
CFLAGS  := -pipe  -DNEW_RAND 
LDFLAGS := -lm -lgcc 
OPTS    := -O3 -fomit-frame-pointer -ffast-math 
#OPTS	+= -msse -msse2 -mfpmath=sse 
#OPTS	+= -ftree-vectorize # use only with gcc-4 or above.
#OPTS	+= -march=k8
OPTS	+= -g 
#OPTS	+= -DTSI_DEBUG
DEBUG   := -g -ggdb -DTSI_DEBUG2
#DEBUG	+= -m32
#DEBUG  += -pg
DEBUG	+= -Wall -std=gnu99 
DEBUG	+= -Wextra

ifeq ($(COMPILER), icc)
RELEASE_TARGET := tsi-icc
DEBUG_TARGET   := tsi-iccd
CC      := icc
CPP     := icc
OPTS    += -ipo
DEBUG   := -std=gnu99 -Wp64 -DTSI_DEBUG2 -DNEW_RAND -Wall
endif


ifeq ($(COMPILER), mpicc)
RELEASE_TARGET := tsi-mpi
DEBUG_TARGET   := tsi-mpid
CC       := mpicc
CPP      := mpiCC
# use gcc settings
OPTS     += -DTSI_MPI
DEBUG    += -DTSI_MPI
endif

#OS=$(shell uname -s)
#ifeq ($(OS), Linux)
#	CFLAGS += -DLINUX
#endif


# Include files
INCLUDE_DIR := include
INCLUDE_FILES := $(wildcard $(INCLUDE_DIR)/*.h)

# Source files
SOURCE_DIR := src
SOURCE_FILES := $(wildcard $(SOURCE_DIR)/*.c)
DSS_DIR := $(SOURCE_DIR)/dss
DSS_FILES := $(wildcard $(DSS_DIR)/*.c)
SI_DIR  := $(SOURCE_DIR)/si
SI_FILES  := $(wildcard $(SI_DIR)/*.c)
#COMP_DIR  := $(SOURCE_DIR)/compare
#COMP_FILES  := $(wildcard $(COMP_DIR)/*.c)

# Build locations
RELEASE_DIR := obj/release
SRC_DEPS := $(SOURCE_FILES:.c=.o)
SRC_OBJS := $(subst $(SOURCE_DIR),$(RELEASE_DIR),$(SRC_DEPS))
SI_DEPS := $(SI_FILES:.c=.o)
SI_OBJS := $(subst $(SI_DIR),$(RELEASE_DIR),$(SI_DEPS))
DSS_DEPS := $(DSS_FILES:.c=.o)
DSS_OBJS := $(subst $(DSS_DIR),$(RELEASE_DIR),$(DSS_DEPS))
#COMP_DEPS := $(COMP_FILES:.c=.o)
#COMP_OBJS := $(subst $(COMP_DIR),$(RELEASE_DIR),$(COMP_DEPS))

DEBUG_DIR := obj/debug
SRC_DEPS_DEBUG := $(SOURCE_FILES:.c=.o)
SRC_OBJS_DEBUG := $(subst $(SOURCE_DIR),$(DEBUG_DIR),$(SRC_DEPS_DEBUG))
SI_DEPS_DEBUG := $(SI_FILES:.c=.o)
SI_OBJS_DEBUG := $(subst $(SI_DIR),$(DEBUG_DIR),$(SI_DEPS_DEBUG))
DSS_DEPS_DEBUG := $(DSS_FILES:.c=.o)
DSS_OBJS_DEBUG := $(subst $(DSS_DIR),$(DEBUG_DIR),$(DSS_DEPS_DEBUG))
#COMP_DEPS_DEBUG := $(COMP_FILES:.c=.o)
#COMP_OBJS_DEBUG := $(subst $(COMP_DIR),$(DEBUG_DIR),$(COMP_DEPS_DEBUG))


# RELEASE RULES
release: $(RELEASE_TARGET)

$(RELEASE_TARGET): $(SRC_OBJS) $(SI_OBJS) $(DSS_OBJS)
	$(CC) $(OPTS) $(SRC_OBJS) $(SI_OBJS) $(DSS_OBJS) $(LDFLAGS) -o $@

$(SRC_OBJS): $(RELEASE_DIR)/%.o: $(SOURCE_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)

$(DSS_OBJS): $(RELEASE_DIR)/%.o: $(DSS_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)

$(SI_OBJS): $(RELEASE_DIR)/%.o: $(SI_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(OPTS) -c -o $@ $< -I$(INCLUDE_DIR)


# DEBUG RULES
debug: $(DEBUG_TARGET)

$(DEBUG_TARGET): $(SRC_OBJS_DEBUG) $(SI_OBJS_DEBUG) $(DSS_OBJS_DEBUG)
	$(CC) $(DEBUG) $(SRC_OBJS_DEBUG) $(SI_OBJS_DEBUG) $(DSS_OBJS_DEBUG) $(LDFLAGS) -o $@ 

$(SRC_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(SOURCE_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

$(DSS_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(DSS_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)

$(SI_OBJS_DEBUG): $(DEBUG_DIR)/%.o: $(SI_DIR)/%.c $(INCLUDE_FILES) Makefile
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< -I$(INCLUDE_DIR)


# OTHER
wc:
	wc -l $(SOURCE_FILES) $(DSS_FILES) $(SI_FILES) $(INCLUDE_FILES)

clean: FORCE
	rm -f $(RELEASE_DIR)/*.* $(DEBUG_DIR)/*.* $(DEBUG_TARGET) $(RELEASE_TARGET)
	rm -f *~ $(SOURCE_DIR)/*~ $(DSS_DIR)/*~ $(SI_DIR)/*~ $(INCLUDE_DIR)/*~ gmon.out
	rm -f core

all: release debug

rebuild: clean all

FORCE:
#end of Makefile
