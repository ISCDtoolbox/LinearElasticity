# Set compiler's options depending on the architecture
ifeq ($(ARCHI),i86)
  CC     = cc  
  CFLAGS = -c -O3 -ftree-vectorize -Wuninitialized -Wunused -Winline -Wshadow
  LDFLAGS=
endif
ifeq ($(ARCHI),ppc)
  CC    = gcc
  CFLAGS= -O3 -c\
          -Wuninitialized -Wunused -Winline -Wshadow \
          -fexpensive-optimizations -funroll-loops
  LDFLAGS= -s
endif
ifeq ($(ARCHI), gccub)
  CC      = gcc
  CFLAGS  = -c -O3 -fast -isysroot /Developer/SDKs/MacOSX10.5.sdk -arch i386 -arch ppc -mcpu=G4
  LDFLAGS = -isysroot /Developer/SDKs/MacOSX10.5.sdk -arch i386 -arch ppc -mcpu=G4
endif
ifeq ($(ARCHI),win32)
  CC    = gcc
  CFLAGS= -c -O3 -mno-cygwin -Wuninitialized -Wunused -Winline -Wshadow
  LDFLAGS=-v -s -mno-cygwin
endif

# working dirs
EXEDIR = $(HOME)/bin/$(ARCHI)
SRCDIR = sources
OBJDIR = objects/$(ARCHI)
ARCDIR = archives
DIRDIR = objects $(EXEDIR) $(OBJDIR) $(ARCDIR)
INCDIR = 
LDLDIR =
VPATH  = $(SRCDIR)

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
prog   = elastic

#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(INCDIR) $(CFLAGS) $< -o $@

$(EXEDIR)/$(prog):$(DIRDIR) $(objs) $(LDLDIR)
	echo "#define COMPIL " '"' `date` '"' > $(SRCDIR)/compil.date
	$(CC) $(CFLAGS) $(INCDIR) $(SRCDIR)/elastic.c -o $(OBJDIR)/elastic.o
	$(CC) $(LDFLAGS) $(LDLDIR) $(objs) -o $@ -lm

$(objs):$(header)

$(DIRDIR):
	@[ -d $@ ] || mkdir $@

clean:
	-rm $(objs) $(EXEDIR)/$(prog)

tar:$(DIRDIR)
	tar czfL $(ARCDIR)/$(prog).`date +"%Y.%m.%d"`.tgz sources makefile

target: $(EXEDIR)/$(prog)
