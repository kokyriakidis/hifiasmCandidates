CXX=		g++
CC=			gcc
# -ffunction-sections/-fdata-sections + -Wl,--gc-sections lets the linker drop
# the unused error-correction / consensus / phasing code that lives alongside
# the base-level alignment path in Correct.cpp and ecovlp.cpp, so only the
# candidate-detection + alignment-filter functions are kept in the binary.
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall -ffunction-sections -fdata-sections
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
LDFLAGS=	-Wl,--gc-sections
OBJS=		CommandLines.o Process_Read.o Hash_Table.o \
			kthread.o htab.o hist.o sketch.o anchor.o extract.o sys.o \
			kalloc.o Correct.o ecovlp.o candidates.o
EXE=		hifiasm
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

CommandLines.o: CommandLines.h ketopt.h
Hash_Table.o: Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h kdq.h
Hash_Table.o: CommandLines.h ksort.h
Process_Read.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
anchor.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
anchor.o: ksort.h Hash_Table.h
extract.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h khashl.h
extract.o: kseq.h
hist.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
htab.o: kthread.h khashl.h kseq.h ksort.h htab.h Process_Read.h Overlaps.h
htab.o: kvec.h kdq.h CommandLines.h
kthread.o: kthread.h
main.o: CommandLines.h Process_Read.h htab.h
sketch.o: kvec.h htab.h Process_Read.h Overlaps.h kdq.h CommandLines.h
sys.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
kalloc.o: kalloc.h
candidates.o: CommandLines.h Process_Read.h Hash_Table.h Overlaps.h htab.h
candidates.o: kthread.h ksort.h