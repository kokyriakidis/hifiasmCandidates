CXX=		g++
CC=			gcc
# -ffunction-sections/-fdata-sections + -Wl,--gc-sections lets the linker drop
# the unused error-correction / consensus / phasing code that lives alongside
# the base-level alignment path in Correct.cpp and ecovlp.cpp, so only the
# candidate-detection + alignment-filter functions are kept in the binary.
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall -fPIC -ffunction-sections -fdata-sections
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
LDFLAGS=	-Wl,--gc-sections
OBJS=		CommandLines.o Process_Read.o Hash_Table.o \
			kthread.o htab.o hist.o sketch.o anchor.o extract.o sys.o \
			kalloc.o Correct.o ecovlp.o candidates.o hifiasm_overlaps.o
EXE=		hifiasm
LIB=		libhifiasm_overlaps.a
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend test

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

# Static library exposing the hifiasm_overlaps.h C API (no CLI main), for
# linking into other programs such as dinara. $(OBJS) already includes
# hifiasm_overlaps.o and excludes main.o.
lib:$(LIB)

$(LIB):$(OBJS)
		$(AR) rcs $@ $^

# Standalone bridge tests/benchmarks. Each links the C API static library.
# test_sketch_minimizers : unfiltered no-HPC minimizer bridge
# test_sketch_filter      : overlap-parity filter (hf + subsampling), incl. a
#                           direct parity check against mz1_ha_sketch
# bench_sketch_filter     : marker-count / timing comparison filtered vs not
TEST_CXXFLAGS=-std=c++11 -O2
TEST_LIBS=-lz -lpthread -lm

test_sketch_minimizers:test_sketch_minimizers.cpp $(LIB)
		$(CXX) $(TEST_CXXFLAGS) -I. $< $(LIB) $(TEST_LIBS) -o $@

test_sketch_filter:test_sketch_filter.cpp $(LIB)
		$(CXX) $(TEST_CXXFLAGS) -I. $< $(LIB) $(TEST_LIBS) -o $@

bench_sketch_filter:bench_sketch_filter.cpp $(LIB)
		$(CXX) $(TEST_CXXFLAGS) -I. $< $(LIB) $(TEST_LIBS) -o $@

# test_store_overlaps     : shared-read-store path parity vs the file-based path
test_store_overlaps:test_store_overlaps.cpp $(LIB)
		$(CXX) $(TEST_CXXFLAGS) -I. $< $(LIB) $(TEST_LIBS) -o $@

# test_cli_parity         : file-mem AND store paths vs the original hifiasm CLI
#                           (needs the $(EXE) binary; builds it first)
test_cli_parity:test_cli_parity.cpp $(LIB) $(EXE)
		$(CXX) $(TEST_CXXFLAGS) -I. $< $(LIB) $(TEST_LIBS) -o $@

# Build and run the correctness tests.
test:test_sketch_minimizers test_sketch_filter test_store_overlaps test_cli_parity
		./test_sketch_minimizers
		./test_sketch_filter
		./test_store_overlaps
		./test_cli_parity ./$(EXE)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM \
			test_sketch_minimizers test_sketch_filter bench_sketch_filter \
			test_store_overlaps test_cli_parity

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