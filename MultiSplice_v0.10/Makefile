all: SeparateChrFromGTF SeparateGeneFromGTF SeparateChrFromSAM fragment Estimate process MultiSplice

OPTFLAGS = -O3

#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)

SeparateChrFromGTF:
	g++ $(CFLAGS) -o bin/SeparateChrFromGTF src/SeparateChrFromGTF.cpp 2>log

SeparateGeneFromGTF:
	g++ $(CFLAGS) -o bin/SeparateGeneFromGTF src/SeparateGeneFromGTF.cpp 2>log

SeparateChrFromSAM:
	g++ $(CFLAGS) -o bin/SeparateChrFromSAM src/SeparateChrFromSAM.cpp 2>log

fragment:
	g++ $(CFLAGS) -o bin/fragment src/parseFragment.cpp 2>log

Estimate:
	g++ $(CFLAGS) -o bin/Estimate src/Estimate.cpp 2>log

process:
	g++ $(CFLAGS) -o bin/process src/MergeOverlapGTF.cpp 2>log

MultiSplice:
	g++ $(CFLAGS) -o MultiSplice src/MultiSplice.cpp 2>log

clean: