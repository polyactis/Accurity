
CXXTARGETS	=
#CXXLDFLAGS 	= $(BoostLib)
CTARGETS	=
SharedLibTargets	=
StaticLibTargets =


SRCS	= infer.cpp read_para.cpp BaseGADA.cc GADA.cc

ExtraTargets = infer GADA

infer:	%:	%.o read_para.o prob.o BaseGADA.o format.o
	$(CXXCOMPILER) $< read_para.o prob.o BaseGADA.o format.o $(CXXFLAGS) -o $@ $(CXXLDFLAGS) -lgsl -lgslcblas $(BoostLib)

GADA:   %:   %.o BaseGADA.o BaseGADA.h read_para.o format.o
	$(CXXCOMPILER) $< BaseGADA.o read_para.o format.o $(CXXFLAGS) -o $@ -lm $(CXXLDFLAGS) $(BoostLib)

recall_precision:	%:	%.o
	$(CXXCOMPILER) $< $(CXXFLAGS) -o $@ $(CXXLDFLAGS)


#:= is different from =. The latter will cause the function evaluation every time the make variable is invoked.
#currentUnixTime:=$(shell echo "import time; print str(time.time()).replace('.', '_')"|python)
currentTime:=$(shell echo "from datetime import datetime; print str(datetime.now()).replace(' ', '_').replace(':', '')"|python)


debug: all ../src/main.rs
	-rm -rf debug
	-mkdir debug
	-mkdir -p ../target/debug/
	cargo build
	git checkout -- ../src/main.rs
	cp -r __init__.py ../LICENSE GADA ../target/debug/accurity main.py configure infer plotCPandMCP.py plot_autocor_diff.py plot_coverage_after_normalization.py plot_tre.py plot.tre.autocor.R plot_snp_maf_exp.py plot_snp_maf_peak.py debug/
	tar -cavf debug.$(currentTime).tar.gz debug/

release: all ../src/main.rs
	-rm -rf release
	-mkdir release
	-mkdir -p ../target/release/
	cargo build --release
	git checkout -- ../src/main.rs
	cp -r __init__.py ../LICENSE GADA ../target/release/accurity main.py configure infer plotCPandMCP.py plot_autocor_diff.py plot_coverage_after_normalization.py plot_tre.py plot.tre.autocor.R plot_snp_maf_exp.py plot_snp_maf_peak.py release/
	tar -cavf release.$(currentTime).tar.gz release/


include ../Makefile.common

# leave it at last to avoid being overridden by Makefile.common
clean: common_clean
	-rm -r release*
	-rm -r debug*
