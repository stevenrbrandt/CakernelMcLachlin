# -*-Makefile-*-

all: McLachlan.out McLachlanW.out WaveToy.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan.out: McLachlan.m
	rm -rf ML_ADM ML_BSSN
	./runmath.sh $^
	for thorn in ML_ADM ML_BSSN; do ./copy-if-changed.sh $$thorn ../$$thorn; done

McLachlanW.out: McLachlanW.m
	rm -rf ML_BSSNW
	./runmath.sh $^
	for thorn in ML_BSSNW; do ./copy-if-changed.sh $$thorn ../$$thorn; done

McLachlanUp.out: McLachlanUp.m
	rm -rf ML_BSSNUp
	./runmath.sh $^
	for thorn in ML_BSSNUp; do ./copy-if-changed.sh $$thorn ../$$thorn; done

McLachlan6.out: McLachlan6.m
	rm -rf ML_BSSN6
	./runmath.sh $^
	for thorn in ML_BSSN6; do ./copy-if-changed.sh $$thorn ../$$thorn; done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy ML_FOWaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy ML_FOWaveToy; do ./copy-if-changed.sh $$thorn ../$$thorn; done

clean:
	rm -rf McLachlan.out McLachlan.err ML_ADM ML_BSSN
	rm -rf McLachlanW.out McLachlanW.err ML_BSSNW
	rm -rf WaveToy.out WaveToy.err ML_WaveToy ML_FOWaveToy

.PHONY: all clean
