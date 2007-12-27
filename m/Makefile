# -*-Makefile-*-

all: McLachlan.out WaveToy.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan.out: McLachlan.m
	rm -rf ML_ADM ML_BSSN
	./runmath.sh $^
	for thorn in ML_ADM ML_BSSN; do ./copy-if-changed.sh $$thorn ../$$thorn; done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy ML_FOWaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy ML_FOWaveToy; do ./copy-if-changed.sh $$thorn ../$$thorn; done

clean:
	rm -rf McLachlan.out McLachlan.err ML_ADM ML_BSSN
	rm -rf WaveToy.out WaveToy.err ML_WaveToy ML_FOWaveToy

.PHONY: all clean
