# -*-Makefile-*-

all: McLachlan_ADM.out McLachlan_BSSN.out WaveToy.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan_ADM.out: McLachlan_ADM.m
#	rm -rf ML_ADM*
#	./runmath.sh $^
#	for thorn in ML_ADM*; do ./copy-if-changed.sh $$thorn ../$$thorn; done

McLachlan_BSSN.out: McLachlan_BSSN.m
	rm -rf ML_BSSN*
	./runmath.sh $^
	for thorn in ML_BSSN*; do ./copy-if-changed.sh $$thorn ../$$thorn; done
WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy ML_FOWaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy ML_FOWaveToy; do ./copy-if-changed.sh $$thorn ../$$thorn; done

clean:
	rm -rf McLachlan_ADM.out McLachlan_ADM.err ML_ADM*
	rm -rf McLachlan_BSSN.out McLachlan_BSSN.err ML_BSSN*
	rm -rf WaveToy.out WaveToy.err ML_WaveToy ML_FOWaveToy

.PHONY: all clean
