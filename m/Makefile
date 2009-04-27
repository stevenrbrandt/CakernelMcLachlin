# -*-Makefile-*-

all: McLachlan_ADM.out McLachlan_BSSN.out McLachlanW.out WaveToy.out hydro.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan_ADM.out: McLachlan_ADM.m
	rm -rf ML_ADM*
	./runmath.sh $^
	for thorn in ML_ADM*; do ./copy-if-changed.sh $$thorn ../$$thorn; done

McLachlan_BSSN.out: McLachlan_BSSN.m
	rm -rf ML_BSSN*
	./runmath.sh $^
	for thorn in ML_BSSN*; do ./copy-if-changed.sh $$thorn ../$$thorn; done

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

McLachlan_Psilon.out: McLachlan_Psilon.m
	rm -rf ML_Psilon*
	./runmath.sh $^
	for thorn in ML_Psilon*; do ./copy-if-changed.sh $$thorn ../$$thorn; done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy ML_FOWaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy ML_FOWaveToy; do ./copy-if-changed.sh $$thorn ../$$thorn; done

hydro.out: hydro.m
	rm -rf ML_hydro
	./runmath.sh $^
	for thorn in ML_hydro; do ./copy-if-changed.sh $$thorn ../$$thorn; done

clean:
	rm -rf McLachlan_ADM.out McLachlan_ADM.err ML_ADM*
	rm -rf McLachlan_BSSN.out McLachlan_BSSN.err ML_BSSN*
	rm -rf McLachlanW.out McLachlanW.err 
	rm -rf McLachlanUp.out McLachlanUp.err 
	rm -rf McLachlan6.out McLachlan6.err 
	rm -rf WaveToy.out WaveToy.err ML_WaveToy ML_FOWaveToy

.PHONY: all clean
