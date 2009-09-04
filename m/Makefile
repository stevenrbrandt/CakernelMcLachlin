# -*-Makefile-*-

all: McLachlan_ADM.out McLachlan_BSSN.out McLachlanW.out McLachlanUp.out McLachlan6.out WaveToy.out WaveToyFO.out hydro.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan_ADM.out: McLachlan_ADM.m
	rm -rf ML_ADM*
	./runmath.sh $^
	for thorn in ML_ADM*; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
	done

McLachlan_BSSN.out: McLachlan_BSSN.m
	rm -rf ML_BSSN*
	./runmath.sh $^
	for thorn in ML_BSSN*; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
		./create-helper-thorn.sh $$thorn; \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

McLachlanW.out: McLachlanW.m
	rm -rf ML_BSSNW
	./runmath.sh $^
	for thorn in ML_BSSNW; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
		./create-helper-thorn.sh $$thorn; \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

McLachlanUp.out: McLachlanUp.m
	rm -rf ML_BSSNUp
	./runmath.sh $^
	for thorn in ML_BSSNUp; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
		./create-helper-thorn.sh $$thorn; \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

McLachlan6.out: McLachlan6.m
	rm -rf ML_BSSN6
	./runmath.sh $^
	for thorn in ML_BSSN6; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
		./create-helper-thorn.sh $$thorn; \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
	done

WaveToyFO.out: WaveToyFO.m
	rm -rf ML_FOWaveToy
	./runmath.sh $^
	for thorn in ML_FOWaveToy; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
	done

hydro.out: hydro.m
	rm -rf ML_hydro
	./runmath.sh $^
	for thorn in ML_hydro; do \
		./copy-if-changed.sh $$thorn ../$$thorn; \
	done

clean:
	rm -rf ML_ADM* ML_BSSN*
	rm -rf ML_WaveToy ML_FOWaveToy
	rm -rf ML_hydro
	rm -rf McLachlan_ADM.out McLachlan_ADM.err
	rm -rf McLachlan_BSSN.out McLachlan_BSSN.err
	rm -rf McLachlanW.out McLachlanW.err
	rm -rf McLachlanUp.out McLachlanUp.err
	rm -rf McLachlan6.out McLachlan6.err
	rm -rf WaveToy.out WaveToy.err
	rm -rf hydro.out hydro.err

.PHONY: all clean
