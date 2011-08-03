# -*-Makefile-*-

all: McLachlan_ADM.out McLachlan_BSSN.out McLachlan_ADMConstraints.out McLachlan_ADMQuantities.out WaveToy.out WaveToyFO.out hydro.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

McLachlan_ADM.out: McLachlan_ADM.m
	rm -rf ML_ADM*
	./runmath.sh $^
	for thorn in ML_ADM*; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

McLachlan_BSSN.out: McLachlan_BSSN.m
	rm -rf ML_BSSN*
	./runmath.sh $^
	for thorn in ML_BSSN*; do					   \
		./copy-if-changed.sh $$thorn ../$$thorn &&		   \
		./create-helper-thorn.sh $$thorn &&			   \
		./copy-if-changed.sh $${thorn}_Helper ../$${thorn}_Helper; \
	done

McLachlan_ADMConstraints.out: McLachlan_ADMConstraints.m
	rm -rf ML_ADMConstraints*
	./runmath.sh $^
	for thorn in ML_ADMConstraints*; do			\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

McLachlan_ADMQuantities.out: McLachlan_ADMQuantities.m
	rm -rf ML_ADMQuantities*
	./runmath.sh $^
	for thorn in ML_ADMQuantities*; do			\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

WaveToy.out: WaveToy.m
	rm -rf ML_WaveToy
	./runmath.sh $^
	for thorn in ML_WaveToy; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

WaveToyFO.out: WaveToyFO.m
	rm -rf ML_WaveToyFO
	./runmath.sh $^
	for thorn in ML_WaveToyFO; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

hydro.out: hydro.m
	rm -rf ML_hydro
	./runmath.sh $^
	for thorn in ML_hydro; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

clean:
	rm -rf ML_ADM* ML_BSSN* ML_ADMConstraints*
	rm -rf ML_WaveToy ML_WaveToyFO
	rm -rf ML_hydro
	rm -f McLachlan_ADM.out McLachlan_ADM.err
	rm -f McLachlan_BSSN.out McLachlan_BSSN.err
	rm -f McLachlan_ADMConstraints.out McLachlan_ADMConstraints.err
	rm -f McLachlan_ADMQuantities.out McLachlan_ADMQuantities.err
	rm -f WaveToy.out WaveToy.err
	rm -f hydro.out hydro.err

.PHONY: all clean
