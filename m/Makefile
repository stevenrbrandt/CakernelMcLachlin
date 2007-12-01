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
	rm -rf ML_WaveToy MLFO_Wavetoy
	./runmath.sh $^
	for thorn in ML_WaveToy MLFO_Wavetoy; do ./copy-if-changed.sh $$thorn ../$$thorn; done

clean:
	rm -rf McLachlan.out McLachlan.err ML_ADM ML_BSSN
	rm -rf Wavetoy.out Wavetoy.err ML_WaveToy MLFO_Wavetoy

.PHONY: all clean
