default: TestRun

TestRun: Execute
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MCMerged.root \
			     --EvtSelEffName EvtSelEff.root --MakeEvtSelEffCorrFactor true
	./Execute --Input /data/ALEPH/MC/LEP1MC/LEP1MCMerged.root \
			     --EvtSelEffName EvtSelEff.root

Execute: evtSelEffCorr.cpp
	g++ evtSelEffCorr.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o