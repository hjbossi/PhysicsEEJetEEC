Execute2DDoubleLog: RooUnfoldSimple_DoubleLogBins_ZTheta.cpp
	g++ RooUnfoldSimple_DoubleLogBins_ZTheta.cpp -o Execute2DDoubleLog \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs`

ExecuteReweighting: createReweightingHisto.cpp
	g++ createReweightingHisto.cpp -o ExecuteReweighting \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


ExecutePlotReweighting: plotReweighting.cpp
	g++ plotReweighting.cpp -o ExecutePlotReweighting \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o