Execute2DDoubleLog: RooUnfoldSimple_DoubleLogBins_ZTheta.cpp
	g++ RooUnfoldSimple_DoubleLogBins_ZTheta.cpp -o Execute2DDoubleLog \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs` -I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


Execute: fillPairTree.cpp
	g++ fillPairTree.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExecutePlotUnfoldingStab: plotUnfoldingStab.cpp
	g++ plotUnfoldingStab.cpp -o ExecutePlotUnfoldingStab \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


ExecuteClosureCheck: FullAnalysisChainClosureCheck.cpp
	g++ FullAnalysisChainClosureCheck.cpp -o ExecuteClosureCheck \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs` -I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o