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
		
		
ExecuteCov: preUnfoldingCovariance.cpp
	g++ preUnfoldingCovariance.cpp -o ExecuteCov \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

ExecuteCovUnfolding: postUnfoldingCovariance.cpp
	g++ postUnfoldingCovariance.cpp -o ExecuteCovUnfolding \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
		
ExecuteUnfoldE2C: unfoldE2C_flatten_fullcov.C
	g++ unfoldE2C_flatten_fullcov.C -o ExecuteUnfoldE2C \
		-O2 -march=native -std=c++17 \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs` -I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
		
		
ExecuteCov1D: preUnfolding1DCovariance.cpp
	g++ preUnfolding1DCovariance.cpp -o ExecuteCov1D \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
		
		
		
ExecuteUnfoldRW: RooSimpleDoubleLog_WithPoissonRW.cpp
	g++ RooSimpleDoubleLog_WithPoissonRW.cpp -o ExecuteUnfoldRW \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs` -I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


ExecuteUnfoldProjectedCovariance: RooEEC_ProjectedCovariance.cpp
	g++ RooEEC_ProjectedCovariance.cpp -o ExecuteUnfoldProjectedCovariance \
		-I$(RooUnfoldBase)/src -L$(RooUnfoldBase)/ -lRooUnfold -lCore  -Wl,-rpath=$(RooUnfoldBase)/ \
		`root-config --cflags --glibs --libs` -I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o

