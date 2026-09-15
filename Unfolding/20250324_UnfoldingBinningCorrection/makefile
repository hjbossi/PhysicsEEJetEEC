default: TestRun

TestRun: Execute
	./Execute --Input  /home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/unfoldingE2C_DataUnfolding_03242025.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_z.root \
			     --UnfoldingBinCorrArgName z --MakeUnfoldingBinCorrFactor true
	./Execute --Input  /home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/unfoldingE2C_DataUnfolding_03242025.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_z.root \
			     --UnfoldingBinCorrArgName z
	./Execute --Input  /home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/unfoldingE2C_DataUnfolding_03242025.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_theta.root \
			     --UnfoldingBinCorrArgName theta --MakeUnfoldingBinCorrFactor true
	./Execute --Input  /home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/unfoldingE2C_DataUnfolding_03242025.root \
			     --UnfoldingBinCorrName UnfoldingBinCorr_with_theta.root \
			     --UnfoldingBinCorrArgName theta

Execute: unfoldingBinningCorrection.cpp
	g++ unfoldingBinningCorrection.cpp -o Execute \
		`root-config --glibs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o