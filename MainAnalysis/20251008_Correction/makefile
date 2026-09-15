

default: TestRun

TestRun: Execute
	./Execute --InputData /home/hbossi/PhysicsEEJetEEC/Unfolding/20251002_UnfoldingStatErr/unfoldingE2C_DataUnfolding_StatErrCheck_10022025.root

StatUpdate: Execute
	./Execute --InputData /home/hbossi/PhysicsEEJetEEC/Unfolding/20251002_UnfoldingStatErr/unfoldingE2C_ProjectedCovariance_09122026.root



Execute: Correction.cpp
	g++ Correction.cpp -o Execute `root-config --cflags --libs` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o
