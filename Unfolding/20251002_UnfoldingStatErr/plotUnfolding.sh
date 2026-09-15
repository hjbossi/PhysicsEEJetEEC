
echo "-------> Starting to plot unfolding closure"
./ExecutePlotUnfoldingStab --Input "/home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/unfoldingE2C_DataUnfolding_03242025.root" \
        --Output UnfoldingStab \
        --Label "Iteration 4 (Nominal)","Iteration 1","Iteration 2","Iteration 3","Iteration 4"\
        --Prefix "03242025" \
        --DoRatio true \
        --Iter 4 \
        --DoWeight false