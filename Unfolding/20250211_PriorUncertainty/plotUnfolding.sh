inputFileReweighted=unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption1_Reweighted02122025.root
inputFileNonReweighted=unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption1_NOTReweighted02122025.root


echo "-------> Starting to plot unfolding closure"
./ExecutePlotReweighting --Input "unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption1_Reweighted_FittingMethod02122025.root","unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption1_NOTReweighted_NewErrorTreatment02122025.root" \
        --Output ReweightingUncertiantyFirstLook \
        --Label "Reweighted","Nominal"\
        --Prefix "02192025" \
        --DoRatio true \
        --Iter 4 \
        --DoWeight false