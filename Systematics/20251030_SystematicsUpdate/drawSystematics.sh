# inputFile=SystematicsE2C_Theta_10302025.root
# ./Execute --Input $inputFile \
#         --Variable "Theta"\
#         --Label "TPCHits","Matching","Reweighting","MCBinning","Regularization","MCStat","Total"\
#         --Prefix "10302025"

inputFile=SystematicsE2C_Z_10302025.root
./Execute --Input $inputFile \
        --Variable "Z"\
        --Label "TPCHits","Matching","Reweighting","MCBinning","Regularization","MCStat","Total"\
        --Prefix "10302025"