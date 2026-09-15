# EEC measurement ALEPH


Reminder: When you checkout this repository, you must also update the path for the base directory of RooUnfold in SetupAnalysis.sh. This is essential for being able to run unfolding code associated with this analysis.


## Analysis Steps
1. Create the right tree structure for the input data for the unfolding using fillPair
2. Perform the matching, which forms what will fill in the response matrix, using matching definition 2. 
3. Perform the Unfolding
4. Compute the necessary corrections
5. Compute the Systematic Uncertainties
6. Tie it all together and plot the final result
7. Crosscheck of the full framework
