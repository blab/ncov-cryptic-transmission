# augmentTransTree

augmentTransTree adds additional attributes to an already-simulated transmission tree, formatted as a treeClusterSim linelist with who-infected-whom data

addFactorAttributes.m is used to add categorical attributes like sex, age, vaccination status, etc.  For each factor, you specify the probability of each level.  REMEMBER that the probabilities and levels describe people who have been infected (and are thus in the linelist), and not a random person from the population!

addImportationsField.m finds the root node of each connected component of the transmission tree, from which each infection descends, and adds that to the linelist as the importation label. 

addSpatialAttributes.m overlays a spatial tranmission model on the transmission tree, specified as a discrete hopping model with a weighted connectivity graph. 

(If ever implemented in this toolset) simulateSequencesOnTree.m will simulate a sequence alignment on the tree, either neutral or informed by attributes in the linelist (like beta or infectiousDuration). 


