% demo augmentTransTree functions

addpath(genpath('../utils'))

%% add factor attributes

infection = addFactorAttributes('example/configFactors.yaml');


%% add spatial attributes

infection = addSpatialAttributes('example/configSpatial.yaml');

%% add importation attribute

infection = addImportationsField('example/configSpatial.yaml');

%% add genetic sequences
% long-term to-do