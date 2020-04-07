clc; close all

%% example phytree
b = [1 2 ; 3 4 ];
d = [1; 2; 1.5; 1; 0];
n = {'A','B','C','D','E'};

% b = [1 2 ; 4 5; 3 7; 6 8];
% d = [1; 2; 1.5; 1; 1.5; 2; 1; 1; 0];
% n = {'A','B','C','D','E','F','G','H','I'};

tre=phytree(b,d,n);
plot(tre)

treStr = getnewickstr(tre,'distances',true,'branchnames',true)


%% equivalent digraph
B = [ 5 3; 5 4; 4 1; 4 2];
D = [ 1.5;   1;   1;   2];

% B = [7 4; 7 5; 8 7; 8 3; 9 8; 9 6; 6 1; 6 2];
% D = [ 1; 1.5; 1; 1.5; 1; 2; 1; 2];

dg = digraph(B(:,1),B(:,2),D,n)
plot(dg)


%% newick from digraph

[dgStr, subDgStr] = bifurcatingDigraphToNewick(dg);

all(treStr(1:end-1) == dgStr(1:end-3))

%% equal?

dgTre = phytreeread(dgStr);

plot(dgTre)