clear; clc; close all; format compact;
%% Matrix math speed test
nIterations = 5; % number of each test to perform
mSize = 5000;   % size of NxN matrix

% anon function to display results
dispRes = @(timeS, nIterations, maxtime, mintime) ...
fprintf('\n meantime \t%2.6f\n maxtime \t%2.6f\n mintime \t%2.6f\n', ...
    timeS/nIterations, maxtime, mintime );

%%
fprintf('\nElement by Element Multiplication')
% init counters
timeS = 0;
maxtime = 0;
mintime = 10000;

for n = 1:nIterations
    x = rand(mSize,mSize);
    y = rand(mSize,mSize);
    Z = zeros(mSize,mSize);
    tic
    Z = x.*y;
    tmptime = toc;
    timeS = timeS + tmptime;
    maxtime = max(maxtime, tmptime);
    mintime = min(mintime, tmptime);
    fprintf('.')
    if mod(n,20) == 0
        fprintf('\n')
    end
    clear x y Z
end
dispRes(timeS, nIterations, maxtime, mintime)

%%
fprintf('\nTrue Matrix Multiplication')
% init counters
timeS = 0;
maxtime = 0;
mintime = 10000;

for n = 1:nIterations
    x = rand(mSize,mSize);
    y = rand(mSize,mSize);
    Z = zeros(mSize,mSize);
    tic
    Z = x*y;
    tmptime = toc;
    timeS = timeS + tmptime;
    maxtime = max(maxtime, tmptime);
    mintime = min(mintime, tmptime);
    fprintf('.')
    if mod(n,20) == 0
        fprintf('\n')
    end
    clear x y Z
end
dispRes(timeS, nIterations, maxtime, mintime)

%%
fprintf('\nLower Triangle Multiplication')
% init counters
timeS = 0;
maxtime = 0;
mintime = 10000;

for n = 1:nIterations
    x = tril(rand(mSize,mSize));
    y = tril(rand(mSize,mSize));
    Z = zeros(mSize,mSize);
    tic
    Z = x*y;
    tmptime = toc;
    timeS = timeS + tmptime;
    maxtime = max(maxtime, tmptime);
    mintime = min(mintime, tmptime);
    fprintf('.')
    if mod(n,20) == 0
        fprintf('\n')
    end
    clear x y Z
end
dispRes(timeS, nIterations, maxtime, mintime)

%%
fprintf('\nDiagonal Matrix Multiplication')
% init counters
timeS = 0;
maxtime = 0;
mintime = 10000;

for n = 1:nIterations
    x = diag(diag(rand(mSize,mSize)));
    y = diag(diag(rand(mSize,mSize)));
    Z = zeros(mSize,mSize);
    tic
    Z = x*y;
    tmptime = toc;
    timeS = timeS + tmptime;
    maxtime = max(maxtime, tmptime);
    mintime = min(mintime, tmptime);
    fprintf('.')
    if mod(n,20) == 0
        fprintf('\n')
    end
    clear x y Z
end
dispRes(timeS, nIterations, maxtime, mintime)

%% Results
%{ 
MATLAB R2017a
Element by Element Multiplication.....
 meantime 	0.035887
 maxtime 	0.042510
 mintime 	0.033445

True Matrix Multiplication.....
 meantime 	1.801634
 maxtime 	1.857192
 mintime 	1.749742

Lower Triangle Multiplication.....
 meantime 	1.722980
 maxtime 	1.782160
 mintime 	1.668384

Diagonal Matrix Multiplication.....
 meantime 	1.774029
 maxtime 	1.821005
 mintime 	1.671652

OCTAVE 5.2.0
Element by Element Multiplication.....
 meantime       0.091087
 maxtime        0.093216
 mintime        0.089758

True Matrix Multiplication.....
 meantime       1.330038
 maxtime        1.354000
 mintime        1.308069

Lower Triangle Multiplication.....
 meantime       1.389301
 maxtime        1.426544
 mintime        1.342429

Diagonal Matrix Multiplication.....
 meantime       0.017534
 maxtime        0.019450
 mintime        0.015193

Discussion:
Octave must have some optimization for diagnoal matricies as results are
faster than element by element multiplication.
Multiple order of magnitude speed up comparted to MATLAB computation for
diagonal matricies

TODO:
Explore sparse matrix operations and possible speed ups.
Investigate possiblity of diagonalizing... 
    'Easily' possible with real and distinct eighen vectors in NxN matrix.
    Speed up vs vectorization?
    Most useful in generators...
%}