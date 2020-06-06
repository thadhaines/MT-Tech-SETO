clear; close all; clc; format compact
%% make some global with data that will be updated
global g
g.data.list = zeros(5,4);
g.data.list(3,:) = 3;

%% Use create test class
b = TestClass();
b.dataNdx = 3;

b.getAllData()
b.setData(1:3,[7,8,9])
b.getkData(3)

g.data.list

%% Test class 2
b = TestClass2();
b.dataNdx = 3;
b.getData()
%% Test class 3
clear; close all; clc; format compact
% make some global with data that will be updated
global g
g.data.list = zeros(5,4);
g.data.list(3,:) = 3;

% create test class object and initialize with ndx
a = TestClass3();
a.dataNdx = 3;


a.getData()
a.setData(1:4,[2,3,2,4]);
a.getData(4)
b = a.getData()
