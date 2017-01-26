% Check if dataset forms a single cluster

clc
clear all
close all

FM         = load(strcat(pwd,'/NETWORK/INT_FDMatrix/33/33_69_red.txt'));
x          = FM(FM~=0);
hist(x)
[x2, n, b] = compute_xpdf(x);




[dip,p_value] = HartigansDipSignifTest(x2,500)

