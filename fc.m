function [y,dy]=fc(c,s)
% function that calculate rho_n(c) and its derivative
% Inputs: c- nx1 vector c
%         s- nxn marix s
%         For details see Agarwal R, Chen Z, Sarma SV, A Novel Nonparametric Maximum Likelihood
%         Estimator for Probability Density Functions. IEEE TPAMI 2016. 
% Outputs:
%         y- rho_n(c) value
%        dy- rho_n(s) derivative value
% author: Rahul Agarwal. rahul.jhu@gmail.com
y=s*c-1./c;
dy=s+diag(1./(c.^2));