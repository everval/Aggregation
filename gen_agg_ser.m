function [agg_ser] = gen_agg_ser(T,N,params,t)
%Creates cross-sectional aggregated series for long memory generation.
%Syntax:
%[agg_ser] = gen_agg_ser(T,N,params,t)
%
%Inputs:
%   T - time dimension
%   N - cross-sectional dimension
%   params = [p q] - parameters of the beta distribution
%   t - burn-in initial sample. Default t = 100
%
%This version requires the parfor command from the parallel toolbox for
%speed
%
%J. Eduardo Vera-Valdés
%eduardo@cimat.mx
%This version: February 2016

if nargin < 4
    t = 100;
end
a = params(1,1);
b = params(2,1);
series = zeros(T+t,N);

coefs = betarnd(a,b,N,1);

parfor i = 1:N
    series(:,i) = 1/sqrt(N)*filter(1,[ 1, -sqrt(coefs(i,1)) ] , randn(T+t,1));
end

agg_ser = sum( series, 2);