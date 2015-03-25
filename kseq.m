function [ out ] = kseq( sigma  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s_0 = -0.65;

if(sigma<s_0)
    out = abs(sigma-s_0)/(1+abs(sigma-s_0));
else
out = 0;

end

