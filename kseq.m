function [ out ] = kseq( sigma  )
%function accompanies supercoilingGFPRFP.m

s_0 = -0.65;

if(sigma<s_0)
    out = abs(sigma-s_0)/(1+abs(sigma-s_0));
else
out = 0;

end

