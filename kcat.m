function [ out ] = kcat( sigma, max, TL  )
%function accompanies supercoilingGFPRFP.m

out = (max/TL)/(abs(sigma+0.65)+1);

end

