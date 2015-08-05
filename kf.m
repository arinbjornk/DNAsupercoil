function [ out ] = kf( sigma, max )
%function accompanies supercoilingGFPRFP.m

out = max/(abs(sigma+0.65)+1);

end

