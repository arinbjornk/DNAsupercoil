function [ out ] = kcat( sigma, max, TL  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out = (max/TL)/(abs(sigma+0.65)+1);

end

