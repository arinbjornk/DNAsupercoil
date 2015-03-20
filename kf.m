function [ out ] = kf( sigma, max )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out = max/(abs(sigma+0.65)+1);

end

