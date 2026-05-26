function [COG,Amp,Width,rtrckr_bin]=Retracker1_OCOG(W,tn,l)

% this function retracks a waveform according to OCOG retracking algorithm

% input 
%       W: waveform n x 1 or 1 x n
%      tn: the number of aliased sample for instance 8 for TOPEX
%       l: length of the waveform
% output
%     COG: center of gravity of a waveform
%     Amp: Amplitude of a waveform
%   Width: Width of the OCOG box
%   retrckr_bin: retracking bin

% Mohammad J. Tourian, November 2008
% tourian@gis.uni-stuttgart.de

C=0;
G=0;
A=0;

for i=tn:l
    c=i*W(i)^2;
    C=C+c;
    
    g=W(i)^2;
    G=G+g;
    
    a=W(i)^4;
    A=A+a;
end

COG=C/G;
Amp=sqrt(A/G);
Width=G^2/A;

rtrckr_bin=COG-Width/2;