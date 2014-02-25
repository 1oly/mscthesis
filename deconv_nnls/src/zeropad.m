function Bpad = zeropad(B, m, n, type)
%ZEROPAD Pad array with zeros to avoid wrap-around effects in decovolution
%
%function Bpad = zeropad(B, m, n, type)
%
%      Bpad = zeropad(B, m, n, type);
%
%  Pad B with zeros to make it an m-by-n array. 
%  Requires Image Processing Toolbox.
%
%  Input:
%     B     Array to be zero-padded
%     m, n  Desired dimension of padded array.
%     Type  String 'pad' or 'unpad'.
%
%  Output:
%        Bpad  Padded m-by-n array.
%
% Reference: See Chapter 4,
%            "Deblurring Images - Matrices, Spectra, and Filtering"
%            by P. C. Hansen, J. G. Nagy, and D. P. O'Leary,
%            SIAM, Philadelphia, 2006.
%
% Oliver Lylloff
% Date: 6/11/2013
% Revision: 28/11/2013
%

[M,N] = size(B);
switch type
    case {'pad'}
        s = [ceil((m-M)/2),ceil((n-N)/2)];
        Bpad = padarray(B,s);
    case {'unpad'}
        s = [floor(m/2),floor(n/2)];
        Bpad = B(1+s(1):end-s(1),1+s(2):end-s(2));
end
end