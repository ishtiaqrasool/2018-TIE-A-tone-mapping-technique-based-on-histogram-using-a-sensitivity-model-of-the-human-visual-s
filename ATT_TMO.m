function ldr = ATT_TMO(hdr)

%------------------------------------
% ATT stands for Adaptive TVI Tone-mapping. The function tone-maps a high
% dynamic range image (hdr) to a low dynamic range image (ldr).

% This is implementation of the algorithm presented in
% I. R. Khan, S. Rahardja, M. M. Khan, M. Movania, F. Abed, “A new
% tone-mapping technique based on histogram using a sensitivity model of
% the human visual system,” IEEE Transactions on Industrial Electronics,
% 65(4), 3469-3479, 2018.
% Copyright (C) 2018 Ishtiaq Rasool Khan (irkhan@uj.edu.sa)
% Use of this TMO and the following implementation software is permitted
% for non-commercial research purposes free of charge. For commercial use
% contact irkhan@uj.edu.sa
%-------------------------------------


[x, y] = design_tmo(hdr);
ldr = apply_tmo(hdr, x, y, 1/1.5);

end

%-----------------------------------
function [x, y] = design_tmo(hdr)
%------------------------------------
% Input parameter in this function is 3-channel HDR image in RGB
% format. Output parameters [x, y] give an LUT that defines the TMO.

nbins = [70; 90];
scale_factor = 179;
hdr = scale_factor * hdr;
hdr_lum = 1.0*(0.2126.*hdr(:,:,1)+0.7152.*hdr(:,:,2)+0.0722.*hdr(:,:,3));
[~,midx] = min(hdr_lum); hdr_lum(midx) = min(hdr(:));
[~,midx] = max(hdr_lum); hdr_lum(midx) = max(hdr(:));
hdr_lum_log = log10(hdr_lum);

%
% calculate bin centers equi-spaced in terms of JND
%
n = 16;
while(1)
    jnd = [];
    bins = min(hdr_lum(:));
    while (bins(end) < max(hdr_lum(:)))
        jnd = [jnd; 10 ^ (tvi (log10(bins(end))))];
        bins = [bins; bins(end) + n * jnd(end)];
    end
    bins(end) = max(hdr_lum(:));
    if(numel(bins) < nbins(1))
        n = n - 2;
        if n==0; n = 0.1; end
    elseif(numel(bins) > nbins(2))
        n = n + 1.5;
    else break;
    end
end

%
% jnd measured above is based on bin edge. here we adjust it based on bin
% centers for additional accuracy.
%
centers = 0.5*(bins(2:end)+bins(1:end-1));
centers = log10(centers);
jnd = 10.^tvi(centers);
counts = hist(hdr_lum_log(:), centers)';

%
% refined counts ignoring indistinguishable pixels in each bin
%
rcounts = zeros(size(centers));
for i = 1 : numel(rcounts)
    j = find(hdr_lum >= bins(i) & hdr_lum < bins(i+1));
    pix = sort(hdr_lum(j));
    
    while(1)
        diff = pix(2:2:end) - pix(1:2:end-1);
        k = find(diff < jnd(i));
        if(isempty(k))break; end
        pix(2*k) = [];
    end
    
    rcounts(i) = max(size(pix));
end

%
% Refine the histogram
%
w1 = 0.8;
counts = counts / sum(counts);
rcounts = rcounts / sum (rcounts);
counts = w1*counts + (1-w1)*rcounts;
%figure, plot(counts);

%
% generate cumulative histogram
%
w = [0; counts/max(counts)];
for i = 2 : max(size(w))
    w(i) = w(i) + w(i-1);
end

%
% The lookup table (LUT)
%
x = double(bins)/scale_factor;
y = double(w) / max(w);
%figure, plot(x, y);

end

%------------------------------------------
function out = apply_tmo(in, x, y, gamma)
%------------------------------------------
% input parameter "in" is HDR image in RGB format. [x, y] are the LUT that
% define the TMO. The last optional parameter is for gamma correction. If
% it is defined, then HDR luminance is tone-mapped and then a ratio image
% is used to transform RGB images and gamma corretion is also applied. If
% the last parameter is not defined, then RGB values of each channel are
% tone-mapped independently. We recoomend the first option and gamma =
% 1/1.5.


if (nargin<4)
    % Each channel tone-mapped independently
    out = interp1(x, y ,in, 'linear', 'extrap');
else
    % Luminance channel tone-mapped first and then same scaling applied to each
    % of RGB channels.
    in_lum = 1.0*(0.2126.*in(:,:,1)+0.7152.*in(:,:,2)+0.0722.*in(:,:,3));
    out_lum = interp1(x, y, in_lum, 'linear', 'extrap');
    out = zeros(size(in));
    out(:, :, 1) = (in(:, :, 1) ./ (in_lum)) .^ gamma .* (out_lum) ;
    out(:, :, 2) = (in(:, :, 2) ./ (in_lum)) .^ gamma .* (out_lum) ;
    out(:, :, 3) = (in(:, :, 3) ./ (in_lum)) .^ gamma .* (out_lum) ;
    out(out > 1) = 1;
end
end

%---------------------------------------
function threshold = tvi(intensity)
%---------------------------------------
% The output is the sensitivity of HVS measured as just noticeable
% difference in the ambient intensity defined by the input parameter.

threshold = zeros(size(intensity));

idx = find(intensity < -3.94);
threshold (idx) = -2.86;

idx = find(intensity >= -3.94 & intensity < -1.44);
threshold (idx) = (0.405 * intensity(idx) + 1.6) .^ 2.18 - 2.86;

idx = find(intensity >= -1.44 & intensity < -0.0184);
threshold (idx) = intensity(idx) - 0.395;

idx = find(intensity >= -0.0184 & intensity < 1.9);
threshold(idx) = (0.249 * intensity(idx) + 0.65) .^ 2.7 - 0.72;

idx = find(intensity >= 1.9);
threshold (idx) = intensity(idx) - 1.255;

threshold = threshold - 0.95;
end