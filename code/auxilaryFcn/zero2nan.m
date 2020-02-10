function mask = zero2nan(mask)

mask = 1*mask;
mask(~mask) = NaN;