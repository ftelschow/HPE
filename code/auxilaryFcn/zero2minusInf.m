function mask = zero2minusInf(mask)

mask = 1*mask;

mask(~mask(:)) = -Inf;