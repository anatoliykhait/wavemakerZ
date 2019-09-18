function V012 = V012(w0, w1, w2, k0, k1)
global grav;
V012 = 1/(8*pi) * sqrt(grav*w2/(2*w0*w1)) * (k0*k1 + (w0*w1/grav)^2);