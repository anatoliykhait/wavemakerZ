function V11 = V11(w0, w1, w2, k0, k1, k2)
V11 = -2*V012(w0,w1,w2,-k0,k1) + V012(w1,w2,w0,k1,k2);