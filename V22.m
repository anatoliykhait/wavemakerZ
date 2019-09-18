function V22 = V22(w0, w1, w2, k0, k1, k2)
V22 = 2 * (V012(w0,w1,w2,k0,k1) - V012(w0,w2,w1,-k0,k2) - V012(w1,w2,w0,-k1,k2));