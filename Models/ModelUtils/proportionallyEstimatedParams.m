function [l, lf, mf, ank, linkmass, b, c] = proportionallyEstimatedParams(m, h)

    l = h*0.575;
    lf = 0.152 * h;
    mf = 2*0.0145*m;
    ank = 0.19*lf;

    b = 0.039*h;
    c = 0.5*lf-ank;
    
    linkmass = m - mf;

end