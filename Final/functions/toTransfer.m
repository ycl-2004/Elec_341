function [f] = toTransfer(FV, ZETA, Wn)
    num = FV*Wn^2;
    den = [1 2*ZETA*Wn Wn^2];
    f = tf(num,den);
end