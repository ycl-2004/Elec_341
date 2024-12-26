%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heurRCGTune()
% Heuristic tuning when given RCG values we must satisfy
%
% Syntax:
%   [Kp_best, Ki_best, Kd_best] = heurRCGTune(KIN, Ts, Osu, Ess, p, G, H)
% 
% Input:
%   KIN  :  struct of Kp, Ki, Kd to be tuned
%   Ts  : Rise time target
%   Osu : OSu target
%   Ess : Ess target
%   G     :  fwd path
%   H     :  fb path
%   p    :  derivative pole (Must include)
% Output:
%   List of best Ks : [Kp_best, Ki_best, Kd_best]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KOUT = heurRCGTune(KIN, Ts, Osu, Ess, p, G, H)

    warning('off');

    s = tf('s');

    function [ts_, os_, ess_, tr_, valid] = getStep(kp, ki, kd)
        try
            D = kp + ki/s + kd * -p*s/(s-p);
            cltf = G * D / (1 + G * D * H);
            stp = stepinfo(cltf);
            ts_ = stp.SettlingTime;
            tr_ = stp.RiseTime;
            os_ = stp.Overshoot;
            ess_ = abs(1 - dcgain(cltf));
            valid = true;
        catch
            ts_ = Inf; os_ = Inf; ess_ = Inf; tr_ = Inf; % edge case where we have bad stuff
            valid = false;
        end
    end

    function cost = costFunction(params)
        kp = params(1); ki = params(2); kd = params(3);
        [ts_, os_, ess_, ~, valid] = getStep(kp, ki, kd);
        if ~valid || isnan(ts_) || isnan(os_) || isnan(ess_)
            cost = Inf; % make invalid outputs seem terrible to the algo
        else
            cost = (Ts - ts_)^2 / max(1, Ts^2) + ...
                   (Osu - os_)^2 / max(1, Osu^2) + ...
                   (Ess - ess_)^2 / max(1, Ess^2);
        end
    end

    % init params
    Kp = KIN.Kp * KIN.K;
    Ki = KIN.Ki * KIN.K;
    Kd = KIN.Kd * KIN.K;
    lb = [0, 0, 0]; % lower bounds so no 0 ks
    ub = [10, 10, 10]; % upp bounds so the ks dont blow up

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 1000);
    [best_ks, ~] = fmincon(@costFunction, [Kp, Ki, Kd], [], [], [], [], lb, ub, [], options);

    Kp_best = best_ks(1);
    Ki_best = best_ks(2);
    Kd_best = best_ks(3);
    KOUT = [Kp_best, Ki_best, Kd_best];
    warning('on');

    fprintf('==========================================\n');
    fprintf('Found PID Parameters:\n');
    fprintf('Kp = %.4f\nKi = %.4f\nKd = %.4f\n', Kp_best, Ki_best, Kd_best);
    [ts, os, ess, tr, ~] = getStep(Kp_best, Ki_best, Kd_best);
    fprintf('Resulting in step response: Ts = %.4f, Osu = %.4f, Ess = %.4f | Tr = %.4f\n', ts, os, ess, tr);
end