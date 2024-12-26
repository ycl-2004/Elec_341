function [zero, PM] = optimizeSingleZero(start,pdy, res)
    if ~exist('res','var')
      res = 0.1; % default resolution
    end
    currPM = getPM(start, pdy);
    z = start;
    done = false;
    haveChecked = [z];
    while ~done
        disp([num2str(z),' PM: ', num2str(currPM)])
        inc = z + res;
        if ~isin(haveChecked, inc)
            haveChecked(end+1) = inc;
            incPM = getPM(inc, pdy);
            if (incPM > currPM)
                z = inc;
                currPM = incPM;
                continue
            end % if
        end % if
        dec = z - res;
        if ~isin(haveChecked, dec)
            haveChecked(end+1) = dec;
            decPM = getPM(dec, pdy);
            if (decPM > currPM)
                z = dec;
                currPM = decPM;
                continue
            end % if
        end % if
        done = true;
    end
    zero = z;
    PM = getPM(zero, pdy);
end

function pm = getPM(zero,pdy)
    s = tf('s');
    Dz = (s-zero)/(-zero);
    openloop = Dz * pdy;
    [bbb,~] = phsMargin(openloop);
    pm = bbb;
end

function [pm wcp] = phsMargin(KDGH)

  % Find pm x-over freq
  [x pm x wcp] = margin(KDGH);            % x is a dummy variable

  % Use bode to re-calculate pm WRT -180 (not -540)
  if ~isnan(wcp)
    [x phase] = bode(KDGH, wcp);          % x is a dummy variable
    pm = phase + 180;
  end

end % function

function out = isin(array, element)
    arr_ismember = ismembertol(array,element);
    out = sum(arr_ismember) > 0;
end
