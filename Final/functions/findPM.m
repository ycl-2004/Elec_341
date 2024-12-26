function PM = findPM(K, G, H, target_PM)
    % Define the closed-loop transfer function with gain K
  
    % Calculate OSu in percentage
    %PM = margin;
    [~,PM,~,~] = margin(K*G*H);
    % Return the difference from the target OSu
   PM = PM - target_PM;
end