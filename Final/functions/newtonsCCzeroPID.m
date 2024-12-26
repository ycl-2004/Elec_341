function [Zret, PMret] = newtonsCCzeroPID(Wn, Zeta, GH)
    % 主程序
    warning('off', 'MATLAB:colon:nonIntegerIndex');

    % 初始化步长
    WnRes = 0.1; % rad/s
    ZetaRes = 0.01; % pure
    max_iter = 5000; % 最大迭代次数
    no_improvement_count = 0; % 无改进计数
    maxWnRes = 1.0;  % 最大步长
    minWnRes = 0.01; % 最小步长
    maxZetaRes = 0.1;
    minZetaRes = 0.001;

    % 初始零点
    p1 = -Zeta * Wn + 1j * Wn * sqrt(1 - Zeta^2); % 复数极点
    Z0.orig = p1;

    % 初始化相位裕度
    PM.orig = 0;
    Dz = CCzero(Z0.orig);
    if (~test_margin(GH * Dz))
        [~, PM.orig] = margin(GH);
    else
        error('Invalid initial zero.');
    end

    PM_prev = PM.orig; % 前一次的相位裕度
    spiral = 1; % 迭代次数

    while (spiral <= max_iter)
        % 单方向搜索
        [Z0.up, PM.up] = searchTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'up');
        [Z0.down, PM.down] = searchTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'down');
        [Z0.right, PM.right] = searchTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'right');
        [Z0.left, PM.left] = searchTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'left');

        % 组合方向搜索
        [Z0.upRight, PM.upRight] = combinedTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'up', 'right');
        [Z0.upLeft, PM.upLeft] = combinedTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'up', 'left');
        [Z0.downRight, PM.downRight] = combinedTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'down', 'right');
        [Z0.downLeft, PM.downLeft] = combinedTrace(Z0.orig, GH, PM.orig, WnRes, ZetaRes, 'down', 'left');

        % 汇总所有方向
        PM_arr = [PM.orig, PM.up, PM.down, PM.right, PM.left, ...
                  PM.upRight, PM.upLeft, PM.downRight, PM.downLeft];
        [~, dir] = max(PM_arr(~isinf(PM_arr))); % 找到最优方向

        % 更新零点位置
        switch dir
            case 2, Z0.orig = Z0.up; PM.orig = PM.up;
            case 3, Z0.orig = Z0.down; PM.orig = PM.down;
            case 4, Z0.orig = Z0.right; PM.orig = PM.right;
            case 5, Z0.orig = Z0.left; PM.orig = PM.left;
            case 6, Z0.orig = Z0.upRight; PM.orig = PM.upRight;
            case 7, Z0.orig = Z0.upLeft; PM.orig = PM.upLeft;
            case 8, Z0.orig = Z0.downRight; PM.orig = PM.downRight;
            case 9, Z0.orig = Z0.downLeft; PM.orig = PM.downLeft;
        end

        Zret = Z0.orig;
        PMret = PM.orig;

        % 动态调整步长和计数
        if dir == 1 % 没有改进
            no_improvement_count = no_improvement_count + 1;
            if no_improvement_count > 10
                disp(['No significant improvement after 50 attempts. Stopping with Zero: ', num2str(Zret)]);
                break;
            end
        else
            no_improvement_count = 0; % 重置计数
        end
        
        % 动态调整步长

        if dir == 1
            WnRes = min(maxWnRes, WnRes * 2); % 增大步长
            ZetaRes = min(maxZetaRes, ZetaRes * 2);
            no_improvement_count = no_improvement_count + 1; % 增加无改进计数
        else
            WnRes = max(minWnRes, WnRes / 2); % 减小步长
            ZetaRes = max(minZetaRes, ZetaRes / 2);
            no_improvement_count = 0; % 重置无改进计数
        end


        disp(['Iteration ', num2str(spiral), ': Best Zero: ', num2str(Zret), ', Best Phase Margin: ', num2str(PMret)]);
        spiral = spiral + 1;
    end

    % 输出最终结果
    disp(['Optimization completed. Final Zero: ', num2str(Zret)]);
    disp(['Final Phase Margin: ', num2str(PMret)]);
    warning('on', 'MATLAB:colon:nonIntegerIndex');
end

%% 辅助函数：方向搜索
function [bestZ0, bestPM] = searchTrace(Z0, GH, PM0, WnRes, ZetaRes, dir)
    bestZ0 = Z0; % 初始化最佳零点
    bestPM = PM0; % 初始化最佳相位裕度
    stp = 0; % 用于计数成功计算 margin 的次数
    c = 1; % 用于计数连续失败的次数

    while (stp < 2) % 至少计算 2 次 margin，确保结果可靠
        % 按固定步长调整零点
        switch dir
            case 'up', Z0_new = Z0 + WnRes * 1j;
            case 'down', Z0_new = Z0 - WnRes * 1j;
            case 'left', Z0_new = Z0 - ZetaRes;
            case 'right', Z0_new = Z0 + ZetaRes;
            otherwise, error('Invalid direction');
        end

        % 计算新的零点模块传递函数
        Dz = CCzero(Z0_new);

        % 测试 margin 是否可以计算
        if (~test_margin(GH * Dz))
            % margin 可以计算，更新最佳值
            stp = stp + 1; % 成功次数 +1
            [~, Phm] = margin(GH * Dz); % 计算相位裕度

            if (Phm > bestPM) % 如果相位裕度更好，更新最佳值
                bestPM = Phm;
                bestZ0 = Z0_new;
            end
        else
            % margin 失败，增加失败计数
            c = c + 1;
            if (c > 50) % 连续失败超过 50 次，跳出循环并警告
                warning(['Warning: margin failed in direction "', dir, '". Skipping this direction.']);
                return; % 退出当前方向搜索，保留初始值
            end
        end
    end
end

%% 辅助函数：组合方向搜索
function [bestZ0, bestPM] = combinedTrace(Z0, GH, PM0, WnRes, ZetaRes, dir1, dir2)
    bestZ0 = Z0; % 初始化最佳零点
    bestPM = PM0; % 初始化最佳相位裕度
    stp = 0; c = 1; % 成功次数和失败计数

    while (stp < 2) % 至少尝试两次有效计算
        Z0_new = Z0;

        % 按组合方向调整零点
        switch dir1
            case 'up', Z0_new = Z0_new + WnRes * 1j;
            case 'down', Z0_new = Z0_new - WnRes * 1j;
        end

        switch dir2
            case 'right', Z0_new = Z0_new + ZetaRes;
            case 'left', Z0_new = Z0_new - ZetaRes;
        end

        % 测试 margin
        Dz = CCzero(Z0_new);
        if (~test_margin(GH * Dz))
            stp = stp + 1; % 成功计数
            [~, Phm] = margin(GH * Dz);
            if (Phm > bestPM) % 更新最佳值
                bestPM = Phm;
                bestZ0 = Z0_new;
            end
        else
            c = c + 1; % 增加失败计数
            if (c > 50)
                warning(['Warning: margin failed in combination direction "', dir1, '-', dir2, '". Skipping.']);
                return; % 跳过当前方向
            end
        end
    end
end

%% 辅助函数：测试 margin 是否可用
function C = test_margin(sys)
    C = 0;
    try
        allmargin(sys);
    catch
        C = 1;
    end
end

%% 辅助函数：生成复数零点传递函数
function Dzero = CCzero(Z)
    s = tf('s');
    Z_re = real(Z);
    Z_im = imag(Z);
    Z1 = Z_re + Z_im * 1j;
    Z2 = Z_re - Z_im * 1j;
    Dzero = 1 / (Z1 * Z2) * (s - Z1) * (s - Z2);
end
