
function processed_data = QuantumOsc(current, Samplewidth, Samplelength, thickness, filePath,symmetric_magnetic_field, VxxColumn, VyxColumn)
    % 如果VyxColumn 没有 则设置VyxColumn ='';
    % 读取数据
    data = readtable(filePath, 'Delimiter', '\t');

    electron = 1.602e-19;
    hbar = 1.05457e-34 ;
    vK = hbar*2*pi/(electron^2);% 25812.807

    % 计算xx电阻
    data.Rxx = data.(VxxColumn) / current;

    % 计算sheet resistance 方块电阻
    data.Rsq = data.Rxx * Samplewidth / Samplelength;

% %     PPMS Oe单位转化为特斯拉
%     data.Magnet = data.Magnet/10000;

    % 去除重复磁场值
    [unique_magnet, unique_indices] = unique(data.Magnet);
    unique_Rsq = data.Rsq(unique_indices);

% % 保存原始数据备份
% unique_Rsq_raw = unique_Rsq;
% 
% % ---- 去除跳点 ----
% diff_Rsq = [0; abs(diff(unique_Rsq))];
% threshold = 5 * median(diff_Rsq);
% valid_idx = diff_Rsq < threshold;
% unique_Rsq(~valid_idx) = interp1(unique_magnet(valid_idx), unique_Rsq(valid_idx), unique_magnet(~valid_idx), 'linear', 'extrap');
% 
% % ---- 可选：再轻度平滑一次 ----
% unique_Rsq = medfilt1(unique_Rsq, 5);
% figure;
% plot(unique_magnet, unique_Rsq_raw, 'r--'); hold on;
% plot(unique_magnet, unique_Rsq, 'b-');
% legend('Raw', 'Filtered');
% xlabel('Magnet (T)');
% ylabel('R_{sq} (\Omega)');
% title('Raw vs Filtered Rsq');
% % 设置不希望被平滑的范围
% mask_keep_raw = (unique_magnet >= -0.1) & (unique_magnet <= 0.1);
% mask_smooth = ~mask_keep_raw;
% 
% % 分开磁场和电阻数据
% magnet_raw = unique_magnet(mask_keep_raw);
% Rsq_raw = unique_Rsq(mask_keep_raw);
% 
% magnet_smooth = unique_magnet(mask_smooth);
% Rsq_smooth = unique_Rsq(mask_smooth);
% 
% % 平滑外部数据
% span = 0.005;
% Rsq_smooth_smoothed = smooth(magnet_smooth, Rsq_smooth, span, 'loess');
% 
% % 合并平滑后的和未平滑的数据
% magnet_combined = [magnet_raw; magnet_smooth];
% Rsq_combined = [Rsq_raw; Rsq_smooth_smoothed];
% 
% % 排序（因为拼接顺序可能乱了）
% [unique_magnet_sorted, sort_idx] = sort(magnet_combined);
% unique_Rsq_sorted = Rsq_combined(sort_idx);
% 
% % 使用 sorted 版本继续后续处理
% 选择在 -0.1T 到 0.1T 范围内的磁场和 Rsq 数据
%     valid_indices = (unique_magnet >= -0.05) & (unique_magnet <= 0.05);
%     valid_magnet = unique_magnet(valid_indices);
%     valid_Rsq = unique_Rsq(valid_indices);
%     [~, min_index] = min(valid_Rsq);
%     zero_field_offset = valid_magnet(min_index); % 作为校准磁场
%     unique_magnet = unique_magnet - zero_field_offset; % 校准磁场数据
%     disp(['用于校准的磁场值: ', num2str(zero_field_offset)]);
    % --- 插值设置 ---
    interp_resolution = 1e-3;  % 1 mT
    num_interp_points = round(2 * symmetric_magnetic_field / interp_resolution) + 1;

    % --- 构建对称插值磁场向量 ---
    interp_magnet = linspace(-symmetric_magnetic_field, symmetric_magnetic_field, num_interp_points);
    
%     % --- 检查对称性 ---
%     tolerance = 1e-10;
%     if any(abs(interp_magnet + flip(interp_magnet)) > tolerance)
%         error('interp_magnet 不对称！请检查 symmetric_magnetic_field 与 interp_resolution。');
%     end
    
    % --- 排序保证插值单调 ---
    [unique_magnet, sortIdx] = sort(unique_magnet);
    unique_Rsq = unique_Rsq(sortIdx);
    
    % --- 插值（pchip 方法 + 允许外推）---
    interp_Rsq = interp1(unique_magnet, unique_Rsq, interp_magnet, 'pchip', 'extrap');
    
    % --- 对称化处理 ---
    symmetrized_Rsq = (interp_Rsq + flip(interp_Rsq)) / 2;
    
    % --- 最终结果赋值 ---
    Rsq = symmetrized_Rsq;

% 计算MR 
    zero_field_index = abs(interp_magnet) == min(abs(interp_magnet)); % 找到interp_magnet最接近0的索引
    
    % 提取symmetrized_Rsq在磁场为0时的值
    Rsq_at_zero_field = min(abs(Rsq(zero_field_index)));

    %计算MR (%)
    MR = (Rsq - Rsq_at_zero_field)/Rsq_at_zero_field *100;
% 返回处理后的数据
    % Rxx
    rhoxx = Rsq * thickness;
    processed_data.interp_magnet = interp_magnet;
    processed_data.Rsq = Rsq;
    processed_data.rhoxx = rhoxx;
    processed_data.MR = MR;

% corrected_Rsq_invB 相关 这一部分需要去Origin确定拟合参数再回来输出；
    % 选择 1 到 13T 的数据
    selected_indices = (interp_magnet >= 1) & (interp_magnet <= symmetric_magnetic_field);% symmetric_magnetic_field
    select_magnet = interp_magnet(selected_indices);
    selected_Rsq = Rsq(selected_indices);

    % 使用 loess 进行平滑 减去背底
    span = 0.3; % 可调整的平滑参数 这个参数能把强的振荡基本平滑掉
    selected_Rsq_loess = smooth(select_magnet, selected_Rsq, span, 'loess');
    selected_Rsq_loess = selected_Rsq_loess';
    % 计算差值
    corrected_Rsq = selected_Rsq - selected_Rsq_loess;

%     % 多项式拟合
%     poly_degree = 4;  % 你可以调整多项式的阶数
%     poly_coefficients = polyfit(half_magnet, selected_smoothed_Rsq, poly_degree);
%     % 计算多项式拟合结果
%     fit_result = polyval(poly_coefficients, half_magnet);  
%     % 减去拟合结果
%     corrected_Rsq = selected_Rsq - fit_result;
%     corrected_smoothed_Rsq = selected_smoothed_Rsq - fit_result;
    % 1/B 的计算
    invB = 1 ./ select_magnet;
    
% 线性插值 与 FFT
    num_interp_points = 100000;
    invB_interp = linspace(min(invB), max(invB), num_interp_points);
    corrected_Rsq_invB = interp1(invB, corrected_Rsq, invB_interp, 'linear');
    
    % FFT
    spacing = (max(invB) - min(invB)) / (num_interp_points - 1);  % 插值点间隔
    N = length(corrected_Rsq_invB);  % 数据长度
    y_fft = fft(corrected_Rsq_invB);  % FFT 计算
    
    % 计算频率轴 (非整数频率)
    freq_axis = (0:N-1) / (N * spacing);  % 频率轴，包含正负频率
    
    % 只保留正频率部分（FFT 对称性）
    half_index = ceil(N / 2);
    x_fft = freq_axis(1:half_index);
    y_fft = y_fft(1:half_index);  % 取正频率部分的复数结果
    
    % 归一化处理
    y_fft = abs(y_fft) / N;  % 对振幅进行归一化
    y_fft(2:end) = y_fft(2:end) * 2;  % 非直流分量乘以 2
    
    % 截取 0 ~ 350 Hz 的频率范围
    freq_limit = 250;  % 设置频率限制
    valid_indices = (x_fft >= 0) & (x_fft <= freq_limit);  % 筛选有效频率范围
    
    % 使用筛选后的频率和幅度
    x_fft_limited = x_fft(valid_indices);
    y_fft_limited = y_fft(valid_indices);  % 对应的幅度

    % === 自动寻峰部分 ===
    % 寻找峰值（只考虑明显的波峰）
    [min_prominence, max_num_peaks] = deal(0.01, 9); % 设置最小显著性和最多峰数·
    [peak_values, peak_locs] = findpeaks(y_fft_limited, x_fft_limited, ...
        'MinPeakProminence', min_prominence, ...
        'SortStr', 'descend', ...
        'NPeaks', max_num_peaks);
    
    % 保存主峰信息
    processed_data.PeakFrequencies = peak_locs;
    processed_data.PeakAmplitudes = peak_values;
        
%     % 绘制图形
%     figure;
%     plot(x_fft_limited, y_fft_limited, '.-');  % 使用点线图展示结果
%     xlabel('Frequency (Hz)');
%     ylabel('Normalized Amplitude');
%     title('Normalized FFT of the Corrected Signal (0 to 350 Hz)');
    
    % 返回结果
    processed_data.Frequency = x_fft_limited;
    processed_data.Amplitude = y_fft_limited;

    processed_data.select_magnet = select_magnet;
    processed_data.corrected_Rsq = corrected_Rsq;

    processed_data.invB = invB_interp;
    processed_data.corrected_Rsq_invB = corrected_Rsq_invB;

    % -----------------------------
    % Hall 相关分析部分
    % -----------------------------
    if ~isempty(VyxColumn)
        % 计算 Hall 电阻
        data.Ryx = data.(VyxColumn) / current;
        unique_Ryx = data.Ryx(unique_indices);
    
        % 平滑与插值
        span = 0.01; % loess 平滑参数
        unique_Ryx = smooth(unique_magnet, unique_Ryx, span, 'loess');
        interp_Ryx = interp1(unique_magnet, unique_Ryx, interp_magnet);
        interp_Ryx = fillmissing(interp_Ryx, 'next');
        interp_Ryx = fillmissing(interp_Ryx, 'previous');
    
        % 反对称化（去除接触不对称影响）
        Antisymmetrized_Ryx = 0.5 * (interp_Ryx - flip(interp_Ryx));
        Ryx = Antisymmetrized_Ryx;
    
        % -----------------------------
        % 三个区间：0~0.1T, 1~3T, 7~13T
        % -----------------------------
        field_ranges = [0, 0.1; 1, 3; 7, 13];
        range_names = {'low', 'mid', 'high'};
    
        electron = 1.602176634e-19; % 元电荷 (C)
    
        for i = 1:3
            range_min = field_ranges(i,1);
            range_max = field_ranges(i,2);
            range_label = range_names{i};
    
            % 选择该区间的数据
            sel_idx = (interp_magnet >= range_min) & (interp_magnet <= range_max);
            B_sel = interp_magnet(sel_idx);
            Ryx_sel = Ryx(sel_idx);
    
            if numel(B_sel) >= 5  % 至少 5 个点再拟合
                % 线性拟合 Ryx = RH_2D * B + offset
                p = polyfit(B_sel, Ryx_sel, 1);
                RH_2D = p(1);
                RH_3D = RH_2D * thickness;
    
                % 载流子浓度
                CarrierDensity_2D = 1 / (RH_2D * electron);
                CarrierDensity_3D = 1 / (RH_3D * electron);
    
                % 迁移率
                mu = abs(RH_2D / Rsq_at_zero_field);
    
                % 保存结果（带区间标签）
                processed_data.(['RH_2D_' range_label]) = RH_2D;
                processed_data.(['RH_3D_' range_label]) = RH_3D;
                processed_data.(['CarrierDensity_2D_' range_label]) = CarrierDensity_2D;
                processed_data.(['CarrierDensity_3D_' range_label]) = CarrierDensity_3D;
                processed_data.(['mu_' range_label]) = mu;
            else
                warning('磁场区间 %.2f~%.2f T 数据点过少，跳过拟合。', range_min, range_max);
            end
        end
    
        % -----------------------------
        % 电阻率与电导张量部分
        % -----------------------------
        rhoyx = Ryx * thickness;
        processed_data.Ryx = Ryx;
        processed_data.rhoyx = rhoyx;
    
        % 电导计算
        Gxx = Rsq ./ (Rsq.^2 + Ryx.^2);
        Gyx = Ryx ./ (Rsq.^2 + Ryx.^2);
        processed_data.Gxx = Gxx;
        processed_data.Gyx = Gyx * vK; % 若不需要量子化，去掉 vK
    
        sigmaxx = rhoxx ./ (rhoxx.^2 + rhoyx.^2);
        sigmayx = rhoyx ./ (rhoxx.^2 + rhoyx.^2);
        processed_data.sigmaxx = sigmaxx;
        processed_data.sigmayx = sigmayx;
    end

    % 处理电导率
    Gxx = Rsq ./(Rsq.^2+ Ryx.^2);
    Gyx = Ryx ./(Rsq.^2+ Ryx.^2);
    processed_data.Gxx = Gxx;
    processed_data.Gyx = Gyx*vK; %这里我计算了量子电导，不需要的把vK去掉即可
    sigmaxx = rhoxx ./(rhoxx.^2+ rhoyx.^2);
    sigmayx = rhoyx ./(rhoxx.^2+ rhoyx.^2);
    processed_data.sigmaxx = sigmaxx;
    processed_data.sigmayx = sigmayx;

end
