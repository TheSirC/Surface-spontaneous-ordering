clear; clc; close all;
working_dictionary = '.';
margins_ = [-2, -3, -0.5, 0.1]; papersize_ = [18, 18 / 2^0.5]; font_size_ = 20; marker_size_ = 6; line_width_ = 0.5;
total_average_num = 2;
total_pulse_num = 2;
all_arrange = [1, 2]; %%% how to arrange the subplot
k_range = 5; %%% k0, the range of wave number to see
E_ab_slice_fft_all_all = [];
rough_height_fft_all_all = [];
for i = 1:total_average_num
    exci = ['cold\Dz_excitation\for_average_run_', num2str(i), '\'];
    E_ab_slice_fft_all = importdata([working_dictionary, exci, 'E_ab_slice_fft_all.mat']);
    rough_height_fft_all = importdata([working_dictionary, exci, 'rough_height_fft_all.mat']);
    E_ab_slice_fft_all_all = [E_ab_slice_fft_all_all; E_ab_slice_fft_all]; %#ok<AGROW>
    rough_height_fft_all_all = [rough_height_fft_all_all; rough_height_fft_all]; %#ok<AGROW>
    F_x = importdata([working_dictionary, exci, 'F_x.mat']);
    F_y = importdata([working_dictionary, exci, 'F_y.mat']);
end
F_xy_mesh = meshgrid(F_x, F_y);
%%% use to defind the colormap for the fft image
F_colormap = 1; %%% k_0
colormap_F = find(abs(F_xy_mesh) > F_colormap);
color_factor = 0.5;
%%% use to defind the colormap for the fft image
%%% do the average %%%
rough_height_fft_all_av = [];
for ii = 1:total_pulse_num;
    temp = 0;
    for i = 1:total_average_num
        temp = temp + abs(rough_height_fft_all_all{i, ii});
    end
    temp = temp / total_average_num;
    rough_height_fft_all_av{ii} = temp; %#ok<SAGROW>
end

E_ab_slice_fft_all_av = [];
for ii = 1:total_pulse_num;
    temp = 0;
    for i = 1:total_average_num
        temp = temp + abs(E_ab_slice_fft_all_all{i, ii});
    end
    temp = temp / total_average_num;
    E_ab_slice_fft_all_av{ii} = temp; %#ok<SAGROW>
end

for i = 1:total_average_num
    figure('name', 'rough_height_fft_all');
    for ii = 1:total_pulse_num
        subaxis(all_arrange(1), all_arrange(2), ii, 'spacing', 0, 'Margin', 0);
        pcolor(F_y, F_x, (abs(rough_height_fft_all_all{i, ii}))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
        caxis([0, color_factor * max(max(abs((rough_height_fft_all_all{i, ii}(colormap_F)))))]);
        ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); xlim([-k_range, k_range]); ylim([-k_range, k_range]);
        set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -k_range:1:k_range, 'YTick', -k_range:1:k_range); axis off;
    end
end

figure('name', 'rough_height_fft_all_av');
for ii = 1:total_pulse_num
    subaxis(all_arrange(1), all_arrange(2), ii, 'spacing', 0, 'Margin', 0);
    pcolor(F_y, F_x, (abs(rough_height_fft_all_av{ii}))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
    caxis([0, color_factor * max(max(abs((rough_height_fft_all_av{ii}(colormap_F)))))]);
    ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); xlim([-k_range, k_range]); ylim([-k_range, k_range]);
    set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -k_range:1:k_range, 'YTick', -k_range:1:k_range); axis off;
end

for i = 1:total_average_num
    figure('name', 'E_ab_slice_fft_all');
    for ii = 1:total_pulse_num
        subaxis(all_arrange(1), all_arrange(2), ii, 'spacing', 0, 'Margin', 0);
        pcolor(F_y, F_x, (abs(E_ab_slice_fft_all_all{i, ii}))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
        caxis([0, color_factor * max(max(abs((E_ab_slice_fft_all_all{i, ii}(colormap_F)))))]);
        ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); xlim([-k_range, k_range]); ylim([-k_range, k_range]);
        set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -k_range:1:k_range, 'YTick', -k_range:1:k_range); axis off;
    end
end

figure('name', 'E_ab_slice_fft_all_av');
for ii = 1:total_pulse_num
    subaxis(all_arrange(1), all_arrange(2), ii, 'spacing', 0, 'Margin', 0);
    pcolor(F_y, F_x, (abs(E_ab_slice_fft_all_av{ii}))); shading interp; colormap(jet(2^10)); xlabel('k_x (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_);
    caxis([0, color_factor * max(max(abs((E_ab_slice_fft_all_av{ii}(colormap_F)))))]);
    ylabel('k_y (k_0)', 'Interpreter', 'tex', 'FontSize', font_size_); xlim([-k_range, k_range]); ylim([-k_range, k_range]);
    set(gca, 'fontsize', font_size_, 'FontAngle', 'italic', 'XTick', -k_range:1:k_range, 'YTick', -k_range:1:k_range); axis off;
end

