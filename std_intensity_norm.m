function std_intensity_norm(input_dir)
    nii_files = dir(fullfile(input_dir, '**', 'w_realigned.nii'));

    pons_path = 'reference_regions/wfu_pons.nii';
    if exist(pons_path, 'file') ~= 2
        wfu_path = 'reference_regions/TD_lobe.nii';
        wfu_header = spm_vol(wfu_path);
        wfu = spm_read_vols(wfu_header);
        pons = (wfu == 7);
        pons_header = wfu_header;
        pons_header.fname = pons_path;
        pons_header.dt(1) = 2;
        spm_write_vol(pons_header, pons);
    end
    gm_path = 'reference_regions/GM.nii';

    output_pdf = fullfile(input_dir, 'intensity_QC.pdf');
    if exist(output_pdf, 'file')
        delete(output_pdf);
    end

    for i = 1:length(nii_files)
        file_path = fullfile(nii_files(i).folder, nii_files(i).name);

        if ~exist(file_path, 'file')
            warning('File not found: %s. Skipping.', file_path);
            continue;
        end

        fprintf('📈 Intensity Normalization... %s\n', file_path);

        pet_hdr = spm_vol(file_path);
        pet_vol = spm_read_vols(pet_hdr);

        copied_pons_path = fullfile(nii_files(i).folder, 'wfu_pons.nii');
        copyfile(pons_path, copied_pons_path);
        copied_pons_hdr = spm_vol(copied_pons_path);
        P = char(pet_hdr.fname, copied_pons_hdr.fname);
        flags = struct('interp', 0, 'wrap', [0 0 0], 'mask', 0, 'which', 1, 'mean', 0);
        spm_reslice(P, flags);
        delete(copied_pons_path);

        resliced_pons_path = fullfile(nii_files(i).folder, 'rwfu_pons.nii');
        pons_header = spm_vol(resliced_pons_path);
        pons = spm_read_vols(pons_header);
        pons_mask = (pons == 1);
        
	pons_mask_pre = (pons == 1);
	se = strel('sphere', 2);
	pons_mask = imerode(pons_mask_pre, se);
	fprintf('Pons voxels before: %d, after: %d\n', nnz(pons_mask_pre), nnz(pons_mask));

        pons_mean = mean(pet_vol(pons_mask), 'omitnan');

        if pons_mean == 0 || isnan(pons_mean)
            warning('⚠️ Pons mean is zero or NaN in %s. Skipping pons normalization.', file_path);
        else
            norm_pons = pet_vol / pons_mean;
            pons_hdr = pet_hdr;
            [~, name, ext] = fileparts(nii_files(i).name);
            pons_hdr.fname = fullfile(nii_files(i).folder, ['pons_' name ext]);
            spm_write_vol(pons_hdr, norm_pons);
        end

        copied_gm_path = fullfile(nii_files(i).folder, 'GM.nii');
        copyfile(gm_path, copied_gm_path);
        copied_gm_hdr = spm_vol(copied_gm_path);
        P = char(pet_hdr.fname, copied_gm_hdr.fname);
        flags = struct('interp', 0, 'wrap', [0 0 0], 'mask', 0, 'which', 1, 'mean', 0);
        spm_reslice(P, flags);
        delete(copied_gm_path);

        resliced_gm_path = fullfile(nii_files(i).folder, 'rGM.nii');
        gm_header = spm_vol(resliced_gm_path);
        gm = spm_read_vols(gm_header);
        gm_mask = (gm == 1);
        gm_mean = mean(pet_vol(gm_mask), 'omitnan');

        if gm_mean == 0 || isnan(gm_mean)
            warning('⚠️ GM mean is zero or NaN in %s. Skipping GM normalization.', file_path);
        else
            norm_gm = pet_vol / gm_mean;
            gm_hdr = pet_hdr;
            [~, name, ext] = fileparts(nii_files(i).name);
            gm_hdr.fname = fullfile(nii_files(i).folder, ['gm_' name ext]);
            spm_write_vol(gm_hdr, norm_gm);
        end

	% Save images for QC
	% Axial slice through the middle of the pons mask
	[~, ~, pons_z] = ind2sub(size(pons_mask), find(pons_mask));
	pons_mid_z = round(mean(pons_z));     % central Z of the pons
	orig_ax_pons = flipud(squeeze(pet_vol(:, :, pons_mid_z))');
	pons_ax      = flipud(squeeze(pons_mask(:, :, pons_mid_z))');

	% Sagittal slice through the middle of the pons mask
	[pons_x, ~, ~] = ind2sub(size(pons_mask), find(pons_mask));
	pons_mid_x = round(mean(pons_x));     % central X of the pons
	orig_sag_pons = flipud(squeeze(pet_vol(pons_mid_x, :, :))');
	pons_sag      = flipud(squeeze(pons_mask(pons_mid_x, :, :))');

	% Axial slice through the middle of the GM mask
	[~, ~, gm_z] = ind2sub(size(gm_mask), find(gm_mask));
	gm_mid_z = round(mean(gm_z));         % central Z of the GM
	orig_ax_gm = flipud(squeeze(pet_vol(:, :, gm_mid_z))');
	gm_ax      = flipud(squeeze(gm_mask(:, :, gm_mid_z))');

	fig = figure('Name', 'Intensity Norm QC', 'NumberTitle', 'off');
	colormap gray;

	% --- Pons (top row) ---
	subplot(2, 3, 1);
	imagesc(orig_ax_pons); axis image off;
	title('PET slice (pons)');

	subplot(2, 3, 2);
	imagesc(orig_ax_pons); axis image off; hold on;
	red = cat(3, ones(size(pons_ax)), zeros(size(pons_ax)), zeros(size(pons_ax)));
	h = imagesc(red);
	set(h, 'AlphaData', double(pons_ax) * 0.5);
	title('PET + eroded pons (50%)');

	subplot(2, 3, 3);
	imagesc(orig_sag_pons); axis image off; hold on;
	red = cat(3, ones(size(pons_sag)), zeros(size(pons_sag)), zeros(size(pons_sag)));
	h = imagesc(red);
	set(h, 'AlphaData', double(pons_sag) * 0.5);
	title('PET sagittal + pons (50%)');

	% --- GM (bottom row) ---
	subplot(2, 3, 4);
	imagesc(orig_ax_gm); axis image off;
	title('PET slice (GM)');

	subplot(2, 3, 5);
	imagesc(orig_ax_gm); axis image off; hold on;
	red = cat(3, ones(size(gm_ax)), zeros(size(gm_ax)), zeros(size(gm_ax)));
	h = imagesc(red);
	set(h, 'AlphaData', double(gm_ax) * 0.5);
	title('PET + GM (50%)');

	[filepath_parent, ~, ~] = fileparts(file_path);
	[filepath_gdparent, parent_folder] = fileparts(filepath_parent);
	[~, gdparent_folder] = fileparts(filepath_gdparent);
	full_title = sprintf('Intensity Norm QC for: %s/%s/%s', gdparent_folder, parent_folder, nii_files(i).name);
	sgtitle(full_title, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 10);
	exportgraphics(fig, output_pdf, 'Append', true, 'ContentType', 'image');
	close(fig);
    end
end
