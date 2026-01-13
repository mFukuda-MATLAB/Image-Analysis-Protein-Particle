%%=========================================================================
%
% PDV mask image:
%           The TRITC image of the red signals (nanoparticles fluorescence) in a microparticle
%           is subtracted from FITC image of the green signals (AF488-NHS fluorescence).
%           The segmented region of interest (ROI), containing an isolated cell with a single
%           engulfed microparticle, is manually selected, and the total fluorescence
%           intensity of the PDV mask image within the cell in the ROI is measured.
%
% File inputs: (decomposed TIFF images into separate channels using ImageJ)
%           TRITC images (8 bit tifs with format OVA_##.tif) 
%           FITC images (8 bit tifs with format OVA_##.tif)
%           DAPI images (8 bit tifs with format OVA_##.tif)
%           Phase (Bright Field) images (8 bit tifs with format OVA_##.tif)
%       
% Manual selection of FITC image:
%           The segmented ROI images containing a cell with a single engulfed microparticle
%           and corresponding PDV mask images are stored in a designated folder.
%           The PDV mask images are selected from the designated folder.
%
% File Outputs:
%           1) The corresponding images of overlaid images (bright field and PDV mask images)
%           with the segmented FITC images by PDV mask images are stored in a designated
%           folder.
%           2) Total fluorescence intensity is measured from manually selected images
%           and output as Excel file in a designated folder.
%
%

clc; clear; close all; warning('off','all')

file = dir(['D:/Image Process/Images/OVA/Red' '/*.tif']); % Example path
    
n = numel(file);

ImageCount_S = 0;
ImageCount_A = 0;

for z = 1
    % Reading images (TRITC, FITC, DAPI, and Bright field) from example path
    Img_1 = ['D:/Image Process/Images/OVA/Phase/OVA_' num2str(z,'%d') '.tif'];
    Img_2 = ['D:/Image Process/Images/OVA/Red/OVA_' num2str(z,'%d') '.tif'];
    Img_3 = ['D:/Image Process/Images/OVA/Green/OVA_' num2str(z,'%d') '.tif'];
    Img_4 = ['D:/Image Process/Images/OVA/Blue/OVA_' num2str(z,'%d') '.tif'];

    BF = imread(Img_1);     % Bright field image
    Red = imread(Img_2);        % TRITC signal image 
    Green = imread(Img_3);      % FITC signal image
    Blue = imread(Img_4);       % DAPI signal image

% ======================================================================== 
% 
%                   Cell segmentation
%
% ========================================================================

% Nuclear mask
    Nucleus_BW = im2bw(Blue,0.2);   % threshold for RAW264.7, BMDM, BV2
%     Nucleus_BW = im2bw(Blue,0.1);     % threshold for bEnd.3 cells and Pericytes
    Nucleus_BW_close = imopen(Nucleus_BW, ones(10,10));
    Nucleus_mask = imfill(Nucleus_BW_close,"holes");

% Cell edge mask
    Bright_Field_eq = adapthisteq(BF,'clipLimit',0.01,'Distribution','rayleigh');
    Bright_Field_grad = imgradient(Bright_Field_eq,'prewitt');
    Bright_Field_grad_inv = imcomplement(Bright_Field_grad);

    Segmentation = imimposemin(Bright_Field_grad_inv, Nucleus_mask);
    Cell_segmentation = watershed(Segmentation);
    ws = Cell_segmentation == 0;
    
    Invert_Img = imcomplement(ws);
    label = bwlabel(Invert_Img);    % labeling each ROI
    
    cent_num = regionprops(label,'Centroid');

    % Segmentation Figure (Overlay bright field w/ watershed image)
    figure()
    imshow(imoverlay(BF,ws));
    hold on;
    for i = 1:numel(cent_num)
        text(cent_num(i).Centroid(1),cent_num(i).Centroid(2),num2str(i),'color','r')
    end

% ======================================================================== 
% 
%                   Separation of ROI
%
% ========================================================================

% Red fluorescence mask
    Red_BW = im2bw(Red,0.2);
    Red_BW_filter = medfilt2(Red_BW,[5 5]);     % filtering out small objects (noises)
    Red_BW_mask = imclearborder(Red_BW_filter);

% Green fluorescence mask
    Green_BW = im2bw(Green,0.2);
    Green_BW_filter = medfilt2(Green_BW,[5 5]);     % filtering out small objects (noises)
    Green_BW_mask = imclearborder(Green_BW_filter);

% PDV mask
    PDV_BW_mask = Green_BW_mask - Red_BW_mask;
        % Subtraction (red fluorescence signals from green fluorescence signals)
     
    Green_mask = {};
    Red_mask = {};
    PDV_mask = {};
    Org_Green_Img = {};
    Org_BF_Img = {};
    Select_images = {};
    PDV_img = {};
 
% ROI identification / isolation
    for i = 1:max(max(label))
        [row, col] = find(label==i);
        len = max(row) -min(row);
        breadth = max(col) -min(col);
        target1 = zeros([len breadth]);
        target2 = zeros([len breadth]);
        target3 = zeros([len breadth]);
        target4 = uint8(zeros([len breadth]));
        target5 = uint8(zeros([len breadth]));
        sy = min(col)-1;
        sx = min(row)-1;
    
        for j = 1:size(row,1)
            x = row(j,1)-sx;
            y = col(j,1)-sy;
            target1(x,y) = Red_BW_mask(row(j,1),col(j,1));
            target2(x,y) = Green_BW_mask(row(j,1),col(j,1));
            target3(x,y) = PDV_BW_mask(row(j,1),col(j,1));
            target4(x,y) = Green(row(j,1),col(j,1));
            target5(x,y) = BF(row(j,1),col(j,1));
        end
    
        % Classification of ROI image (presence of red and green fluorescence in cells)
        if nnz(target1) > 0 && nnz(target2) > 0
            Red_mask{end+1} = target1;  % Red fluorescence mask (BW)
            Green_mask{end+1} = target2;  % Green fluorescence mask (BW)
            PDV_mask{end+1} = target3;  % PDV mask (BW)
            Org_Green_Img{end+1} = target4;  % Green fluorescence mask (Original)
            Org_BF_Img{end+1} = target5;  % Bright field mask
        end
    
        % Identification of PDV image (from original FITC signal image)
        for k = 1:length(Green_mask)
            [Row, Col] = size(Green_mask{k});
            for l = 1:Row
                for m = 1:Col
                    if PDV_mask{k}(l,m) > 0
                        PDV_img{k}(l,m) = Org_Green_Img{k}(l,m);
                    else
                        PDV_img{k}(l,m) = 0;
                    end
                end
            end
        end
    end

    single_particle_cell = {};
    single_particle_cell_BF = {};
    image_samples = {};

    % Classification of cells with a single engulfed microparticle in ROI
    for o = 1:length(Red_mask)
        Current_Red = Red_mask{o};
        Current_Red_close = bwareaopen(Current_Red,50);
        Current_Red_dilate = imdilate(Current_Red_close,[strel('line',1,0), strel('line',1,90)]);
        Current_Red_fill = imfill(Current_Red_dilate,'holes');
        Current_Red_erode = imerode(Current_Red_fill,strel('diamond',2));
        label_R_par = bwlabel(Current_Red_erode);
        if max(max(label_R_par)) == 1  % for RAW, BV2, and BMDM
%         if max(max(label_R_par)) >= 1  % for bEnd.3 cells and Pericytes
            single_particle_cell{end+1} = PDV_img{o};
            single_particle_cell_BF{end+1} = Org_BF_Img{o};
        end
    end

    % Creation of overlay images (isolated bright field and isolated PDV images)
    for t = 1:length(single_particle_cell)
        image_samples{end+1} = imoverlay(single_particle_cell_BF{t},single_particle_cell{t});
    end

    % Organization and storing of isolated and classified images (FITC images and Overlaid images)
    for u = 1:length(image_samples)
            imwrite(image_samples{u},fullfile('Image Selection','OVA','Overlay Images',...
                sprintf('img_%03d.tif', u+ImageCount_S)),'WriteMode','append'); % Example path
            imwrite(single_particle_cell{u},fullfile('Image Selection','OVA','FITC images',...
                sprintf('img_%03d.tif', u+ImageCount_S)),'WriteMode','append'); % Example path
    end

    ImageCount_S = ImageCount_S + length(image_samples);
    ImageCount_A = ImageCount_A + length(single_particle_cell);

end

% ======================================================================== 
% 
%                   Intensity Measurement
%
% ========================================================================

% Manual selection criteria
%           Selected ROI includes:
%               1) an isolated single cell with a single engulfed particle
%               2) a clear and visible entire cell edge 

% Selection (multiple selections) of images for analysis (total fluorescence intensity)
[files, path] = uigetfile({'*.tif', 'Image Files (*.tif)'},...
    'Select Image Files', 'MultiSelect', 'on');

Intensity_Single_wPDV = [];

% Calculation of total fluorescence intensity of an isolated cell with a single engulfed microparticle
for v = 1:length(files)
   Images{v} = imread(fullfile(path, files{v}));
   Intensity_Single_wPDV(v) = sum(sum(Images{v}./max(max(Images{v}))));
end

% Organization and storing of overlaid images correspond to the manually selected images into a folder
files_S = dir(fullfile('D:\Image Process\Image Selection\OVA\Overlay Images','*.tif')); % Example path

for a = 1:length(files_S)
    Source_file = files_S(a).name;
    for b = 1:length(files)
        Target_file = files(b);
        if strcmp(Source_file, Target_file)
            copyfile(fullfile('D:\Image Process\Image Selection\OVA\Overlay Images',Source_file),...
                fullfile('D:\Image Process\Image Selection\OVA\Selected',Source_file)); % Example path
        end 
    end
end

% Output fluorescence intensity of an isolated cell with a single engulfed microparticle (from selected images) (Excel file)
T_int = table(Intensity_Single_wPDV');
writetable(T_int,'D:/Protein degradation paper/MATLAB/Intensity_Selected/OVA_PLL_intensity.xlsx','WriteVariableName',false); % Example path


