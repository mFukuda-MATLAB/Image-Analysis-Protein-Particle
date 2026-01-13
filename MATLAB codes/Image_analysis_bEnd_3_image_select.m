%% Image Analysis
clc; clear; close all; warning('off','all')

file = dir(['D:/Image Process/Images/bEnd.3/Red' '/*.tif']);
n = numel(file);

ImageCount_S = 0;
ImageCount_A = 0;

for z = 1:n
%% Reading images

    Img_1 = ['D:/Image Process/Images/bEnd.3/Phase/bEnd.3_' num2str(z,'%d') '.tif'];
    Img_2 = ['D:/Image Process/Images/bEnd.3/Red/bEnd.3_' num2str(z,'%d') '.tif'];
    Img_3 = ['D:/Image Process/Images/bEnd.3/Green/bEnd.3_' num2str(z,'%d') '.tif'];
    Img_4 = ['D:/Image Process/Images/bEnd.3/Blue/bEnd.3_' num2str(z,'%d') '.tif'];

    BF = imread(Img_1);
    Red = imread(Img_2);
    Green = imread(Img_3);
    Blue = imread(Img_4);

%% Cell segmentation

    Nucleus_BW = im2bw(Blue,0.1);
    Nucleus_BW_close = imopen(Nucleus_BW, ones(10,10));
    Nucleus_mask = imfill(Nucleus_BW_close,"holes");
    
    Bright_Field_eq = adapthisteq(BF,'clipLimit',0.01,'Distribution','rayleigh');
    Bright_Field_grad = imgradient(Bright_Field_eq,'prewitt');
    Bright_Field_grad_inv = imcomplement(Bright_Field_grad);
    
    Segmentation = imimposemin(Bright_Field_grad_inv, Nucleus_mask);
    Cell_segmentation = watershed(Segmentation);
    ws = Cell_segmentation == 0;
    
    Invert_Img = imcomplement(ws);
    label = bwlabel(Invert_Img);
    
%     cent_num = regionprops(label,'Centroid');
%     
%     figure()
%     imshow(imoverlay(BF,ws));
%     hold on;
%     for i = 1:numel(cent_num)
%         text(cent_num(i).Centroid(1),cent_num(i).Centroid(2),num2str(i),'color','r')
%     end

%% Separation of individual cell

    Red_BW = im2bw(Red,0.2);
    Red_BW_filter = medfilt2(Red_BW,[5 5]);
    Red_BW_mask = imclearborder(Red_BW_filter);
    
    Green_BW = im2bw(Green,0.2);
    Green_BW_filter = medfilt2(Green_BW,[5 5]);
    Green_BW_mask = imclearborder(Green_BW_filter);
    
    PDV_BW_mask = Green_BW_mask - Red_BW_mask;
     
    Green_mask = {};
    Red_mask = {};
    PDV_mask = {};
    Org_Green_Img = {};
    Org_BF_Img = {};
    Select_images = {};
    PDV_img = {};
    
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
    
        if nnz(target1) > 0 && nnz(target2) > 0
            Red_mask{end+1} = target1;
            Green_mask{end+1} = target2;
            PDV_mask{end+1} = target3;
            Org_Green_Img{end+1} = target4;
            Org_BF_Img{end+1} = target5;
        end
    
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

    for o = 1:length(Red_mask)
        Current_Red = Red_mask{o};
        Current_Red_close = bwareaopen(Current_Red,50);
        Current_Red_dilate = imdilate(Current_Red_close,[strel('line',1,0), strel('line',1,90)]);
        Current_Red_fill = imfill(Current_Red_dilate,'holes');
        Current_Red_erode = imerode(Current_Red_fill,strel('diamond',2));
        label_R_par = bwlabel(Current_Red_erode);
        if max(max(label_R_par)) >= 1
            single_particle_cell{end+1} = PDV_img{o};
            single_particle_cell_BF{end+1} = Org_BF_Img{o};
        end
    end

    for t = 1:length(single_particle_cell)
        image_samples{end+1} = imoverlay(single_particle_cell_BF{t},single_particle_cell{t});
    end

    for u = 1:length(image_samples)
            imwrite(image_samples{u},fullfile('Image Selection','bEnd.3','Overlay Images',sprintf('img_%03d.tif', u+ImageCount_S)),'WriteMode','append');
            imwrite(single_particle_cell{u},fullfile('Image Selection','bEnd.3','FITC images',sprintf('img_%03d.tif', u+ImageCount_S)),'WriteMode','append');
    end

    ImageCount_S = ImageCount_S + length(image_samples);
    ImageCount_A = ImageCount_A + length(single_particle_cell);

end

%% Intensity Measurement

[files, path] = uigetfile({'*.tif', 'Image Files (*.tif)'}, 'Select Image Files', 'MultiSelect', 'on');

Intensity_Single_wPDV = [];

for v = 1:length(files)
   Images{v} = imread(fullfile(path, files{v}));
   Intensity_Single_wPDV(v) = sum(sum(Images{v}./max(max(Images{v}))));
end

files_S = dir(fullfile('D:\Image Process\Image Selection\bEnd.3\Overlay Images','*.tif'));

for a = 1:length(files_S)
    Source_file = files_S(a).name;
    for b = 1:length(files)
        Target_file = files(b);
        if strcmp(Source_file, Target_file)
            copyfile(fullfile('D:\Image Process\Image Selection\bEnd.3\Overlay Images',Source_file),...
                fullfile('D:\Image Process\Image Selection\bEnd.3\Selected',Source_file));
        end
    end
end

T_int = table(Intensity_Single_wPDV');
writetable(T_int,'D:/Protein degradation paper/MATLAB/Intensity_Selected/OVA_PLL+bEnd.3_intensity.xlsx','WriteVariableName',false);


