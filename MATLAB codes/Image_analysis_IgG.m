%% Image Analysis
clc; clear; close all; warning('off','all')

file = dir(['D:/Image Process/Images/IgG/Red' '/*.tif']);
n = numel(file);

currentRow = 1;

for z = 1:n
%% Reading images

    Img_1 = ['D:/Image Process/Images/IgG/Phase/IgG_' num2str(z,'%d') '.tif'];
    Img_2 = ['D:/Image Process/Images/IgG/Red/IgG_' num2str(z,'%d') '.tif'];
    Img_3 = ['D:/Image Process/Images/IgG/Green/IgG_' num2str(z,'%d') '.tif'];
    Img_4 = ['D:/Image Process/Images/IgG/Blue/IgG_' num2str(z,'%d') '.tif'];

    BF = imread(Img_1);
    Red = imread(Img_2);
    Green = imread(Img_3);
    Blue = imread(Img_4);

%% Cell segmentation

    Nucleus_BW = im2bw(Blue,0.2);
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
    
    % cent_num = regionprops(label,'Centroid');
    
    % figure()
    % imshow(imoverlay(Test_Img_BF,ws));
    % hold on;
    % for i = 1:numel(cent_num)
    %     text(cent_num(i).Centroid(1),cent_num(i).Centroid(2),num2str(i),'color','r')
    % end

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
    
    for i = 1:max(max(label))
        [row, col] = find(label==i);
        len = max(row) -min(row);
        breadth = max(col) -min(col);
        target1 = zeros([len breadth]);
        target2 = zeros([len breadth]);
        target3 = zeros([len breadth]);
        target4 = uint8(zeros([len breadth]));
        sy = min(col)-1;
        sx = min(row)-1;
    
        for j = 1:size(row,1)
            x = row(j,1)-sx;
            y = col(j,1)-sy;
            target1(x,y) = Red_BW_mask(row(j,1),col(j,1));
            target2(x,y) = Green_BW_mask(row(j,1),col(j,1));
            target3(x,y) = PDV_BW_mask(row(j,1),col(j,1));
            target4(x,y) = Green(row(j,1),col(j,1));
        end
    
        if nnz(target1) > 0 && nnz(target2) > 0
            Red_mask{end+1} = target1;
            Green_mask{end+1} = target2;
            PDV_mask{end+1} = target3;
            Org_Green_Img{end+1} = target4;
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

    for o = 1:length(Red_mask)
        Current_Red = Red_mask{o};
        Current_Red_close = bwareaopen(Current_Red,50);
        Current_Red_dilate = imdilate(Current_Red_close,[strel('line',1,0), strel('line',1,90)]);
        Current_Red_fill = imfill(Current_Red_dilate,'holes');
        Current_Red_erode = imerode(Current_Red_fill,strel('diamond',2));
        label_R_par = bwlabel(Current_Red_erode);
        if max(max(label_R_par)) == 1
            single_particle_cell{end+1} = PDV_img{o};
        end
    end

%% Intensity Measurement
    Single_particle_intensity = [];
    etpy_wPDV = isempty(single_particle_cell);
    
    for v = 1:length(single_particle_cell)
        if etpy_wPDV == 0
            if nnz(single_particle_cell{v}) ~= 0
                Intensity_Single_wPDV = sum(sum(single_particle_cell{v}./max(max(single_particle_cell{v}))));
                Single_particle_intensity(end+1) = Intensity_Single_wPDV;
            end
        end
    end

    range = sprintf('A%d', currentRow);
    
    T_int = table(Single_particle_intensity);

    writetable(T_int,'D:/Protein degradation paper/MATLAB/Intensity/IgG_PDAC_intensity.xlsx','Range',range,'WriteVariableName',false);

    currentRow = currentRow + size(Single_particle_intensity,1);
end