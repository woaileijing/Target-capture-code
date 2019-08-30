close all;

H=imread('.\9.jpg');
I=rgb2gray(H);
%subplot(1,2,1);
figure('NumberTitle','off','Name','原图像')
imshow(H);
title('原图像');
figure('NumberTitle','off','Name','灰度化图')
imshow(I);
imwrite(I,'./灰度化图.jpg');
title('Grayscale');

%小波变化消噪
X=imread('9.jpg');            
 
X=rgb2gray(X);
 
subplot(221);          
 
imshow(X);             
 
title('Grayscale');                  
 
% 生成含噪图像并图示
 
init=2055615866;       
 
randn('seed',init);      
 
X=double(X);
 
% 添加随机噪声
 
XX=X+8*randn(size(X));  
 
subplot(222);             
 
imshow(uint8(XX));              
 
title(' Noisy image ');       
 
%用小波函数coif2对图像XX进行2层
 
% 分解
 
[c,l]=wavedec2(XX,2,'coif3'); 
 
% 设置尺度向量
 
n=[1,2];                  
 
% 设置阈值向量 , 对高频小波系数进行阈值处理
 
p=[10.28,24.08]; 
 
nc = wthcoef2('h',c,l,n,p,'s');

nc = wthcoef2('v',nc,l,n,p,'s');

nc = wthcoef2('d',nc,l,n,p,'s');

 
% 图像的二维小波重构
 
X1=waverec2(nc,l,'coif3');
 
subplot(223);              
 
imshow(uint8(X1));                
 
%colormap(map);            
 
title(' First elimination of noise image '); 
 
%再次对高频小波系数进行阈值处理
 
mc=wthcoef2('v',nc,l,n,p,'s');
 
% 图像的二维小波重构
 
X2=waverec2(mc,l,'coif3');  
 
subplot(224);             
 
imshow(uint8(X2));               
 imwrite(uint8(X2),'./消噪图.jpg');
title('Second elimination of noise image');

%维纳滤波
k=0.01;
blur2=uint8(X2);
bf = fftshift(fft2(blur2));
[M,N]=size(blur2);
H=zeros(M,N);
rad= 5;
length = 2*rad+1;
psf=zeros(length,length);
for i=1:length
    for j=1:length
        
        if((  (i-rad-1)*(i-rad-1) + (j-rad-1)*(j-rad-1) ) < rad*rad )
					H(i,j) =1.0/(rad*rad*3.14159);
                    psf(i,j) =1.0/(rad*rad*3.14159);
        end
        
        
    end
end
 
 
 
H= fftshift(fft2(H));
Hg = conj(H);
spec = Hg.*H;
spectt = spec +k;
can = spec./spectt;
 
jk = can./H;
 Gg=bf.*jk;
 
 
Gg=ifftshift(Gg);
pp=im2uint8(mat2gray(real(ifft2(Gg))));%退化后的图像

figure('NumberTitle','off','Name','维纳滤波');
subplot(2,2,1);
imshow(blur2);
subplot(2,2,2);
imshow(pp);
imwrite(pp,'./维纳滤波图.jpg');

%sober滤波
soberI=edge(pp,'Canny',0.1);

figure('NumberTitle','off','Name','自适应阈值二值化');
level =graythresh(pp);            %得到大津阈值法阈值
J = imbinarize(pp,0.43);              %实现图像二值化
subplot(121),imshow(pp);title('原图像');
subplot(122),imshow(J);title('大津法阈值分割后图像');
imwrite(J,'./阈值分割图.jpg');




%霍夫变换
J_reverse = imcomplement(J);

se=strel('disk',2);
J_reverse=imerode(J_reverse,se);
figure('NumberTitle','off','Name','腐蚀结果');
imshow(J_reverse);

[x, y] = mylineextraction(J_reverse);
% Plot the line in the image
figure('NumberTitle','off','Name','直线检测结果');
 imshow(uint8(J_reverse), [min(min(uint8(J_reverse))) max(max(uint8(J_reverse)))]), hold on
plot([x(1) y(1)], [x(2) y(2)],'LineWidth',2,'Color','blue');
plot(x(1),x(2),'x','LineWidth',2,'Color','red');
plot(y(1),y(2),'x','LineWidth',2,'Color','red');
hold off
fprintf('飞行器质心坐标=(%.1f,%.1f)',((x(1)+x(2))/2),((y(1)+y(2))/2));

function [bp, ep] = mylineextraction(BW)
%   The function extracts the longest line segment from the given binary image
%       Input parameter:
%       BW = A binary image.
%
%       Output parameters:
%       [bp, ep] = beginning and end points of the longest line found
%       in the image.
[n, m] = size(BW);
[H,T,R] = hough(BW);
P  = houghpeaks(H,20,'threshold',ceil(0.2*max(H(:))));
lines= houghlines(BW,T,R,P,'FillGap',3,'MinLength',10);
minLength = 10000;
threshLenth=5;
for k = 1:length(lines)  
    xy = [lines(k).point1; lines(k).point2];
    %if ((((xy(1,1) - xy(2,1))^2 + (xy(1,2) - xy(2,2))^2) > threshLenth)&&(((xy(1,1) - xy(2,1))^2 + (xy(1,2) - xy(2,2))^2) < minLength ))
    if ((((xy(1,1) - xy(2,1))^2 + (xy(1,2) - xy(2,2))^2) < minLength ))
        minLength  = (xy(1,1) - xy(2,1))^2 + (xy(1,2) - xy(2,2))^2;
        bp = xy(1,:);
        ep = xy(2,:);
    end
end 
x=[278.0,310.0,314.0,320.0,339.0,366.5,427.5,439.0,467.0,506.5,523.5]; y=[283.0,314.0,320.0,322.0,344.0,369.5,432.5,442.5,471.5,512,526.5];
values = spcrv([[x(1) x x(end)];[y(1) y y(end)]],3);
plot(values(1,:),values(2,:), 'g')；
end

