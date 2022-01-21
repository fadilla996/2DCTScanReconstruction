%% Program Rekonstruksi Data Hasil Proyeksi CT Scan%%
%%Fadilla Sofa Amatullah%%
%%10217012%%
clc
clear all
 
%% Membaca Data Hasil Proyeksi CT Scan
 
% Baca direktori file
N = 361;
img_dir = 'G:\Polyspace\R2020b\bin\Data'
strfile = 'image-000001.dcm'; 
img = dicomread(fullfile(img_dir, strfile));
 
% Membaca ukuran pixel gambar
siz_img = size(img);
 
% Membaca Nilai Pixel/Gray Level dari File Direktori
ct3d = NaN([siz_img N]);
ct3d(:,:,1) = img;    
for ii=2:N
    strfile = sprintf('image-00%04d.dcm',ii);
    ct3d(:,:,ii)= dicomread(fullfile(img_dir, strfile));
end
 
nx = siz_img(1,1);                       
ny = siz_img(1,2);                       
ps = 1;
x = zeros(1,nx);
y = zeros(ny,1);
xx = zeros(nx,ny);
yy = zeros(nx,ny);
x(1,1) = -nx/2;
y(1,1) = ny/2;
 
for i = 2 : nx
    x(1,i)= x(1,1)+(i-1)*ps;
    y(i,1)= y(1,1)-(i-1)*ps;
end
 
for i = 1 : nx
    for j = 1 : ny
        xx(i,j) = x(1,j);
        yy(i,j) = y(i,1);
    end
end
 
% Memilih contoh slice yang akan direkonstruksi
slice = input('Masukkan slice yang akan direkonstruksi : ');
 
%Kalibrasi nilai gray level
m = max(max(ct3d(:,:,slice)))+abs(min(min(ct3d(:,:,slice))));
ct3d(:,:,slice) = ((ct3d(:,:,slice) + abs(min(min(ct3d(:,:,slice)))))/m)*4095;
 
figure (1)
imagesc(x,y,ct3d(:,:,slice))
colormap(gray)
title('Slice 50')
xlabel('Position')
ylabel('Position')
colorbar
%% Akuisisi Data Proyeksi
 
disp(sprintf('Akuisisi data proyeksi...'));
 
% Menentukan Jumlah Berkas Sinar
[M N] = size(ct3d(:,:,slice));
rhomax = round(sqrt(M^2 + N^2));
if mod(rhomax,2)==0;
    nrho = rhomax + 1;
    rho = zeros(1,nrho);
    rho(1,1)=-rhomax/2*ps;
else
    nrho = rhomax + 2;
    rho = zeros(1,nrho);
    rho(1,1)=-(rhomax+1)/2*ps;
end
 
for i=2:nrho
    rho(1,i)=rho(1,1)+((i-1)*ps);
end
 
% Menentukan jumlah sudur proyeksi
ntheta = 180;
maxtheta = 180;
theta = zeros(ntheta,1);
theta(1,1)=maxtheta/ntheta;
for i = 2:ntheta
    theta(i,1)=theta(1,1)+((i-1)*theta(1,1));
end
 
res = zeros(nrho, ntheta);
post = zeros(nx,ny);
 
% Proses Akuisisi Data Proyeksi
for i=1:44
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+ct3d(j,k,slice)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end
 
for i=45
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=((sqrt(ps^2+ps^2))/2);
                    res(l,i)=res(l,i)+(ct3d(j,k,slice)*2*abs(post(j,k)-rho(1,l)));
                end
            end
        end
    end
end
for i=46:134
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+ct3d(j,k,slice)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end
 
for i=135
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=((sqrt(ps^2+ps^2))/2);
                    res(l,i)=res(l,i)+(ct3d(j,k,slice)*2*abs(post(j,k)-rho(1,l)));
                end
            end
        end
    end
end
 
for i=136:180
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+ct3d(j,k,slice)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end
 
% Menampilkan Gambar Hasil Sinogram
figure(2)
imagesc(theta,rho,res);
colormap(gray);
title('Sinogram Slice 50')
xlabel('Angle')
ylabel('Ray Position')
colorbar
 
%% Rekonstruksi Citra dengan Proyeksi Balik Tanpa Filter
 
disp(sprintf('\nProyeksi Balik Citra...'));
% Proyeksi balik data sinogram
sinogram = res(108:619,:);
lamin = zeros(nx,ny);  
  for ia = 1:180
    projection_ia=sinogram(:,ia); 
    projection_smear=repmat (projection_ia,1,512);
    rot= imrotate(projection_smear', ia, 'bicubic','crop');
    lamin=lamin+rot';   
  end
 
%Kalibrasi nilai gray level
m = max(max(lamin))+abs(min(min(lamin)));
lamin = ((lamin + abs(min(min(lamin))))/m)*4095;
 
% Menampilkan citra hasil proyeksi balik
figure(3)
imagesc(x, y, lamin'); colormap('gray'); 
axis('image');
title('Simple Backprojection Image');
xlabel('mm');  
ylabel('mm');
colorbar;

%Hitung Nilai PSNR
input = ct3d(:,:,slice);
output = lamin';
M = 512;
N = 512;
peakval = 4095;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Ram-Lak

disp(sprintf('\nFilter Ram-Lak...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Ram-Lak
sinogramfiltered=sinogram;     
a = length(sinogram);
freq=linspace(-1, 1, a).';
Filter = abs(freq);
Filtersp = ifftshift(ifft(ifftshift(Filter)));
Filtersp = real(Filtersp);
sinogramfilt = zeros(512,180);
for i = 1:180
    sinogramfilt(:,i)=imfilter(sinogramfiltered(:,i),Filtersp,'same','conv');    
end  

%Proyeksi balik data sinogram
bpf_recon = zeros(nx,ny);
  for ia = 1:180
    bpf_ia=sinogramfilt(:,ia);
    bpf_smear=repmat(bpf_ia,1,512);
    rot1= imrotate(bpf_smear', ia, 'bicubic','crop');   
    bpf_recon=bpf_recon+(rot1'/(180));
  end
  
%Kalibrasi nilai gray level
m = max(max(bpf_recon))+abs(min(min(bpf_recon)));
bpf_recon = ((bpf_recon + abs(min(min(bpf_recon))))/m)*4095;

 figure(4)
    imagesc(theta,rho, sinogramfilt)   
    colormap('gray')                 
    title('Sinogram Filtered Ram-Lak')
    xlabel('Angle')
    ylabel('Ray Position')
    colorbar
    
% Plot Frekuensi Filter Ram-Lak
figure(5)
plot(freq, Filtersp);  
title('Filter Ram-Lak')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Ram-Lak
figure(6)
imagesc(x, y, bpf_recon'); colormap('gray'); 
axis('image')  
title('Filter Ram-Lak')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = ct3d(:,:,slice);
output = bpf_recon';
M = 512;
N = 512;
peakval = 4095;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Hamming  

disp(sprintf('\nFilter Hamming...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hamming
sinogramfiltered2=sinogram;
b = length(sinogram);
hamm = zeros(b,1);
freq=linspace(-1, 1, b).';
absfreq = abs(freq);
for i=1:b
    hamm(i,1)=0.54+0.46*cos(pi*freq(i,1));
end
hamm = hamm.*absfreq; 
hammsp = real(ifftshift(ifft(ifftshift(hamm))));
sinogramfit2 = zeros(512,180);
for i = 1:180
    sinogramfilt2(:,i)=imfilter(sinogramfiltered2(:,i),hammsp,'same','conv');    
end

% Proyeksi balik data sinogram
bpf_recon2 = zeros(nx,ny);
  for ia = 1:180
    bpf_ia2=sinogramfilt2(:,ia);
    bpf_smear2=repmat(bpf_ia2,1,512);
    rot2= imrotate(bpf_smear2', ia, 'bicubic','crop');   
    bpf_recon2=bpf_recon2+(rot2'/(180));
  end
  
%Kalibrasi nilai gray level
m = max(max(bpf_recon2))+abs(min(min(bpf_recon2)));
bpf_recon2 = ((bpf_recon2 + abs(min(min(bpf_recon2))))/m)*4095;

% Respon Frekuensi Filter Hamming
figure(7)
plot(freq, hammsp);  
title('Filter Hamming')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Hamming
figure(8)
imagesc(x, y, bpf_recon2'); colormap('gray'); 
axis('image')  
title('Filter Hamming')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = ct3d(:,:,slice);
output = bpf_recon2';
M = 512;
N = 512;
peakval = 4095;
Fungsi(input,output,M,N,peakval);

%% Rekonstruksi Citra dengan Penerapan Filter Hanning

disp(sprintf('\nFilter Hanning...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hanning
sinogramfiltered3=sinogram;
c = length(sinogram);
hann = zeros(c,1);
freq=linspace(-1, 1, c).';
absfreq = abs(freq);
for i=1:c
    hann(i,1)=0.5*(1+cos(pi*freq(i,1)));
end
hann = hann.*absfreq;
hannsp = real(ifftshift(ifft(ifftshift(hann))));
sinogramfilt3 = zeros(512,180);
for i = 1:180
    sinogramfilt3=imfilter(sinogramfiltered3,hannsp,'same','conv');  
end   

%Proyeksi balik data sinogram
bpf_recon3 = zeros(nx,ny);
  for ia = 1:180
    bpf_ia3=sinogramfilt3(:,ia);
    bpf_smear3=repmat(bpf_ia3,1,512);
    rot3= imrotate(bpf_smear3', ia, 'bicubic','crop');   
    bpf_recon3=bpf_recon3+(rot3'/(180));
  end

%Kalibrasi nilai gray level
m = max(max(bpf_recon3))+abs(min(min(bpf_recon3)));
bpf_recon3 = ((bpf_recon3 + abs(min(min(bpf_recon3))))/m)*4095;

% Respon Frekuensi Filter Hanning
figure(9)
plot(freq, hannsp);  
title('Filter Hanning')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Hanning
figure(10)
imagesc(x, y, bpf_recon3'); colormap('gray'); 
axis('image')  
title('Filter Hanning')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = ct3d(:,:,slice);
output = bpf_recon3';
M = 512;
N = 512;
peakval = 4095;
Fungsi(input,output,M,N,peakval);

function [MSE,PSNR] = Fungsi(input,output,M,N,peakval)
    MSE = sum(sum((input-output).^2))/(M*N);
    PSNR = 10*log10((peakval^2)/MSE);
    fprintf('MSE : %7.2f ', MSE);
    fprintf('\nPSNR : %9.7f dB \n', PSNR);
end