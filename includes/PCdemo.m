function [row, col] = PCdemo(I, I2)
%I,I2 are reference and target images
%[row col] are row, column shifts
%SCd 4/2010
%

%Fourier transform both images
fi = fft2(double(I));
fr = fft2(double(I2));

%Perform phase correlation (amplitude is normalized)
fc = fi .* conj(fr);
fcn = fc ./abs(fc);

%Inverse fourier of peak correlation matrix and max location
peak_correlation_matrix = abs(ifft2(fcn));
[peak, idx] = max(peak_correlation_matrix(:));

%Calculate actual translation
[row, col] = ind2sub(size(peak_correlation_matrix),idx);
if row < size(peak_correlation_matrix,1)/2
    row = -(row - 1);
else
    row = size(peak_correlation_matrix,1) - (row - 1);
end;
if col < size(peak_correlation_matrix,2)/2
    col = -(col - 1);
else
    col = size(peak_correlation_matrix,2) - (col - 1);
end



if 0
    %%
    I = imread('cameraman.tif'); %Read in cameraman image
    I2 = imtranslate(I,[-11.3 17.8]); %Translate 11.3 down 17.8 over
    
    figure;
    imshow(I);
    title('Cameraman Original');
    figure;
    imshow(I2);
    title('Cameraman Shifted');
    
    [row col] = PCdemo(I, I2);
    
    disp('The shifts are:')
    disp([row col])
end

