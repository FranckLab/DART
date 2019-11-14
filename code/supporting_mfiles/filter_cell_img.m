function [] = filter_cell_img(filename)
%Find cell mask and cell surface
tic 

% Load images
load(filename, 'vol');
vol = double(vol);
vol = vol/(2^12);

% Filter image
I = vol;
I = imgaussfilt(medfilt3(vol), 2.5);
clear vol
vol = uint16(I*2^12); 

% Save image back to the same file
save(filename, 'vol', '-append')
disp(strcat('Found cell in files:  ', filename, ...
            sprintf('  Time elapsed: %0.2f seconds', toc)))

end

