function [I] = process_images(I)
% Process images for histogram equalization and cleaning image for particle
% detection

y = double(bpass3dMB(I, 0.5*[1,1,1], 4*[1,1,1], [0,0]));

% intensity line
y_int = squeeze(sum(sum(y, 1), 2));

% Fit a power function to the intensity line
[~, idx] = max(y_int); idx = idx-4;
if idx<1
    idx=1;
end
fit_x = (idx:length(y_int))';
fit_y = y_int(idx:end);

% Perform fit
f_int = fit(fit_x,fit_y,'power2');

% Intensity
norm_int = f_int(1:length(y_int));
norm_int(1:idx+5) = norm_int(idx+5);
norm_int = norm_int/max(norm_int);
norm_int_factor = ones(1,1,length(norm_int));
norm_int_factor(:) = 1./norm_int;

% Normalize image intensity
I = bsxfun(@times, I, norm_int_factor);
I = I/max(I(:));

% Bandpass filter on image intensity normalized images.
I = double(bpass3dMB(I, 0.5*[1,1,1], 4*[1,1,1], [0,0]));
I = I/max(I(:));

end

