function [] = nd2mat(file, folder, save_prefix, varargin)
%Convert nd2files to .mat files.
% Individual xyz stack files are saved for each series (multi-point),
% channel, and time

% File info
r = bfGetReader(fullfile(folder, file));

% Parse Inputs
p = inputParser;
% default_multipt_name = [];
% default_channel_name = [];
addParamValue(p,'multipt', []);
addParamValue(p,'channel', []);
parse(p, varargin{:});
multipt = p.Results.multipt;
channel = p.Results.channel;

% get image sizes
sizeX = r.getSizeX();
sizeY = r.getSizeY();
sizeC = r.getSizeC();
sizeT = r.getSizeT();
sizeZ = r.getSizeZ();
sizeM = r.getSeriesCount;

% Check if the names for multipt and channel match given size in the image
if ~isempty(multipt)
    if length(multipt) ~= sizeM
        error('Size of multipt input names does not match number of multipts in the image')
    end
end
if ~isempty(channel)
    if length(channel) ~= sizeC
        error('Size of Channel input names does not match number of channels in the image')
    end
end



tic
for m = 1:(sizeM)
    % Set current multi-time-point
    r.setSeries(m-1);
    
    nImages = r.getImageCount();
    
    for t = 1:sizeT%1:(sizeT)
        for c = 1:(sizeC)
            
            %Initialize
            vol = zeros([sizeX,sizeY,sizeZ],'uint16');
            counter = 0;
                    
            for i = 1:nImages
                zct = r.getZCTCoords(i-1) + 1;
                
                if c == zct(2) & t == zct(3)
                    xy = r.openBytes(i-1);
                    xy = typecast(xy, 'uint16');
                    vol(:,:,zct(1)) = reshape(xy, [sizeX sizeY]);
                end
            end
            
            % Multi_pt name
            if isempty(multipt)
                mp_name = sprintf('_M%0.2d', m);
            else
                mp_name = strcat('_', multipt{m});
            end
            % Channel name
            if isempty(channel)
                c_name = sprintf('_C%0.1d', c);
            else
                c_name = strcat('_', channel{c});
            end
            % Time name
            t_name = sprintf('_T%0.3d', t);
            
            % Save file
            filename = strcat(save_prefix, mp_name, c_name, t_name, '.mat');    
            save(fullfile(folder, filename), 'vol', '-v7.3');
            
            display(strcat('Saved file: ', filename, '    Time elapsed: ',...
                sprintf('%0.2f', toc), ' sec'))
        end
    end
end
display('File conversion complete')

end

