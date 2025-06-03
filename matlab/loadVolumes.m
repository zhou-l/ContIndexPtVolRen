function [Vlist, nVols, volInfoList] = loadVolumes(fileList)
% load the nhdr/nrrd file
nVols = length(fileList);
Vlist = cell(nVols, 1);
volInfoList = cell(nVols,1);
for i = 1:nVols
    filename = fileList{i};
    headerInfo = nhdr_nrrd_read(filename, true);
    volInfoList{i} = headerInfo;
    [filepath, name, ext] = fileparts(filename);
    V = headerInfo.data;
    dimX = size(V,1);
    dimY = size(V,2);
    dimZ = size(V,3);
    for z = 1:dimZ
        for y = 1:dimY
            for x = 1:dimX
               if isnan(V(x,y,z))
                   V(x,y,z) = -6;
               end
            end
        end
    end
%     volumeViewer(V);
    Vlist{i} = V;
end 