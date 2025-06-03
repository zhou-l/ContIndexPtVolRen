function [outFileNameList,writeOkList] = writeVolumes(vols, fileList)
outFileNameList = fileList;
nVols = length(vols);
writeOkList = zeros(nVols,1);

for i = 1:nVols
    fileName = fileList{i};
    ok = nrrdWriter(fileName, single(vols{i}), [1 1 1], [0 0 0], 'raw');
    outFileNameList{i} = fileName;
    writeOkList(i) = ok;
end
writeOkList = logical(writeOkList);
end