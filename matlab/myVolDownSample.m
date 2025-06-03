function dsList = myVolDownSample(Vlist, sampleRatePerAxis)
nVols = length(Vlist);
if nVols < 1
    dsList = [];
    return;
end

dim = size(Vlist{1});
% stepSize = ceil(1 / sampleRatePerAxis);
dsList = cell(nVols, 1);
dsDim = floor(dim .* sampleRatePerAxis);%floor((dim - [1 1 1]) ./ stepSize);
% cntZ = 1;
for d = 1:nVols
    dsList{d} = zeros(dsDim(1),dsDim(2),dsDim(3));
    dsList{d} = imresize3(Vlist{d}, [dsDim(1) dsDim(2) dsDim(3)]);
%     cntX = 1; cntY = 1; cntZ = 1;
%     for z = 1:stepSize:dim(3)
%         cntY = 1;
%         for y = 1:stepSize:dim(2)
%             cntX = 1;
%             for x = 1:stepSize:dim(1)
%                 dsList{d}(cntX,cntY,cntZ) = Vlist{d}(x,y,z);
%                 cntX = cntX + 1;
%             end
%             cntY = cntY + 1;
%         end
%         cntZ = cntZ + 1;
%     end
end

end