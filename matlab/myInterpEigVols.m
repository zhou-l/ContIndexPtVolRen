function fullSzEigVols = myInterpEigVols(nVols, fullDim, dsDim, dsMajEigVols)
%    3---2
%   /|  /|
%  0---1 |
%  | 4-|-5
%  |/  |/
%  7---6
checkPosList =...
   [0 0 0;
    1 0 0;
    1 1 0;
    0 1 0;
    0 1 1;
    1 1 1;
    1 0 1;
    0 0 1];

% interpolate eigen information?
fullSzEigVols = cell(nVols,1);
for i = 1:nVols
    fullSzEigVols{i} = zeros(fullDim(1),fullDim(2),fullDim(3));
end
numNeighbors = 8;
for z = 1:fullDim(3)
    for y = 1:fullDim(2)
        for x = 1:fullDim(1)
            pos = [x y z];
            fpos = pos ./ fullDim;
            % get its eight neighbors
            dsPos = fpos .* dsDim + [1 1 1];
            idsPos = (floor(dsPos));
            
            posList = idsPos + checkPosList;
            % value list for interpolations
            valList = zeros(numNeighbors,nVols);
 
            % blend factor creation
            aList = createWeights(dsPos, idsPos);
            
            for nn = 1:numNeighbors
                npos = posList(nn,:);
%                 idx = find(checkPosList(nn,:)>0);
                
                if npos(1) >= dsDim(1) || npos(2) >= dsDim(2) || npos(3) >= dsDim(3)
                    continue;
                end
                for d = 1:nVols
                    valList(nn,d,:) = dsMajEigVols{d}(npos(1),npos(2),npos(3));
                end
            end
            % interpolation
            eigVec = eigInterp3D(valList, aList);
%             eigVec = eigInterpPerElement3D(valList, aList);

            for d = 1:nVols
                fullSzEigVols{d}(pos(1),pos(2),pos(3)) = eigVec(d); % the major eigenvector
            end
        end
    end
end
end

function alphaList = createWeights(dsPos, idsPos)
% checkPosList =...
%    [0 0 0;
%     1 0 0;
%     0 1 0;
%     0 0 1;
%     1 1 0;
%     0 1 1;
%     1 0 1;
%     1 1 1];
alphaList = zeros(6,1);
alphaList(2,:) = dsPos(1) - idsPos(1);% weights on x axis 
alphaList(1,:) = 1.0 - alphaList(2);
alphaList(4,:) = dsPos(2) - idsPos(2);% y axis
alphaList(3,:) = 1.0 - alphaList(4);
alphaList(6,:) = dsPos(3) - idsPos(3);% z axis
alphaList(5,:) = 1.0 - alphaList(6);
end