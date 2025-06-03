function mu = localMu_spaceDom(pos, vols, R)
nVols = length(vols);
D = zeros((2*R+1)*(2*R+1)*(2*R+1),nVols);
cnt = 1;
dim = size(vols{1});
% compute with Monte Carlo sampling!!!
s = rng;
for z = -R:R
    for y = -R:R
        for x = -R:R
            % get random position
            offset = 2.*R.*rand(1,3) - [R,R,R];
            % with integral position
%             npos = pos + [x y z];
            % with random position
            npos = pos + offset;
            if npos(1) < 1 || npos(2) < 1 || npos(3) < 1 ||...
                    npos(1) > dim(1) || npos(2) > dim(2) || npos(3) > dim(3)
                continue;
            end
            val = zeros(1,nVols);
            for d = 1:nVols
                % with value fetching
%                val(1,d) = vols{d}(npos(1),npos(2),npos(3));
                % with interpolation
                val(1,d) = volInterp(vols{d}, npos);
            end
            D(cnt,:)  = val;
            cnt = cnt + 1;
        end
    end
end
D(cnt:end,:) = [];
% the columns of the input data are variables, rows are data samples
DM = D;% D';
mu = mean(DM);
end