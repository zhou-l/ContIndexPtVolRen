
function testInterpTriangle(imgDataList, m)
%     pts = [253, 261; 257 261; 255 258];
    pts = [452, 412; 448 420; 455 418; 455 412];
    nEV = 3; % the number of eigenvectors 
    Vlist = zeros(nEV,3*m);
    xmin = min(pts(:,1));
    xmax = max(pts(:,1));
    ymin = min(pts(:,2));
    ymax = max(pts(:,2));
    R = 3;
    pca_test_gold = zeros(ymax-ymin+1, xmax-xmin+1, 3*m );
    pca_test_interp = zeros(ymax-ymin+1, xmax-xmin+1, 3*m );

    % find the gold solution
    for x = xmin:xmax
        for y = ymin:ymax
              pca_gold = findSamplesInBallGold([y,x], imgDataList, R);
              offpos = [y,x] - [ymin, xmin] + 1;
              for k = 1:nEV
                pca_test_gold(offpos(1),offpos(2),(k-1)*m+1:k * m) = pca_gold(:,k);
%                 pca_test_gold(offpos(1),offpos(2),m+1:2*m) = pca_gold(:,2);
%                 pca_test_gold(offpos(1),offpos(2),2*m+1:3*m) = pca_gold(:,3);
              end
        end
    end

    % interpolation scheme
    pp = pts - [xmin ymin] + 1;
    for i = 1:size(pp,1)
        Vlist(i,:) = pca_test_gold(pp(i,2),pp(i,1),:);
    end
    meanV = mean(Vlist, 1);
    C = [1 2 3; 1 3 4];
    TR = triangulation(C,pts);
    triplot(TR);
    for x = xmin:xmax
        for y = ymin:ymax
            Pq = [x y];
            [ti,bc] = pointLocation(TR,Pq);
            offpos = [y,x] - [ymin,xmin] + 1;
            if isnan(ti)
                % outside of the triangle
                pca_test_interp(offpos(1),offpos(2),:) = meanV;
                continue;
            end
            % find the correct correspondence of vectors
            corresEvs = sortAssocEvs(Vlist, m, bc);

            % and then perform the interpolation
            for nn = 1:nEV
                EvInterp(1,(nn-1)*m + 1:nn*m) = eigInterp2D(corresEvs(:,(nn-1)*m+1:nn*m), bc);
                pca_test_interp(offpos(1),offpos(2),(nn-1)*m + 1:nn*m) = EvInterp(1,(nn-1)*m + 1:nn*m) ;
                % interpolate the gold computation          
            end
        end
    end
   figure;
   subplot(2,3,1);imshow(pca_test_interp(:,:,1),[]);
   subplot(2,3,2);imshow(pca_test_interp(:,:,2),[]);
   subplot(2,3,3);imshow(pca_test_interp(:,:,3),[]);

   subplot(2,3,4);imshow(pca_test_gold(:,:,1),[]);
   subplot(2,3,5);imshow(pca_test_gold(:,:,2),[]);
   subplot(2,3,6);imshow(pca_test_gold(:,:,3),[]);
end