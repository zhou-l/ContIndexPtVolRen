function hit = aabbIntersect(AminBounds, AmaxBounds, BminBounds, BmaxBounds)
   hit = AminBounds(1) <= BmaxBounds(1) && AmaxBounds(1) >= BminBounds(1)...
      && AminBounds(2) <= BmaxBounds(2) && AmaxBounds(2) >= BminBounds(2)...
      && AminBounds(3) <= BmaxBounds(3) && AmaxBounds(3) >= BminBounds(3);
%     function intersect(a, b) {
%   return (
%     a.minX <= b.maxX &&
%     a.maxX >= b.minX &&
%     a.minY <= b.maxY &&
%     a.maxY >= b.minY &&
%     a.minZ <= b.maxZ &&
%     a.maxZ >= b.minZ<
%   );
% }
end