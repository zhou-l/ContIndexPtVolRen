function common_to_use = commonString(Scell)
    Schar = char(Scell(:));
    diffChar = diff(Schar);
    anyNoneZero = any(diffChar,1);
    common_to_use = Scell{1}(1:find(anyNoneZero,1,'first')-1);   
    if isempty(common_to_use)
        common_to_use='?'
    end
end