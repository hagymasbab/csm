function cc = filterList2Components(gestalts,nullComp,Dx)
    k = size(gestalts,1);
    cc = {};
    for kk=1:k
        cc{kk} = zeros(Dx);
        for i = 1:size(gestalts,2)
            cc{kk}(gestalts(kk,i),gestalts(kk,i)) = 1;
            for j =i+1:size(gestalts,2)
               cc{kk}(gestalts(kk,i),gestalts(kk,j)) = 0.5; 
               cc{kk}(gestalts(kk,j),gestalts(kk,i)) = cc{kk}(gestalts(kk,i),gestalts(kk,j));
            end
        end
        if ~nullComp
            cc{kk} = cc{kk} + 0.1*eye(Dx);
        end
    end
    if nullComp
        cc{end+1} = eye(Dx);
    end
end