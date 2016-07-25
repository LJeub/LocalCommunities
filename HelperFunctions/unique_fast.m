function u=unique_fast(a)
% fast computation of unique elments of array without type-check overhead
if ~isempty(a)
a=sort(a(:));
d=diff(a)~=0;
u=[a(d);a(end)];
else
    u=[];
end

end
