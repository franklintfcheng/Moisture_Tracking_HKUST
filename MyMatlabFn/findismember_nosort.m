function idx = findismember_nosort(A, B)

[tf,loc] = ismember(A, B);
[~,p] = sort(loc(tf));
idx = find(tf);
idx = idx(p);
end