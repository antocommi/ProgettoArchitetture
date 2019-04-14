function save_sparse(G,filename)
    if (size(G,2)~=2)
        error('Not a valid sparse graph');
    end
    n = max(max(G));
    m = size(G,1);
    fname = sprintf('%s.graph',filename);
    fid = fopen(fname,'w');
    fwrite(fid,n,'uint32');
    fwrite(fid,m,'uint32');
    for i = 1:m
        fwrite(fid,uint32(G(i,1:2)),'uint32');
    end
    fclose(fid);
end
