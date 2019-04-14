function save_dense(P,filename)
    [n,m] = size(P);
    fname = sprintf('%s.matrix',filename);
    fid = fopen(fname,'w');
    fwrite(fid,n,'uint32');
    fwrite(fid,m,'uint32');
    for i = 1:n
        fwrite(fid,P(i,:),'double');
    end
    fclose(fid);
end
