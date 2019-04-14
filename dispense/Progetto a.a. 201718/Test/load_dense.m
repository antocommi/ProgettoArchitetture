function P = load_dense(filename)
    fname = sprintf('%s.matrix',filename);
    fid = fopen(fname,'r');
    n = fread(fid,1,'uint32');
    m = fread(fid,1,'uint32');
    P = zeros(n,m);
    for i = 1:n
        P(i,:) = fread(fid,m,'double');
    end
    fclose(fid);
end
