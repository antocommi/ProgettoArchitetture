function G = load_sparse(filename)
    fname = sprintf('%s.graph',filename);
    fid = fopen(fname,'r');
    n = fread(fid,1,'uint32');
    m = fread(fid,1,'uint32');
    G = fread(fid,[2 m],'uint32');
    fclose(fid);
    G = double(G');
end
