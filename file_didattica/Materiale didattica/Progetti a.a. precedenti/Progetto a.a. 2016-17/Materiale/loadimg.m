function f = loadimg(fname)

    if (~exist(fname,'file'))
        error('File %s not found.', fname);
    end
    
    file=fopen(fname,'r');
    dim=fread(file,[1 2],'integer*4');
    f=fread(file,[dim(1) dim(2)],'float');
    fclose(file);
    f=f';

    imshow(mat2gray(f));
    
end
