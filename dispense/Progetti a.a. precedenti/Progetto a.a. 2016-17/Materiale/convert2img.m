function f = convert2img(sorg,dest)
    f = imread(sorg);
    f = rgb2gray(f);
    f = im2double(f);
    f = single(f);
    saveimg(dest,f);
end


function saveimg(fname,x)
	x=x';
	dim=size(x);
	file=fopen(fname,'w');
	fwrite(file,[dim(1) dim(2)],'integer*4');
	fwrite(file,x,'float');
	fclose(file);
end
