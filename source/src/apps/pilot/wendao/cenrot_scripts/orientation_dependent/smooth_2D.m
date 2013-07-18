function smoothed = smooth_2D(data, param1, param2)
data = squeeze(data);
[a, b] = size(data);
wr = param1(1);
ar = param1(2);
flagr = param1(3);
wa = param2(1);
aa = param2(2);
flaga = param2(3);

%filter have half size of data
[xf,yf] = meshgrid(linspace(-wr/4,wr/4,a/2), linspace(-wa/4,wa/4,b/2));
Gaussian = exp(-(xf.^2/ar^2+yf.^2/aa^2));
filter = Gaussian / sum(sum(Gaussian));

bda = int8(a/4);
bdb = int8(b/4);

%col
if flaga== 1
	extended = [data(:,(b-bdb+1):b), data, data(:,(1:bdb))];
else 
	symm = fliplr(data);
	extended = [symm(:,(b-bdb+1):b), data, symm(:,1:bdb)];
end
%row
if flagr == 1
	extended = [extended(a-bda+1:a,:); extended; extended(1:bda,:)];
else 
	symm = fliplr(extended')';
	extended = [symm(a-bda+1:a,:); extended; symm(1:bda,:)];
end

smoothed_ext = conv2(extended, filter,'same');
smoothed = smoothed_ext( (bda+1):(bda+a), (bdb+1):(bdb+b));
end
