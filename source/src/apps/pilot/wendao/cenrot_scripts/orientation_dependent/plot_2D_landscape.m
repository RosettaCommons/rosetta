function plot_2D_landscape(a,b,ai,bi, flag)
dat1 = squeeze(a(ai,bi,:,:));
dat2 = squeeze(b(bi,ai,:,:));

for kk = 1:30
	dat1(kk,:) = dat1(kk,:);
	dat2(kk,:) = dat2(kk,:);
end

data = (dat1+dat2);
[a, b] = size(dat1);
[xx, yy] = meshgrid(1:b, 1:a);
data = smooth_2D(data, [15, 1, 2], [2, 0.1, flag]);

f = figure();
set(f, "visible", "off")
subplot (2, 2, 1);
contourf(xx,yy,dat1);
subplot (2, 2, 2);
contourf(xx,yy,dat2);
subplot (2, 2, 3);
contourf(xx,yy,dat1+dat2);
subplot (2, 2, 4);
contourf(xx,yy,data);
print("out.png", "-dpng")

