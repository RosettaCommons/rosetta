function plot_2D_landscape_sym(c, ai, bi)
dat1 = squeeze(c(ai,bi,:,:));
dat2 = squeeze(c(bi,ai,:,:));

[a, b] = size(dat1);
[xx, yy] = meshgrid(1:b, 1:a);
data = smooth_2D(dat1+dat2, [15, 1.0, 2], [360, 18, 1]);

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

