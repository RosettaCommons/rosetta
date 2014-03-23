function [w] = gausswin(L, a)
	x = linspace(-(L-1)/L, (L-1)/L, L);
	w = exp ( -(a*x).^2/2 );
