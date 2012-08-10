function y=airfoil(x, xa,yaul)
y=interp1(xa,yaul,x,'spline');
end