function res=nextnode(dx,x,d,xa,yaul)
    dx=abs(dx);
    %if(x<0.25)
       % if x<0.12
       %     d=d/1.5;
       % else
       %     d=d/1.3;
       % end
   % end
    res=sqrt(dx^2+(airfoil(x+dx, xa,yaul)-airfoil(x, xa,yaul))^2)-d;
end
