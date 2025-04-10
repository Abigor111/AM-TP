function y = funcODE(f,a,b,n,y0)
   
h = (b-a)/n;
    t = a:h:b;
    [~,yODE45] = ode45(f,t,y0);
    y = yODE45.';


end