% bowl parameters
zmin = 80;
ep = 0.01;

if (PlotType == 1)
    hold on;
    
    % Patches
    showpatchborders(1:3);
    % showgridlines(1:2);    
    setpatchborderprops('linewidth',2);
    set(gca,'zlim',[0,10]);
    % set(gca,'zlim',[-20 1]);   % Need so that all patchborders show up
    
    colormap(parula);
    colorbar;
    tol = -0.8;
    c1 = -0.14;
    c2 = 0.15;
    caxis([c1,c2]);
    
    % Contour lines (some are hidden by patch boundaries)
    cv = linspace(c1,c2,21);
    drawcontourlines(cv);
    
    % Plot boundary of true solution
    a = 1;
    sigma = 0.5;
    h0 = 0.1;
    grav = 9.81;
    omega = sqrt(2*grav*h0) / a;

    xe = linspace(-2,2,200);
    ye = linspace(-2,2,200);
    [xem,yem] = meshgrid(xe,ye);
    
    y = 0;
    B = h0*(xem.^2 + yem.^2)/a^2 - h0;
    eta1 = sigma*h0/a.^2 * (2*xem*cos(omega*t) + 2.*yem*sin(omega*t) -sigma);
    
    % contour(xem,yem,B-eta1,[0 0],'k','linewidth',5);
   
    hold off;

    fprintf('%10s %g\n','qmin',qmin)
    fprintf('%10s %g\n','qmax',qmax)    
    
    % Axes
    axis([-2 2 -2 2])
    daspect([1 1 1]);
    set(gca,'fontsize',16);
else
    hold on;
    r = linspace(0,100,500);
    r2 = r.^2;
    b = ep*r.^2 - zmin;
    b(b >= 0) = 0;
    plot(r,b,'k','linewidth',2);  
    xlim([0,100]);
    ylim([-zmin,20]);    
    daspect([1 1 1])    
    [h_amr, labels_amr] = getlegendinfo;
    h = legend(h_amr,labels_amr,'location','southeast');
    hold off;
end

title(sprintf('slosh (amrclaw) at time %g',t),'fontsize',18);

plot_surf = true;
if plot_surf
    figure(2);
    a = 1;
    h0 = 0.1;
    s = 0;
    
    [xm,ym] = meshgrid(xcenter,ycenter);
    h = surf(xm,ym,q'+s);
    set(h,'edgecolor','none');
    hold on;
    
    if (exist('hp','var'))
        delete(hp);
    end
    xf = linspace(-2,2,500);
    yf = xf;
    [xfm,yfm] = meshgrid(xf,yf);
    P = h0*(xfm.^2 + yfm.^2)/a.^2 - h0 + s;
    m = P > 0;
    % P(m) = 0;
    hp = surf(xfm,yfm,P);
    set(hp,'facecolor',[1,1,1]*0.2);
    set(hp,'edgecolor','none','facealpha',0.6);
    c1 = -0.14 + s;
    c2 = 0.15 + s;
    caxis([c1,c2]);
    daspect([1,1,0.25])
    axis off;
    set(gca,'zlim',[2*c1,c2]);
    camlight;
    
    axis([-2,2,-2,2]);
    % view(3);
    view([-58,22]);
    hold off;
    figure(1);
end




NoQuery = 0;
MaxFrames = 16;
prt = true;
if (prt)
    figure(2);
    filename = sprintf('parabola%04d.png',Frame);
    fprintf('Print file %s\n',filename);
    print('-dpng',filename);
    figure(1);
end

shg
