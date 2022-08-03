%%%%%%%INPUTS
%If set true, only plots a single layer in the three plane options below
singlelayer=0;
%layeraverage=0;

xyplot=0; %If set true, plots a top-down view of the field in a given layer.
xzplot=1;
yzplot=0;

%%%%%%%%

%%% Vertical layer index (1 starting from the surface) to use for top-down xy plots
xylayer = 5;

zlaymin = 30;
zlaymax = 58;%Nr-2;

%%% Layer to plot in the y/z plane
yzlayer = 40;

xlaymin = 50;
xlaymax = 78; %Nx-1;

%%% Layer to plot in the x/z plane
xzlayer = 40;

ylaymin = 60;%1;
ylaymax = 98;%Ny-2;


%%% Mesh grids for plotting
hFac = hFacC;
kmax = zeros(Nx,Ny);
if (xyplot)
  [YY,XX] = meshgrid(yy/1000,xx/1000);
  kmax = sum(ceil(hFac),3);
  kmax(kmax==0) = 1;
else  
  %%% Create mesh grid with vertical positions adjusted to sit on the bottom
  %%% topography and at the surface
  [ZZ,YY] = meshgrid(zz,yy/1000);  
  for j=1:Ny
    if (yzavg)
      hFacC_col = squeeze(hFacC(:,j,:));    
      hFacC_col = max(hFacC_col,[],1);    
    else
      hFacC_col = squeeze(hFacC(yzlayer,j,:))';
    end
    kmax = length(hFacC_col(hFacC_col>0));  
    zz_botface = -sum(hFacC_col.*delR);
    ZZ(j,1) = 0;
    if (kmax>0)
      ZZ(j,kmax) = zz_botface;
    end
  end
  
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.11 0.15 0.76 0.73];
  framepos = [327    80   941   885];
end

%%% Set up the figure
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');

Amean = [];
tdays = [];
%%% retrieve the last frame uvel, vvel, wvel (the 3D time-averaged fields over the
%%% last day at day 10)


n=length(dumpIters)
% for n=50:length(dumpIters)
  
  t = dumpIters(n)*deltaT/86400;
  
  if (n > 1)
    Aprev = A;
  end
%   A = rdmdsWrapper(fullfile(exppath,'results',outfname),0);          
  uvel = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));   %most important line    
vvel = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));   %most important line   
wvel = rdmdsWrapper(fullfile(exppath,'results','WVEL'),dumpIters(n));   %most important line   


%calculate the individual advection terms

%first term in advection operator
delX3D=repmat(delX',[1,Ny,Nr]);
dvdx = (vvel(2:end,:,:)-vvel(1:end-1,:,:))./delX3D(1:end-1,:,:);
u = uvel(1:end-1,:,:)+uvel(2:end,:,:)./2;
dvdxav = dvdx(1:end-1,1:end-1,:)./4+dvdx(2:end,1:end-1,:)./4+dvdx(1:end-1,2:end,:)./4+dvdx(2:end,2:end,:)./4;

u_dvdx = u(1:end-1,1:end-1,:).*dvdxav;


  %second term in advection operator
delY3D=repmat(delY,[Nx,1,Nr]);
dvdy = (vvel(:,2:end,:)-vvel(:,1:end-1,:))./delY3D(:,1:end-1,:);
v = vvel(:,1:end-1,:)+vvel(:,2:end,:)./2;
dvdyav = dvdy(1:end-1,1:end-1,:)./4+dvdy(2:end,1:end-1,:)./4+dvdy(1:end-1,2:end,:)./4+dvdy(2:end,2:end,:)./4;

v_dvdy = v(1:end-1,1:end-1,:).*dvdyav;


   %third term in advection operator
newdelR = reshape(delR,[1,1,60]);
delR3D=repmat(newdelR,[Nx,Ny,1]);
dvdz = (vvel(:,:,2:end)-vvel(:,:,1:end-1))./delR3D(:,:,1:end-1);
w = wvel(:,:,1:end-1)+wvel(:,:,2:end)./2;
dvdzav = dvdz(:,1:end-1,1:end-1)./4+dvdz(:,2:end,1:end-1)./4+dvdz(:,1:end-1,2:end)./4+dvdz(:,2:end,2:end)./4;

w_dvdz = w(:,1:end-1,1:end-1).*dvdzav;

%sum the advection terms

adv = u_dvdx(:,1:end-1,1:end-2) + v_dvdy(1:end-1,:,1:end-2) + w_dvdz(1:end-2,1:end-1,:);


  if (isempty(adv))
    error(['Ran out of data at t=,',num2str(tdays(n)),' years']);
  end      
  
  max(max(max(abs(adv))))
 
  tdays(n) = t;

  %{
  %%% x/y plot
    
    if (botplot)      
      FF = zeros(Nx,Ny);
      for i=1:Nx
        for j=1:Ny
          FF(i,j) = adv(i,j,kmax(i,j));
%           FF(i,j) = A(i,j,kmax(i,j))*hFac(i,j,kmax(i,j));
        end
      end
    else
      FF = squeeze(adv(:,:,xylayer,outfidx));        
    end
    
    %FF(hFacC(:,:,xylayer)==0) = NaN;
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX(1:end-1,1:end-1),YY(1:end-1, 1:end-1),FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');
  %}


   %{ 
  %%% y/z zonally-averaged plot


    if (yzavg)
      Ayz = squeeze(nanmean(adv(:,:,:)));    
    else
      Ayz = squeeze(adv(yzlayer,:,:,outfidx));
    end
    
      
    pcolor(YY(1:end-1,1:end-1),ZZ(1:end-1,1:end-1)/1000,Ayz);
    shading interp;


    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex'); 
   %}


% Converting axes from m to km 
  xxnew = xx/1000;
  yynew = yy/1000;
  zznew = zz/1000;


% xy plot of advection
  if xyplot==1

      if singlelayer==1
          advlayer = squeeze(adv(:,:,xylayer,outfidx));
      else%if layeraverage==1
          advlayer = squeeze(mean(adv(:,:,zlaymin:zlaymax),3,'omitnan'));
      end
      
      pcolor(xxnew(1:end-2),yynew(1:end-2),advlayer')
      shading interp; 
      xlabel('x (km)', 'interpreter','latex');
      ylabel('y (km)', 'interpreter','latex');
      caxis([-5e-7 5e-7])
      cmocean('balance')
      handle = colorbar;
      set(handle,'FontSize',fontsize);
      set(gca,'FontSize',fontsize);
      title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
      set(gca,'Position',plotloc);


% xz plot of advection 
  elseif xzplot==1

      if singlelayer==1
          advlayer = squeeze(adv(:,xzlayer,:,outfidx));
      else%if layeraverage==1
         advlayer = squeeze(mean(adv(:,ylaymin:ylaymax,:),2,'omitnan'));
         advlayer_u = squeeze(mean(u_dvdx(:,ylaymin:ylaymax,:),2,'omitnan'));
         advlayer_v = squeeze(mean(v_dvdy(:,ylaymin:ylaymax,:),2,'omitnan'));
         advlayer_w = squeeze(mean(w_dvdz(:,ylaymin:ylaymax,:),2,'omitnan'));
      end

      
      pcolor(xxnew(1:end-2),zznew(1:end-2),advlayer')
      shading interp;
      handle = colorbar;
      xlabel('x (km)', 'interpreter','latex')
      ylabel('Height $z$ (km)','interpreter','latex'); 
      cmocean('balance')
      caxis([-5e-7 5e-7])
      set(handle,'FontSize',fontsize);
      set(gca,'FontSize',fontsize);
      title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
      %set(gca,'Position',plotloc);

      % subplot of the u,v,w components 
      figure
      subplot(2,2,1)
      pcolor(xxnew(1:end-2),zznew(1:end-2),advlayer_u(1:end,1:end-2)')
      shading interp; 
      cmocean('balance')
      handle = colorbar; 
      xlabel('x (km)', 'interpreter','latex')
      ylabel('Height $z$ (km)','interpreter','latex')
      title('Averaged u advection')
      caxis([-5e-7 5e-7])
      handle=colorbar;
      set(handle,'FontSize',13);
      set(gca,'FontSize',13);
  

      subplot(2,2,2)
      pcolor(xxnew(1:end-2),zznew(1:end-2),advlayer_v(1:end-1,1:end-2)')
      shading interp; 
      cmocean('balance')
      handle = colorbar; 
      xlabel('x (km)','interpreter','latex')
      ylabel('Height $z$ (km)','interpreter','latex') 
      title('Averaged v advection')
      caxis([-5e-7 5e-7])
      handle=colorbar;
      set(handle,'FontSize',13);
      set(gca,'FontSize',13);
    

     subplot(2,2,3)
     pcolor(xxnew(1:end-2),zznew(1:end-2),advlayer_w(1:end-2,1:end)')
     shading interp; 
     cmocean('balance')
     handle = colorbar; 
     xlabel('x (km)', 'interpreter','latex')
     ylabel('Height $z$ (km)','interpreter','latex')
     title('Averaged w advection')
     caxis([-5e-7 5e-7])
     handle=colorbar;
     set(handle,'FontSize',13);
     set(gca,'FontSize',13);

     sgtitle(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');
  
   
      %{
      figure
      pcolor(xx(1:end-2),zz(1:end-2),advlayer_u(1:end,1:end-2)')
      shading interp; 
      cmocean('balance')
      handle = colorbar; 
      xlabel('x (km)')
      ylabel('Height $z$ (km)','interpreter','latex'); 
      %caxis([])
      figure
      pcolor(xx(1:end-2),zz(1:end-2),advlayer_v(1:end-1,1:end-2)')
      shading interp; 
      figure
      pcolor(xx(1:end-2),zz(1:end-2),advlayer_w(1:end-2,1:end)')
      shading interp; 
      %}
  
% yz plot of advection 

  elseif yzplot==1

      if singlelayer==1
            advlayer = squeeze(adv(yzlayer,:,:,outfidx));
      else%if layeraverage==1
            advlayer = squeeze(mean(adv(xlaymin:xlaymax,:,:),1,'omitnan'));
            advlayer_v = squeeze(mean(v_dvdy(xlaymin:xlaymax,:,:),1,'omitnan'));
      end
      
      pcolor(yynew(1:end-2),zznew(1:end-2),advlayer')
      shading interp;
      cmocean('balance')
      handle=colorbar;
      xlabel('Offshore $y$ (km)','interpreter','latex')
      ylabel('Height $z$ (km)','interpreter','latex'); 
      set(handle,'FontSize',fontsize);
      set(gca,'FontSize',fontsize);
      title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
      %set(gca,'Position',plotloc);
      caxis([-5e-7 5e-7])
      

      figure
      pcolor(yynew(1:end-2),zznew(1:end-2),advlayer_v(1:end,1:end-2)')
      shading interp;
      cmocean('balance')
      xlabel('Offshore $y$ (km)','interpreter','latex');
      ylabel('Height $z$ (km)','interpreter','latex'); 
      handle=colorbar;
      set(handle,'FontSize',13);
      set(gca,'FontSize',13);
      subtitle('Averaged v advection', 'Fontsize', 13)
      title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
      set(gca,'Position',plotloc);
      caxis([-5e-7 5e-7])

     
  



  end

     

    %cmocean('balance')
    %caxis([-5e-7 5e-7])
  

    
  %handle=colorbar;
  %set(handle,'FontSize',fontsize);
  %set(gca,'FontSize',fontsize);
  %title(['$t=',num2str(tdays(n),'%.1f'),'$ days'],'interpreter','latex');  
  %set(gca,'Position',plotloc);



  

  

