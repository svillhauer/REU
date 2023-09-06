%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.,
%%%
function nTimeSteps = setParams (inputpath,codepath,listterm,Nx,Ny,Nr)  

  %%% Load EOS utilities
 addpath /Users/sarahvillhauer/Desktop/MITgcm-master/MITgcm_SC/GSW;
  addpath /Users/sarahvillhauer/Desktop/MITgcm-master/MITgcm_SC/GSW/html;
  addpath /Users/sarahvillhauer/Desktop/MITgcm-master/MITgcm_SC/GSW/library;
  addpath /Users/sarahvillhauer/Desktop/MITgcm-master/MITgcm_SC/GSW/pdf;
    
  %%% TODO ADD SEA ICE AND ATMOSPHERIC FORCING; REMOVE SW, LW, SENS HEAT FLUXES      
  %%% TODO add sponges for waves?
  
  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%      
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;      
  fignum = 1;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% hours in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000; 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SIMULATION CONTROL PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nonHydrostatic = false; %%% Whether to run in nonhydrostatic mode
  use_seaIce = false; %%% Whether to run with sea ice (not yet implemented)
  use_3D = true; %%% Whether to run a 3D vs 2D simulation
  Ypoly = 0; %%% Latitudinal location 
  Wpoly = 5*m1km; %%% Latitudinal width 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  simTime = 10*t1day; %%% Simulation time  
  nIter0 = 0; %%% Initial iteration 
  Lx = 3*m1km; %%% Domain size in x 
  Ly = 3*m1km; %%% Domain size in y   
  H = 101.25; %33.75; %29.85; %29.75; %29.6500; %41.5; %%% Domain size in z 
  g = 9.81; %%% Gravity
  Omega = 2*pi*366/365/86400;  
  lat0 = 77; %%% Latitude at southern boundary
  f0 = -1.4e-4; %0; %2*Omega*sind(lat0); % f0 = 2*Omega*sind(lat0); %%% Coriolis parameter      
  rho0 = 1000;  
  
  %%% Diffusion parameters
  viscAh = 0; %%% Horizontal viscosity    
  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
%   viscA4Grid = 0.2; %%% Grid-dependent biharmonic viscosity    
%   viscC4smag = 0; %%% Smagorinsky biharmonic viscosity
%   difK4Tgrid = 0.2; %%% Grid-dependent biharmonic diffusivity
  viscA4Grid = 0; %%% Grid-dependent biharmonic viscosity    
  viscC4smag = 4.00000000e+00; %%% Smagorinsky biharmonic viscosity
  diffK4Tgrid = 0.0; %%% Grid-dependent biharmonic diffusivity
  viscAr = 0; %%% Vertical viscosity
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKrT = 1.3e-7; %%% Vertical temp diffusion     
  diffKhS = 0; %%% Horizontal salt diffusion
  diffKrS = 7.2e-10; %%% Vertical salt diffusion     
 smag3D_coeff=2e-02; %3D smag coefficient

  %%% Parameters related to periodic forcing
  %{ 
periodicExternalForcing = true;
  externForcingCycle = 2*simTime;  
  if (~periodicExternalForcing)
    externForcingPeriod = externForcingCycle;
    nForcingPeriods = 1;
  else
    externForcingPeriod = 0.01*t1day;
    nForcingPeriods = externForcingCycle/externForcingPeriod;
  end
  %}

  
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('useSmag3D',true,PARM_BOOL);
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);      
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL); 
  parm01.addParm('viscC2smag',0,PARM_REAL); 
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL); 
  parm01.addParm('smag3D_coeff',smag3D_coeff,PARM_REAL);

  %%% diffusivity
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);  
  parm01.addParm('diffKrS',diffKrS,PARM_REAL);
  parm01.addParm('diffKhS',diffKhS,PARM_REAL);  
  %%% advection schemes
  parm01.addParm('tempAdvScheme',33,PARM_INT);
  parm01.addParm('saltAdvScheme',33,PARM_INT);
%   parm01.addParm('tempAdvScheme',7,PARM_INT);
%   parm01.addParm('saltAdvScheme',7,PARM_INT);
  parm01.addParm('multiDimAdvection',true,PARM_BOOL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
%   parm01.addParm('convertFW2Salt',-1,PARM_REAL);
  %%% equation of state
  parm01.addParm('eosType','JMD95Z',PARM_STR);   
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',0,PARM_REAL); %2e-3
  %parm01.addParm('sideDragFactor',2,PARM_REAL);
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  parm01.addParm('rhonil',rho0,PARM_REAL);
  parm01.addParm('rhoConst',rho0,PARM_REAL);
 
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL);
  parm01.addParm('implicitDiffusion',~nonHydrostatic,PARM_BOOL);
  parm01.addParm('implicitViscosity',~nonHydrostatic,PARM_BOOL);
  parm01.addParm('nonHydrostatic',nonHydrostatic,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);  

  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);  
  parm01.addParm('debugLevel',-1,PARM_INT);
%   parm01.addParm('debugLevel',5,PARM_INT);

  %%%%% Vertical advection
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
%   parm01.addParm('momImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('tempImplVertAdv',true,PARM_BOOL);
%   parm01.addParm('saltImplVertAdv',true,PARM_BOOL);


  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);

  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
  parm02.addParm('cg3dMaxIters',300,PARM_INT);  
  parm02.addParm('cg3dTargetResidual',1e-7,PARM_REAL);  
 
  %%% PARM03
  %parm03.addParm('alph_AB',1/2,PARM_REAL);
  %parm03.addParm('beta_AB',5/12,PARM_REAL);
  parm03.addParm('momDissip_In_AB',false,PARM_BOOL);
  
  %%%%%%% Added to make SEAICE work
  %parm03.addParm('momForcingOutAB',1,PARM_INT);
  parm03.addParm('tracForcingOutAB',1,PARM_INT);


  %%%%%%%
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',0.05*t1day,PARM_REAL);
  parm03.addParm('pChkptFreq',0.05*t1day,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',0*t1year,PARM_REAL);
  parm03.addParm('monitorFreq',1*t1day,PARM_REAL);
  parm03.addParm('cAdjFreq',0,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);
 %{
  if (periodicExternalForcing)
    parm03.addParm('periodicExternalForcing',periodicExternalForcing,PARM_BOOL);
    parm03.addParm('externForcingPeriod',externForcingPeriod,PARM_REAL);
    parm03.addParm('externForcingCycle',externForcingCycle,PARM_REAL);
  end

 %}
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
 
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    
    

  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);
  yy = cumsum(dy);
  yy = yy-mean(yy);
  
  %%% Zonal grid
%{
  refinement = 1.5; %# controls spacing near surface (higher means finer spaced)
  stretching = 10;  %# controls rate of stretching at bottom
  %# "Warped" height coordinate
  for k=1:Nx/2
      h(k) = (Nx/2 + 1 - k) / Nx/2;
      %# Linear near-surface generator
      zeta(k) = 1 + (h(k) - 1) / refinement;
      %# Bottom-intensified stretching function
      sigma(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching));
  end
  for k=1:Nx/2
      %# Generating function
      xx(Nx/2-k+1) =  Lx/2 * (zeta(Nx/2-k+1) * sigma(Nx/2-k+1) -1) ;
  end
  half_x = -flip(diff(xx));
  dx_1 =[half_x(1) half_x]; 
  dx = [dx_1 flip(dx_1)];

  plot(dx)
%} 
  dxmin=8; %m, smallest dx
  dxmax = 24; %m, largest dx
halfwidth= 125/2;
xloc_stretch = (halfwidth+dxmax*9)/dxmin;
   xind=[1:1:Nx/2];
  dx = dxmin* tanh(0.05* (xind-xloc_stretch)) + (dxmin+dxmax)/2;
  dx = [flip(dx) dx];
xx = [-flip(cumsum(dx)) cumsum(dx)];
  
%{  
  if (use_3D)
    dx = Lx/Nx*ones(1,Nx);  
    xx = cumsum(dx);
    xx = xx-mean(xx);
  else
    dx = dy(1);
    xx = 0;
  end
%}
  
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
 %dz = H/Nr*ones(1,Nr);
%dz = [ones(1,round(0.7*Nr))*.3 linspace(.3,2,Nr-round(0.7*Nr))]; 
%dz = [ones(1,45)*.3 linspace(.3,2,15)];
%resolution in the top 70% of Nr gridopints is 0.3 m and resolution in
%bottom 30% of NR gridpts linearly telescopes from 0.3 to 2m vertical res.
%dz = [ones(1,55)*.3 linspace(.3,2,15)]; %[ones(1,53)*.3 linspace(.3,2,20)];  %[ones(1,53)*.3 linspace(.3,1.4,17)];
%dz = [ones(1,200)*.3]; 
%dz=[ones(1,100)]* H/Nr;

%dz = [ones(1,45)*.3 linspace(.3,2,15)];

 % zz = -cumsum((dz+[0 dz(1:end-1)])/2);


  refinement = 8; %# controls spacing near surface (higher means finer spaced)
  stretching = 2;  %# controls rate of stretching at bottom
  %# "Warped" height coordinate
  for k=1:Nr
      h(k) = (Nr + 1 - k) / Nr;
      %# Linear near-surface generator
      zeta(k) = 1 + (h(k) - 1) / refinement;
      %# Bottom-intensified stretching function
      sigma(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching));
  end
  for k=1:Nr
      %# Generating function
      zz(Nr-k+1) =  H* (zeta(Nr-k+1) * sigma(Nr-k+1) -1) ;
  end
  dz=[ -zz(2) -diff(zz)];

  %%% Store grid spacings
  parm04.addParm('delX',dx,PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);          
  

  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   

  %%% Flat bottom  
  h = -H*ones(Nx,Ny);
  
 %h(:,1)=zeros;
 %h(:,end)=zeros;

  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
  
  
  
  
  
  
  
  

   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% REFERENCE TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
  %%% Quasi-tanh-shaped T/S profiles
 %%%%South BC
 shelfthickness= 50; %idea here is to lower T and S profiles by shelfthickness and to make the surface shelfthickness layer relatively unstratified
  Zpyc = -shelfthickness; % -10-shelfthickness; %southern/inflow boundary pycnocline mid-depth (depth scale)
  Wpyc = 10; %5 %pycnocline width scale
  %Smin = 34.0350; %34.2; %34.2;
  %Smax = 34.4175; %34.22; %34.22; %34.7
  %Tmin = -0.6825; %-0.6; %-0.6 ;% 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
  %Tmax= -0.2576; %-0.15; %-0.15;  
  gam_h = 0.01;
  %tRef = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); %repmat(Tmin,1,60); %Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
  %saltmx = Smin:(Smax-Smin)/(length(zz)-1):Smax;

  %sRef =  Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); %interp1(depthmx,saltmx,zz,'linear'); %Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); 
 ctd26 = readtable('PetermannIceShelf_26_1.dat');
  ctd16 = readtable('PetermannIceShelf_03_1.dat'); 
  ctd26new = table2array(ctd26);
  ctd16new = table2array(ctd16);
ctd16new = sortrows(ctd16new);
ctd26new = sortrows(ctd26new);
ctd16new = ctd16new(4813:8658,:);
ctd26new = ctd26new(3512:4507,:);
[~,u1] = unique(ctd16new(:,1),'stable');
[~,u2] = unique(ctd26new(:,1),'stable');
ctd16new = ctd16new(u1,:);
ctd26new = ctd26new(u2,:);

new16x = ctd16new(1,1):min(diff(ctd16new)):ctd16new(end,1);
salt16 = interp1(ctd16new(:,1), ctd16new(:,4), new16x);
temp16 = interp1(ctd16new(:,1), ctd16new(:,2), new16x);
new26x = ctd26new(1,1):min(diff(ctd26new)):ctd26new(end,1);
salt26 = interp1(ctd26new(:,1), ctd26new(:,4), new26x);
temp26 = interp1(ctd26new(:,1), ctd26new(:,2), new26x);

depth16mx = 0:(-H)/(length(new16x)-1):-H;
depth26mx = 0:(-H)/(length(new26x)-1):-H;



  tRef = interp1(depth16mx,temp16,zz,'spline'); 
  sRef =  interp1(depth16mx,salt16,zz,'spline');

  %{
  % finding gridpoint halfway down domain 
  for i = 1:Nr 
      sum_dz = sum(dz(1:i));
      if sum_dz > H/2
          break;
      end
      z_half_index = i;
  end

   tRef_1 = tRef(1:z_half_index);
   tRef_2 = tRef(z_half_index+1:length(tRef));
   tRef_1 = smooth(tRef_1);
   tRef_2 = smooth(tRef_2);

    sRef_1 = sRef(1:z_half_index);
   sRef_2 = sRef(z_half_index+1:length(sRef));
   sRef_1 = smooth(sRef_1);
   sRef_2 = smooth(sRef_2);

   zz_1 = zz(1:z_half_index);
   zz_2 = zz(z_half_index+1:length(zz));
  %}
 

  tRef_1 = tRef(1:Nr/2);
   tRef_2 = tRef(Nr/2+1:length(tRef));
   tRef_1 = smooth(tRef_1);
   tRef_2 = smooth(tRef_2);

    sRef_1 = sRef(1:Nr/2);
   sRef_2 = sRef(Nr/2+1:length(sRef));
   sRef_1 = smooth(sRef_1);
   sRef_2 = smooth(sRef_2);

   zz_1 = zz(1:Nr/2);
   zz_2 = zz(Nr/2+1:length(zz));


%% first half of profile
   % Create a polynomial fit (you can adjust the degree as needed)
degree = 3; % Adjust the degree as needed
p = polyfit(zz_1, tRef_1, degree);
p_3 = polyfit(zz_1, sRef_1, degree);



% Generate a set of points for the fitted curve
%tFit = linspace(min(tRef_1), max(tRef_1), Nr/2); % Adjust the number of points as needed
tFit1 = polyval(p, zz_1);
sFit1 = polyval(p_3,zz_1);


%% second half of profile 
p_2 = polyfit(zz_2, tRef_2, degree);
p_4 = polyfit(zz_2,sRef_2,degree);

% Generate a set of points for the fitted curve
%tFit = linspace(min(tRef_1), max(tRef_1), Nr/2); % Adjust the number of points as needed
tFit2 = polyval(p_2, zz_2);
sFit2 = polyval(p_4,zz_2);




tRef = [tFit1 tFit2];
tRef = smooth(tRef); 
tRef = tRef';
sRef = [sFit1 sFit2];
sRef = smooth(sRef);
sRef = sRef';




  %tRef = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
% sRef =  Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));



%%%%North BC
   Zpyc = -shelfthickness; % -10-shelfthickness; %%northern/outflow boundary pycnocline mid-depth
  Wpyc = 10;
  %{
  Smin = 33.9911; %34.2-.2;%33.95;
  Smax = 34.4091; %34.02; %34.7 %34.02
  Tmin = -0.6764; %-.61;%-0.65 ;%0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
  Tmax=  -0.2586; %-0.15;  
  %}
  gam_h = 0.01;
  %tRefout = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); %repmat(Tmin,1,60); %Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
  %saltmx = Smin:(Smax-Smin)/(length(zz)-1):Smax;
  depthmx = 0:(-length(ctd26new(:,1)))/(length(dz)-1):-H;
  %sRefout = Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); %interp1(depthmx,saltmx,zz,'linear'); %Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); 


  tRefout = interp1(depth26mx,temp26,zz,'spline'); 
  sRefout =  interp1(depth26mx,salt26,zz,'spline');


diffsalt = (sRefout - sRef)./10;
difftemp = (tRefout - tRef)./10;

sRefout = sRef +diffsalt; 
tRefout = tRef +difftemp;

 tRefout_1 = tRefout(1:Nr/2);
   tRefout_2 = tRefout(Nr/2+1:length(tRefout));
   tRefout_1 = smooth(tRefout_1);
   tRefout_2 = smooth(tRefout_2);

    sRefout_1 = sRefout(1:Nr/2);
   sRefout_2 = sRefout(Nr/2+1:length(sRefout));
   sRefout_1 = smooth(sRefout_1);
   sRefout_2 = smooth(sRefout_2);

   zz_1 = zz(1:Nr/2);
   zz_2 = zz(Nr/2+1:length(zz));

%% first half of profile
   % Create a polynomial fit (you can adjust the degree as needed)
degree = 3; % Adjust the degree as needed
p_5 = polyfit(zz_1, tRefout_1, degree);
p_7 = polyfit(zz_1, sRefout_1, degree);



% Generate a set of points for the fitted curve
%tFit = linspace(min(tRef_1), max(tRef_1), Nr/2); % Adjust the number of points as needed
toutFit1 = polyval(p_5, zz_1);
soutFit1 = polyval(p_7,zz_1);


%% second half of profile 
p_6 = polyfit(zz_2, tRefout_2, degree);
p_8 = polyfit(zz_2,sRefout_2,degree);

% Generate a set of points for the fitted curve
%tFit = linspace(min(tRef_1), max(tRef_1), Nr/2); % Adjust the number of points as needed
toutFit2 = polyval(p_6, zz_2);
soutFit2 = polyval(p_8,zz_2);




tRefout = [toutFit1 toutFit2];
tRefout = smooth(tRefout); 
tRefout = tRefout';
sRefout = [soutFit1 soutFit2];
sRefout = smooth(sRefout);
sRefout = sRefout';

s_avg = (sRef + sRefout)./2;  
t_avg = (tRef + tRefout)./2;

sRef = s_avg;
sRefout = s_avg; 
tRef = t_avg;
tRefout = t_avg; 

  

%{
Tmin = min(tRefout);
Tmax = max(tRefout);
Smin = min(sRefout);
Smax = max(sRefout);
%}
%{
 Tmax = -0.5 + mean(difftemp);
 Tmin = Tmax - 0.4; 
 Smax = 34.5 + mean(diffsalt);
 Smin = Smax- 0.4; 
  

tRefout = Tmin + (Tmax-Tmin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2)); 
sRefout = Smin + (Smax-Smin)*0.5*(1+0.5*(sqrt((1-(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2)-sqrt((1+(zz-Zpyc)/Wpyc).^2 + 4*gam_h*((zz-Zpyc)/Wpyc).^2))/(1+gam_h)^(1/2));
%}

%[a, b, c, d] = fitTanhFunction(tRefout, zz);
%[a, b, c, d] = fitTanhFunction(sRefout, zz);

% interp1(depthmx,saltmx,zz,'linear');
 %sRefout(round(Nr/4):end)=linspace(sRefout(round(Nr/4)),sRefout(end),Nr-round(Nr/4)+1);
%sRefout(round(Nr/4):end)=sRef(round(Nr/4):end);
 %sRefout=smoothdata(sRefout);%

 %tRefout(round(Nr/4):end)=linspace(tRefout(round(Nr/4)),tRefout(end),Nr-round(Nr/4)+1);
%tRefout(round(Nr/4):end)=tRef(round(Nr/4):end);
% tRefout=smoothdata(tRefout);

  %%% Quasi-linear near-surface stratification
%   Zpyc = 0;
%   Hpyc = 300;  
%   Smin = 34.3;
%   Smax = 34.7;
%   Tmin = 0.0901-0.0575*Smin; %%% MITgcm surface freezing temperature
%   Tmax= 1;  
%   gam_h = 0.01;
%   tRef = Tmin + (Tmax-Tmin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);
%   sRef = Smin + (Smax-Smin)*0.5*(sqrt((1-zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2)-sqrt((1+zz/Hpyc).^2 + 4*gam_h*(zz/Hpyc).^2))/(1+gam_h)^(1/2);
  
  %%% Plot the reference temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tRef,zz); hold on
    plot(tRefout,zz);
    xlabel('\theta_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference temperature');
  end
  
  %%% Plot the reference salinity
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(sRef,zz); hold on
     plot(sRefout,zz);
    xlabel('S_r_e_f');
    ylabel('z','Rotation',0);
    title('Reference salinity');
  end
  
  
  
  
  
  



  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Check Brunt-Vaisala frequency using full EOS
  pp = - zz;
%   SA = gsw_SA_from_SP(sRef,pp,-50,-64);  
%   CT = gsw_CT_from_pt(sRef,tRef);
%   [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
 [N2 pp_mid] = gsw_Nsquared(sRef,tRef,pp);
  dzData = zz(1:end-1)-zz(2:end);

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(max(N2,0)); 
  Cig = sum(N.*dzData);
  Rd = Cig./(pi*abs(f0));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2,-pp_mid/1000);
    xlabel('N^2 (km)');
    ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
  end
  
  if (showplots)
    ['First baroclinic Rossby deformation radius: ' num2str(Rd)]
  end
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 0.225*.84/1;  
  Wmax = 0.045*1.42857143/8.00*.84/1; %%% NEEDS TO BE INCREASED FOR DEEP PYCNOCLINES
  %%% Max gravity wave speed 
  cmax = 0;%max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Advective CFL
  deltaT_adv= min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% Vertical advective CFL
  deltaT_vadv = min(0.5*dz/Wmax);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  if (nonHydrostatic)
    deltaT = min([deltaT deltaT_vadv]);
  end
  deltaT = round(deltaT);
  deltaT= deltaT/3; %ad hoc: we found that the normal dT wasn't working (approx. dT=14s)

  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% TRACER DIFFUSION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %diffK4T = diffK4Tgrid * max([dx dy])^4 / (32*deltaT);
  %diffK4S = diffK4T;
  %parm01.addParm('diffK4T',diffK4T,PARM_REAL); 
  %parm01.addParm('diffK4S',diffK4S,PARM_REAL); 
 
  
  
  
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %{
  %%% Fixed salt flux in the BW formation region  
  salt_flux = zeros(Nx,Ny,nForcingPeriods);    
  tt = zeros(1,nForcingPeriods);
  for n=1:nForcingPeriods   
    t = (n-1)*externForcingPeriod;
    tt(n) = t;
    for j=1:Ny      
      for i=1:Nx             
        if (abs(yy(j)-Ypoly)<Wpoly/2)
          
          %%% Linear decrease of salt flux from salt_amp to zero after Tpoly days
%           salt_amp = -1e-1;
%           Tpoly = 3*t1day;
%           salt_flux(i,j,n) = salt_amp - salt_amp*min(1,t/Tpoly);             
          
          %%% Empirical ice growth equation from Anderson (1969)
          C1 = 1e4; %%% Constants from Anderson (1969) 
          C2 = 510;
          C3 = 6.7/86400;
          C4 = C2/(2*C1);
          C5 = C3/C1;
          Theta = 40; %%% Difference between atmospheric temperature and ocean freezing temperature 
          rho_i = 920; %%% Density of ice
          S_i = 5; %%% Salinity of ice
          e = -C4 + sqrt(C4^2+C5*t*Theta); %%% Ice thickness (m) as a function of time
          dedt = C3*Theta/(2*C1*e+C2); %%% Ice growth rate (m/s) as a function of time
          salt_flux(i,j,n) = -dedt*rho_i*(Smin-S_i); %%% Equivalent salt flux (g/m^2/s) as a function of time
          
        end
      end
    end         
  end  
  
 
  
  %%% Plot the surface salt flux
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
%     contourf(X,Y,squeeze(salt_flux(:,:,1)));
    colorbar;
%     plot(yy,squeeze(salt_flux(1,:,1)));    
%     xlabel('y');
    plot(tt,squeeze(salt_flux(1,Ny/2,:)));    
    xlabel('t');
    ylabel('salt flux');
    title('Surface salt flux');
  end  
  
  %%% Save as parameters  
  writeDataset(salt_flux,fullfile(inputpath,'saltFluxFile.bin'),ieee,prec);
  parm05.addParm('saltFluxFile','saltFluxFile.bin',PARM_STR);
  











  %}
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = true;
  useRBCsalt = true;
  useRBCuVel = true;
  useRBCvVel = true;
  tauRelaxT = 0.1*t1day;
  tauRelaxS = 0.1*t1day;
  tauRelaxU = -1/f0;%0.05*t1day;
  tauRelaxV = -1/f0;%0.05*t1day;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION TEMPERATURE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %tRefout (outflowing/north BC T/S) defined above
%{
  for i=1:size(tRef,2)
Tmat(i,:)=linspace(tRef(1,i), tRefout(1,i),Ny);
Smat(i,:)=linspace(sRef(1,i), sRefout(1,i),Ny);
  end
  
 Tmatfull=repmat(Tmat,[1 1 Nx]);
Tmatfinal=permute(Tmatfull,[3,2,1]);
Smatfull=repmat(Smat,[1 1 Nx]);
Smatfinal=permute(Smatfull,[3,2,1]);

%}
%%% Setting lower quarter of domain to be nudging region 


%{
lqrt = zz(end)-(zz(end)/4);

saltnorth = sRefout(end-10:end);
saltsouth = sRef(end-10:end);
tempnorth = tRefout(end-10:end);
tempsouth = tRef(end-10:end); 

tempfinal = -2*ones(1,11);
saltfinal = 35*ones(1,11);
%}
%{
  for i=1:size(tempsouth,2)
Tmat(i,:)=linspace(tempsouth(1,i), tempnorth(1,i),Ny);
Smat(i,:)=linspace(saltsouth(1,i), saltnorth(1,i),Ny);
  end
%}
  
  %{
%northern conditions
  for i=1:size(tempnorth,2)
Tmatnorth(i,:)=linspace(tempnorth(1,i), tempfinal(1,i),Ny);
Smatnorth(i,:)=linspace(saltnorth(1,i), saltfinal(1,i),Ny);
  end
%}
 %final 

 tRef=linspace(-0.7, -.25,Nr);
 sRef=linspace(34,34.15,Nr);

 % for i=1:size(tRef,2)
%Tmat(i,:)=linspace(tRef(1,i), tRefout(1,i),Ny);
%Smat(i,:)=linspace(sRef(1,i), sRefout(1,i),Ny);
 % end

  for i=1:size(tRef,2)
Tmat(i,:)=linspace(tRef(1,i), tRef(1,i),Ny);
Smat(i,:)=linspace(sRef(1,i), sRef(1,i),Ny);
  end

 Tmatfull = repmat(Tmat,[1 1 Nx]);
 Tmatfinal = permute(Tmatfull,[3,2,1]);
 

Smatfull = repmat(Smat,[1 1 Nx]);
Smatfinal = permute(Smatfull,[3,2,1]);

vvelfinal = 0*ones(Nx,Ny,Nr);
uvelfinal = 0.1*ones(Nx,Ny,Nr);
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Random noise amplitude
  tNoise = 0;%0.01;  
  sNoise = 0;%0.001;
      
  %%% Align initial temp with background
  %changed initial data to be the relaxed 3D fields (from RBCS)
  hydroTh = Tmatfinal;%ones(Nx,Ny,Nr);
  hydroSa = Smatfinal;%ones(Nx,Ny,Nr);
  %for k=1:1:Nr
  %  hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tRef(k);
  %  hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sRef(k);
  %end
  
  %%% Add some random noise
  hydroTh = hydroTh + tNoise*(0*rand(Nx,Ny,Nr)-0);
  hydroSa = hydroSa + sNoise*(0*rand(Nx,Ny,Nr)-0);


  %%% Adiing initial velocity 
  
  
  vvelinitial=zeros(Nx,Ny,Nr);
  for i = 1:Nr
      vvelinitial(:,:,i) = vvelfinal(:,:,i) - (1/2)*(vvelfinal(:,:,i)/Nr)*(Nr-i);
  end
  
  %vvelinitial = rdmdsWrapper(fullfile('/Users/sarahvillhauer/Desktop/MITgcm-master/MITgcm_SC/experiments/channel/icetopotest','results','VVEL_inst'),249231); 

  %vvelinitial = vvelfinal;

  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);

  
  writeDataset(vvelinitial,fullfile(inputpath,'vVelInitFile.bin'),ieee,prec); 
  parm05.addParm('vVelInitFile','vVelInitFile.bin',PARM_STR);

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  


  
  %%% Set relaxation temp equal to surface freezing temperature
  % temp_relax = Tmin*ones(Nx,Ny,Nr);  

  %%% Save as parameters
  writeDataset(Tmatfinal,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);  

  writeDataset(Smatfinal,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR); 

  writeDataset(uvelfinal,fullfile(inputpath,'sponge_uvel.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxUFile','sponge_uvel.bin',PARM_STR); 

  %writeDataset(vvelfinal,fullfile(inputpath,'sponge_vvel.bin'),ieee,prec); 
  %rbcs_parm01.addParm('relaxVFile','sponge_vvel.bin',PARM_STR); 
  
  %%%%%%%%%%%%%%%%%%%%%  
  %%%%% RBCS MASK %%%%%
  %%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Mask is zero everywhere by default, i.e. no relaxation
  
msku=zeros(Nx,Ny,Nr);
mskv=zeros(Nx,Ny,Nr);
  mskt=zeros(Nx,Ny,Nr);
  msks=zeros(Nx,Ny,Nr);

    
  %%% Mask only nonzero at surface in the polynya
 %{
  for j=1:Ny      
    for i=1:Nx             
      if (yy(j)<-0.4*Ly) || (yy(j)>0.4*Ly) 
        msk(:,j,:) = 1;
      end
    end
  end         
 %}  

      for i = 1:Nr
          if zz(i)>-H % -0.9*H
              msku(:,:,i) = 1;
          end
      end

       for i = 1:Nr
          if zz(i)< -0.9*H % -0.9*H
              mskt(:,:,i) = 1;
              msks(:,:,i) = 1;
          end
      end
  %%% Save as an input parameter
  writeDataset(mskt,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR); 

  writeDataset(msks,fullfile(inputpath,'rbcs_salt_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_salt_mask.bin',PARM_STR); 

  writeDataset(msku,fullfile(inputpath,'rbcs_uvel_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskUFile','rbcs_uvel_mask.bin',PARM_STR); 

  writeDataset(mskv,fullfile(inputpath,'rbcs_vvel_mask.bin'),ieee,prec); 
 rbcs_parm01.addParm('relaxMaskVFile','rbcs_vvel_mask.bin',PARM_STR); 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  
  
  

 
  
  

  
  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%% SEA ICE $ %%%%%%%%%%%
% %   
% %   
% %   % to store parameter names and values
%   seaice_parm01 = parmlist;
%   SEAICE_PARM = {seaice_parm01};
% % %   
% % %   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%% SEA ICE  %%%%%%%%%%%%
%     %%%%%%%% PARAMETERS %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Oringinal albedos from llc_1/48th  other values from llc_2160
% % 
%   SEAICEwriteState   = true;
%   SEAICEuseDYNAMICS  = true;
%   SEAICE_multDim     = 7;
% %   SEAICE_dryIceAlb   = 0.8783;
%   SEAICE_dryIceAlb   = 0.8509;
% %   SEAICE_wetIceAlb   = 0.7869;
%   SEAICE_wetIceAlb   = 0.7284;
% %   SEAICE_drySnowAlb  = 0.9482;
%   SEAICE_drySnowAlb  = 0.7754;
% %   SEAICE_wetSnowAlb  = 0.8216;
%   SEAICE_wetSnowAlb  = 0.7753;
%   SEAICE_waterDrag   = 5.5399;
%   SEAICE_drag        = 0.002;
% %   HO                 = 0.1;
%   HO                 = .05;
% 
% 
%   SEAICE_no_slip          = false;
% %   SEAICE_no_slip          = true;
% 
%   SEAICEadvScheme         = 7;
% %   SEAICEadvScheme         = 33;
% 
% 
%   %%%SOSEdoesn't have a seaice dataset for salinity, they used this value
%   %%%in their estimate
%   
%   LSR_ERROR               = 1.0e-4;
% %   LSR_ERROR               = 2.0e-4; 
%   MIN_ATEMP               = -40;
%   MIN_TICE                = -40;
%   SEAICE_area_reg         = 0.15;
%   SEAICE_hice_reg         = 0.1;
%   IMAX_TICE               = 6;
%   SEAICE_EPS		      = 1.0e-8;
% %   SEAICE_EPS              = 2.0e-9;
%   SEAICE_doOpenWaterMelt  = true;
%   SEAICE_areaLossFormula  = 1;
%   SEAICE_wetAlbTemp       = 0.0;
%   SEAICE_saltFrac         = 0.3;
% %   SEAICE_frazilFrac       = 0.003;
%  SEAICE_frazilFrac       = 0.01;
% %   SEAICE_frazilFrac       = 1.0;
% 
%   
%   seaice_parm01.addParm('LSR_ERROR',LSR_ERROR,PARM_REAL);
%   seaice_parm01.addParm('SEAICEwriteState',SEAICEwriteState,PARM_BOOL);
%   seaice_parm01.addParm('SEAICEuseDYNAMICS',SEAICEuseDYNAMICS,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_multDim',SEAICE_multDim,PARM_INT);
%   seaice_parm01.addParm('SEAICE_dryIceAlb',SEAICE_dryIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetIceAlb',SEAICE_wetIceAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drySnowAlb',SEAICE_drySnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_wetSnowAlb',SEAICE_wetSnowAlb,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_waterDrag',SEAICE_waterDrag,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_drag',SEAICE_drag,PARM_REAL);
%   seaice_parm01.addParm('HO',HO,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_dryIceAlb_south',SEAICE_dryIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetIceAlb_south',SEAICE_wetIceAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drySnowAlb_south',SEAICE_drySnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_wetSnowAlb_south',SEAICE_wetSnowAlb_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_waterDrag_south',SEAICE_waterDrag_south,PARM_REAL);
% %   seaice_parm01.addParm('SEAICE_drag_south',SEAICE_drag_south,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_no_slip',SEAICE_no_slip,PARM_BOOL);
% %   seaice_parm01.addParm('SEAICE_salinity',SEAICE_salinity,PARM_REAL);
%   seaice_parm01.addParm('SEAICEadvScheme',SEAICEadvScheme,PARM_INT);
%   seaice_parm01.addParm('MIN_ATEMP',MIN_ATEMP,PARM_REAL);
%   seaice_parm01.addParm('MIN_TICE',MIN_TICE,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_area_reg',SEAICE_area_reg,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_hice_reg',SEAICE_hice_reg,PARM_REAL);
%   seaice_parm01.addParm('IMAX_TICE',IMAX_TICE,PARM_INT);
%   seaice_parm01.addParm('SEAICE_EPS',SEAICE_EPS,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_doOpenWaterMelt',SEAICE_doOpenWaterMelt,PARM_BOOL);
%   seaice_parm01.addParm('SEAICE_areaLossFormula',SEAICE_areaLossFormula,PARM_INT);
%   seaice_parm01.addParm('SEAICE_wetAlbTemp',SEAICE_wetAlbTemp,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_saltFrac',SEAICE_saltFrac,PARM_REAL);
%   seaice_parm01.addParm('SEAICE_frazilFrac',SEAICE_frazilFrac,PARM_REAL);
% 
%   
%   seaice_parm01.addParm('HeffFile',HeffFile,PARM_STR);
%   seaice_parm01.addParm('AreaFile',AreaFile,PARM_STR);
%   seaice_parm01.addParm('HsnowFile',HsnowFile,PARM_STR);
%   seaice_parm01.addParm('HsaltFile',HsaltFile,PARM_STR);
%   seaice_parm01.addParm('uIceFile',uIceFile,PARM_STR);
%   seaice_parm01.addParm('vIceFile',vIceFile,PARM_STR);
%   
%   
%      
% %   
% %   
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%% WRITE THE 'data.seaice' FILE %%%
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% %   
%   write_data_seaice(inputpath,SEAICE_PARM,listterm,realfmt);  
%   
%   
%   
%   
%   
%   
%   
%   
%   
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%EXF PKG%%%%%
%     %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%
%     
%     
%   %%% To store parameter names and values
%     
%      EXF_NML_01 = parmlist;
%      EXF_NML_02 = parmlist;
%      EXF_NML_03 = parmlist;
%      EXF_NML_04 = parmlist;
%      EXF_NML_OBCS = parmlist;
% 
%      EXF_PARM = {EXF_NML_01,EXF_NML_02,EXF_NML_03,EXF_NML_04,EXF_NML_OBCS};  
%     
%     
% 
% 
% % %     EXF_NML_01
% 
%   	exf_albedo        = 0.15;
%  	exf_scal_BulkCdn  = 1.015;
%  	exf_iprec         = 64;  
%  	useExfYearlyFields= false;
%  	useExfCheckRange  = false;
%  	useRelativeWind   = true;
% %  	useRelativeWind   = false;
%     repeatPeriod      = 283996800.0;
% %     exf_offset_atemp =  273.16;
%     
%     
% %%%runoff from ERA is in hours, need to convert to seconds
% %     exf_inscal_runoff = 1.14e-04;
%     
%     
%     apressurefile     = 'apressurefile.bin';
%     atempfile         = 'atempfile.bin';
%     aqhfile           = 'aqhfile.bin';
%     uwindfile         = 'uwindfile.bin';
%     vwindfile         = 'vwindfile.bin';
%     precipfile        = 'precipfile.bin';
%     swdownfile        = 'swdownfile.bin';
%     lwdownfile        = 'lwdownfile.bin';
% %     runofffile        = 'runofffile.bin';
%     
%     
% % "*period=-12" specifies monthly-mean forcing
%    apressurestartdate1 = 20070101;
%    apressurestartdate2 = 000000;
% %    apressureperiod     = 86400.0;
%    apressureperiod     = 10800.0;
%     
%     aqhstartdate1 = 20070101;
%     aqhstartdate2 = 000000;
% %     aqhperiod = 86400.0;
%     aqhperiod           = 10800.0;
% 
%     atempstartdate1 = 20070101;
%     atempstartdate2 = 000000;
% %     atempperiod = 86400.0;
%     atempperiod         = 10800.0;
%  
% 
%     uwindstartdate1 = 20070101;
%     uwindstartdate2 = 000000;
% %    uwindperiod = 86400.0;
%     uwindperiod         = 10800.0;
%  
%     vwindstartdate1 = 20070101;
%     vwindstartdate2 = 000000;
% %     vwindperiod = 86400.0;
%     vwindperiod         = 10800.0; 
%  
%     precipstartdate1 = 20070101;
%     precipstartdate2 = 000000;
% %     precipperiod = 86400.0;
%     precipperiod        = 10800.0; 
%  
% 
%     swdownstartdate1 = 20070101;
%     swdownstartdate2 = 000000;
% %     swdownperiod = 86400.0;
%     swdownperiod        = 10800.0;
% % 
% 
%     lwdownstartdate1 = 20070101;
%     lwdownstartdate2 = 000000;
% %     lwdownperiod = 86400.0;
%     lwdownperiod        = 10800.0;
% 
% 
% %     runoffstartdate1 = 20070101;
% %     runoffstartdate2 = 000000;
% %     runoffperiod = 2592000.0;
% %    runoffperiod        = 10800.0;
%  
%    precip_lon0 = xmin;
%    precip_lon_inc = dmxg(1);
%    precip_lat0 = ymin;
%    precip_lat_inc = dmyg;
%    precip_nlon = Nx;
%    precip_nlat = Ny;
%    
%    atemp_lon0 = xmin;
%    atemp_lon_inc = dmxg(1);
%    atemp_lat0 = ymin;
%    atemp_lat_inc = dmyg;
%    atemp_nlon = Nx;
%    atemp_nlat = Ny;
%    
%    apressure_lon0 = xmin;
%    apressure_lon_inc = dmxg(1);
%    apressure_lat0 = ymin;
%    apressure_lat_inc = dmyg;
%    apressure_nlon = Nx;
%    apressure_nlat = Ny;
%     
%    aqh_lon0 = xmin;
%    aqh_lon_inc = dmxg(1);
%    aqh_lat0 = ymin;
%    aqh_lat_inc = dmyg;
%    aqh_nlon = Nx;
%    aqh_nlat = Ny;
%    
%    uwind_lon0 = xmin;
%    uwind_lon_inc = dmxg(1);
%    uwind_lat0 = ymin;
%    uwind_lat_inc = dmyg;
%    uwind_nlon = Nx;
%    uwind_nlat = Ny;
%    
%    vwind_lon0 = xmin;
%    vwind_lon_inc = dmxg(1);
%    vwind_lat0 = ymin;
%    vwind_lat_inc = dmyg;
%    vwind_nlon = Nx;
%    vwind_nlat = Ny;
%    
%    swdown_lon0 = xmin;
%    swdown_lon_inc = dmxg(1);
%    swdown_lat0 = ymin;
%    swdown_lat_inc = dmyg;
%    swdown_nlon = Nx;
%    swdown_nlat = Ny;
%    
%    lwdown_lon0 = xmin;
%    lwdown_lon_inc = dmxg(1);
%    lwdown_lat0 = ymin;
%    lwdown_lat_inc = dmyg;
%    lwdown_nlon = Nx;
%    lwdown_nlat = Ny;
% 
% %    runoff_lon0 = xmin;
% %    runoff_lon_inc = dmxg(1);
% %    runoff_lat0 = ymin;
% %    runoff_lat_inc = dmyg;
% %    runoff_nlon = Nx;
% %    runoff_nlat = Ny;
%   
%  

%      
% %%%%%%%%%%%%%%%%%%%% ADD PARAMETERS
% 
%   EXF_NML_01.addParm('exf_albedo',exf_albedo,PARM_INT);
%   EXF_NML_01.addParm('exf_scal_BulkCdn',exf_scal_BulkCdn,PARM_REAL);
%   EXF_NML_01.addParm('exf_iprec',exf_iprec,PARM_INT);
%   EXF_NML_01.addParm('useExfYearlyFields',useExfYearlyFields,PARM_BOOL);
%   EXF_NML_01.addParm('useExfCheckRange',useExfCheckRange,PARM_BOOL);
%   EXF_NML_01.addParm('useRelativeWind',useRelativeWind,PARM_BOOL);
%   EXF_NML_01.addParm('repeatPeriod',repeatPeriod,PARM_REAL);
% %   EXF_NML_03.addParm('exf_offset_atemp',exf_offset_atemp,PARM_REAL);
% %   EXF_NML_03.addParm('exf_inscal_runoff',exf_inscal_runoff,PARM_REAL);
%   
%   
%   EXF_NML_02.addParm('apressurefile',apressurefile,PARM_STR);
%   EXF_NML_02.addParm('atempfile',atempfile,PARM_STR);
%   EXF_NML_02.addParm('aqhfile',aqhfile,PARM_STR);
%   EXF_NML_02.addParm('uwindfile',uwindfile,PARM_STR);
%   EXF_NML_02.addParm('vwindfile',vwindfile,PARM_STR);
%   EXF_NML_02.addParm('precipfile',precipfile,PARM_STR);
%   EXF_NML_02.addParm('swdownfile',swdownfile,PARM_STR);
%   EXF_NML_02.addParm('lwdownfile',lwdownfile,PARM_STR);
% %   EXF_NML_02.addParm('runofffile',runofffile,PARM_STR);
% 
% 
%   EXF_NML_02.addParm('apressurestartdate1',apressurestartdate1,PARM_INT);
%   EXF_NML_02.addParm('apressurestartdate2',apressurestartdate2,PARM_INT);
%   EXF_NML_02.addParm('apressureperiod',apressureperiod,PARM_REAL);
% 
%   EXF_NML_02.addParm('atempstartdate1',atempstartdate1,PARM_INT);
%   EXF_NML_02.addParm('atempstartdate2',atempstartdate2,PARM_INT);
%   EXF_NML_02.addParm('atempperiod',atempperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('aqhstartdate1',aqhstartdate1,PARM_INT);
%   EXF_NML_02.addParm('aqhstartdate2',aqhstartdate2,PARM_INT);
%   EXF_NML_02.addParm('aqhperiod',aqhperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('uwindstartdate1',uwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('uwindstartdate2',uwindstartdate2,PARM_INT);
%   EXF_NML_02.addParm('uwindperiod',uwindperiod,PARM_REAL);
%   
%   EXF_NML_02.addParm('vwindperiod',vwindperiod,PARM_REAL);
%   EXF_NML_02.addParm('vwindstartdate1',vwindstartdate1,PARM_INT);
%   EXF_NML_02.addParm('vwindstartdate2',vwindstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('precipperiod',precipperiod,PARM_REAL);
%   EXF_NML_02.addParm('precipstartdate1',precipstartdate1,PARM_INT);
%   EXF_NML_02.addParm('precipstartdate2',precipstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('swdownperiod',swdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('swdownstartdate1',swdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('swdownstartdate2',swdownstartdate2,PARM_INT);
%   
%   EXF_NML_02.addParm('lwdownperiod',lwdownperiod,PARM_REAL);
%   EXF_NML_02.addParm('lwdownstartdate1',lwdownstartdate1,PARM_INT);
%   EXF_NML_02.addParm('lwdownstartdate2',lwdownstartdate2,PARM_INT);
%   
% %   EXF_NML_02.addParm('runoffperiod',runoffperiod,PARM_REAL);
% 
%   EXF_NML_04.addParm('precip_lon0',precip_lon0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lon_inc',precip_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat0',precip_lat0,PARM_REAL);
%   EXF_NML_04.addParm('precip_lat_inc',precip_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('precip_nlon',precip_nlon,PARM_INT);
%   EXF_NML_04.addParm('precip_nlat',precip_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('atemp_lon0',atemp_lon0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lon_inc',atemp_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat0',atemp_lat0,PARM_REAL);
%   EXF_NML_04.addParm('atemp_lat_inc',atemp_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('atemp_nlon',atemp_nlon,PARM_INT);
%   EXF_NML_04.addParm('atemp_nlat',atemp_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('apressure_lon0',apressure_lon0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lon_inc',apressure_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat0',apressure_lat0,PARM_REAL);
%   EXF_NML_04.addParm('apressure_lat_inc',apressure_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('apressure_nlon',apressure_nlon,PARM_INT);
%   EXF_NML_04.addParm('apressure_nlat',apressure_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('aqh_lon0',aqh_lon0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lon_inc',aqh_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat0',aqh_lat0,PARM_REAL);
%   EXF_NML_04.addParm('aqh_lat_inc',aqh_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('aqh_nlon',aqh_nlon,PARM_INT);
%   EXF_NML_04.addParm('aqh_nlat',aqh_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('uwind_lon0',uwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lon_inc',uwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat0',uwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('uwind_lat_inc',uwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('uwind_nlon',uwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('uwind_nlat',uwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('vwind_lon0',vwind_lon0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lon_inc',vwind_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat0',vwind_lat0,PARM_REAL);
%   EXF_NML_04.addParm('vwind_lat_inc',vwind_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('vwind_nlon',vwind_nlon,PARM_INT);
%   EXF_NML_04.addParm('vwind_nlat',vwind_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('swdown_lon0',swdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lon_inc',swdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat0',swdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('swdown_lat_inc',swdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('swdown_nlon',swdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('swdown_nlat',swdown_nlat,PARM_INT);
% 
%   EXF_NML_04.addParm('lwdown_lon0',lwdown_lon0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lon_inc',lwdown_lon_inc,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat0',lwdown_lat0,PARM_REAL);
%   EXF_NML_04.addParm('lwdown_lat_inc',lwdown_lat_inc,PARM_REALS);
%   EXF_NML_04.addParm('lwdown_nlon',lwdown_nlon,PARM_INT);
%   EXF_NML_04.addParm('lwdown_nlat',lwdown_nlat,PARM_INT);
% 
% %   EXF_NML_04.addParm('runoff_lon0',runoff_lon0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lon_inc',runoff_lon_inc,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat0',runoff_lat0,PARM_REAL);
% %   EXF_NML_04.addParm('runoff_lat_inc',runoff_lat_inc,PARM_REALS);
% %   EXF_NML_04.addParm('runoff_nlon',runoff_nlon,PARM_INT);
% %   EXF_NML_04.addParm('runoff_nlat',runoff_nlat,PARM_INT);
% 
% 
%   %%z% Create the data.exf file
%   write_data_exf(inputpath,EXF_PARM,listterm,realfmt);





 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% SHELF ICE  %%%%%%%%%%
  %%%%     PARAMETERS       %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi0surf=zeros(Nx,Ny);
fid=fopen(fullfile(inputpath,'SHELFICEloadAnomalyFile.bin'), 'w','b'); 
fwrite(fid,phi0surf,prec);fclose(fid);

shelfthickness= 50; %15.75;

depth= - shelfthickness; %default -15 deep channel
icetopo=depth*ones(Nx,Ny);
%halfwidth= 500/2; % defined in the dx section %150; %30 %10 for default half channel width
%icetopo(round(Nx/2)-round(halfwidth/dx(1)):round(Nx/2)+round(halfwidth/dx(1)),: )= - sum(dz(1:3));



% finding index of side of channel 
for i = 1:Nx/2
     sum_dx = sum(dx(i:Nx/2));
if sum_dx < halfwidth
    break;
end
 side_index = i;
end

side_diff_index = round(Nx/2)-side_index;
side_index2 = round(Nx/2)+side_diff_index;

%{
% tanh function
for i = side_index:Nx/2
    icetopo(i,:) = depth + (shelfthickness/2)*tanh(0.007*mean(dx)*(i+1-side_index)-2)+(shelfthickness/2);
    icetopo_stored(i+1-side_index,:) = depth + (shelfthickness/2)*tanh(0.007*mean(dx)*(i+1-side_index)-2)+(shelfthickness/2); 
end
icetopo_stored = flip(icetopo_stored);

for i = 0:Nx/2-side_index
    icetopo(Nx/2 + i,: ) = icetopo_stored(i+1);
end
%}


% making ~1/3 of the channel width a flat top 

for i = round(Nx/2):Nx
    sum_top = sum(dx(round(Nx/2):i));
    if sum_top > halfwidth/3 
        break;
    end
    top_index2 = i; % finding the index to the right of the domain / in the positive direction
end

% find the index to the left of the domain / in the negative direction 
diff_top_index = top_index2 - round(Nx/2); 
top_index1 = round(Nx/2) - diff_top_index; 

% making the middle 1/3 of the channel start 3 gridpoints down 
for i = top_index1:top_index2 
    icetopo(i,:) = -sum(dz(1:3));
end 

% first 1/3 of the channel (left/negative side of x domain)
% arbitrarily defining first, middle, and last points and fitting cubic function
x_initial = 0; 
middle_index1 = round((top_index1-side_index)/4); 
middle_index2 = round((top_index1-side_index)/2);
middle_index3 = round(3*(top_index1-side_index)/4); 
x_middle1 = sum(dx(side_index:side_index+middle_index1));
x_middle2 = sum(dx(side_index:side_index+middle_index2));
x_middle3 = sum(dx(side_index:side_index+middle_index3));
x_last = sum(dx(side_index:top_index1));

y_initial = -shelfthickness; 
y_middle1 = -(9/10)*shelfthickness;
y_middle2 = -(1/3)*shelfthickness;
y_middle3 = -(1/10)*shelfthickness;
y_last = -sum(dz(1:3));

x_domain = [x_initial x_middle1 x_middle2 x_middle3 x_last];
y_domain = [y_initial y_middle1 y_middle2 y_middle3 y_last]; 

p_icetopo = polyfit(x_domain, y_domain, 3);
Fit_icetopo = polyval(p_icetopo,x_domain);


% interpolating to full domain 
%channel_domain1 = zeros(1,21); 
for i = side_index:top_index1
channel_domain1(i-side_index+1) = sum(dx(side_index:i)); 
end

channel_interp1 = interp1(x_domain,Fit_icetopo,channel_domain1,'spline'); 

for i = 1:length(channel_interp1)
    if channel_interp1(i) > -sum(dz(1:3))
        channel_interp1(i) = -sum(dz(1:3));
        break;
    end
    channel_interp1_index = i; 
end

for i = channel_interp1_index+1:length(channel_interp1)
channel_interp1(i) = -sum(dz(1:3));
end

for i = 1:length(channel_interp1)
    if channel_interp1(i) < -shelfthickness
        channel_interp1(i) = -shelfthickness;
    end
end



for i = side_index:top_index1
    icetopo(i,:) = channel_interp1(i-side_index+1);
end


% last 1/3 of the channel (right/positive side of x domain) 
channel_interp2 = flip(channel_interp1); 

for i = top_index2:side_index2
    icetopo(i,:) = channel_interp2(i-top_index2+1);
end




%{
% creating a tanh function as a function of x
%  defining x domain
initial_xpoint = sum(dx(1:top_index1));
final_xpoint = sum(dx(1:top_index2));
x_length = final_xpoint - initial_xpoint; % finding the distance between the left side of the channel and the start of the middle of the channel
x_domain = 0:1/x_length:x_length; 
% defining z domain
initial_zpoint = -shelfthickness;
final_zpoint = - sum(dz(1:3)); 
z_length = final_zpoint - initial_zpoint; 
z_domain = 0:1/z_length:z_length; 

% creating a z = tanh(x) function that has a z range equal to z_length 
% creating arbitrary x domain 
x = -1:1/top_index1:1; 
z_domain = (tanh(x) + 1) * (z_length / 2);

% interpolating to dx 
for i = 1:top_index1
channel_domain1(i) = sum(dx(1:i)); 
end

channel_interp1 = interp1(x,z_domain,channel_domain1,'spline'); 
%}










%{
%rounder channel 
% finding index of side of channel 
for i = 1:Nx/2
     sum_dx = sum(dx(i:Nx/2));
if sum_dx < halfwidth
    break;
end
 side_index1 = i;
end

index_diff = round(Nx/2)- side_index1; 
side_index2 = round(Nx/2) + index_diff; 

% finding index for bottom of the channel 
for i = 1:Nr 
    sum_dz = sum(dz(1:i));
    if sum_dz > shelfthickness 
        break;
    end
    bot_index = i;
end

% creating gradient 

icetopo_stored = [];
for i = side_index1:round(Nx/2) 
    icetopo(i,:) = depth + (shelfthickness/index_diff)*(i-side_index1+1);

    if icetopo(i,:) > -sum(dz(1:3))
        icetopo(i,:) = -sum(dz(1:3));
    end
    icetopo_stored(i-side_index1+1,:) = icetopo(i,:);
end

icetopo_stored = flip(icetopo_stored); 

% mirroring icetopo
for i = round(Nx/2):side_index2 
    icetopo(i,:) = icetopo_stored(i-round(Nx/2)+1,:);
end

%{
top_length = halfwidth/3; 

for i = round(Nx/2):1
sum_top = sum(dz((Nx/2):i)
end 

incr_chan = diff_index/5; 
for i = round(Nx/2)-5:round(Nx/2)
icetopo(i,:) = -sum(dz(1:3))
end





%}






%icetopo(round(Nx/2),:) = -dz(1);
%}

    
%{

dx_og = 12.5; 

for i = 0:round(halfwidth/dx_og) 
    icetopo(round(Nx/2)-round(halfwidth/dx_og)+i,:) = depth + (shelfthickness/2)*tanh(0.02*dx_og*i-2)+(shelfthickness/2); 
    icetopo_stored(i+1,:) = depth + (shelfthickness/2)*tanh(0.02*dx_og*i-2)+(shelfthickness/2); 
end

icetopo_stored = flip(icetopo_stored);

for i = 0:round(halfwidth/dx_og) 
    icetopo(round(Nx/2) + i,: ) = icetopo_stored(i+1); %depth - (shelfthickness/2)*tanh(0.02*dx(1)*i-2)+(shelfthickness/2); 
end 

icetopo(round(Nx/2),:) = -dz(1);
%}

%{
%rounding bottom of channel
icetopo(round(Nx/2)-round(halfwidth/dx(1)):round(Nx/2)+round(halfwidth/dx(1)),: )= - 0.9;
icetopo(round(Nx/2)-round(halfwidth/dx(1))-1,:) = depth+2*0.9; 
icetopo(round(Nx/2)+round(halfwidth/dx(1))+1,:) = depth+2*0.9;
icetopo(round(Nx/2)-round(halfwidth/dx(1))-2,:) = depth+2*0.6;
icetopo(round(Nx/2)+round(halfwidth/dx(1))+2,:) = depth+2*0.6;
icetopo(round(Nx/2)-round(halfwidth/dx(1))-3,:) = depth+2*0.3;
icetopo(round(Nx/2)+round(halfwidth/dx(1))+3,:) = depth+2*0.3;
icetopo(round(Nx/2)-round(halfwidth/dx(1))-4,:) = depth+0.3;
icetopo(round(Nx/2)+round(halfwidth/dx(1))+4,:) = depth+0.3;


%rounding top of channel 
icetopo(round(Nx/2)-round(halfwidth/dx(1)),:) = -0.9-2*0.9; 
icetopo(round(Nx/2)+round(halfwidth/dx(1)),:) = -0.9-2*0.9;
icetopo(round(Nx/2)-round(halfwidth/dx(1))+1,:) = -0.9-2*0.6; 
icetopo(round(Nx/2)+round(halfwidth/dx(1))-1,:) = -0.9-2*0.6;
icetopo(round(Nx/2)-round(halfwidth/dx(1))+2,:) = -0.9-2*0.3; 
icetopo(round(Nx/2)+round(halfwidth/dx(1))-2,:) = -0.9-2*0.3;
icetopo(round(Nx/2)-round(halfwidth/dx(1))+3,:) = -0.9-0.3; 
icetopo(round(Nx/2)+round(halfwidth/dx(1))-3,:) = -0.9-0.3;
%}





fid=fopen(fullfile(inputpath,'SHELFICEtopoFile.bin'), 'w','b'); 
fwrite(fid,icetopo,prec);fclose(fid);

% to store parameter names and values
shelfice_parm01 = parmlist;
SHELFICE_PARM = {shelfice_parm01};

SHELFICEloadAnomalyFile = 'SHELFICEloadAnomalyFile.bin';
SHELFICEtopoFile = 'SHELFICEtopoFile.bin';
SHELFICEuseGammaFrict = true;
SHELFICEboundaryLayer = true;
SHELFICEconserve = true;
%   SHELFICEheatTransCoeff = .0005;
%   SHELFICEheatTransCoeff = .0001;
SHELFICEheatTransCoeff = 0;
SHELFICEwriteState = true;

SHELFICEDragQuadratic = 2e-3;
SHELFICEselectDragQuadr = 0;


shelfice_parm01.addParm('SHELFICEloadAnomalyFile',SHELFICEloadAnomalyFile,PARM_STR);
shelfice_parm01.addParm('SHELFICEtopoFile',SHELFICEtopoFile,PARM_STR);
shelfice_parm01.addParm('SHELFICEuseGammaFrict',SHELFICEuseGammaFrict,PARM_BOOL);
shelfice_parm01.addParm('SHELFICEboundaryLayer',SHELFICEboundaryLayer,PARM_BOOL);
shelfice_parm01.addParm('SHELFICEconserve',SHELFICEconserve,PARM_BOOL);
shelfice_parm01.addParm('SHELFICEheatTransCoeff',SHELFICEheatTransCoeff,PARM_REAL);
shelfice_parm01.addParm('SHELFICEwritestate',SHELFICEwriteState,PARM_BOOL);

shelfice_parm01.addParm('SHELFICEDragQuadratic',SHELFICEDragQuadratic,PARM_REAL);
shelfice_parm01.addParm('SHELFICEselectDragQuadr',SHELFICEselectDragQuadr,PARM_INT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WRITE THE 'data.shelfice' FILE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_data_shelfice(inputpath,SHELFICE_PARM,listterm,realfmt);

  

  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%



%{
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%


  %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  obcs_parm03 = parmlist;
  obcs_parm04 = parmlist;
  %obcs_parm05 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02,obcs_parm03,obcs_parm04};


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

  %%%%%% tides ^%%%%%%%% 
  useOBCStides = use_tides;

  tidalPeriod=[44714, 43200, 45570, 43082, 86164, 92950, 86637, 96726,1180300,2380706];
%   tidalPeriod=[44714, 43200, 45570, 43082, 86164, 92950, 86637];

  OBNamFile= 'OBNamFile.bin';
  OBNphFile= 'OBNphFile.bin';
  OBEamFile= 'OBEamFile.bin';
  OBEphFile= 'OBEphFile.bin';



  obcs_parm01.addParm('useOBCStides',useOBCStides,PARM_BOOL);

  obcs_parm01.addParm('tidalPeriod',tidalPeriod,PARM_INTS);

  obcs_parm01.addParm('OBNamFile',OBNamFile,PARM_STR);
  obcs_parm01.addParm('OBNphFile',OBNphFile,PARM_STR);
  obcs_parm01.addParm('OBEamFile',OBEamFile,PARM_STR);
  obcs_parm01.addParm('OBEphFile',OBEphFile,PARM_STR);

%}


  %%% Enables an Orlanski radiation condition at the northern boundary


    useOrlanskiNorth = false;

    OB_Jsouth = 1*ones(1,Nx);
%     OB_Iwest = 1*ones(1,Ny);
    OB_Jnorth = -1*ones(1,Nx);
  
  %OB_Jnorth(1:idx_obcs_n) = 0;
  %OB_Ieast(1:idx_obcs_e) = 0;


 % obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS);
  obcs_parm01.addParm('OB_Jsouth',OB_Jsouth,PARM_INTS);
%   obcs_parm01.addParm('OB_Iwest',OB_Iwest,PARM_INTS);    



    useOBCSsponge = true;
   

%     useOBCSsponge = false;

    useOBCSprescribe = true;


    OBNtFile = 'NBCt.bin';
    OBNsFile = 'NBCs.bin';

%     OBWtFile = 'OBWtFile.bin';


    OBStFile = 'SBCt.bin';
    OBSsFile = 'SBCs.bin';

%     OBWsFile = 'OBWsFile.bin';

%{
  fid = fopen(fullfile(inputpath,OBNtFile),'r','b');
  OBNa = fread(fid,[Nx 12],'real*8');
  fclose(fid);


    OBNhFile = 'OBNhFile.bin';
    OBEhFile = 'OBEhFile.bin';
%     OBWhFile = 'OBWhFile.bin';

    OBNsnFile = 'OBNsnFile.bin';
    OBEsnFile = 'OBEsnFile.bin';
%     OBWsnFile = 'OBWsnFile.bin';

    OBNuiceFile = 'OBNuiceFile.bin';
    OBEuiceFile = 'OBEuiceFile.bin';
% %     OBWuiceFile = 'OBWuiceFile.bin';
% 
    OBNviceFile = 'OBNviceFile.bin';
    OBEviceFile = 'OBEviceFile.bin';
%     OBWviceFile = 'OBWviceFile.bin';
%}

  obcs_parm01.addParm('useOBCSsponge',useOBCSsponge,PARM_BOOL);
  obcs_parm01.addParm('useOBCSprescribe',useOBCSprescribe,PARM_BOOL);

  obcs_parm01.addParm('OBNtFile',OBNtFile,PARM_STR);
  obcs_parm01.addParm('OBNsFile',OBNsFile,PARM_STR);
%   obcs_parm01.addParm('OBWtFile',OBWtFile,PARM_STR);  

  obcs_parm01.addParm('OBStFile',OBStFile,PARM_STR);
  obcs_parm01.addParm('OBSsFile',OBSsFile,PARM_STR);
%   obcs_parm01.addParm('OBWsFile',OBWsFile,PARM_STR);



shelfz=Nr;
ny=Ny;nx=Nx;
clear EBCu EBCs EBCt  WBCu WBCs WBCt NBCt NBCu NBCs 
EBCu = zeros(ny,shelfz);
EBCs = zeros(ny,shelfz);
EBCt = zeros(ny,shelfz);
EBCv = zeros(ny,shelfz);

WBCu = zeros(ny,shelfz);
WBCs = zeros(ny,shelfz);
WBCt = zeros(ny,shelfz);
WBCv = zeros(ny,shelfz);

NBCu = zeros(nx,shelfz);
NBCs = zeros(nx,shelfz);
NBCt = zeros(nx,shelfz);
NBCv = zeros(nx,shelfz);

SBCu = zeros(nx,shelfz);
SBCs = zeros(nx,shelfz);
SBCt = zeros(nx,shelfz);
SBCv = zeros(nx,shelfz);


%repmat(linspace(saltini(round(h2loc/deltaX)-2,1,i),saltini(round(h2loc/deltaX)+50,1,i), 103),ny,1)'


        NBCt = repmat(tRefout,[Nx,1]);
        NBCs = repmat(sRefout,[Nx,1]);
       
        SBCt=repmat(tRef,[Nx,1]);
        SBCs = repmat(sRef,[Nx,1]);



% Apply barotropic velocity to balance input of runoff


fid=fopen(fullfile(inputpath,'NBCt.bin'), 'w','b');  fwrite(fid,NBCt,prec);fclose(fid);
fid=fopen(fullfile(inputpath,'NBCs.bin'), 'w','b');  fwrite(fid,NBCs,prec);fclose(fid);
fid=fopen(fullfile(inputpath,'SBCt.bin'), 'w','b');  fwrite(fid,SBCt,prec);fclose(fid);
fid=fopen(fullfile(inputpath,'SBCs.bin'), 'w','b');  fwrite(fid,SBCs,prec);fclose(fid);

  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow
  %useOBCSbalance = true;
  %OBCS_balanceFacN = 1;
  %OBCS_balanceFacE = 0;
%   OBCS_balanceFacS = 0;
%   OBCS_balanceFacW = -1;
  %obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);
  %obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL);
%  obcs_parm01.addParm('OBCS_balanceFacE',OBCS_balanceFacE,PARM_REAL);
%   obcs_parm01.addParm('OBCS_balanceFacS',OBCS_balanceFacS,PARM_REAL);  
%   obcs_parm01.addParm('OBCS_balanceFacW',OBCS_balanceFacW,PARM_REAL);  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Velocity averaging time scale - must be larger than deltaT.
  %%% The Orlanski radiation condition computes the characteristic velocity
  %%% at the boundary by averaging the spatial derivative normal to the 
  %%% boundary divided by the time step over this period.
  %%% At the moment we're using the magic engineering factor of 3.
%   cvelTimeScale = 3*deltaT;




  %%% Max dimensionless CFL for Adams-Basthforth 2nd-order method
%   CMAX = 0.45; 
%   
%   obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
%   obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Sponge Layer Parms (OBCS_PARM03) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Values taken from SOSE.
%%% Urelaxobcsinner = relaxation time scale at the innermost sponge layer point of a meridional OB
%%% Vrelaxobcsinner = relaxation time scale at the innermost sponge layer point of a zonal OB
%%% Urelaxobcsbound = relaxation time scale at the outermost sponge layer point of a meridional OB
%%% Vrelaxobcsbound = relaxation time scale at the outermost sponge layer point of a zonal OB



    Urelaxobcsinner = 2160;  %%% 10 days
    Urelaxobcsbound = 216;  %%% half a day
    Vrelaxobcsinner = 2160;
    Vrelaxobcsbound = 216;

%%%%%% sponge thickness - finer in high-res simulation due to placement of
%%%%%% eastern boundary, increased number of gridpoints
  %if (res_fac == 24)
    spongethickness =round(Ny/10);
  %else
  %  spongethickness = round(10*res_fac/3);
%     spongethickness = 5;
  %end


  obcs_parm03.addParm('Urelaxobcsinner',Urelaxobcsinner,PARM_REAL);
  obcs_parm03.addParm('Urelaxobcsbound',Urelaxobcsbound,PARM_REAL);
  obcs_parm03.addParm('Vrelaxobcsinner',Vrelaxobcsinner,PARM_REAL);
  obcs_parm03.addParm('Vrelaxobcsbound',Vrelaxobcsbound,PARM_REAL);

  obcs_parm03.addParm('spongethickness',spongethickness,PARM_INT);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Sea ice Sponge Parms (OBCS_PARM05) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
  seaiceSpongeThickness = spongethickness;
  Arelaxobcsinner = Urelaxobcsinner;
  Arelaxobcsbound = Urelaxobcsbound;
  Hrelaxobcsinner = Urelaxobcsinner;
  Hrelaxobcsbound = Urelaxobcsbound;
  SLrelaxobcsinner = Urelaxobcsinner;
  SLrelaxobcsbound = Urelaxobcsbound;
  SNrelaxobcsinner = Urelaxobcsinner;
  SNrelaxobcsbound = Urelaxobcsbound;

  obcs_parm05.addParm('Arelaxobcsinner',Arelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('Arelaxobcsbound',Arelaxobcsbound,PARM_REAL);
  obcs_parm05.addParm('Hrelaxobcsinner',Hrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('Hrelaxobcsbound',Hrelaxobcsbound,PARM_REAL);
  obcs_parm05.addParm('SLrelaxobcsinner',SLrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('SLrelaxobcsbound',SLrelaxobcsbound,PARM_REAL);
  obcs_parm05.addParm('SNrelaxobcsinner',SNrelaxobcsinner,PARM_REAL);
  obcs_parm05.addParm('SNrelaxobcsbound',SNrelaxobcsbound,PARM_REAL);
  obcs_parm05.addParm('seaiceSpongeThickness',seaiceSpongeThickness,PARM_INT);
%}


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);

%}















  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  diag_fields_avg = ...
  {...
     'UV_VEL_Z','WU_VEL','WV_VEL', ... %%% Momentum fluxes
     'UVELSLT','VVELSLT','WVELSLT', ... %%% Salt fluxes
      'UVEL','VVEL','WVEL', ... %%% Velocities
      'THETA', ... %%% Temperature
      'SALT', ... %%% Salinity
      'PHIHYD','PHI_NH', 'Vm_dPhiY', 'Um_dPhiX', ... %%% Pressure
      'SHIgammT','SHIgammS','SHIuStar','SHI_mass', ... %%%%% Ice shelf melt
      'SHIfwFlx','SHIhtFlx','SHIForcT','SHIForcS', 'SHI_TauY', 'SHI_TauX', ... 
      'UVELTH','VVELTH','WVELTH', ... %%% Temperature fluxes
      'Um_Advec', 'Vm_Advec', 'Wm_Advec', ... %%%%% Advective fluxes
      'VSidDrag', ... %%% Momentum tendency from side drag 
      'Vm_Diss', 'Vm_ImplD', ... %%% Momentum tendency from dissipation 
      'Vm_Cori', ... %%% Momentum tendency from coriolis
      'Vm_Ext', ... %%% Momentum tendency from external forcing 
      'Vm_AdvZ3', 'Vm_AdvRe', ... %%% Momentum tendency from vorticity and vertical advection
      'TOTVTEND', ... %%% V tendency


  };
  % 'SHIVDrag', ...
  %     'PHIHYD', ... %%% Pressure
  %'SHIVDrag'
%     'LaUH1TH','LaVH1TH','LaHw1TH','LaHs1TH' ... %%% LAYERS fluxes
%     'ADVr_TH','ADVx_TH','ADVy_TH', ... %%%%% Advective fluxes of H.B.
%     'Um_Advec', 'Vm_Advec', 'Wm_Advec', %%%%% Advective fluxes 
%     'KPPg_TH', ... %%%% KPP non-local flux of P.T.
%     'DFrE_TH','DFrI_TH','DFxE_TH','DFyE_TH', ... %%%% Diffusive Fluxes
%     'oceQsw','TOTTTEND','WTHMASS','TFLUX'...
%      'UVELSQ','VVELSQ','WVELSQ',... %%% For Kinetic Energy
   % 'SHIgammT','SHIgammS','SHIuStar','SHI_mass', ... %%%%% Ice shelf melt
%      'SIarea','SIheff','SIhsnow','SIhsalt','SIuice','SIvice'...%%% Sea ice state

  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 1*t1day;
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',true,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
  diag_fields_inst = ...
  {...  
    'ETAN', ...  %%% SSH
    'UVEL','VVEL','WVEL'...%%%velocities
    'THETA', ... %%% Temperature
    'SALT', ... %%% Salinity
     %'KPPhbl' ... %%% BL depth
    }; 
  
  numdiags_inst = length(diag_fields_inst);  
  if (use_3D)
    diag_freq_inst = 0.05*t1day; % 0.05
  else
    diag_freq_inst = 0.05*t1day;
  end
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,Nr);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);
  packages.addParm('useSHELFICE',true,PARM_BOOL);
  %packages.addParm('useKPP',~nonHydrostatic,PARM_BOOL);
  packages.addParm('useRBCS',~use_seaIce,PARM_BOOL);      
  %packages.addParm('useOBCS',true,PARM_BOOL);     
  packages.addParm('useSEAICE',false,PARM_BOOL);  
  packages.addParm('useEXF',false,PARM_BOOL);  
  
  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
  ALL_PARMS = [PARM PACKAGE_PARM DIAG_MATLAB_PARM];
  if (~nonHydrostatic)
    ALL_PARMS = [ALL_PARMS];
  end
  if (use_seaIce)
    ALL_PARMS = [ALL_PARMS SEAICE_PARM EXF_PARM];
  %else
    %ALL_PARMS = [ALL_PARMS OBCS_PARM];
  end    
  write_matlab_params(inputpath,ALL_PARMS,realfmt);
  
end


