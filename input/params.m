%%% 
%%% params.m
%%% 
%%% Automatically-generated Matlab file containing parameters
%%% used in this MITgcm experiment.
%%% 

vectorInvariantMomentum=1;
viscAr=0.00000000e+00;
viscA4=0.00000000e+00;
viscAh=0.00000000e+00;
viscA4Grid=0.00000000e+00;
viscAhGrid=0.00000000e+00;
viscA4GridMax=5.00000000e-01;
viscAhGridMax=1.00000000e+00;
useAreaViscLength=0;
useFullLeith=1;
useSmag3D=1;
viscC4leith=0.00000000e+00;
viscC4leithD=0.00000000e+00;
viscC2leith=0.00000000e+00;
viscC2leithD=0.00000000e+00;
viscC2smag=0.00000000e+00;
viscC4smag=4.00000000e+00;
smag3D_coeff=2.00000000e-02;
diffKrT=1.30000000e-07;
diffKhT=0.00000000e+00;
diffKrS=7.20000000e-10;
diffKhS=0.00000000e+00;
tempAdvScheme=33;
saltAdvScheme=33;
multiDimAdvection=1;
tempStepping=1;
saltStepping=1;
staggerTimeStep=1;
eosType='JMD95Z';
no_slip_sides=0;
no_slip_bottom=0;
bottomDragLinear=0.00000000e+00;
bottomDragQuadratic=0.00000000e+00;
f0=-1.40000000e-04;
gravity=9.81000000e+00;
rhonil=1.00000000e+03;
rhoConst=1.00000000e+03;
ivdc_kappa=0.00000000e+00;
implicitDiffusion=1;
implicitViscosity=1;
nonHydrostatic=0;
exactConserv=1;
useCDscheme=0;
readBinaryPrec=64;
useSingleCpuIO=1;
debugLevel=-1;
vectorInvariantMomentum=1;
useJamartWetPoints=1;
useJamartMomAdv=1;

useSRCGSolver=1;
cg2dMaxIters=1000;
cg2dTargetResidual=1.00000000e-12;
cg3dMaxIters=300;
cg3dTargetResidual=1.00000000e-07;

momDissip_In_AB=0;
tracForcingOutAB=1;
nIter0=0;
abEps=1.00000000e-01;
chkptFreq=4.32000000e+03;
pChkptFreq=4.32000000e+03;
taveFreq=0.00000000e+00;
dumpFreq=0.00000000e+00;
monitorFreq=8.64000000e+04;
cAdjFreq=0.00000000e+00;
dumpInitAndLast=1;
pickupStrictlyMatch=0;
endTime=864002.6667;
deltaT=7.66666667e+00;

usingCartesianGrid=1;
usingSphericalPolarGrid=0;
delX=[ 2.39764266e+01 2.39739514e+01 2.39712168e+01 2.39681956e+01 2.39648581e+01 2.39611711e+01 2.39570984e+01 2.39525998e+01 2.39476310e+01 2.39421432e+01 2.39360826e+01 2.39293900e+01 2.39220001e+01 2.39138410e+01 2.39048335e+01 2.38948905e+01 2.38839162e+01 2.38718054e+01 2.38584423e+01 2.38437000e+01 2.38274391e+01 2.38095068e+01 2.37897357e+01 2.37679427e+01 2.37439276e+01 2.37174718e+01 2.36883369e+01 2.36562632e+01 2.36209685e+01 2.35821464e+01 2.35394653e+01 2.34925665e+01 2.34410636e+01 2.33845409e+01 2.33225532e+01 2.32546246e+01 2.31802491e+01 2.30988907e+01 2.30099842e+01 2.29129375e+01 2.28071337e+01 2.26919353e+01 2.25666887e+01 2.24307308e+01 2.22833963e+01 2.21240271e+01 2.19519836e+01 2.17666571e+01 2.15674841e+01 2.13539627e+01 2.11256694e+01 2.08822781e+01 2.06235786e+01 2.03494960e+01 2.00601085e+01 1.97556647e+01 1.94365975e+01 1.91035350e+01 1.87573072e+01 1.83989472e+01 1.80296868e+01 1.76509456e+01 1.72643134e+01 1.68715274e+01 1.64744426e+01 1.60749978e+01 1.56751787e+01 1.52769783e+01 1.48823576e+01 1.44932072e+01 1.41113120e+01 1.37383203e+01 1.33757180e+01 1.30248087e+01 1.26867010e+01 1.23623012e+01 1.20523126e+01 1.17572403e+01 1.14774003e+01 1.12129331e+01 1.09638190e+01 1.07298964e+01 1.05108803e+01 1.03063816e+01 1.01159258e+01 9.93897040e+00 9.77492198e+00 9.62315083e+00 9.48300446e+00 9.35381916e+00 9.23492996e+00 9.12567889e+00 9.02542174e+00 8.93353345e+00 8.84941222e+00 8.77248256e+00 8.70219734e+00 8.63803916e+00 8.57952089e+00 8.52618584e+00 8.52618584e+00 8.57952089e+00 8.63803916e+00 8.70219734e+00 8.77248256e+00 8.84941222e+00 8.93353345e+00 9.02542174e+00 9.12567889e+00 9.23492996e+00 9.35381916e+00 9.48300446e+00 9.62315083e+00 9.77492198e+00 9.93897040e+00 1.01159258e+01 1.03063816e+01 1.05108803e+01 1.07298964e+01 1.09638190e+01 1.12129331e+01 1.14774003e+01 1.17572403e+01 1.20523126e+01 1.23623012e+01 1.26867010e+01 1.30248087e+01 1.33757180e+01 1.37383203e+01 1.41113120e+01 1.44932072e+01 1.48823576e+01 1.52769783e+01 1.56751787e+01 1.60749978e+01 1.64744426e+01 1.68715274e+01 1.72643134e+01 1.76509456e+01 1.80296868e+01 1.83989472e+01 1.87573072e+01 1.91035350e+01 1.94365975e+01 1.97556647e+01 2.00601085e+01 2.03494960e+01 2.06235786e+01 2.08822781e+01 2.11256694e+01 2.13539627e+01 2.15674841e+01 2.17666571e+01 2.19519836e+01 2.21240271e+01 2.22833963e+01 2.24307308e+01 2.25666887e+01 2.26919353e+01 2.28071337e+01 2.29129375e+01 2.30099842e+01 2.30988907e+01 2.31802491e+01 2.32546246e+01 2.33225532e+01 2.33845409e+01 2.34410636e+01 2.34925665e+01 2.35394653e+01 2.35821464e+01 2.36209685e+01 2.36562632e+01 2.36883369e+01 2.37174718e+01 2.37439276e+01 2.37679427e+01 2.37897357e+01 2.38095068e+01 2.38274391e+01 2.38437000e+01 2.38584423e+01 2.38718054e+01 2.38839162e+01 2.38948905e+01 2.39048335e+01 2.39138410e+01 2.39220001e+01 2.39293900e+01 2.39360826e+01 2.39421432e+01 2.39476310e+01 2.39525998e+01 2.39570984e+01 2.39611711e+01 2.39648581e+01 2.39681956e+01 2.39712168e+01 2.39739514e+01 2.39764266e+01 ];
delY=[ 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 6.00000000e+01 ];
delR=[ 6.39296926e-01 6.39296926e-01 6.50921494e-01 6.62858521e-01 6.75116354e-01 6.87703563e-01 7.00628942e-01 7.13901521e-01 7.27530569e-01 7.41525598e-01 7.55896376e-01 7.70652925e-01 7.85805537e-01 8.01364774e-01 8.17341478e-01 8.33746777e-01 8.50592096e-01 8.67889158e-01 8.85650001e-01 9.03886976e-01 9.22612766e-01 9.41840384e-01 9.61583190e-01 9.81854896e-01 1.00266958e+00 1.02404167e+00 1.04598602e+00 1.06851783e+00 1.09165272e+00 1.11540672e+00 1.13979630e+00 1.16483832e+00 1.19055014e+00 1.21694953e+00 1.24405476e+00 1.27188457e+00 1.30045818e+00 1.32979533e+00 1.35991628e+00 1.39084182e+00 1.42259330e+00 1.45519259e+00 1.48866219e+00 1.52302514e+00 1.55830513e+00 1.59452645e+00 1.63171401e+00 1.66989341e+00 1.70909090e+00 1.74933343e+00 1.79064863e+00 1.83306490e+00 1.87661134e+00 1.92131784e+00 1.96721506e+00 2.01433447e+00 2.06270836e+00 2.11236987e+00 2.16335299e+00 2.21569263e+00 2.26942459e+00 2.32458560e+00 2.38121337e+00 2.43934659e+00 2.49902494e+00 2.56028918e+00 2.62318108e+00 2.68774353e+00 2.75402055e+00 2.82205728e+00 ];

bathyFile='bathyFile.bin';
hydrogThetaFile='hydrogThetaFile.bin';
hydrogSaltFile='hydrogSaltFile.bin';
vVelInitFile='vVelInitFile.bin';

useDiagnostics=1;
useSHELFICE=1;
useRBCS=1;
useSEAICE=0;
useEXF=0;

diag_fields{1,1}='UV_VEL_Z';
diag_fileNames{1}='UV_VEL_Z';
diag_frequency(1)=8.64000000e+04;
diag_timePhase(1)=0.00000000e+00;
diag_fields{1,2}='WU_VEL';
diag_fileNames{2}='WU_VEL';
diag_frequency(2)=8.64000000e+04;
diag_timePhase(2)=0.00000000e+00;
diag_fields{1,3}='WV_VEL';
diag_fileNames{3}='WV_VEL';
diag_frequency(3)=8.64000000e+04;
diag_timePhase(3)=0.00000000e+00;
diag_fields{1,4}='UVELSLT';
diag_fileNames{4}='UVELSLT';
diag_frequency(4)=8.64000000e+04;
diag_timePhase(4)=0.00000000e+00;
diag_fields{1,5}='VVELSLT';
diag_fileNames{5}='VVELSLT';
diag_frequency(5)=8.64000000e+04;
diag_timePhase(5)=0.00000000e+00;
diag_fields{1,6}='WVELSLT';
diag_fileNames{6}='WVELSLT';
diag_frequency(6)=8.64000000e+04;
diag_timePhase(6)=0.00000000e+00;
diag_fields{1,7}='UVEL';
diag_fileNames{7}='UVEL';
diag_frequency(7)=8.64000000e+04;
diag_timePhase(7)=0.00000000e+00;
diag_fields{1,8}='VVEL';
diag_fileNames{8}='VVEL';
diag_frequency(8)=8.64000000e+04;
diag_timePhase(8)=0.00000000e+00;
diag_fields{1,9}='WVEL';
diag_fileNames{9}='WVEL';
diag_frequency(9)=8.64000000e+04;
diag_timePhase(9)=0.00000000e+00;
diag_fields{1,10}='THETA';
diag_fileNames{10}='THETA';
diag_frequency(10)=8.64000000e+04;
diag_timePhase(10)=0.00000000e+00;
diag_fields{1,11}='SALT';
diag_fileNames{11}='SALT';
diag_frequency(11)=8.64000000e+04;
diag_timePhase(11)=0.00000000e+00;
diag_fields{1,12}='PHIHYD';
diag_fileNames{12}='PHIHYD';
diag_frequency(12)=8.64000000e+04;
diag_timePhase(12)=0.00000000e+00;
diag_fields{1,13}='PHI_NH';
diag_fileNames{13}='PHI_NH';
diag_frequency(13)=8.64000000e+04;
diag_timePhase(13)=0.00000000e+00;
diag_fields{1,14}='Vm_dPhiY';
diag_fileNames{14}='Vm_dPhiY';
diag_frequency(14)=8.64000000e+04;
diag_timePhase(14)=0.00000000e+00;
diag_fields{1,15}='Um_dPhiX';
diag_fileNames{15}='Um_dPhiX';
diag_frequency(15)=8.64000000e+04;
diag_timePhase(15)=0.00000000e+00;
diag_fields{1,16}='SHIgammT';
diag_fileNames{16}='SHIgammT';
diag_frequency(16)=8.64000000e+04;
diag_timePhase(16)=0.00000000e+00;
diag_fields{1,17}='SHIgammS';
diag_fileNames{17}='SHIgammS';
diag_frequency(17)=8.64000000e+04;
diag_timePhase(17)=0.00000000e+00;
diag_fields{1,18}='SHIuStar';
diag_fileNames{18}='SHIuStar';
diag_frequency(18)=8.64000000e+04;
diag_timePhase(18)=0.00000000e+00;
diag_fields{1,19}='SHI_mass';
diag_fileNames{19}='SHI_mass';
diag_frequency(19)=8.64000000e+04;
diag_timePhase(19)=0.00000000e+00;
diag_fields{1,20}='SHIfwFlx';
diag_fileNames{20}='SHIfwFlx';
diag_frequency(20)=8.64000000e+04;
diag_timePhase(20)=0.00000000e+00;
diag_fields{1,21}='SHIhtFlx';
diag_fileNames{21}='SHIhtFlx';
diag_frequency(21)=8.64000000e+04;
diag_timePhase(21)=0.00000000e+00;
diag_fields{1,22}='SHIForcT';
diag_fileNames{22}='SHIForcT';
diag_frequency(22)=8.64000000e+04;
diag_timePhase(22)=0.00000000e+00;
diag_fields{1,23}='SHIForcS';
diag_fileNames{23}='SHIForcS';
diag_frequency(23)=8.64000000e+04;
diag_timePhase(23)=0.00000000e+00;
diag_fields{1,24}='SHI_TauY';
diag_fileNames{24}='SHI_TauY';
diag_frequency(24)=8.64000000e+04;
diag_timePhase(24)=0.00000000e+00;
diag_fields{1,25}='SHI_TauX';
diag_fileNames{25}='SHI_TauX';
diag_frequency(25)=8.64000000e+04;
diag_timePhase(25)=0.00000000e+00;
diag_fields{1,26}='UVELTH';
diag_fileNames{26}='UVELTH';
diag_frequency(26)=8.64000000e+04;
diag_timePhase(26)=0.00000000e+00;
diag_fields{1,27}='VVELTH';
diag_fileNames{27}='VVELTH';
diag_frequency(27)=8.64000000e+04;
diag_timePhase(27)=0.00000000e+00;
diag_fields{1,28}='WVELTH';
diag_fileNames{28}='WVELTH';
diag_frequency(28)=8.64000000e+04;
diag_timePhase(28)=0.00000000e+00;
diag_fields{1,29}='Um_Advec';
diag_fileNames{29}='Um_Advec';
diag_frequency(29)=8.64000000e+04;
diag_timePhase(29)=0.00000000e+00;
diag_fields{1,30}='Vm_Advec';
diag_fileNames{30}='Vm_Advec';
diag_frequency(30)=8.64000000e+04;
diag_timePhase(30)=0.00000000e+00;
diag_fields{1,31}='Wm_Advec';
diag_fileNames{31}='Wm_Advec';
diag_frequency(31)=8.64000000e+04;
diag_timePhase(31)=0.00000000e+00;
diag_fields{1,32}='VSidDrag';
diag_fileNames{32}='VSidDrag';
diag_frequency(32)=8.64000000e+04;
diag_timePhase(32)=0.00000000e+00;
diag_fields{1,33}='Vm_Diss';
diag_fileNames{33}='Vm_Diss';
diag_frequency(33)=8.64000000e+04;
diag_timePhase(33)=0.00000000e+00;
diag_fields{1,34}='Vm_ImplD';
diag_fileNames{34}='Vm_ImplD';
diag_frequency(34)=8.64000000e+04;
diag_timePhase(34)=0.00000000e+00;
diag_fields{1,35}='Vm_Cori';
diag_fileNames{35}='Vm_Cori';
diag_frequency(35)=8.64000000e+04;
diag_timePhase(35)=0.00000000e+00;
diag_fields{1,36}='Vm_Ext';
diag_fileNames{36}='Vm_Ext';
diag_frequency(36)=8.64000000e+04;
diag_timePhase(36)=0.00000000e+00;
diag_fields{1,37}='Vm_AdvZ3';
diag_fileNames{37}='Vm_AdvZ3';
diag_frequency(37)=8.64000000e+04;
diag_timePhase(37)=0.00000000e+00;
diag_fields{1,38}='Vm_AdvRe';
diag_fileNames{38}='Vm_AdvRe';
diag_frequency(38)=8.64000000e+04;
diag_timePhase(38)=0.00000000e+00;
diag_fields{1,39}='TOTVTEND';
diag_fileNames{39}='TOTVTEND';
diag_frequency(39)=8.64000000e+04;
diag_timePhase(39)=0.00000000e+00;
diag_fields{1,40}='ETAN';
diag_fileNames{40}='ETAN_inst';
diag_frequency(40)=-4.32000000e+03;
diag_timePhase(40)=0.00000000e+00;
diag_fields{1,41}='UVEL';
diag_fileNames{41}='UVEL_inst';
diag_frequency(41)=-4.32000000e+03;
diag_timePhase(41)=0.00000000e+00;
diag_fields{1,42}='VVEL';
diag_fileNames{42}='VVEL_inst';
diag_frequency(42)=-4.32000000e+03;
diag_timePhase(42)=0.00000000e+00;
diag_fields{1,43}='WVEL';
diag_fileNames{43}='WVEL_inst';
diag_frequency(43)=-4.32000000e+03;
diag_timePhase(43)=0.00000000e+00;
diag_fields{1,44}='THETA';
diag_fileNames{44}='THETA_inst';
diag_frequency(44)=-4.32000000e+03;
diag_timePhase(44)=0.00000000e+00;
diag_fields{1,45}='SALT';
diag_fileNames{45}='SALT_inst';
diag_frequency(45)=-4.32000000e+03;
diag_timePhase(45)=0.00000000e+00;

