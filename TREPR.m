clear%function[PP3]=triplet_TREPR_Sp(D,E,pX,pY,pZ,dc,df,T2)
tic;
% This program is used for the TREPR lineshape analysis of the excited triplet states.
% Microwave power effect is not included.

ll=150; %%%number of points to simulate

time=0.2e-6; % time
T = 300;
ks=12e6;kss=1e1;kt=1e4;kd=1e1; % kS: 1TT->S0 kSS: 1(T+T)->S0 kt=2e4;Triplet lifetime 
%%%%% microwave frequency (MHz),        center magnetic field (mT), sweep width (mT)
	   u0=9591.0 ;u0=u0*2*pi*1e6;             Bav=342.20;                 Bw=160;

       dB=Bw/ll;         Bs=Bav-Bw/2;
% D value and E value : cm-1  
%dc=d_theta df=d_fai : increment in the field angles (degrees)
%T2 = phase memory time (s)
D=1495*1/29979.3; E=-65*1/29979.3; dc=3; df=3; T22=7e-9;T2=6e-9; % differences in field angle % spectrum width parameters
DDT=-155/29979.3; % Interspin D in MHz   TT
DDT2=-45/29979.3; % Interspin D  in MHz interspin T+T 

%polar angles with respect to the Z and X in the ZFS interaction (degrees)
%Not in use
deltam=0;gammam=0;
deltam=deltam/360*2*pi;gammam=gammam/360*2*pi;

D=D*2*pi*29979.3*1e6;    E=E*2*pi*29979.3*1e6;  %(rad/s)
HD(1,1)=-D/3+E; HD(2,2)=-D/3-E; HD(3,3)=D*2/3;

%inter-triplet dipolar tensor
DDT=DDT*2*pi*29979.3*1e6; 
HDD(1,1)=-DDT/3; HDD(2,2)=-DDT/3; HDD(3,3)=DDT*2/3;

% Exchange Couplings / gauss
			J2=-35000;   J2=J2*2*pi*28.0249*0.1*1e6;%(rad/s)TT1
            J3=-3500;   J3=J3*2*pi*28.0249*0.1*1e6;%(rad/s)TT2
            J4=0; J4=J4*2*pi*28.0249*0.1*1e6;%(rad/s) 1.5G(0.4us), -0.5G(1us), 3G(2us) 3G T+T
% correlation time        T+T dissociation rate
            tau=3.5e-12; % molecular motion to transition state
            kJT=1.7e6; % Rate constant of singlet T+T dissociation 
            kJTQ=0.8e6; % Rate constant of quinitet T+T dissociation 
            deltaE=40; %cm-1  activation energy in the TT -> T...T
            
           kJ = 1/tau; kJQ= kJ*exp(deltaE/0.6950/T);  % molecular motion back to large J 
           kJTT= kJT*exp(-130/0.6950/T); %% T+T -> TT back-dissociation in quintet TT  (-110 cm-1 for parafin)
           kJTTS= kJT*exp(-130/0.6950/T);% back-dissociation reactions in Singlet TT
           %modulation frequency in cm-1
           vib=(kJ+kJQ)/2/pi/3e10
           %  dissociation-quintet equilibrium time constant
    time_eq=1/(kJT+kJTTS);

kab=1/1e-1; % spin lattice relaxation time
kbc=1/1e-1; % spin lattice relaxation time

kba=kab*exp(-(Bav*2*pi*28.0249*1e6)*0.3/6e10/0.6950/T);
kcb=kbc*exp(-(Bav*2*pi*28.0249*1e6)*0.3/6e10/0.6950/T);
            
%B0= 8000; %External Magnetic Field (gauss)
%B1= 100;   %Micromave Field Strength (gauss)

%Initial density matirix condition
%roh(1,1)=pX; roh(2,2)=pY; roh(3,3)=pZ; 
%roh(1,2)=0;roh(1,3)=0;roh(2,1)=0;roh(2,3)=0;roh(3,1)=0;roh(3,2)=0;
% Sublevel populations
%roh00(1,1)=0.247-0.0;roh00(2,2)=0.00-0.0;roh00(3,3)=0.419+0.0;roh00(4,4)=0.021+0.0;roh00(5,5)=0.312;
%roh00(1,1)=0.065-0.0;roh00(2,2)=0.00-0.0;roh00(3,3)=0.42+0.0;roh00(4,4)=0.049+0.0;roh00(5,5)=0.465

      PP3=[];PP4=[];PP03=[];
   for j=1:ll
      PP3=[PP3,0];PP4=[PP4,0];PP03=[PP03,0];
   end

%%%Euler angles (deg) of the D2 tensor principal axis along the D1 (X,Y,Z) 
%% Conformation in TT1 (T1 vs T2)

		f0=0; f0=f0/360*2*pi; % rotation around the Z axis 
		g0=0; g0=g0/360*2*pi; % rotation around the X' axis  %2 @ 0.4 us  42 deg at 2.5 us
		h0=-0; h0=h0/360*2*pi; % rotation around the Z'' axis

		f1=90; f1=f1/360*2*pi; % rotation around the Z axis 
		g1=0; g1=g1/360*2*pi; % rotation around the X' axis  %2 @ 0.4 us  42 deg at 2.5 us
		h1=-90; h1=h1/360*2*pi; % rotation around the Z'' axis

      %% Conformation in TT2 (T1 and T2)
		f01=10; f01=f01/360*2*pi; % rotation around the Z axis 
		g01=30; g01=g01/360*2*pi; % rotation around the X' axis  %2 @ 0.4 us  42 deg at 2.5 us
		h01=-10; h01=h01/360*2*pi; % rotation around the Z'' axis
      
        f2=90; f2=f2/360*2*pi; % rotation around the Z axis 
		g2=-30; g2=g2/360*2*pi; % rotation around the X' axis  %12 deg for TcD-BiPh
		h2=-90; h2=h2/360*2*pi; % rotation around the Z'' axis
  
        %% C63 T+T state
		f02=0; f02=f02/360*2*pi; % rotation around the Z axis 
		g02=0; g02=g02/360*2*pi; % rotation around the X' axis  %2 @ 0.4 us  42 deg at 2.5 us
		h02=0; h02=h02/360*2*pi; % rotation around the Z'' axis

		f3=-30; f3=f3/360*2*pi; % rotation around the Z axis 
		g3=105; g3=g3/360*2*pi; % rotation around the X' axis  %70 deg for PcD-BiPh
		h3=60; h3=h3/360*2*pi; % rotation around the Z'' axis

        %% TT dipole direction Z' (polar angles in deg) with respect to X_Y_Z system in Triplet1
		c22=20+0;	% theta
        v22=160-0; % phai       
        c22=c22/360*2*pi; v22=v22/360*2*pi;  %p22=p22/360*2*pi;   

UMT1(1,1)=sin(v22); UMT1(2,1)=-cos(v22); UMT1(3,1)=0;
UMT1(1,2)=cos(c22)*cos(v22); UMT1(2,2)=cos(c22)*sin(v22); UMT1(3,2)=-sin(c22);
UMT1(1,3)=sin(c22)*cos(v22); UMT1(2,3)=sin(c22)*sin(v22); UMT1(3,3)=cos(c22);

%% T+T dipole direction Z'' (polar angles in deg) with resoect to gX_gY_gZ system in FAD
		c11=100.7;	% theta
        v11=180+30; % phai       
        c11=c11/360*2*pi; v11=v11/360*2*pi;  %p22=p22/360*2*pi;   
        
UMT2(1,1)=sin(v11); UMT2(2,1)=-cos(v11); UMT2(3,1)=0;
UMT2(1,2)=cos(c11)*cos(v11); UMT2(2,2)=cos(c11)*sin(v11); UMT2(3,2)=-sin(c11);
UMT2(1,3)=sin(c11)*cos(v11); UMT2(2,3)=sin(c11)*sin(v11); UMT2(3,3)=cos(c11);

% HDD11=UMT1*HDD*UMT1'; % Inter-triplet ZFS tensor
% HDD22=UMT2*HDD*UMT2'; % Inter-triplet ZFS tensor
        
        %%%% rotation matrix C60
Rc0(1,1)=cos(h0)*cos(f0)-cos(g0)*sin(f0)*sin(h0);
Rc0(1,2)=cos(h0)*sin(f0)+cos(g0)*cos(f0)*sin(h0);
Rc0(1,3)=sin(h0)*sin(g0);
Rc0(2,1)=-sin(h0)*cos(f0)-cos(g0)*sin(f0)*cos(h0);
Rc0(2,2)=-sin(h0)*sin(f0)+cos(g0)*cos(f0)*cos(h0);
Rc0(2,3)=cos(h0)*sin(g0);
Rc0(3,1)=sin(g0)*sin(f0);
Rc0(3,2)=-sin(g0)*cos(f0);
Rc0(3,3)=cos(g0);
HD0=(Rc0)\HD*(Rc0); % ZFS tensor in the first triplet in the T1 frame


        %%%% rotation matrix C60
Rc60(1,1)=cos(h1)*cos(f1)-cos(g1)*sin(f1)*sin(h1);
Rc60(1,2)=cos(h1)*sin(f1)+cos(g1)*cos(f1)*sin(h1);
Rc60(1,3)=sin(h1)*sin(g1);
Rc60(2,1)=-sin(h1)*cos(f1)-cos(g1)*sin(f1)*cos(h1);
Rc60(2,2)=-sin(h1)*sin(f1)+cos(g1)*cos(f1)*cos(h1);
Rc60(2,3)=cos(h1)*sin(g1);
Rc60(3,1)=sin(g1)*sin(f1);
Rc60(3,2)=-sin(g1)*cos(f1);
Rc60(3,3)=cos(g1);
HD2=(Rc60)\HD*(Rc60); % ZFS tensor in the second triplet in the T1 frame

        %%%% rotation matrix C60
Rc01(1,1)=cos(h01)*cos(f01)-cos(g01)*sin(f01)*sin(h01);
Rc01(1,2)=cos(h01)*sin(f01)+cos(g01)*cos(f01)*sin(h01);
Rc01(1,3)=sin(h01)*sin(g01);
Rc01(2,1)=-sin(h01)*cos(f01)-cos(g01)*sin(f01)*cos(h01);
Rc01(2,2)=-sin(h01)*sin(f01)+cos(g01)*cos(f01)*cos(h01);
Rc01(2,3)=cos(h01)*sin(g01);
Rc01(3,1)=sin(g01)*sin(f01);
Rc01(3,2)=-sin(g01)*cos(f01);
Rc01(3,3)=cos(g01);
HD03=(Rc01)\HD*(Rc01); % ZFS tensor in the first triplet in the T1 frame

Rc61(1,1)=cos(h2)*cos(f2)-cos(g2)*sin(f2)*sin(h2);
Rc61(1,2)=cos(h2)*sin(f2)+cos(g2)*cos(f2)*sin(h2);
Rc61(1,3)=sin(h2)*sin(g2);
Rc61(2,1)=-sin(h2)*cos(f2)-cos(g2)*sin(f2)*cos(h2);
Rc61(2,2)=-sin(h2)*sin(f2)+cos(g2)*cos(f2)*cos(h2);
Rc61(2,3)=cos(h2)*sin(g2);
Rc61(3,1)=sin(g2)*sin(f2);
Rc61(3,2)=-sin(g2)*cos(f2);
Rc61(3,3)=cos(g2);
HD3=(Rc61)\HD*(Rc61); % ZFS tensor in the second triplet in the T1 frame

       %%%% rotation matrix C63 T+T
Rc02(1,1)=cos(h02)*cos(f02)-cos(g02)*sin(f02)*sin(h02);
Rc02(1,2)=cos(h02)*sin(f02)+cos(g02)*cos(f02)*sin(h02);
Rc02(1,3)=sin(h02)*sin(g02);
Rc02(2,1)=-sin(h02)*cos(f02)-cos(g02)*sin(f02)*cos(h02);
Rc02(2,2)=-sin(h02)*sin(f02)+cos(g02)*cos(f02)*cos(h02);
Rc02(2,3)=cos(h02)*sin(g02);
Rc02(3,1)=sin(g02)*sin(f02);
Rc02(3,2)=-sin(g02)*cos(f02);
Rc02(3,3)=cos(g02);
HD04=(Rc02)\HD*(Rc02); % ZFS tensor in the first triplet in the T1 frame


Rc63(1,1)=cos(h3)*cos(f3)-cos(g3)*sin(f3)*sin(h3);
Rc63(1,2)=cos(h3)*sin(f3)+cos(g3)*cos(f3)*sin(h3);
Rc63(1,3)=sin(h3)*sin(g3);
Rc63(2,1)=-sin(h3)*cos(f3)-cos(g3)*sin(f3)*cos(h3);
Rc63(2,2)=-sin(h3)*sin(f3)+cos(g3)*cos(f3)*cos(h3);
Rc63(2,3)=cos(h3)*sin(g3);
Rc63(3,1)=sin(g3)*sin(f3);
Rc63(3,2)=-sin(g3)*cos(f3);
Rc63(3,3)=cos(g3);
HD4=(Rc63)\HD*(Rc63); % ZFS tensor in the second triplet in the T1 frame

% vetor of the inter-spin1 in the (X,Y,Z)
vectd1(1,1)=sin(c22)*cos(v22);
vectd1(1,2)=sin(c22)*sin(v22);
vectd1(1,3)=cos(c22);

vectd2(1,1)=sin(c11)*cos(v11);
vectd2(1,2)=sin(c11)*sin(v11);
vectd2(1,3)=cos(c11);

%%%% Spin Opearators in TT pairs %%%%%%% basis=(++,+0,+-,0+,00,0-,-+,-0,--)

S1X=[0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0;1,0,0,0,0,0,1,0,0;0,1,0,0,0,0,0,1,0;0,0,1,0,0,0,0,0,1,;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0]/sqrt(2);
S1Y=[0,0,0,-1,0,0,0,0,0;0,0,0,0,-1,0,0,0,0;0,0,0,0,0,-1,0,0,0;1,0,0,0,0,0,-1,0,0;0,1,0,0,0,0,0,-1,0;0,0,1,0,0,0,0,0,-1,;0,0,0,1,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0,0]*1i/sqrt(2);
S1Z=[1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,;0,0,0,0,0,0,-1,0,0;0,0,0,0,0,0,0,-1,0;0,0,0,0,0,0,0,0,-1];
S2X=[0,1,0,0,0,0,0,0,0;1,0,1,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0;0,0,0,1,0,1,0,0,0;0,0,0,0,1,0,0,0,0,;0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,1,0,1;0,0,0,0,0,0,0,1,0]/sqrt(2);
S2Y=[0,-1,0,0,0,0,0,0,0;1,0,-1,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0;0,0,0,1,0,-1,0,0,0;0,0,0,0,1,0,0,0,0,;0,0,0,0,0,0,0,-1,0;0,0,0,0,0,0,1,0,-1;0,0,0,0,0,0,0,1,0]*1i/sqrt(2);
S2Z=[1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,-1,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,-1];   

% STTX=S1X+S2X;
% STTY=S1Y+S2Y;
% STTZ=S1Z+S2Z;

%%%% Spin Opearators in T %%%%%%% basis=(+,0,-Åj
STX=[0,1,0;1,0,1;0,1,0]/sqrt(2);
STY=[0,-1i,0;1i,0,-1i;0,1i,0]/sqrt(2);
STZ=[1,0,0;0,0,0;0,0,-1];



%%%%% conversion matrix from the basis system of
%%%%% (++,+0,+-,0+,00,0-,-+,-0,--) to the system in the high field limit
%%%%% (S,T-,T0,T+,Q-2,Q-,Q0,Q+1,Q+2)

UU=[0,0,-1/sqrt(3),0,1/sqrt(3),0,-1/sqrt(3),0,0;0,0,0,0,0,-1/sqrt(2),0,1/sqrt(2),0; 0,0,1/sqrt(2),0,0,0,-1/sqrt(2),0,0; 0,1/sqrt(2),0,-1/sqrt(2),0,0,0,0,0; 0,0,0,0,0,0,0,0,1; 0,0,0,0,0,1/sqrt(2),0,1/sqrt(2),0; 0,0,1/sqrt(6),0,2/sqrt(6),0,1/sqrt(6),0,0; 0,1/sqrt(2),0,1/sqrt(2),0,0,0,0,0; 1,0,0,0,0,0,0,0,0];
UU1=[1,0,0,0,0,0,0,0,0;0,1/sqrt(2),0,1/sqrt(2),0,0,0,0,0;0,0,1/sqrt(6),0,2/sqrt(6),0,1/sqrt(6),0,0;0,0,0,0,0,1/sqrt(2),0,1/sqrt(2),0;0,0,0,0,0,0,0,0,1;0,1/sqrt(2),0,-1/sqrt(2),0,0,0,0,0;0,0,1/sqrt(2),0,0,0,-1/sqrt(2),0,0;0,0,0,0,0,-1/sqrt(2),0,1/sqrt(2),0;0,0,-1/sqrt(3),0,1/sqrt(3),0,-1/sqrt(3),0,0];

%UU1*(S1Y+S2Y)*UU1';
%%%%%%%% Angles %%%%%%%%%%
% difference of the angle 'theta' (deg)
	 m=90/dc; dc=dc/360*2*pi;
% difference of the angle 'phai'  (deg)
if abs(E) < 0.000001*2*pi*39979.3*1e6
    mm=0;
else
     mm=m; df=dc;
end
%m=1;

Bav1=Bav*2*pi*28.0249*1e6;%(rad/s)

ZFSmag1=[0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0];
ZFSmag51=0;ZFSmag61=0;ZFSmag71=0;ZFSmag81=0;ZFSmag91=0;
     PP104=0;PP105=0;PP1055=0;PP107=0;
     
     PP1056=[0;0;0;0;0;0;0;0;0;0;0;0];PP1057=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
     PP1058=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
     PP1059=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];PP1060=[0;0;0;0;0;0;0;0;0;0;0;0];


% Angle averaging
    for kk=2:2*m+1
        for l=1:2*mm+1

            y=(l-1)*df; h=(kk-1)*dc; q=0.5*pi;

            
      PP0=[];
      PP=[];
      PPT=[];TTXX=[];TTZZ=[];TTZZZ=[];
% Orientation of the tripet state with respect to X,Y,Z laboratory frame

%c22=44;	v22=24; p22=0; c22=c22/360*2*pi; v22=v22/360*2*pi;  p22=p22/360*2*pi;   % Z' orientation of RP with respect to XYZ in C60 triplet



UM(1,1)=sin(y); UM(2,1)=-cos(y); UM(3,1)=0;
UM(1,2)=cos(h)*cos(y); UM(2,2)=cos(h)*sin(y); UM(3,2)=-sin(h);
UM(1,3)=sin(h)*cos(y); UM(2,3)=sin(h)*sin(y); UM(3,3)=cos(h);

% Direction of the magnetic field:
vectB(1,1)=sin(h)*cos(y);
vectB(2,1)=sin(h)*sin(y);
vectB(3,1)=cos(h);

% direction-cosine between the transition dipolar moment and the magnetic
% field
cosinMB1=vectd1(1,1)*vectB(1,1)+vectd1(1,2)*vectB(2,1)+vectd1(1,3)*vectB(3,1);
cosinMB2=vectd2(1,1)*vectB(1,1)+vectd2(1,2)*vectB(2,1)+vectd2(1,3)*vectB(3,1);


%%% ZFS tensors in the lab frames
HDL1=(UM)\HD*(UM);
HDL0=(UM)\HD0*(UM); HDL2=(UM)\HD2*(UM);  % original Q
HDL03=(UM)\HD03*(UM); HDL3=(UM)\HD3*(UM);  % Transition State 
HDL04=(UM)\HD04*(UM); HDL4=(UM)\HD4*(UM);  % T + T

% HDDL1=(UM)\HDD11*(UM); % interspin Q
% HDDL2=(UM)\HDD11*(UM); % interspin transition Q
% HDDL3=(UM)\HDD22*(UM); % interspin T + T


% Eigenvalue Calc%%%%%%%%%%%%%%

HTTZ = Bav1*(S1Z+S2Z);%Zeemann term in the multiexciton

HTTZFS1=HDL1(1,1)*S1X*S1X+HDL1(1,2)*S1X*S1Y+HDL1(1,3)*S1X*S1Z+HDL1(2,1)*S1Y*S1X+HDL1(2,2)*S1Y*S1Y+HDL1(2,3)*S1Y*S1Z+HDL1(3,1)*S1Z*S1X+HDL1(3,2)*S1Z*S1Y+HDL1(3,3)*S1Z*S1Z;
HTTZFS0=HDL0(1,1)*S1X*S1X+HDL0(1,2)*S1X*S1Y+HDL0(1,3)*S1X*S1Z+HDL0(2,1)*S1Y*S1X+HDL0(2,2)*S1Y*S1Y+HDL0(2,3)*S1Y*S1Z+HDL0(3,1)*S1Z*S1X+HDL0(3,2)*S1Z*S1Y+HDL0(3,3)*S1Z*S1Z;
HTTZFS2=HDL2(1,1)*S2X*S2X+HDL2(1,2)*S2X*S2Y+HDL2(1,3)*S2X*S2Z+HDL2(2,1)*S2Y*S2X+HDL2(2,2)*S2Y*S2Y+HDL2(2,3)*S2Y*S2Z+HDL2(3,1)*S2Z*S2X+HDL2(3,2)*S2Z*S2Y+HDL2(3,3)*S2Z*S2Z;
HTTZFS03=HDL03(1,1)*S1X*S1X+HDL03(1,2)*S1X*S1Y+HDL03(1,3)*S1X*S1Z+HDL03(2,1)*S1Y*S1X+HDL03(2,2)*S1Y*S1Y+HDL03(2,3)*S1Y*S1Z+HDL03(3,1)*S1Z*S1X+HDL03(3,2)*S1Z*S1Y+HDL03(3,3)*S1Z*S1Z;
HTTZFS3=HDL3(1,1)*S2X*S2X+HDL3(1,2)*S2X*S2Y+HDL3(1,3)*S2X*S2Z+HDL3(2,1)*S2Y*S2X+HDL3(2,2)*S2Y*S2Y+HDL3(2,3)*S2Y*S2Z+HDL3(3,1)*S2Z*S2X+HDL3(3,2)*S2Z*S2Y+HDL3(3,3)*S2Z*S2Z;
HTTZFS04=HDL04(1,1)*S1X*S1X+HDL04(1,2)*S1X*S1Y+HDL04(1,3)*S1X*S1Z+HDL04(2,1)*S1Y*S1X+HDL04(2,2)*S1Y*S1Y+HDL04(2,3)*S1Y*S1Z+HDL04(3,1)*S1Z*S1X+HDL04(3,2)*S1Z*S1Y+HDL04(3,3)*S1Z*S1Z;
HTTZFS4=HDL4(1,1)*S2X*S2X+HDL4(1,2)*S2X*S2Y+HDL4(1,3)*S2X*S2Z+HDL4(2,1)*S2Y*S2X+HDL4(2,2)*S2Y*S2Y+HDL4(2,3)*S2Y*S2Z+HDL4(3,1)*S2Z*S2X+HDL4(3,2)*S2Z*S2Y+HDL4(3,3)*S2Z*S2Z;

% Inter-triplet ZFS
% HTTZFSTT1=HDDL1(1,1)*STTX*STTX+HDDL1(1,2)*STTX*STTY+HDDL1(1,3)*STTX*STTZ+HDDL1(2,1)*STTY*STTX+HDDL1(2,2)*STTY*STTY+HDDL1(2,3)*STTY*STTZ+HDDL1(3,1)*STTZ*STTX+HDDL1(3,2)*STTZ*STTY+HDDL1(3,3)*STTZ*STTZ;
% HTTZFSTT2=HDDL2(1,1)*STTX*STTX+HDDL2(1,2)*STTX*STTY+HDDL2(1,3)*STTX*STTZ+HDDL2(2,1)*STTY*STTX+HDDL2(2,2)*STTY*STTY+HDDL2(2,3)*STTY*STTZ+HDDL2(3,1)*STTZ*STTX+HDDL2(3,2)*STTZ*STTY+HDDL2(3,3)*STTZ*STTZ;
% HTTZFSTT3=HDDL3(1,1)*STTX*STTX+HDDL3(1,2)*STTX*STTY+HDDL3(1,3)*STTX*STTZ+HDDL3(2,1)*STTY*STTX+HDDL3(2,2)*STTY*STTY+HDDL3(2,3)*STTY*STTZ+HDDL3(3,1)*STTZ*STTX+HDDL3(3,2)*STTZ*STTY+HDDL3(3,3)*STTZ*STTZ;

HTTZFST1=(3*S1Z*S2Z-(S1X*S2X+S1Y*S2Y+S1Z*S2Z))*(cosinMB1^2-1/3)*DDT;
HTTZFST2=(3*S1Z*S2Z-(S1X*S2X+S1Y*S2Y+S1Z*S2Z))*(cosinMB2^2-1/3)*DDT2;


HTTZFS=(HTTZFS0+HTTZFS2)*0.95;% Strongly Coupled TT
HTTZFS22=(HTTZFS03+HTTZFS3)*0.95;% Transition State 
HTTZFS44=(HTTZFS04+HTTZFS4)*0.97;% T+T


%Spin-spin exchange

HTTEX= -2*J2*(S1X*S2X+S1Y*S2Y+S1Z*S2Z);
HTTEX2= -2*J3*(S1X*S2X+S1Y*S2Y+S1Z*S2Z);
HTTEX3= -2*J4*(S1X*S2X+S1Y*S2Y+S1Z*S2Z);


H00=HTTZ+HTTZFS+HTTEX+HTTZFST1;  % original Quntet TT1
H01=HTTZ+HTTZFS22+HTTEX2+HTTZFST1; % Trasition State  TT2
H02=HTTZ+HTTZFS44+HTTEX3+HTTZFST2;  % T + T

H00=UU1*H00*UU1';
H00Q=UU1*H01*UU1';  % Transition State Spin Hamiltonian in the Q2,Q1,Q0,....
H00TT=UU1*H02*UU1';

[QQ0,lm0]=eig2(H00);  % quintet in the original state

H01Q=QQ0'*H00Q*QQ0; % transition state in the origonal basis 


H0TT=QQ0'*H00TT*QQ0; % dissociated T+T in the origonal basis
[QQb,lmb]=eig2(H0TT);  % nine diagonalized functions in the T+T state
%QQb: eigenvectors from the Q+2, Q+1, .... basis


% non-correlated TT
%H00T=HTTZ+HTTZFS;


%[Qa,lmbb]=eigs(H01,9,1e16);  % quintet

% HTZ=Bav1*(STZ);%Zeemann term
 HTZFS=HDL1(1,1)*STX*STX+HDL1(1,2)*STX*STY+HDL1(1,3)*STX*STZ+HDL1(2,1)*STY*STX+HDL1(2,2)*STY*STY+HDL1(2,3)*STY*STZ+HDL1(3,1)*STZ*STX+HDL1(3,2)*STZ*STY+HDL1(3,3)*STZ*STZ;
% HT0=HTZ+HTZFS;

%s(H0,3,1e12);  %Diagonalization of Hamiltonian

%roh00=Q'*roh*Q;  % Initial density matrix of each eigenstate T+, T0, T-
% charge recombination to S -> T0 state
% roh00(1,1)=0;roh00(2,2)=1;roh00(3,3)=0; 
ZFSmag=abs(UU*HTTZFS*UU');
ZFSmag1=ZFSmag1+ZFSmag;
ZFSmag51=ZFSmag51+ZFSmag(5,1)^2;
ZFSmag61=ZFSmag61+ZFSmag(6,1)^2;
ZFSmag71=ZFSmag71+ZFSmag(7,1)^2;
ZFSmag81=ZFSmag81+ZFSmag(8,1)^2;
ZFSmag91=ZFSmag91+ZFSmag(9,1)^2;

%%%%%%%%%%%%%% Spin Decoherent Matrix in (1,2,3,4,5,6) system
%%% [Q2,Q1,Q0,Q-1,Q-2,T1,T0,T-, S]




% k16=ZFSmag(1,9)^2*tau/(1+((2*Bav*2*pi*28.0249*1e6-6*J2)*tau)^2); k61=k16*exp(-(2*Bav*2*pi*28.0249*1e6-6*J2)*0.3/6e10/0.6950/T);
% k26=ZFSmag(1,8)^2*tau/(1+((Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);   k62=k26*exp(-(Bav*2*pi*28.0249*1e6-6*J2)*0.3/6e10/0.6950/T);
% k36=ZFSmag(1,7)^2*tau/(1+((-6*J2)*tau)^2);                       k63=k36*exp(-(-6*J2)*0.3/6e10/0.6950/T);
% k46=ZFSmag(1,6)^2*tau/(1+((-Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);  k64=k46*exp(-(-Bav*2*pi*28.0249*1e6-6*J2)*0.3/6e10/0.6950/T);
% k56=ZFSmag(1,5)^2*tau/(1+((-2*Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);k65=k56*exp(-(-2*Bav*2*pi*28.0249*1e6-6*J2)*0.3/6e10/0.6950/T);
% 
% k1a=ZFSmag(4,9)^2*tauT/(1+((1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); ka1=k1a*exp(-(1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k1b=ZFSmag(3,9)^2*tauT/(1+((2*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kb1=k1b*exp(-(2*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k1c=ZFSmag(2,9)^2*tauT/(1+((3*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kc1=k1c*exp(-(3*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% 
% k2a=ZFSmag(4,8)^2*tauT/(1+((0*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); ka2=k2a*exp(-(0*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k2b=ZFSmag(3,8)^2*tauT/(1+((1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kb2=k2b*exp(-(1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k2c=ZFSmag(2,8)^2*tauT/(1+((2*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kc2=k2c*exp(-(2*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% 
% k3a=ZFSmag(4,7)^2*tauT/(1+((-1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);ka3=k3a*exp(-(-1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k3b=ZFSmag(3,7)^2*tauT/(1+((0*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kb3=k3b*exp(-(0*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k3c=ZFSmag(2,7)^2*tauT/(1+((1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kc3=k3c*exp(-(1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% 
% k4a=ZFSmag(4,6)^2*tauT/(1+((-2*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);ka4=k4a*exp(-(-2*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k4b=ZFSmag(3,6)^2*tauT/(1+((-1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);kb4=k4b*exp(-(-1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k4c=ZFSmag(2,6)^2*tauT/(1+((0*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2); kc4=k4c*exp(-(0*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% 
% k5a=ZFSmag(4,5)^2*tauT/(1+((-3*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);ka5=k5a*exp(-(-3*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k5b=ZFSmag(3,5)^2*tauT/(1+((-2*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);kb5=k5b*exp(-(-2*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);
% k5c=ZFSmag(2,5)^2*tauT/(1+((-1*Bav*2*pi*28.0249*1e6-4*J2)*tauT)^2);kc5=k5c*exp(-(-1*Bav*2*pi*28.0249*1e6-4*J2)*0.3/6e10/0.6950/T);

%kJ:k12, kJT:kDISS, kJQ:k21   kJTT:kTTback       ks:k1  
%H01Q(1,9)=3e9;H01Q(9,1)=3e9;
%      QQ(largeJ)                                         QS(largeJ)                                                              SQ(largeJ)                                    SS(largeJ)                                QQ(smallJ)                                          QS(smallJ)                                                                                        SQ(smallJ)                                                                           SS(smallJ)                      QQ(T+T)                                                                 QS(T+T)                                                                       SQ(T+T)                                                                                  SS(T+T)
RR1=[-kJ-kJTQ 0 0 0 kJQ 0 0 0 kJTT 0 0 0; 0 -kJ-ks/2-1i*(lm0(1,1)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 +kJTT 0 0; 0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(1,1))-kJT 0 0 0 kJQ 0 0 0 +kJTT 0; 0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 kJTTS; kJ 0 0 0 -kJQ 1i*H01Q(9,1) -1i*H01Q(1,9) 0 0 0 0 0; 0 kJ 0 0 1i*H01Q(1,9) -1i*(H01Q(1,1)-H01Q(9,9))-kJQ 0 -1i*H01Q(1,9) 0 0 0 0; 0 0 kJ 0 -1i*H01Q(9,1) 0 -1i*(H01Q(9,9)-H01Q(1,1))-kJQ 1i*H01Q(9,1) 0 0 0 0; 0 0 0 kJ 0 -1i*H01Q(9,1) 1i*H01Q(1,9) -kJQ 0 0 0 0; kJTQ 0 0 0 0 0 0 0 -kJTT 1i*H0TT(9,1) -1i*H0TT(1,9) 0; 0 kJT 0 0 0 0 0 0 1i*H0TT(1,9) -1i*(H0TT(1,1)-H0TT(9,9))-kss/2-kJTT 0 -1i*H0TT(1,9); 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,1) 0 -1i*(H0TT(9,9)-H0TT(1,1))-kss/2-kJTT 1i*H0TT(9,1);0 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,1) 1i*H0TT(1,9) -kJTTS-kss];
% RR2=[-kJ-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0 0; 0 -kJ-ks/2-1i*(lm0(2,2)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0; 0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(2,2))-kJT 0 0 0 kJQ 0 0 0 kJTT 0; 0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 kJTT; kJ 0 0 0 -kJQ 1i*H01Q(9,2) -1i*H01Q(2,9) 0 0 0 0 0; 0 kJ 0 0 1i*H01Q(2,9) -1i*(H01Q(2,2)-H01Q(9,9))-kJQ 0 -1i*H01Q(2,9) 0 0 0 0; 0 0 kJ 0 -1i*H01Q(9,2) 0 -1i*(H01Q(9,9)-H01Q(2,2))-kJQ 1i*H01Q(9,2) 0 0 0 0; 0 0 0 kJ 0 -1i*H01Q(9,2) 1i*H01Q(2,9) -kJQ 0 0 0 0; kJT 0 0 0 0 0 0 0 -kJTT 1i*H0TT(9,2) -1i*H0TT(2,9) 0; 0 kJT 0 0 0 0 0 0 1i*H0TT(2,9) -1i*(H0TT(2,2)-H0TT(9,9))-kJTT-kss/2 0 -1i*H0TT(2,9); 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,2) 0 -1i*(H0TT(9,9)-H0TT(2,2))-kJTT-kss/2 1i*H0TT(9,2);0 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,2) 1i*H0TT(2,9) -kJTT-kss];
% RR3=[-kJ-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0 0; 0 -kJ-ks/2-1i*(lm0(3,3)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0; 0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(3,3))-kJT 0 0 0 kJQ 0 0 0 kJTT 0; 0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 kJTT; kJ 0 0 0 -kJQ 1i*H01Q(9,3) -1i*H01Q(3,9) 0 0 0 0 0; 0 kJ 0 0 1i*H01Q(3,9) -1i*(H01Q(3,3)-H01Q(9,9))-kJQ 0 -1i*H01Q(3,9) 0 0 0 0; 0 0 kJ 0 -1i*H01Q(9,3) 0 -1i*(H01Q(9,9)-H01Q(3,3))-kJQ 1i*H01Q(9,3) 0 0 0 0; 0 0 0 kJ 0 -1i*H01Q(9,3) 1i*H01Q(3,9) -kJQ 0 0 0 0; kJT 0 0 0 0 0 0 0 -kJTT 1i*H0TT(9,3) -1i*H0TT(3,9) 0; 0 kJT 0 0 0 0 0 0 1i*H0TT(3,9) -1i*(H0TT(3,3)-H0TT(9,9))-kJTT-kss/2 0 -1i*H0TT(3,9); 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,3) 0 -1i*(H0TT(9,9)-H0TT(3,3))-kJTT-kss/2 1i*H0TT(9,3);0 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,3) 1i*H0TT(3,9) -kJTT-kss];
% RR4=[-kJ-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0 0; 0 -kJ-ks/2-1i*(lm0(4,4)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 kJTT 0 0; 0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(4,4))-kJT 0 0 0 kJQ 0 0 0 kJTT 0; 0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 kJTT; kJ 0 0 0 -kJQ 1i*H01Q(9,4) -1i*H01Q(4,9) 0 0 0 0 0; 0 kJ 0 0 1i*H01Q(4,9) -1i*(H01Q(4,4)-H01Q(9,9))-kJQ 0 -1i*H01Q(4,9) 0 0 0 0; 0 0 kJ 0 -1i*H01Q(9,4) 0 -1i*(H01Q(9,9)-H01Q(4,4))-kJQ 1i*H01Q(9,4) 0 0 0 0; 0 0 0 kJ 0 -1i*H01Q(9,4) 1i*H01Q(4,9) -kJQ 0 0 0 0; kJT 0 0 0 0 0 0 0 -kJTT 1i*H0TT(9,4) -1i*H0TT(4,9) 0; 0 kJT 0 0 0 0 0 0 1i*H0TT(4,9) -1i*(H0TT(4,4)-H0TT(9,9))-kJTT-kss/2 0 -1i*H0TT(4,9); 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,4) 0 -1i*(H0TT(9,9)-H0TT(4,4))-kJTT-kss/2 1i*H0TT(9,4);0 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,4) 1i*H0TT(4,9) -kJTT-kss];
RR5=[-kJ-kJTQ 0 0 0 kJQ 0 0 0 kJTT 0 0 0; 0 -kJ-ks/2-1i*(lm0(5,5)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 +kJTT 0 0; 0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(5,5))-kJT 0 0 0 kJQ 0 0 0 +kJTT 0; 0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 kJTTS; kJ 0 0 0 -kJQ 1i*H01Q(9,5) -1i*H01Q(5,9) 0 0 0 0 0; 0 kJ 0 0 1i*H01Q(5,9) -1i*(H01Q(5,5)-H01Q(9,9))-kJQ 0 -1i*H01Q(5,9) 0 0 0 0; 0 0 kJ 0 -1i*H01Q(9,5) 0 -1i*(H01Q(9,9)-H01Q(5,5))-kJQ 1i*H01Q(9,5) 0 0 0 0; 0 0 0 kJ 0 -1i*H01Q(9,5) 1i*H01Q(5,9) -kJQ 0 0 0 0; kJTQ 0 0 0 0 0 0 0 -kJTT 1i*H0TT(9,5) -1i*H0TT(5,9) 0; 0 kJT 0 0 0 0 0 0 1i*H0TT(5,9) -1i*(H0TT(5,5)-H0TT(9,9))-kss/2-kJTT 0 -1i*H0TT(5,9); 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,5) 0 -1i*(H0TT(9,9)-H0TT(5,5))-kss/2-kJTT 1i*H0TT(9,5);0 0 0 kJT 0 0 0 0 0 -1i*H0TT(9,5) 1i*H0TT(5,9) -kJTTS-kss];

%%%%%% initial condition assuming excited singlet state (pIII(7,1)) precursor
%% QQA,QSA,SQA,SSA,QQB,QSB,SQB, SSB
 pIII(1,1)=0;pIII(2,1)=0;pIII(3,1)=0;pIII(4,1)=1;pIII(5,1)=0; pIII(6,1)=0; pIII(7,1)=0; pIII(8,1)=0; pIII(9,1)=0; pIII(10,1)=0; pIII(11,1)=0; pIII(12,1)=0; 
            [VPP1,DPP1] = eig(RR1);
%             [VPP2,DPP2] = eig(RR2);
%             [VPP3,DPP3] = eig(RR3);
%             [VPP4,DPP4] = eig(RR4);
            [VPP5,DPP5] = eig(RR5);
               

            DDEXP1(1,1)=exp(DPP1(1,1)*time); DDEXP1(2,2)=exp(DPP1(2,2)*time); DDEXP1(3,3)=exp(DPP1(3,3)*time); 
            DDEXP1(4,4)=exp(DPP1(4,4)*time); DDEXP1(5,5)=exp(DPP1(5,5)*time); DDEXP1(6,6)=exp(DPP1(6,6)*time); 
            DDEXP1(7,7)=exp(DPP1(7,7)*time); DDEXP1(8,8)=exp(DPP1(8,8)*time); DDEXP1(9,9)=exp(DPP1(9,9)*time); 
            DDEXP1(10,10)=exp(DPP1(10,10)*time); 
            DDEXP1(11,11)=exp(DPP1(11,11)*time); DDEXP1(12,12)=exp(DPP1(12,12)*time); 
            
%             DDEXP2(1,1)=exp(DPP2(1,1)*time); DDEXP2(2,2)=exp(DPP2(2,2)*time); DDEXP2(3,3)=exp(DPP2(3,3)*time); 
%             DDEXP2(4,4)=exp(DPP2(4,4)*time); DDEXP2(5,5)=exp(DPP2(5,5)*time); DDEXP2(6,6)=exp(DPP2(6,6)*time); 
%             DDEXP2(7,7)=exp(DPP2(7,7)*time); DDEXP2(8,8)=exp(DPP2(8,8)*time); DDEXP2(9,9)=exp(DPP2(9,9)*time); 
%             DDEXP2(10,10)=exp(DPP2(10,10)*time); 
%             DDEXP2(11,11)=exp(DPP2(11,11)*time); DDEXP2(12,12)=exp(DPP2(12,12)*time); 
%             
%             DDEXP3(1,1)=exp(DPP3(1,1)*time); DDEXP3(2,2)=exp(DPP3(2,2)*time); DDEXP3(3,3)=exp(DPP3(3,3)*time); 
%             DDEXP3(4,4)=exp(DPP3(4,4)*time); DDEXP3(5,5)=exp(DPP3(5,5)*time); DDEXP3(6,6)=exp(DPP3(6,6)*time); 
%             DDEXP3(7,7)=exp(DPP3(7,7)*time); DDEXP3(8,8)=exp(DPP3(8,8)*time); DDEXP3(9,9)=exp(DPP3(9,9)*time); 
%             DDEXP3(10,10)=exp(DPP3(10,10)*time); 
%             DDEXP3(11,11)=exp(DPP3(11,11)*time); DDEXP3(12,12)=exp(DPP3(12,12)*time); 
%             
%             DDEXP4(1,1)=exp(DPP4(1,1)*time); DDEXP4(2,2)=exp(DPP4(2,2)*time); DDEXP4(3,3)=exp(DPP4(3,3)*time); 
%             DDEXP4(4,4)=exp(DPP4(4,4)*time); DDEXP4(5,5)=exp(DPP4(5,5)*time); DDEXP4(6,6)=exp(DPP4(6,6)*time); 
%             DDEXP4(7,7)=exp(DPP4(7,7)*time); DDEXP4(8,8)=exp(DPP4(8,8)*time); DDEXP4(9,9)=exp(DPP4(9,9)*time); 
%             DDEXP4(10,10)=exp(DPP4(10,10)*time); 
%             DDEXP4(11,11)=exp(DPP4(11,11)*time); DDEXP4(12,12)=exp(DPP4(12,12)*time); 
            
            DDEXP5(1,1)=exp(DPP5(1,1)*time); DDEXP5(2,2)=exp(DPP5(2,2)*time); DDEXP5(3,3)=exp(DPP5(3,3)*time); 
            DDEXP5(4,4)=exp(DPP5(4,4)*time); DDEXP5(5,5)=exp(DPP5(5,5)*time); DDEXP5(6,6)=exp(DPP5(6,6)*time); 
            DDEXP5(7,7)=exp(DPP5(7,7)*time); DDEXP5(8,8)=exp(DPP5(8,8)*time); DDEXP5(9,9)=exp(DPP5(9,9)*time); 
            DDEXP5(10,10)=exp(DPP5(10,10)*time); 
            DDEXP5(11,11)=exp(DPP5(11,11)*time); DDEXP5(12,12)=exp(DPP5(12,12)*time); 
            
     RR3=[-kJ-kJTQ 0 0 0 kJQ 0 0 0 kJTT 0 0 0 0 0 0 0 0 0;...% QQ(largeJ)
       0 -kJ-ks/2-1i*(lm0(3,3)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 0 0 kJTT 0 0 0 0 0 0;... %QS(largeJ)
       0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(3,3))-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 +kJTT 0 0;... %SQ(largeJ)
       0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 kJTTS 0;... %SS(largeJ)
       
       kJ 0 0 0 -kJQ 1i*H01Q(9,3) -1i*H01Q(3,9) 0 0 0 0 0 0 0 0 0 0 0;...  %QQ(smallJ)
       0 kJ 0 0 1i*H01Q(3,9) -1i*(H01Q(3,3)-H01Q(9,9))-kJQ 0 -1i*H01Q(3,9) 0 0 0 0 0 0 0 0 0 0;...  %QS(smallJ)
       0 0 kJ 0 -1i*H01Q(9,3) 0 -1i*(H01Q(9,9)-H01Q(3,3))-kJQ 1i*H01Q(9,3) 0 0 0 0 0 0 0 0 0 0;...  %SQ(smallJ)
       0 0 0 kJ 0 -1i*H01Q(9,3) 1i*H01Q(3,9) -kJQ 0 0 0 0 0 0 0 0 0 0;...  % SS(smallJ)
       
       kJTQ 0 0 0 0 0 0 0  -kJTT 1i*H0TT(7,3) 1i*H0TT(9,3) -1i*H0TT(3,7) 0 0 -1i*H0TT(3,9) 0 0 0;...  % QQ(T+T)
       0 0 0 0 0 0 0 0   1i*H0TT(3,7) -1i*(H0TT(3,3)-H0TT(7,7))-kt/2-kJTT 1i*H0TT(9,7) 0 -1i*H0TT(3,7) 0 0 -1i*H0TT(3,9) 0 0;...  % QT(T+T)
       0 kJT 0 0 0 0 0 0  1i*H0TT(3,9) 1i*H0TT(7,9) -1i*(H0TT(3,3)-H0TT(9,9))-kJTT-kss/2 0 0 -1i*H0TT(3,7) 0 0 -1i*H0TT(3,9) 0;...% QS(T+T)
       0 0 0 0 0 0 0 0   -1i*H0TT(7,3) 0 0 -1i*(H0TT(7,7)-H0TT(3,3))-kt/2-kJTT 1i*H0TT(7,3) 1i*H0TT(9,3) -1i*H0TT(7,9) 0 0 0;...  % TQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(7,3) 0 1i*H0TT(3,7) -kt-kJTT 1i*H0TT(9,7) 0 -1i*H0TT(7,9) 0 0;...  % TT(T+T)
       0 0 0 0 0 0 0 0   0 0 -1i*H0TT(7,3) 1i*H0TT(3,9) 1i*H0TT(7,9) -1i*(H0TT(7,7)-H0TT(9,9))-kt/2-kss/2-kJTT 0 0 -1i*H0TT(7,9) 0;...  % TS(T+T)
       0 0 kJT 0 0 0 0 0 -1i*H0TT(9,3) 0 0 -1i*H0TT(9,7) 0 0 -1i*(H0TT(9,9)-H0TT(3,3))-kJTT-kss/2 1i*H0TT(7,3) 1i*H0TT(9,3) 0;... % SQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(9,3) 0 0 -1i*H0TT(9,7) 0 1i*H0TT(3,7) -1i*(H0TT(9,9)-H0TT(7,7))-kt/2-kss/2-kJTT 1i*H0TT(9,7) 0;...  % ST(T+T)
       0 0 0 kJT 0 0 0 0  0 0 -1i*H0TT(9,3) 0 0 -1i*H0TT(9,7) 1i*H0TT(3,9) 1i*H0TT(7,9) -kJTTS-kss 0;...%SS(T+T)
       
       0 0 0 0 0 0 0 0 0 0 0 0 kJTT 0 0 0 0 -kt];% 3TT0(largeJ) -> 3T0
   
   

    RR2=[-kJ-kJTQ 0 0 0 kJQ 0 0 0 kJTT 0 0 0 0 0 0 0 0 0;...% QQ(largeJ)
       0 -kJ-ks/2-1i*(lm0(2,2)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 0 0 kJTT 0 0 0 0 0 0;... %QS(largeJ)
       0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(2,2))-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 +kJTT 0 0;... %SQ(largeJ)
       0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 kJTTS 0;... %SS(largeJ)
       
       kJ 0 0 0 -kJQ 1i*H01Q(9,2) -1i*H01Q(2,9) 0 0 0 0 0 0 0 0 0 0 0;...  %QQ(smallJ)
       0 kJ 0 0 1i*H01Q(2,9) -1i*(H01Q(2,2)-H01Q(9,9))-kJQ 0 -1i*H01Q(2,9) 0 0 0 0 0 0 0 0 0 0;...  %QS(smallJ)
       0 0 kJ 0 -1i*H01Q(9,2) 0 -1i*(H01Q(9,9)-H01Q(2,2))-kJQ 1i*H01Q(9,2) 0 0 0 0 0 0 0 0 0 0;...  %SQ(smallJ)
       0 0 0 kJ 0 -1i*H01Q(9,2) 1i*H01Q(2,9) -kJQ 0 0 0 0 0 0 0 0 0 0;...  % SS(smallJ)
       
       kJTQ 0 0 0 0 0 0 0  -kJTT 1i*H0TT(6,2) 1i*H0TT(9,2) -1i*H0TT(2,6) 0 0 -1i*H0TT(2,9) 0 0 0;...  % QQ(T+T)
       0 0 0 0 0 0 0 0   1i*H0TT(2,6) -1i*(H0TT(2,2)-H0TT(6,6))-kt/2-kJTT 1i*H0TT(9,6) 0 -1i*H0TT(2,6) 0 0 -1i*H0TT(2,9) 0 0;...  % QT(T+T)
       0 kJT 0 0 0 0 0 0  1i*H0TT(2,9) 1i*H0TT(6,9) -1i*(H0TT(2,2)-H0TT(9,9))-kJTT-kss/2 0 0 -1i*H0TT(2,6) 0 0 -1i*H0TT(2,9) 0;...% QS(T+T)
       0 0 0 0 0 0 0 0   -1i*H0TT(6,2) 0 0 -1i*(H0TT(6,6)-H0TT(2,2))-kt/2-kJTT 1i*H0TT(6,2) 1i*H0TT(9,2) -1i*H0TT(6,9) 0 0 0;...  % TQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(6,2) 0 1i*H0TT(2,6) -kt-kJTT 1i*H0TT(9,6) 0 -1i*H0TT(6,9) 0 0;...  % TT(T+T)
       0 0 0 0 0 0 0 0   0 0 -1i*H0TT(6,2) 1i*H0TT(2,9) 1i*H0TT(6,9) -1i*(H0TT(6,6)-H0TT(9,9))-kt/2-kss/2-kJTT 0 0 -1i*H0TT(6,9) 0;...  % TS(T+T)
       0 0 kJT 0 0 0 0 0 -1i*H0TT(9,2) 0 0 -1i*H0TT(9,6) 0 0 -1i*(H0TT(9,9)-H0TT(2,2))-kJTT-kss/2 1i*H0TT(6,2) 1i*H0TT(9,2) 0;... % SQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(9,2) 0 0 -1i*H0TT(9,6) 0 1i*H0TT(2,6) -1i*(H0TT(9,9)-H0TT(6,6))-kt/2-kss/2-kJTT 1i*H0TT(9,6) 0;...  % ST(T+T)
       0 0 0 kJT 0 0 0 0  0 0 -1i*H0TT(9,2) 0 0 -1i*H0TT(9,6) 1i*H0TT(2,9) 1i*H0TT(6,9) -kJTTS-kss 0;...%SS(T+T)
       
       0 0 0 0 0 0 0 0 0 0 0 0 kJTT 0 0 0 0 -kt];% 3TT(largeJ) -> 3T+
%    
     RR4=[-kJ-kJTQ 0 0 0 kJQ 0 0 0 kJTT 0 0 0 0 0 0 0 0 0;... % QQ(largeJ)
       0 -kJ-ks/2-1i*(lm0(4,4)-lm0(9,9))-kJT 0 0 0 kJQ 0 0 0 0 0 kJTT 0 0 0 0 0 0;... %QS(largeJ)
       0 0 -kJ-ks/2-1i*(lm0(9,9)-lm0(4,4))-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 +kJTT 0 0;... %SQ(largeJ)
       0 0 0 -kJ-ks-kJT 0 0 0 kJQ 0 0 0 0 0 0 0 0 kJTTS 0;... %SS(largeJ)
       
       kJ 0 0 0 -kJQ 1i*H01Q(9,4) -1i*H01Q(4,9) 0 0 0 0 0 0 0 0 0 0 0;...  %QQ(smallJ)
       0 kJ 0 0 1i*H01Q(4,9) -1i*(H01Q(4,4)-H01Q(9,9))-kJQ 0 -1i*H01Q(4,9) 0 0 0 0 0 0 0 0 0 0;...  %QS(smallJ)
       0 0 kJ 0 -1i*H01Q(9,4) 0 -1i*(H01Q(9,9)-H01Q(4,4))-kJQ 1i*H01Q(9,4) 0 0 0 0 0 0 0 0 0 0;...  %SQ(smallJ)
       0 0 0 kJ 0 -1i*H01Q(9,4) 1i*H01Q(4,9) -kJQ 0 0 0 0 0 0 0 0 0 0;...  % SS(smallJ)
       
       kJTQ 0 0 0 0 0 0 0  -kJTT 1i*H0TT(8,4) 1i*H0TT(9,4) -1i*H0TT(4,8) 0 0 -1i*H0TT(4,9) 0 0 0;...  % QQ(T+T)
       0 0 0 0 0 0 0 0   1i*H0TT(4,8) -1i*(H0TT(4,4)-H0TT(8,8))-kt/2-kJTT 1i*H0TT(9,8) 0 -1i*H0TT(4,8) 0 0 -1i*H0TT(4,9) 0 0;...  % QT(T+T)
       0 kJT 0 0 0 0 0 0  1i*H0TT(4,9) 1i*H0TT(8,9) -1i*(H0TT(4,4)-H0TT(9,9))-kJTT-kss/2 0 0 -1i*H0TT(4,8) 0 0 -1i*H0TT(4,9) 0;...% QS(T+T)
       0 0 0 0 0 0 0 0   -1i*H0TT(8,4) 0 0 -1i*(H0TT(8,8)-H0TT(4,4))-kt/2-kJTT 1i*H0TT(8,4) 1i*H0TT(9,4) -1i*H0TT(8,9) 0 0 0;...  % TQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(8,4) 0 1i*H0TT(4,8) -kt-kJTT 1i*H0TT(9,8) 0 -1i*H0TT(8,9) 0 0;...  % TT(T+T)
       0 0 0 0 0 0 0 0   0 0 -1i*H0TT(8,4) 1i*H0TT(4,9) 1i*H0TT(8,9) -1i*(H0TT(8,8)-H0TT(9,9))-kt/2-kss/2-kJTT 0 0 -1i*H0TT(8,9) 0;...  % TS(T+T)
       0 0 kJT 0 0 0 0 0 -1i*H0TT(9,4) 0 0 -1i*H0TT(9,8) 0 0 -1i*(H0TT(9,9)-H0TT(4,4))-kJTT-kss/2 1i*H0TT(8,4) 1i*H0TT(9,4) 0;... % SQ(T+T)
       0 0 0 0 0 0 0 0   0 -1i*H0TT(9,4) 0 0 -1i*H0TT(9,8) 0 1i*H0TT(4,8) -1i*(H0TT(9,9)-H0TT(8,8))-kt/2-kss/2-kJTT 1i*H0TT(9,8) 0;...  % ST(T+T)
       0 0 0 kJT 0 0 0 0  0 0 -1i*H0TT(9,4) 0 0 -1i*H0TT(9,8) 1i*H0TT(4,9) 1i*H0TT(8,9) -kJTTS-kss 0;...%SS(T+T)
       
       0 0 0 0 0 0 0 0 0 0 0 0 kJTT 0 0 0 0 -kt];% 3TT-(largeJ) -> 3T- 
            
   pIIII(1,1)=0;pIIII(2,1)=0;pIIII(3,1)=0;pIIII(4,1)=1;pIIII(5,1)=0; pIIII(6,1)=0; pIIII(7,1)=0; pIIII(8,1)=0; pIIII(9,1)=0; pIIII(10,1)=0; pIIII(11,1)=0; pIIII(12,1)=0; 
    pIIII(13,1)=0; pIIII(14,1)=0; pIIII(15,1)=0; pIIII(16,1)=0; pIIII(17,1)=0; pIIII(18,1)=0; 
  
            [VPP2,DPP2] = eig(RR2);           
            [VPP3,DPP3] = eig(RR3);   
            [VPP4,DPP4] = eig(RR4);
            
            DDEXP2(1,1)=exp(DPP2(1,1)*time); DDEXP2(2,2)=exp(DPP2(2,2)*time); DDEXP2(3,3)=exp(DPP2(3,3)*time); 
            DDEXP2(4,4)=exp(DPP2(4,4)*time); DDEXP2(5,5)=exp(DPP2(5,5)*time); DDEXP2(6,6)=exp(DPP2(6,6)*time); 
            DDEXP2(7,7)=exp(DPP2(7,7)*time); DDEXP2(8,8)=exp(DPP2(8,8)*time); DDEXP2(9,9)=exp(DPP2(9,9)*time); 
            DDEXP2(10,10)=exp(DPP2(10,10)*time); 
            DDEXP2(11,11)=exp(DPP2(11,11)*time); DDEXP2(12,12)=exp(DPP2(12,12)*time); 
            DDEXP2(13,13)=exp(DPP2(13,13)*time); DDEXP2(14,14)=exp(DPP2(14,14)*time); DDEXP2(15,15)=exp(DPP2(15,15)*time); 
            DDEXP2(16,16)=exp(DPP2(16,16)*time); 
            DDEXP2(17,17)=exp(DPP2(17,17)*time); DDEXP2(18,18)=exp(DPP2(18,18)*time); 
%             
            DDEXP3(1,1)=exp(DPP3(1,1)*time); DDEXP3(2,2)=exp(DPP3(2,2)*time); DDEXP3(3,3)=exp(DPP3(3,3)*time); 
            DDEXP3(4,4)=exp(DPP3(4,4)*time); DDEXP3(5,5)=exp(DPP3(5,5)*time); DDEXP3(6,6)=exp(DPP3(6,6)*time); 
            DDEXP3(7,7)=exp(DPP3(7,7)*time); DDEXP3(8,8)=exp(DPP3(8,8)*time); DDEXP3(9,9)=exp(DPP3(9,9)*time); 
            DDEXP3(10,10)=exp(DPP3(10,10)*time); 
            DDEXP3(11,11)=exp(DPP3(11,11)*time); DDEXP3(12,12)=exp(DPP3(12,12)*time); 
            DDEXP3(13,13)=exp(DPP3(13,13)*time); DDEXP3(14,14)=exp(DPP3(14,14)*time); DDEXP3(15,15)=exp(DPP3(15,15)*time); 
            DDEXP3(16,16)=exp(DPP3(16,16)*time); 
            DDEXP3(17,17)=exp(DPP3(17,17)*time);  DDEXP3(18,18)=exp(DPP3(18,18)*time); 
            
 
            DDEXP4(1,1)=exp(DPP4(1,1)*time); DDEXP4(2,2)=exp(DPP4(2,2)*time); DDEXP4(3,3)=exp(DPP4(3,3)*time); 
            DDEXP4(4,4)=exp(DPP4(4,4)*time); DDEXP4(5,5)=exp(DPP4(5,5)*time); DDEXP4(6,6)=exp(DPP4(6,6)*time); 
            DDEXP4(7,7)=exp(DPP4(7,7)*time); DDEXP4(8,8)=exp(DPP4(8,8)*time); DDEXP4(9,9)=exp(DPP4(9,9)*time); 
            DDEXP4(10,10)=exp(DPP4(10,10)*time); 
            DDEXP4(11,11)=exp(DPP4(11,11)*time); DDEXP4(12,12)=exp(DPP4(12,12)*time); 
            DDEXP4(13,13)=exp(DPP4(13,13)*time); DDEXP4(14,14)=exp(DPP4(14,14)*time); DDEXP4(15,15)=exp(DPP4(15,15)*time); 
            DDEXP4(16,16)=exp(DPP4(16,16)*time); 
            DDEXP4(17,17)=exp(DPP4(17,17)*time);  DDEXP4(18,18)=exp(DPP4(18,18)*time); 

  
            
            %%%%% Density matrix at time of t
            %%%% for [p++,pSS,pS0,p0S,p00,p--,pS1S1]
            pIItt1=VPP1*DDEXP1/(VPP1)*pIII;
            pIItt2=VPP2*DDEXP2/(VPP2)*pIIII;
            pIItt3=VPP3*DDEXP3/(VPP3)*pIIII;
            pIItt4=VPP4*DDEXP4/(VPP4)*pIIII;
            pIItt5=VPP5*DDEXP5/(VPP5)*pIII;
            
             pIIdisTT(1,1)=pIItt1(9,1); %roh_Q2_Q2
             
             pIIdisTT(2,2)=pIItt2(9,1); %pIIdisTT(2,9)=0;%pIItt2(11,1); 
%              pIIdisTT(2,9)=pIItt2(10,1); pIIdisTT(9,2)=pIItt2(11,1);
%              %pIIdisTT(9,2)=0;%pIItt2(15,1);% p_Q+1_S

            pIIdisTT(2,6)=pIItt2(10,1); 
            pIIdisTT(6,2)=pIItt2(12,1); % p_T1_Q1
            pIIdisTT(6,9)=pIItt2(14,1); 
            pIIdisTT(9,6)=pIItt2(16,1); % p_T+_S
            pIIdisTT(6,6)=pIItt2(13,1); %T+T+           

            pIIdisTT(3,3)=pIItt3(9,1); pIIdisTT(3,9)=pIItt3(11,1); pIIdisTT(9,3)=pIItt3(15,1);% p_Q0_S
            pIIdisTT(3,7)=pIItt3(10,1); pIIdisTT(7,3)=pIItt3(12,1); % p_T0_Q0
            pIIdisTT(7,9)=pIItt3(14,1); pIIdisTT(9,7)=pIItt3(16,1); % p_T0_S
            pIIdisTT(7,7)=pIItt3(13,1); %T0T0

             pIIdisTT(4,4)=pIItt4(9,1); %pIIdisTT(4,9)=0;%pIItt4(11,1); 
%              %pIIdisTT(9,4)=0;%pIItt4(15,1);% p_Q-1_S
%              %pIIdisTT(2,9)=0;%pIItt2(11,1); 
%              pIIdisTT(4,9)=pIItt4(10,1); pIIdisTT(9,4)=pIItt4(11,1);

           pIIdisTT(4,8)=pIItt4(10,1); 
           pIIdisTT(8,4)=pIItt4(12,1); % p_T-_Q-1
            pIIdisTT(8,9)=pIItt4(14,1); 
            pIIdisTT(9,8)=pIItt4(16,1); % p_T-_S
            pIIdisTT(8,8)=pIItt4(13,1); %T-T-           

            pIIdisTT(5,5)=pIItt5(9,1); pIIdisTT(5,9)=pIItt5(10,1); pIIdisTT(9,5)=pIItt5(11,1);% p_Q-2_S

           pIIdisTT(9,9)=(pIItt1(12,1)+pIItt2(17,1)+pIItt3(17,1)+pIItt4(17,1)+pIItt5(12,1))/5;%  pIItt4(12,1) --> pIItt4(17,1) ?
            
            
            % nine sublevel populations of weakly-coupled TT

pIIdisTT2=QQb'*pIIdisTT*QQb;

%% Depopulations by TTA via S-Q-T mixings
% pIIdisTT2(4,4)=pIIdisTT2(4,4)*exp(-(kss)*QQb(9,4)^2*time);
% pIIdisTT2(5,5)=pIIdisTT2(5,5)*exp(-(kss)*QQb(9,5)^2*time);
% pIIdisTT2(6,6)=pIIdisTT2(6,6)*exp(-(kss)*QQb(9,6)^2*time);
% pIIdisTT2(1,1)=pIIdisTT2(1,1)*exp(-(kt)*QQb(1,1)^2*time);
% pIIdisTT2(2,2)=pIIdisTT2(2,2)*exp(-(kt)*QQb(2,2)^2*time);
% pIIdisTT2(3,3)=pIIdisTT2(3,3)*exp(-(kt)*QQb(9,3)^2*time);
% pIIdisTT2(7,7)=pIIdisTT2(7,7)*exp(-(kt)*QQb(9,7)^2*time);
% pIIdisTT2(8,8)=pIIdisTT2(8,8)*exp(-(kt)*QQb(9,8)^2*time);
% pIIdisTT2(9,9)=pIIdisTT2(9,9)*exp(-(kt)*QQb(9,9)^2*time);

            


% roh00(1,1)=ZFSmag(1,5)^2*tau/(1+((2*Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);
% roh00(2,2)=ZFSmag(1,6)^2*tau/(1+((Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);
% roh00(3,3)=ZFSmag(1,7)^2*tau/(1+((-6*J2)*tau)^2);
% roh00(4,4)=ZFSmag(1,8)^2*tau/(1+((-Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);
% roh00(5,5)=ZFSmag(1,9)^2*tau/(1+((-2*Bav*2*pi*28.0249*1e6-6*J2)*tau)^2);
%       
% Field sweep
            for mmm=1:ll
  		      B0=Bs+mmm*dB;
              B0= B0*2*pi*28.0249*1e6;%(rad/s)
  
            
% Spin Operator (Sx, Sy, Sz) on the basis of ZFS Triplet eigenfunctions (TX,TY,TZ)
% Sx(1,1)=0; Sx(1,2)=-i*sin(h)*sin(q); Sx(1,3)=i*(-sin(y)*cos(q)-cos(h)*sin(q)*cos(y));
% Sx(2,1)=i*sin(h)*sin(q); Sx(2,2)=0; Sx(2,3)=-i*(cos(y)*cos(q)-cos(h)*sin(q)*sin(y));
% Sx(3,1)=-i*(-sin(y)*cos(q)-cos(h)*sin(q)*cos(y)); Sx(3,2)=i*(cos(y)*cos(q)-cos(h)*sin(q)*sin(y)); Sx(3,3)=0;
% 
% Sy(1,1)=0; Sy(1,2)=i*sin(h)*cos(q); Sy(1,3)=i*(-sin(y)*sin(q)+cos(h)*cos(q)*cos(y));
% Sy(2,1)=-i*sin(h)*cos(h); Sy(2,2)=0; Sy(2,3)=-i*(cos(y)*sin(q)+cos(h)*cos(q)*sin(y));
% Sy(3,1)=-i*(-sin(y)*sin(q)+cos(h)*cos(q)*cos(y)); Sy(3,2)=i*(cos(y)*sin(q)+cos(h)*cos(q)*sin(y)); Sy(3,3)=0;

% Sz(1,1)=0; Sz(1,2)=-i*cos(h); Sz(1,3)=i*cos(y)*sin(h);
% Sz(2,1)=i*cos(h); Sz(2,2)=0; Sz(2,3)=-i*sin(y)*sin(h);
% Sz(3,1)=-i*cos(y)*sin(h); Sz(3,2)=i*sin(y)*sin(h); Sz(3,3)=0;


%%%%%%%%% Eigen energies determined by Spin Hamiltonian %%%%%%%%%%%%%%
HTTZ = B0*(S1Z+S2Z);
H011=HTTZ+HTTZFS+HTTEX+HTTZFST1;% quintet
H012=HTTZ+HTTZFS22+HTTEX2+HTTZFST1;% transition quintet
H04=HTTZ+HTTZFS44+HTTEX3+HTTZFST2;% weakly-coupled TT
H0111=UU1*H011*UU1';% quintet
H0112=UU1*H012*UU1';% transition quintet
H00TT=UU1*H04*UU1';% weakly-coupled TT
%H02=HTTZ+HTTZFS;

[QQ,lm]=eig2(H0111);  % quintet
[QQb,lmb]=eig2(H0112);  % transition quintet
[QQ2,lm2]=eig2(H00TT);  % weakly-coupled TT

% Transition probability

%TRANS=QQ2'*UU1*(S1Y+S2Y)*UU1'*QQ2;

 HTZ = B0*(STZ);
 HT1=HTZ+HTZFS;
 
 [QQT,lmT]=eig2(HT1);  % indiviidual triplet


%  H0 = B0*Sz + HD;
%  
%  [Q,lm]=eig(H0);
%s(H0,3,1e12);  %Diagonalization of Hamiltonian



%roh0=Q'*roh*Q;  % Initial density matrix of each eigenstate T+, T0, T-
%C=D/6.0000*(1.0000-3.0000*cos(h)^2)-E/2.0000*(cos(y)^2*sin(h)^2-sin(y)^2*sin(h)^2);

%   CCCC=D/6.000*(1.000-3.000*eee(3,1)^2)-E/2.000*(eee(1,1)^2-eee(2,1)^2);
 
%lm(3,3)=B0+C;lm(2,2)=-2*C;lm(1,1)=-B0+C;

%CC=D/6*(1-3*(cos(h))^2)-E/2*((sin(h))^2*(cos(y))^2-(sin(h))^2*(sin(y))^2);   delta=B0-u0;

%%% Lineshape function

%                            roh21=(roh0(1,1)-roh0(2,2))*(-delta+2*CC+1i/T2)/(1/T2^2+(delta-2*CC)^2 );
%                            roh32=(roh0(2,2)-roh0(3,3))*(-delta-2*CC+1i/T2)/(1/T2^2+(delta+2*CC)^2 );

 %Å@Qiuntet state matrix elements in the liouvlle equations 1,1=site A; 2,2=site B  u0-lm(1,1)+lm(2,2)Å@12-transition
 A12(1,1)=-1i*(-u0+lm(1,1)-lm(2,2))-(kJ+1/T22); A12(1,2)=kJQ;
 A12(2,2)=-1i*(-u0+lmb(1,1)-lmb(2,2))-(kJQ+1/T22); A12(2,1)=kJ;
 B12(1,1)=1i*(pIItt1(1,1)-pIItt2(1,1)); B12(2,1)=1i*(pIItt1(5,1)-pIItt2(5,1));%Q1-Q2(A cite)  Q1-Q2(B cite)
 % solution of the denstity matrix p12
 C12=(A12)\B12;
  
 % for 23 transitions 
% Qiuntet state matrix elements in the liouvlle equations 1,1=site A; 2,2=site B 
 A23(1,1)=-1i*(-u0+lm(2,2)-lm(3,3))-(kJ+1/T22); A23(1,2)=kJQ;
 A23(2,2)=-1i*(-u0+lmb(2,2)-lmb(3,3))-(kJQ+1/T22); A23(2,1)=kJ;
 B23(1,1)=1i*(pIItt2(1,1)-pIItt3(1,1)); B23(2,1)=1i*(pIItt2(5,1)-pIItt3(5,1));%Q2-Q3(A cite)  Q2-Q3(B cite)
 % solution of the denstity matrix p12
 C23=(A23)\B23;
 
% Qiuntet state matrix elements in the liouvlle equations 1,1=site A; 2,2=site B 
 A34(1,1)=-1i*(-u0+lm(3,3)-lm(4,4))-(kJ+1/T22); A34(1,2)=kJQ;
 A34(2,2)=-1i*(-u0+lmb(3,3)-lmb(4,4))-(kJQ+1/T22); A34(2,1)=kJ;
 B34(1,1)=1i*(pIItt3(1,1)-pIItt4(1,1)); B34(2,1)=1i*(pIItt3(5,1)-pIItt4(5,1));
 % solution of the denstity matrix p12
 C34=(A34)\B34;

 % Qiuntet state matrix elements in the liouvlle equations 1,1=site A; 2,2=site B 
 A45(1,1)=-1i*(-u0+lm(4,4)-lm(5,5))-(kJ+1/T22); A45(1,2)=kJQ;
 A45(2,2)=-1i*(-u0+lmb(4,4)-lmb(5,5))-(kJQ+1/T22); A45(2,1)=kJ;
 B45(1,1)=1i*(pIItt4(1,1)-pIItt5(1,1)); B45(2,1)=1i*(pIItt4(5,1)-pIItt5(5,1));
 % solution of the denstity matrix p12
 C45=(A45)\B45;

                        % Triplet State via 3(TT) -> T1S0                            
%                             rohTTT(2,1)=-(pIItt2(18,1)-pIItt3(18,1))/(u0-lmT(1,1)+lmT(2,2)+1i/T22);
%                             rohTTT(3,2)=-(pIItt3(18,1)-pIItt4(18,1))/(u0-lmT(2,2)+lmT(3,3)+1i/T22);
                             
                             % Strongly Coupled Triplet TT
 
                             rohTTT(2,1)=(pIItt2(18,1)-pIItt3(18,1))/(u0-lm(6,6)+lm(7,7)+1i/T22);
                             rohTTT(3,2)=(pIItt3(18,1)-pIItt4(18,1))/(u0-lm(7,7)+lm(8,8)+1i/T22);

%                           weakly-coupled T+T EPR transitions
                            rohTTij(2,1)=-(pIIdisTT2(1,1)-pIIdisTT2(2,2))/(u0-lm2(1,1)+lm2(2,2)+1i/T2);
                            rohTTij(3,1)=-(pIIdisTT2(1,1)-pIIdisTT2(3,3))/(u0-lm2(1,1)+lm2(3,3)+1i/T2);
                            
                            rohTTij(4,2)=-(pIIdisTT2(2,2)-pIIdisTT2(4,4))/(u0-lm2(2,2)+lm2(4,4)+1i/T2);
                            rohTTij(5,2)=-(pIIdisTT2(2,2)-pIIdisTT2(5,5))/(u0-lm2(2,2)+lm2(5,5)+1i/T2);
                            rohTTij(6,2)=-(pIIdisTT2(2,2)-pIIdisTT2(6,6))/(u0-lm2(2,2)+lm2(6,6)+1i/T2);
                            rohTTij(4,3)=-(pIIdisTT2(3,3)-pIIdisTT2(4,4))/(u0-lm2(3,3)+lm2(4,4)+1i/T2);
                            rohTTij(5,3)=-(pIIdisTT2(3,3)-pIIdisTT2(5,5))/(u0-lm2(3,3)+lm2(5,5)+1i/T2);
                            rohTTij(6,3)=-(pIIdisTT2(3,3)-pIIdisTT2(6,6))/(u0-lm2(3,3)+lm2(6,6)+1i/T2);
    
                            rohTTij(7,4)=-(pIIdisTT2(4,4)-pIIdisTT2(7,7))/(u0-lm2(4,4)+lm2(7,7)+1i/T2);
                            rohTTij(8,4)=-(pIIdisTT2(4,4)-pIIdisTT2(8,8))/(u0-lm2(4,4)+lm2(8,8)+1i/T2);
                            rohTTij(7,5)=-(pIIdisTT2(5,5)-pIIdisTT2(7,7))/(u0-lm2(5,5)+lm2(7,7)+1i/T2);
                            rohTTij(8,5)=-(pIIdisTT2(5,5)-pIIdisTT2(8,8))/(u0-lm2(5,5)+lm2(8,8)+1i/T2);
                            rohTTij(7,6)=-(pIIdisTT2(6,6)-pIIdisTT2(7,7))/(u0-lm2(6,6)+lm2(7,7)+1i/T2);
                            rohTTij(8,6)=-(pIIdisTT2(6,6)-pIIdisTT2(8,8))/(u0-lm2(6,6)+lm2(8,8)+1i/T2);
    
                            rohTTij(9,7)=-(pIIdisTT2(7,7)-pIIdisTT2(9,9))/(u0-lm2(7,7)+lm2(9,9)+1i/T2);
                            rohTTij(9,8)=-(pIIdisTT2(8,8)-pIIdisTT2(9,9))/(u0-lm2(8,8)+lm2(9,9)+1i/T2);
                            
                            Iint(2,1)=QQ2(1,1)^2*QQ2(2,2)^2*4 + QQ2(2,1)^2*QQ2(1,2)^2*4 + QQ2(2,1)^2*QQ2(3,2)^2*6+...
                                      QQ2(3,1)^2*QQ2(2,2)^2*6 + QQ2(3,1)^2*QQ2(4,2)^2*6 + QQ2(4,1)^2*QQ2(3,2)^2*6 +...
                                      QQ2(4,1)^2*QQ2(5,2)^2*4 + QQ2(5,1)^2*QQ2(4,2)^2*4 + QQ2(6,1)^2*QQ2(7,2)^2*2+...
                                      QQ2(7,1)^2*QQ2(6,2)^2*2 + QQ2(7,1)^2*QQ2(8,2)^2*2 + QQ2(8,1)^2*QQ2(7,2)^2*2;
                            Iint(3,1)=QQ2(1,1)^2*QQ2(2,3)^2*4 + QQ2(2,1)^2*QQ2(1,3)^2*4 + QQ2(2,1)^2*QQ2(3,3)^2*6+...
                                      QQ2(3,1)^2*QQ2(2,3)^2*6 + QQ2(3,1)^2*QQ2(4,3)^2*6 + QQ2(4,1)^2*QQ2(3,3)^2*6 +...
                                      QQ2(4,1)^2*QQ2(5,3)^2*4 + QQ2(5,1)^2*QQ2(4,3)^2*4 + QQ2(6,1)^2*QQ2(7,3)^2*2+...
                                      QQ2(7,1)^2*QQ2(6,3)^2*2 + QQ2(7,1)^2*QQ2(8,3)^2*2 + QQ2(8,1)^2*QQ2(7,3)^2*2;
                            Iint(4,2)=QQ2(1,2)^2*QQ2(2,4)^2*4 + QQ2(2,2)^2*QQ2(1,4)^2*4 + QQ2(2,2)^2*QQ2(3,4)^2*6+...
                                      QQ2(3,2)^2*QQ2(2,4)^2*6 + QQ2(3,2)^2*QQ2(4,4)^2*6 + QQ2(4,2)^2*QQ2(3,4)^2*6 +...
                                      QQ2(4,2)^2*QQ2(5,4)^2*4 + QQ2(5,2)^2*QQ2(4,4)^2*4 + QQ2(6,2)^2*QQ2(7,4)^2*2+...
                                      QQ2(7,2)^2*QQ2(6,4)^2*2 + QQ2(7,2)^2*QQ2(8,4)^2*2 + QQ2(8,2)^2*QQ2(7,4)^2*2;
                            Iint(5,2)=QQ2(1,2)^2*QQ2(2,5)^2*4 + QQ2(2,2)^2*QQ2(1,5)^2*4 + QQ2(2,2)^2*QQ2(3,5)^2*6+...
                                      QQ2(3,2)^2*QQ2(2,5)^2*6 + QQ2(3,2)^2*QQ2(4,5)^2*6 + QQ2(4,2)^2*QQ2(3,5)^2*6 +...
                                      QQ2(4,2)^2*QQ2(5,5)^2*4 + QQ2(5,2)^2*QQ2(4,5)^2*4 + QQ2(6,2)^2*QQ2(7,5)^2*2+...
                                      QQ2(7,2)^2*QQ2(6,5)^2*2 + QQ2(7,2)^2*QQ2(8,5)^2*2 + QQ2(8,2)^2*QQ2(7,5)^2*2;
                            Iint(6,2)=QQ2(1,2)^2*QQ2(2,6)^2*4 + QQ2(2,2)^2*QQ2(1,6)^2*4 + QQ2(2,2)^2*QQ2(3,6)^2*6+...
                                      QQ2(3,2)^2*QQ2(2,6)^2*6 + QQ2(3,2)^2*QQ2(4,6)^2*6 + QQ2(4,2)^2*QQ2(3,6)^2*6 +...
                                      QQ2(4,2)^2*QQ2(5,6)^2*4 + QQ2(5,2)^2*QQ2(4,6)^2*4 + QQ2(6,2)^2*QQ2(7,6)^2*2+...
                                      QQ2(7,2)^2*QQ2(6,6)^2*2 + QQ2(7,2)^2*QQ2(8,6)^2*2 + QQ2(8,2)^2*QQ2(7,6)^2*2;
                            Iint(4,3)=QQ2(1,3)^2*QQ2(2,4)^2*4 + QQ2(2,3)^2*QQ2(1,4)^2*4 + QQ2(2,3)^2*QQ2(3,4)^2*6+...
                                      QQ2(3,3)^2*QQ2(2,4)^2*6 + QQ2(3,3)^2*QQ2(4,4)^2*6 + QQ2(4,3)^2*QQ2(3,4)^2*6 +...
                                      QQ2(4,3)^2*QQ2(5,4)^2*4 + QQ2(5,3)^2*QQ2(4,4)^2*4 + QQ2(6,3)^2*QQ2(7,4)^2*2+...
                                      QQ2(7,3)^2*QQ2(6,4)^2*2 + QQ2(7,3)^2*QQ2(8,4)^2*2 + QQ2(8,3)^2*QQ2(7,4)^2*2;
                            Iint(5,3)=QQ2(1,3)^2*QQ2(2,5)^2*4 + QQ2(2,3)^2*QQ2(1,5)^2*4 + QQ2(2,3)^2*QQ2(3,5)^2*6+...
                                      QQ2(3,3)^2*QQ2(2,5)^2*6 + QQ2(3,3)^2*QQ2(4,5)^2*6 + QQ2(4,3)^2*QQ2(3,5)^2*6 +...
                                      QQ2(4,3)^2*QQ2(5,5)^2*4 + QQ2(5,3)^2*QQ2(4,5)^2*4 + QQ2(6,3)^2*QQ2(7,5)^2*2+...
                                      QQ2(7,3)^2*QQ2(6,5)^2*2 + QQ2(7,3)^2*QQ2(8,5)^2*2 + QQ2(8,3)^2*QQ2(7,5)^2*2;
                            Iint(6,3)=QQ2(1,3)^2*QQ2(2,6)^2*4 + QQ2(2,3)^2*QQ2(1,6)^2*4 + QQ2(2,3)^2*QQ2(3,6)^2*6+...
                                      QQ2(3,3)^2*QQ2(2,6)^2*6 + QQ2(3,3)^2*QQ2(4,6)^2*6 + QQ2(4,3)^2*QQ2(3,6)^2*6 +...
                                      QQ2(4,3)^2*QQ2(5,6)^2*4 + QQ2(5,3)^2*QQ2(4,6)^2*4 + QQ2(6,3)^2*QQ2(7,6)^2*2+...
                                      QQ2(7,3)^2*QQ2(6,6)^2*2 + QQ2(7,3)^2*QQ2(8,6)^2*2 + QQ2(8,3)^2*QQ2(7,6)^2*2;
                            Iint(7,4)=QQ2(1,4)^2*QQ2(2,7)^2*4 + QQ2(2,4)^2*QQ2(1,7)^2*4 + QQ2(2,4)^2*QQ2(3,7)^2*6+...
                                      QQ2(3,4)^2*QQ2(2,7)^2*6 + QQ2(3,4)^2*QQ2(4,7)^2*6 + QQ2(4,4)^2*QQ2(3,7)^2*6 +...
                                      QQ2(4,4)^2*QQ2(5,7)^2*4 + QQ2(5,4)^2*QQ2(4,7)^2*4 + QQ2(6,4)^2*QQ2(7,7)^2*2+...
                                      QQ2(7,4)^2*QQ2(6,7)^2*2 + QQ2(7,4)^2*QQ2(8,7)^2*2 + QQ2(8,4)^2*QQ2(7,7)^2*2;
                            Iint(8,4)=QQ2(1,4)^2*QQ2(2,8)^2*4 + QQ2(2,4)^2*QQ2(1,8)^2*4 + QQ2(2,4)^2*QQ2(3,8)^2*6+...
                                      QQ2(3,4)^2*QQ2(2,8)^2*6 + QQ2(3,4)^2*QQ2(4,8)^2*6 + QQ2(4,4)^2*QQ2(3,8)^2*6 +...
                                      QQ2(4,4)^2*QQ2(5,8)^2*4 + QQ2(5,4)^2*QQ2(4,8)^2*4 + QQ2(6,4)^2*QQ2(7,8)^2*2+...
                                      QQ2(7,4)^2*QQ2(6,8)^2*2 + QQ2(7,4)^2*QQ2(8,8)^2*2 + QQ2(8,4)^2*QQ2(7,8)^2*2;
                            Iint(7,5)=QQ2(1,5)^2*QQ2(2,7)^2*4 + QQ2(2,5)^2*QQ2(1,7)^2*4 + QQ2(2,5)^2*QQ2(3,7)^2*6+...
                                      QQ2(3,5)^2*QQ2(2,7)^2*6 + QQ2(3,5)^2*QQ2(4,7)^2*6 + QQ2(4,5)^2*QQ2(3,7)^2*6 +...
                                      QQ2(4,5)^2*QQ2(5,7)^2*4 + QQ2(5,5)^2*QQ2(4,7)^2*4 + QQ2(6,5)^2*QQ2(7,7)^2*2+...
                                      QQ2(7,5)^2*QQ2(6,7)^2*2 + QQ2(7,5)^2*QQ2(8,7)^2*2 + QQ2(8,5)^2*QQ2(7,7)^2*2;
                            Iint(8,5)=QQ2(1,5)^2*QQ2(2,8)^2*4 + QQ2(2,5)^2*QQ2(1,8)^2*4 + QQ2(2,5)^2*QQ2(3,8)^2*6+...
                                      QQ2(3,5)^2*QQ2(2,8)^2*6 + QQ2(3,5)^2*QQ2(4,8)^2*6 + QQ2(4,5)^2*QQ2(3,8)^2*6 +...
                                      QQ2(4,5)^2*QQ2(5,8)^2*4 + QQ2(5,5)^2*QQ2(4,8)^2*4 + QQ2(6,5)^2*QQ2(7,8)^2*2+...
                                      QQ2(7,5)^2*QQ2(6,8)^2*2 + QQ2(7,5)^2*QQ2(8,8)^2*2 + QQ2(8,5)^2*QQ2(7,8)^2*2;
                            Iint(7,6)=QQ2(1,6)^2*QQ2(2,7)^2*4 + QQ2(2,6)^2*QQ2(1,7)^2*4 + QQ2(2,6)^2*QQ2(3,7)^2*6+...
                                      QQ2(3,6)^2*QQ2(2,7)^2*6 + QQ2(3,6)^2*QQ2(4,7)^2*6 + QQ2(4,6)^2*QQ2(3,7)^2*6 +...
                                      QQ2(4,6)^2*QQ2(5,7)^2*4 + QQ2(5,6)^2*QQ2(4,7)^2*4 + QQ2(6,6)^2*QQ2(7,7)^2*2+...
                                      QQ2(7,6)^2*QQ2(6,7)^2*2 + QQ2(7,6)^2*QQ2(8,7)^2*2 + QQ2(8,6)^2*QQ2(7,7)^2*2;
                            Iint(8,6)=QQ2(1,6)^2*QQ2(2,8)^2*4 + QQ2(2,6)^2*QQ2(1,8)^2*4 + QQ2(2,6)^2*QQ2(3,8)^2*6+...
                                      QQ2(3,6)^2*QQ2(2,8)^2*6 + QQ2(3,6)^2*QQ2(4,8)^2*6 + QQ2(4,6)^2*QQ2(3,8)^2*6 +...
                                      QQ2(4,6)^2*QQ2(5,8)^2*4 + QQ2(5,6)^2*QQ2(4,8)^2*4 + QQ2(6,6)^2*QQ2(7,8)^2*2+...
                                      QQ2(7,6)^2*QQ2(6,8)^2*2 + QQ2(7,6)^2*QQ2(8,8)^2*2 + QQ2(8,6)^2*QQ2(7,8)^2*2;
                            Iint(9,7)=QQ2(1,7)^2*QQ2(2,9)^2*4 + QQ2(2,7)^2*QQ2(1,9)^2*4 + QQ2(2,7)^2*QQ2(3,9)^2*6+...
                                      QQ2(3,7)^2*QQ2(2,9)^2*6 + QQ2(3,7)^2*QQ2(4,9)^2*6 + QQ2(4,7)^2*QQ2(3,9)^2*6 +...
                                      QQ2(4,7)^2*QQ2(5,9)^2*4 + QQ2(5,7)^2*QQ2(4,9)^2*4 + QQ2(6,7)^2*QQ2(7,9)^2*2+...
                                      QQ2(7,7)^2*QQ2(6,9)^2*2 + QQ2(7,7)^2*QQ2(8,9)^2*2 + QQ2(8,7)^2*QQ2(7,9)^2*2;
                            Iint(9,8)=QQ2(1,8)^2*QQ2(2,9)^2*4 + QQ2(2,8)^2*QQ2(1,9)^2*4 + QQ2(2,8)^2*QQ2(3,9)^2*6+...
                                      QQ2(3,8)^2*QQ2(2,9)^2*6 + QQ2(3,8)^2*QQ2(4,9)^2*6 + QQ2(4,8)^2*QQ2(3,9)^2*6 +...
                                      QQ2(4,8)^2*QQ2(5,9)^2*4 + QQ2(5,8)^2*QQ2(4,9)^2*4 + QQ2(6,8)^2*QQ2(7,9)^2*2+...
                                      QQ2(7,8)^2*QQ2(6,9)^2*2 + QQ2(7,8)^2*QQ2(8,9)^2*2 + QQ2(8,8)^2*QQ2(7,9)^2*2;
                         
                            
%                             Iint(2,1)=QQ2(1,1)*QQ2(2,2)*2+QQ2(2,1)*QQ2(1,2)*2 + QQ2(2,1)*QQ2(3,2)*2.45+QQ2(3,1)*QQ2(2,2)*2.45 + QQ2(3,1)*QQ2(4,2)*2.45+QQ2(4,1)*QQ2(3,2)*2.45 + QQ2(4,1)*QQ2(5,2)*2+QQ2(5,1)*QQ2(4,2)*2 + QQ2(6,1)*QQ2(7,2)*1.41+QQ2(7,1)*QQ2(6,2)*1.41 + QQ2(7,1)*QQ2(8,2)*1.41+QQ2(8,1)*QQ2(7,2)*1.41;
%                             Iint(3,1)=QQ2(1,1)*QQ2(2,3)*2+QQ2(2,1)*QQ2(1,3)*2 + QQ2(2,1)*QQ2(3,3)*2.45+QQ2(3,1)*QQ2(2,3)*2.45 + QQ2(3,1)*QQ2(4,3)*2.45+QQ2(4,1)*QQ2(3,3)*2.45 + QQ2(4,1)*QQ2(5,3)*2+QQ2(5,1)*QQ2(4,3)*2 + QQ2(6,1)*QQ2(7,3)*1.41+QQ2(7,1)*QQ2(6,3)*1.41 + QQ2(7,1)*QQ2(8,3)*1.41+QQ2(8,1)*QQ2(7,3)*1.41;
%                             Iint(4,2)=QQ2(1,2)*QQ2(2,4)*2+QQ2(2,2)*QQ2(1,4)*2 + QQ2(2,2)*QQ2(3,4)*2.45+QQ2(3,2)*QQ2(2,4)*2.45 + QQ2(3,2)*QQ2(4,4)*2.45+QQ2(4,2)*QQ2(3,4)*2.45 + QQ2(4,2)*QQ2(5,4)*2+QQ2(5,2)*QQ2(4,4)*2 + QQ2(6,2)*QQ2(7,4)*1.41+QQ2(7,2)*QQ2(6,4)*1.41 + QQ2(7,2)*QQ2(8,4)*1.41+QQ2(8,2)*QQ2(7,4)*1.41;
%                             Iint(5,2)=QQ2(1,2)*QQ2(2,5)*2+QQ2(2,2)*QQ2(1,5)*2 + QQ2(2,2)*QQ2(3,5)*2.45+QQ2(3,2)*QQ2(2,5)*2.45 + QQ2(3,2)*QQ2(4,5)*2.45+QQ2(4,2)*QQ2(3,5)*2.45 + QQ2(4,2)*QQ2(5,5)*2+QQ2(5,2)*QQ2(4,5)*2 + QQ2(6,2)*QQ2(7,5)*1.41+QQ2(7,2)*QQ2(6,5)*1.41 + QQ2(7,2)*QQ2(8,5)*1.41+QQ2(8,2)*QQ2(7,5)*1.41;
%                             Iint(6,2)=QQ2(1,2)*QQ2(2,6)*2+QQ2(2,2)*QQ2(1,6)*2 + QQ2(2,2)*QQ2(3,6)*2.45+QQ2(3,2)*QQ2(2,6)*2.45 + QQ2(3,2)*QQ2(4,6)*2.45+QQ2(4,2)*QQ2(3,6)*2.45 + QQ2(4,2)*QQ2(5,6)*2+QQ2(5,2)*QQ2(4,6)*2 + QQ2(6,2)*QQ2(7,6)*1.41+QQ2(7,2)*QQ2(6,6)*1.41 + QQ2(7,2)*QQ2(8,6)*1.41+QQ2(8,2)*QQ2(7,6)*1.41;
%                             Iint(4,3)=QQ2(1,3)*QQ2(2,4)*2+QQ2(2,3)*QQ2(1,4)*2 + QQ2(2,3)*QQ2(3,4)*2.45+QQ2(3,3)*QQ2(2,4)*2.45 + QQ2(3,3)*QQ2(4,4)*2.45+QQ2(4,3)*QQ2(3,4)*2.45 + QQ2(4,3)*QQ2(5,4)*2+QQ2(5,3)*QQ2(4,4)*2 + QQ2(6,3)*QQ2(7,4)*1.41+QQ2(7,3)*QQ2(6,4)*1.41 + QQ2(7,3)*QQ2(8,4)*1.41+QQ2(8,3)*QQ2(7,4)*1.41;
%                             Iint(5,3)=QQ2(1,3)*QQ2(2,5)*2+QQ2(2,3)*QQ2(1,5)*2 + QQ2(2,3)*QQ2(3,5)*2.45+QQ2(3,3)*QQ2(2,5)*2.45 + QQ2(3,3)*QQ2(4,5)*2.45+QQ2(4,3)*QQ2(3,5)*2.45 + QQ2(4,3)*QQ2(5,5)*2+QQ2(5,3)*QQ2(4,5)*2 + QQ2(6,3)*QQ2(7,5)*1.41+QQ2(7,3)*QQ2(6,5)*1.41 + QQ2(7,3)*QQ2(8,5)*1.41+QQ2(8,3)*QQ2(7,5)*1.41;
%                             Iint(6,3)=QQ2(1,3)*QQ2(2,6)*2+QQ2(2,3)*QQ2(1,6)*2 + QQ2(2,3)*QQ2(3,6)*2.45+QQ2(3,3)*QQ2(2,6)*2.45 + QQ2(3,3)*QQ2(4,6)*2.45+QQ2(4,3)*QQ2(3,6)*2.45 + QQ2(4,3)*QQ2(5,6)*2+QQ2(5,3)*QQ2(4,6)*2 + QQ2(6,3)*QQ2(7,6)*1.41+QQ2(7,3)*QQ2(6,6)*1.41 + QQ2(7,3)*QQ2(8,6)*1.41+QQ2(8,3)*QQ2(7,6)*1.41;
%                             Iint(7,4)=QQ2(1,4)*QQ2(2,7)*2+QQ2(2,4)*QQ2(1,7)*2 + QQ2(2,4)*QQ2(3,7)*2.45+QQ2(3,4)*QQ2(2,7)*2.45 + QQ2(3,4)*QQ2(4,7)*2.45+QQ2(4,4)*QQ2(3,7)*2.45 + QQ2(4,4)*QQ2(5,7)*2+QQ2(5,4)*QQ2(4,7)*2 + QQ2(6,4)*QQ2(7,7)*1.41+QQ2(7,4)*QQ2(6,7)*1.41 + QQ2(7,4)*QQ2(8,7)*1.41+QQ2(8,4)*QQ2(7,7)*1.41;
%                             Iint(8,4)=QQ2(1,4)*QQ2(2,8)*2+QQ2(2,4)*QQ2(1,8)*2 + QQ2(2,4)*QQ2(3,8)*2.45+QQ2(3,4)*QQ2(2,8)*2.45 + QQ2(3,4)*QQ2(4,8)*2.45+QQ2(4,4)*QQ2(3,8)*2.45 + QQ2(4,4)*QQ2(5,8)*2+QQ2(5,4)*QQ2(4,8)*2 + QQ2(6,4)*QQ2(7,8)*1.41+QQ2(7,4)*QQ2(6,8)*1.41 + QQ2(7,4)*QQ2(8,8)*1.41+QQ2(8,4)*QQ2(7,8)*1.41;
%                             Iint(7,5)=QQ2(1,5)*QQ2(2,7)*2+QQ2(2,5)*QQ2(1,7)*2 + QQ2(2,5)*QQ2(3,7)*2.45+QQ2(3,5)*QQ2(2,7)*2.45 + QQ2(3,5)*QQ2(4,7)*2.45+QQ2(4,5)*QQ2(3,7)*2.45 + QQ2(4,5)*QQ2(5,7)*2+QQ2(5,5)*QQ2(4,7)*2 + QQ2(6,5)*QQ2(7,7)*1.41+QQ2(7,5)*QQ2(6,7)*1.41 + QQ2(7,5)*QQ2(8,7)*1.41+QQ2(8,5)*QQ2(7,7)*1.41;
%                             Iint(8,5)=QQ2(1,5)*QQ2(2,8)*2+QQ2(2,5)*QQ2(1,8)*2 + QQ2(2,5)*QQ2(3,8)*2.45+QQ2(3,5)*QQ2(2,8)*2.45 + QQ2(3,5)*QQ2(4,8)*2.45+QQ2(4,5)*QQ2(3,8)*2.45 + QQ2(4,5)*QQ2(5,8)*2+QQ2(5,5)*QQ2(4,8)*2 + QQ2(6,5)*QQ2(7,8)*1.41+QQ2(7,5)*QQ2(6,8)*1.41 + QQ2(7,5)*QQ2(8,8)*1.41+QQ2(8,5)*QQ2(7,8)*1.41;
%                             Iint(7,6)=QQ2(1,6)*QQ2(2,7)*2+QQ2(2,6)*QQ2(1,7)*2 + QQ2(2,6)*QQ2(3,7)*2.45+QQ2(3,6)*QQ2(2,7)*2.45 + QQ2(3,6)*QQ2(4,7)*2.45+QQ2(4,6)*QQ2(3,7)*2.45 + QQ2(4,6)*QQ2(5,7)*2+QQ2(5,6)*QQ2(4,7)*2 + QQ2(6,6)*QQ2(7,7)*1.41+QQ2(7,5)*QQ2(6,6)*1.41 + QQ2(7,6)*QQ2(8,7)*1.41+QQ2(8,6)*QQ2(7,7)*1.41;
%                             Iint(8,6)=QQ2(1,6)*QQ2(2,8)*2+QQ2(2,6)*QQ2(1,8)*2 + QQ2(2,6)*QQ2(3,8)*2.45+QQ2(3,6)*QQ2(2,8)*2.45 + QQ2(3,6)*QQ2(4,8)*2.45+QQ2(4,6)*QQ2(3,8)*2.45 + QQ2(4,6)*QQ2(5,8)*2+QQ2(5,6)*QQ2(4,8)*2 + QQ2(6,6)*QQ2(7,8)*1.41+QQ2(7,6)*QQ2(6,8)*1.41 + QQ2(7,6)*QQ2(8,8)*1.41+QQ2(8,6)*QQ2(7,8)*1.41;
%                             Iint(9,7)=QQ2(1,7)*QQ2(2,9)*2+QQ2(2,7)*QQ2(1,9)*2 + QQ2(2,7)*QQ2(3,9)*2.45+QQ2(3,7)*QQ2(2,9)*2.45 + QQ2(3,7)*QQ2(4,9)*2.45+QQ2(4,7)*QQ2(3,9)*2.45 + QQ2(4,7)*QQ2(5,9)*2+QQ2(5,7)*QQ2(4,9)*2 + QQ2(6,7)*QQ2(7,9)*1.41+QQ2(7,7)*QQ2(6,9)*1.41 + QQ2(7,7)*QQ2(8,9)*1.41+QQ2(8,7)*QQ2(7,9)*1.41;
%                             Iint(9,8)=QQ2(1,8)*QQ2(2,9)*2+QQ2(2,8)*QQ2(1,9)*2 + QQ2(2,8)*QQ2(3,9)*2.45+QQ2(3,8)*QQ2(2,9)*2.45 + QQ2(3,7)*QQ2(4,9)*2.45+QQ2(4,7)*QQ2(3,9)*2.45 + QQ2(4,7)*QQ2(5,9)*2+QQ2(5,7)*QQ2(4,9)*2 + QQ2(6,7)*QQ2(7,9)*1.41+QQ2(7,7)*QQ2(6,9)*1.41 + QQ2(7,7)*QQ2(8,9)*1.41+QQ2(8,7)*QQ2(7,9)*1.41;
%                              
% (C12(1,1)+C23(1,1)+C12(2,1)+C23(2,1))                            
                            
              % PA=weakly coupled TT             
               PA=-(Iint(2,1)*rohTTij(2,1)+Iint(3,1)*rohTTij(3,1)+Iint(4,2)*rohTTij(4,2)+Iint(5,2)*rohTTij(5,2)+Iint(6,2)*rohTTij(6,2)+Iint(4,3)*rohTTij(4,3)+Iint(5,3)*rohTTij(5,3)+Iint(6,3)*rohTTij(6,3)+Iint(7,4)*rohTTij(7,4)+Iint(8,4)*rohTTij(8,4)+Iint(7,5)*rohTTij(7,5)+Iint(8,5)*rohTTij(8,5)+Iint(7,6)*rohTTij(7,6)+Iint(8,6)*rohTTij(8,6)+Iint(9,7)*rohTTij(9,7)+Iint(9,8)*rohTTij(9,8))*df*(sin((kk-1)*dc))*dc*exp(-kab*time);
%                PA=-(Iint(4,2)*rohTTij(4,2)*1+Iint(4,3)*rohTTij(4,3)*1+Iint(7,4)*rohTTij(7,4)*1+Iint(8,4)*rohTTij(8,4)*1)*1-...
%                (Iint(5,2)*rohTTij(5,2)*1+Iint(5,3)*rohTTij(5,3)*1+Iint(7,5)*rohTTij(7,5)*1+Iint(8,5)*rohTTij(8,5)*1)*1-...
%                (Iint(6,2)*rohTTij(6,2)*1+Iint(6,3)*rohTTij(6,3)*1+Iint(7,6)*rohTTij(7,6)*1+Iint(8,6)*rohTTij(8,6)*1)*1;
               PA=(-2*rohTTT(2,1)-2*rohTTT(3,2))*df*(sin((kk-1)*dc))*dc*exp(-kab*time);
               % PT=Quintet EPR with the A and B sites and T1-S0 state via 3TT -> T1-S0            
             PT=(4*(C12(1,1)+C12(2,1))+6*(C23(1,1)+C23(2,1))+6*(C34(1,1)+C34(2,1))+4*(C45(1,1)+C45(2,1)))*df*dc*(sin((kk-1)*dc));
            % PT=rohTT(2,1)+rohTT(3,2);
                
                               PP=[PP,PA];PPT=[PPT,PT];%PPPT=[,PTT];
             end
             
             PP3=PP3+PP; PP4=PP4+PPT;
             
pIItt1(4,1)=pIItt1(4,1)/5;pIItt1(8,1)=pIItt1(8,1)/5;pIItt1(12,1)=pIItt1(12,1)/5;
pIItt2(4,1)=pIItt2(4,1)/5;pIItt2(8,1)=pIItt2(8,1)/5;pIItt2(17,1)=pIItt2(17,1)/5;
pIItt3(4,1)=pIItt3(4,1)/5;pIItt3(8,1)=pIItt3(8,1)/5;pIItt3(17,1)=pIItt3(17,1)/5;pIItt3(18,1)=pIItt3(18,1)/2; % 3(TT)->T1:S0
pIItt4(4,1)=pIItt4(4,1)/5;pIItt4(8,1)=pIItt4(8,1)/5;pIItt4(17,1)=pIItt4(17,1)/5;
pIItt5(4,1)=pIItt5(4,1)/5;pIItt5(8,1)=pIItt5(8,1)/5;pIItt5(12,1)=pIItt5(12,1)/5;

pIItt_total1 = pIItt1;
pIItt_total2 = pIItt2;
pIItt_total3 = pIItt3;
pIItt_total4 = pIItt4;
pIItt_total5 = pIItt5;
                
TTY=pIItt1(1,1)+pIItt1(4,1)+pIItt1(5,1)+pIItt1(8,1) + pIItt2(1,1)+pIItt2(4,1)+pIItt2(5,1)+pIItt2(8,1) + pIItt3(1,1)+pIItt3(4,1)+pIItt3(5,1)+pIItt3(8,1)+pIItt3(18,1) + pIItt4(1,1)+pIItt4(4,1)+pIItt4(5,1)+pIItt4(8,1) + pIItt5(1,1)+pIItt5(4,1)+pIItt5(5,1)+pIItt5(8,1); % TT yield
TpTY=pIItt1(9,1)+pIItt1(12,1)+pIItt2(9,1)+pIItt2(13,1)+pIItt2(17,1)+...
    pIItt3(9,1)+pIItt3(13,1)+pIItt3(17,1)+...
    pIItt4(9,1)+pIItt4(13,1)+pIItt4(17,1)+...
    pIItt5(9,1)+pIItt5(12,1); %T+T yield

PP104=PP104+TTY*df*(sin((kk-1)*dc))*dc; %TT yield*probablity
PP105=PP105+TpTY*df*(sin((kk-1)*dc))*dc; %T+T yield*probablity
PP1055=PP1055+(TTY+TpTY)*df*(sin((kk-1)*dc))*dc; % total multiexcition

% PP1056=PP1056+pIItt_total1*df*(sin((kk-1)*dc))*dc;  %ms=+2
% PP1057=PP1057+pIItt_total2*df*(sin((kk-1)*dc))*dc;  %ms=+1
% PP1058=PP1058+pIItt_total3*df*(sin((kk-1)*dc))*dc;  %ms=0
% PP1059=PP1059+pIItt_total4*df*(sin((kk-1)*dc))*dc;  %ms=-1
% PP1060=PP1060+pIItt_total5*df*(sin((kk-1)*dc))*dc;  %ms=-2
% 
% 
% %             PP1006=PP1006+TTZZ; PP10006=PP10006+TTZZZ; PP106=PP106+TTXX;  % quantum yield of TTstate (Total, TT, T+T )
%              PP107=PP107+df*(sin((kk-1)*dc))*dc; % quantum yield =1
             
             
         end % of do roop for theta
       				
    end % of do roop for phai

     BBB=[];
 for mmm=1:ll
    BB=Bs+mmm*dB;
    BBB=[BBB,BB];
 end
 
% load tcad3us.dat;
%  NN=301;
% for j=50:250
%    vu(j)=tcad3us(j,1);
%    au(j)=tcad3us(j,2);
% end
% BBU0=[];aau0=[];
% for j=50:250
%    BBU0=[BBU0,vu(j)]; aau0=[aau0,au(j)]; 
% end
%  
%  
% load spsim_betam20_89.dat;
%  NN=200;
% for j=1:NN
%    vu(j)=spsim_betam20_89(j,1);
%    au(j)=spsim_betam20_89(j,2);
% end
% BBU=[];aau=[];
% for j=1:NN
%    BBU=[BBU,vu(j)]; aau=[aau,au(j)]; 
% end
% 
% load spsim_betam89_20.dat;
%  NN=200;
% for j=1:NN
%    vu(j)=spsim_betam89_20(j,1);
%    au(j)=spsim_betam89_20(j,2);
% end
% BBU1=[];aau1=[];
% for j=1:NN
%    BBU1=[BBU1,vu(j)]; aau1=[aau1,au(j)]; 
% end

PP5=imag(PP3)*1;
% effective matrix element of the ZFS in the basis of (S,T-,T0,T+,Q-2,Q-,Q0,Q+1,Q+2) 
ZFSmag2=ZFSmag1/(2*m)/(2*mm+1);
ZFSmagSQ(1,1)=ZFSmag51/(2*m)/(2*mm+1);
ZFSmagSQ(1,2)=ZFSmag61/(2*m)/(2*mm+1);
ZFSmagSQ(1,3)=ZFSmag71/(2*m)/(2*mm+1);
ZFSmagSQ(1,4)=ZFSmag81/(2*m)/(2*mm+1);
ZFSmagSQ(1,5)=ZFSmag91/(2*m)/(2*mm+1);

%roh00*4;
% TTX=TTY*df*(sin((kk-1)*dc))*dc; %TT yield*probablity
% TpTX=TpTY*df*(sin((kk-1)*dc))*dc; %T+T yield*probablity
% TTYY=TTX+TpTX; % total multiexcition

 PP1006=PP104/PP107; PP10006=PP105/PP107;  % TT yield T+T yield
ME_Yield_Total=PP1055/PP107
TpT_Yield=PP10006

pIItt_total1=PP1056/PP107; pIItt_total2=PP1057/PP107; pIItt_total3=PP1058/PP107; pIItt_total4=PP1059/PP107; pIItt_total5=PP1060/PP107;% quantum yield of TTstate (Total, TT, T+T )
%pIItt_total=pIItt_total1+pIItt_total2+pIItt_total3+pIItt_total4+pIItt_total5

pIItt_totalQ=pIItt_total1(1,1)+pIItt_total1(5,1)+pIItt_total1(9,1)+pIItt_total2(1,1)+pIItt_total2(5,1)+pIItt_total2(9,1)+pIItt_total3(1,1)+pIItt_total3(5,1)+pIItt_total3(9,1)+pIItt_total4(1,1)+pIItt_total4(5,1)+pIItt_total4(9,1)+pIItt_total5(1,1)+pIItt_total5(5,1)+pIItt_total5(9,1)
pIItt_totalT=pIItt_total2(13,1)+pIItt_total3(13,1)+pIItt_total4(13,1)+pIItt3(18,1)
pIItt_totalS=pIItt_total1(4,1)+pIItt_total1(8,1)+pIItt_total1(12,1)+pIItt_total2(4,1)+pIItt_total2(8,1)+pIItt_total2(17,1)+pIItt_total3(4,1)+pIItt_total3(8,1)+pIItt_total3(17,1)+pIItt_total4(4,1)+pIItt_total4(8,1)+pIItt_total4(17,1)+pIItt_total5(4,1)+pIItt_total5(8,1)+pIItt_total5(12,1)


pIItt_totalQ1=pIItt_total1(1,1)+pIItt_total1(5,1);%+pIItt_total1(9,1)
pIItt_totalQ2=pIItt_total2(1,1)+pIItt_total2(5,1);%+pIItt_total2(9,1)
pIItt_totalQ3=pIItt_total3(1,1)+pIItt_total3(5,1);%+pIItt_total3(9,1)
pIItt_totalQ4=pIItt_total4(1,1)+pIItt_total4(5,1);%+pIItt_total4(9,1)
pIItt_totalQ5=pIItt_total5(1,1)+pIItt_total5(5,1);%+pIItt_total5(9,1)

S0S0=1-pIItt_totalQ-pIItt_totalS-pIItt_totalT

% pIItt1(4,1)=pIItt1(4,1)/5;pIItt1(8,1)=pIItt1(8,1)/5;pIItt1(12,1)=pIItt1(12,1)/5;
% pIItt2(4,1)=pIItt2(4,1)/5;pIItt2(8,1)=pIItt2(8,1)/5;pIItt2(12,1)=pIItt2(12,1)/5;
% pIItt3(4,1)=pIItt3(4,1)/5;pIItt3(8,1)=pIItt3(8,1)/5;pIItt3(12,1)=pIItt3(12,1)/5;
% pIItt4(4,1)=pIItt4(4,1)/5;pIItt4(8,1)=pIItt4(8,1)/5;pIItt4(12,1)=pIItt4(12,1)/5;
% pIItt5(4,1)=pIItt5(4,1)/5;pIItt5(8,1)=pIItt5(8,1)/5;pIItt5(12,1)=pIItt5(12,1)/5;
% 
% pIItt_total = pIItt1+pIItt2+pIItt3+pIItt4+pIItt5
% TTY=pIItt1(1,1)+pIItt1(4,1)+pIItt1(5,1)+pIItt1(8,1)+pIItt1(9,1)+pIItt1(12,1) + pIItt2(1,1)+pIItt2(4,1)+pIItt2(5,1)+pIItt2(8,1)+pIItt2(9,1)+pIItt2(12,1) + pIItt3(1,1)+pIItt3(4,1)+pIItt3(5,1)+pIItt3(8,1)+pIItt3(9,1)+pIItt3(12,1) + pIItt4(1,1)+pIItt4(4,1)+pIItt4(5,1)+pIItt4(8,1)+pIItt4(9,1)+pIItt4(12,1) + pIItt5(1,1)+pIItt5(4,1)+pIItt5(5,1)+pIItt5(8,1)+pIItt5(9,1)+pIItt5(12,1)

%pIItt2(1,1)
%pIItt3
%pIItt4(1,1)
%pIItt5(1,1)
  
figure(153)
 %plot(BBB,imag(PP3))%
 plot(BBB,imag(PP4)+0.5e-8,BBB,imag(PP3),BBB,imag(PP3)+imag(PP4)-0.5e-8)%,'b-',BBB,PP5-10e-9,'r-',BBB,aau-5e-9,'b--',BBU0,aau0/2e6-13e-9,'b')%,BBB,imag(PP3+0.0*PP4),'r-')
 %axis([290 410.0 -1e-11 1e-11])
BBB=BBB'; PP41=imag(PP3.'+PP4.');PP42=imag(PP3.')-0.4e-9;PP51=imag(PP4.')-0.2e-9;
spsim=[BBB,PP41,PP42,PP51];

save 0.005p_200ns_MPD3_TREPR_flip.dat spsim -ascii;
% nutation(4,3)=abs(QQ2(1,3))*abs(QQ2(2,4))*2 + abs(QQ2(2,3))*abs(QQ2(1,4))*2 + abs(QQ2(2,3))*abs(QQ2(3,4))*sqrt(6)+...
%               abs(QQ2(3,3))*abs(QQ2(2,4))*sqrt(6) + abs(QQ2(3,3))*abs(QQ2(4,4))*sqrt(6) + abs(QQ2(4,3))*abs(QQ2(3,4))*sqrt(6) +...
%               abs(QQ2(4,3))*abs(QQ2(5,4))*2 + abs(QQ2(5,3))*abs(QQ2(4,4))*2 + abs(QQ2(6,3))*abs(QQ2(7,4))*sqrt(2)+...
%               abs(QQ2(7,3))*abs(QQ2(6,4))*sqrt(2) + abs(QQ2(7,3))*abs(QQ2(8,4))*sqrt(2) + abs(QQ2(8,3))*abs(QQ2(7,4))*sqrt(2);
% %sqrt(Iint(4,3));
% nutation(4,2)=abs(QQ2(1,2))*abs(QQ2(2,4))*2 + abs(QQ2(2,2))*abs(QQ2(1,4))*2 + abs(QQ2(2,2))*abs(QQ2(3,4))*sqrt(6)+...
%               abs(QQ2(3,2))*abs(QQ2(2,4))*sqrt(6) + abs(QQ2(3,2))*abs(QQ2(4,4))*sqrt(6) + abs(QQ2(4,2))*abs(QQ2(3,4))*sqrt(6) +...
%               abs(QQ2(4,2))*abs(QQ2(5,4))*2 + abs(QQ2(5,2))*abs(QQ2(4,4))*2 + abs(QQ2(6,2))*abs(QQ2(7,4))*sqrt(2)+...
%               abs(QQ2(7,2))*abs(QQ2(6,4))*sqrt(2) + abs(QQ2(7,2))*abs(QQ2(8,4))*sqrt(2) + abs(QQ2(8,2))*abs(QQ2(7,4))*sqrt(2);
% sqrt(Iint(4,2));
nutation(4,2)=sqrt(Iint(4,2));nutation(4,3)=sqrt(Iint(4,3));
nutation(5,2)=sqrt(Iint(5,2));nutation(5,3)=sqrt(Iint(5,3));
nutation(6,2)=sqrt(Iint(6,2));nutation(6,3)=sqrt(Iint(6,3));
nutation(7,6)=sqrt(Iint(7,6));nutation(8,6)=sqrt(Iint(8,6));
nutation

toc;