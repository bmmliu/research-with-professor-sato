%Echebarria B, Karma A.
%Mechanisms for initiation of cardiac discordant alternans.
%The European Physical Journal Special Topics 2007; 146: 217‐31.
%http://link.springer.com/article/10.1140/epjst/e2007-00182-y
%

clear all;

% initial values
num=400;% the number of cells in a cable
v=zeros(num,1);
tmp=zeros(num,1);%temporary array to solve the PDE
h=1.0*ones(num,1);
f=0.9*ones(num,1);
stim=zeros(num,1);

% constants
tauso=15;
taufi=0.8;
tauh1=4.8;
tauh2=10.0;

tausi=4;
tauf1=100;
tauf2=30;

% constants for the PDE
dfu=0.0005;%diffusion coefficient
dx=0.015;%cell size

dt=0.1;% time step (0.1 ms)

pcl=140;%pacing cycle length (ms)
itr=10;% the number of beats
tmax=pcl*itr;%total time

%convert to integer
tnmax=round(tmax/dt);
pcln=round(pcl/dt);
durn=round(1.0/dt);

%array for results
resultv=[];

% main loop
for tn=0:tnmax

  % stimulation
  if (mod(tn,pcln)<durn)
    stim(1:5)=0.2;
  elseif(mod(tn,pcln)==durn)
    stim(1:5)=0;
  end

  minf=((v/0.2).^6)./(1+((v/0.2).^6));
  hinf=1./(1+((v/0.1).^6));
  dinf=((v/0.4).^4)./(1+((v/0.4).^4));
  finf=1./(1+((v/0.1).^4));

  tauh=tauh1+tauh2*exp(-20*((v-0.1).^2));
  tauf=tauf2+(tauf1-tauf2).*v.^3;

  jfi=1 *h.*minf.*(v-1.3)/taufi;
  jsi=1 * f.*dinf.*(v-1.4)/tausi;
  jso=(1-exp(-4*v))/tauso;
  ion=-(jfi+jsi+jso-stim);

  % update variables
  v=v+ion*dt;
  h=h+dt*(hinf-h)./tauh;
  f=f+dt*(finf-f)./tauf;

  %non-flux boundary condition
  
  v(1)=v(3);
  v(num)=v(num-2);
  %solve diffusion equation
  for c=2:num-1
    tmp(c)=v(c)+(v(c-1)+v(c+1)-2*v(c))*dfu*dt/(dx*dx);
    s=v(c);
    s1=v(c-1);
    s2=v(c+1);
    s4 = tmp(c);
  end

  v=tmp;

  if (mod(tn,10)==0)
    resultv=[resultv v];
  end
end

%space-time plot
surf(resultv);
shading interp
axis([0 inf 0 inf]);
view(3)