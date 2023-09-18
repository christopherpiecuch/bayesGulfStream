%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_main_code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code for producing results in the main text of Piecuch and Beal (2023) 
%   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL105170
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code written by CGP (last changed 18 September 2023) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_main_code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Say hello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(0.1)
disp('Hello.  Things have started.')
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN_burn=20000; % 20,000 burn-in draws
NN_post=10000; % 10,000 post burn-in draws
thin_period=50; % thin and only keep every 50th draw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% name of experiment and file to be saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_name=[date,'_output_file'];  % file name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the seeds of the random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(249*sum(clock))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cable data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('20211210_daily_fs_transport_sv_19820318_20211206.mat','tr','dn')
    cab_dn=dn;
    cab_tr=tr;  

    cab_id=nan*cab_tr;
    cab_id(find(str2num(datestr(cab_dn,10))<=1992))=1;
    cab_id(find(str2num(datestr(cab_dn,10))>=1993&str2num(datestr(cab_dn,10))<=1998))=2;
    cab_id(find(str2num(datestr(cab_dn,10))>=2000&str2num(datestr(cab_dn,10))<=2005))=3;
    cab_id(find(str2num(datestr(cab_dn,10))>=2006))=4;
    cab_er(find(str2num(datestr(cab_dn,10))<=1992))=1.0; 
    cab_er(find(str2num(datestr(cab_dn,10))>=1993&str2num(datestr(cab_dn,10))<=1998))=2.0;
    cab_er(find(str2num(datestr(cab_dn,10))>=2000&str2num(datestr(cab_dn,10))<=2005))=1.5;
    cab_er(find(str2num(datestr(cab_dn,10))>=2006))=1.0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('FC_section_transport.mat')
    sec_dn=datenum(Year,Month,Day);
    sec_tr=Transport;
    sec_id=Type;

    sec_er(find(sec_id==1))=1.0; % dropsonde
    sec_er(find(sec_id==2))=1.5; % ladcp
    sec_er(find(sec_id==3))=1.0; % pegasus
    sec_er(find(sec_id==4))=1.0; % pegasus in dropsonde mode

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Altimetry data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('fc_transport_from_altimetry.mat')
    ii=find(Flag==1&Taltimetry<1e4);
    alt_dn=datenum(Year(ii),Month(ii),Day(ii));
    alt_tr=Taltimetry(ii);
    alt_id=ones(size(alt_tr));
    alt_er=2.0*ones(size(alt_tr));

clearvars -except L M N P *_dn *_tr *_er NN* thin* save_name 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on data, define time span
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=min([min(cab_dn) min(sec_dn) min(alt_dn)])-1; % 1 day before first measurement
tm1=min([min(cab_dn) min(sec_dn) min(alt_dn)])-2; % 2 days before first measurement
tm2=min([min(cab_dn) min(sec_dn) min(alt_dn)])-3; % 3 days before first measurement
t=min([min(cab_dn) min(sec_dn) min(alt_dn)]):(max([max(cab_dn) max(sec_dn) max(alt_dn)]));
t=t';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(numel(t),1);
delta2=Inf*ones(size(x)); 
for nn=1:numel(x)
 ii=[]; ii=find(sec_dn==t(nn));
 if ~isempty(ii)
     x(nn)=mean(sec_tr(ii));
     delta2(nn)=1/numel(ii)*sum(sec_er(ii).^2);
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cable data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=zeros(numel(t),1);
epsilon2=Inf*ones(size(y)); 
for nn=1:numel(y)
 ii=[]; ii=find(cab_dn==t(nn)&~isnan(cab_tr));
 if ~isempty(ii)
	y(nn)=mean(cab_tr(ii));
	epsilon2(nn)=1/numel(ii)*sum(cab_er(ii).^2);
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Altimetry data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=zeros(numel(t),1);
omega2=Inf*ones(size(z)); 
for nn=1:numel(z)
 ii=[]; ii=find(alt_dn==t(nn));
 if ~isempty(ii)
	z(nn)=mean(alt_tr(ii));
	omega2(nn)=1/numel(ii)*sum(alt_er(ii).^2);
 end
end
clear *_dn *_tr *_er ii nn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define time on more natural interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=mean(t);
t0=t0-mm;
tm1=tm1-mm;
tm2=tm2-mm;
t=t-mm; clear mm
tSa=365.25;
tSsa=365.25/2;
% map t onto the interval ~[-1,1]
tSa=tSa/max(t);
tSsa=tSsa/max(t);
t0=t0/max(t);
tm1=tm1/max(t);
tm2=tm2/max(t);
t=t/max(t);
K=numel(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define design matrix [PxK]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=[]; w=[ones(size(t')); t'; cos(2*pi*t'/tSa); sin(2*pi*t'/tSa); cos(2*pi*t'/tSsa); sin(2*pi*t'/tSsa)];
w0=[1; t0; cos(2*pi*t0/tSa); sin(2*pi*t0/tSa); cos(2*pi*t0/tSsa); sin(2*pi*t0/tSsa)];
wm1=[1; tm1; cos(2*pi*tm1/tSa); sin(2*pi*tm1/tSa); cos(2*pi*tm1/tSsa); sin(2*pi*tm1/tSsa)];
wm2=[1; tm2; cos(2*pi*tm2/tSa); sin(2*pi*tm2/tSa); cos(2*pi*tm2/tSsa); sin(2*pi*tm2/tSsa)];
P=numel(w)/K;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define iteration parameters based on input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN_burn_thin=NN_burn/thin_period;    % Total number of burn-in to keep
NN_post_thin=NN_post/thin_period;     % Total number of post-burn-in to keep
NN=NN_burn+NN_post;                       % Total number of draws to take 
NN_thin=NN_burn_thin+NN_post_thin;% Total number of draws to keep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_beta=[30 0 0 0 0 0]'; % mean, trend, Sa cos and sin, Ssa cos and sin; cf. w above
Z_beta=25*eye(P);
zeta2_rho=0.5^2;
zeta2_theta=0.5^2;
zeta2_phi=0.5^2;
xi_sigma2=0.5;
chi_sigma2=0.5*1; % units are Sv^2
xi_tau2=0.5;
chi_tau2=0.5*1; % units are Sv^2
mu_t0=30; % units are Sv
zeta2_t0=(5)^2; % units are Sv^2
epsilon02=epsilon2(find(epsilon2~=Inf,1,'first'));
if isempty(epsilon02)
    epsilon02=Inf; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate space for the sample arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR0=zeros(NN,1);
TRM1=zeros(NN,1);
TRM2=zeros(NN,1);
TR=zeros(NN,K);
E0=zeros(NN,1);
EM1=zeros(NN,1);
E=zeros(NN,K);
S=zeros(NN,K);
G=zeros(NN,K);
F=zeros(NN,K);
D=zeros(NN,K);
U0=zeros(NN,1);
U=zeros(NN,K);
BETA=zeros(NN,P);
SIGMA2=zeros(NN,1);
TAU2=zeros(NN,1);
RHO1=zeros(NN,1);
RHO2=zeros(NN,1);
RHO3=zeros(NN,1);
THETA1=zeros(NN,1);
THETA2=zeros(NN,1);
PHI=zeros(NN,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=[]; beta=zeros(P,1);%mvnrnd(mu_beta,Z_beta)';
beta(1)=nanmean(y(y~=0));
rho1=[]; rho1=mvnrnd(0,zeta2_rho);
rho2=[]; rho2=mvnrnd(0,zeta2_rho);
rho3=[]; rho3=mvnrnd(0,zeta2_rho);
theta1=[]; theta1=mvnrnd(0,zeta2_theta);
theta2=[]; theta2=mvnrnd(0,zeta2_theta);
phi=[]; phi=mvnrnd(0,zeta2_phi);
sigma2=[]; sigma2=min([10 1/randraw('gamma', [0,1/(chi_sigma2),(xi_sigma2)], [1,1])]); 
tau2=[]; tau2=min([10 1/randraw('gamma', [0,1/(chi_tau2),(xi_tau2)], [1,1])]); 
tr0=[]; tr0=mu_t0+sqrt(zeta2_t0)*randn(1);
trm1=[]; trm1=mu_t0+sqrt(zeta2_t0)*randn(1);
trm2=[]; trm2=mu_t0+sqrt(zeta2_t0)*randn(1);
tr=[]; tr=mu_t0*ones(K,1)+sqrt(zeta2_t0)*randn(K,1);
e0=[]; e0=sqrt(epsilon02)*randn(1);
em1=[]; em1=sqrt(epsilon02)*randn(1);
e=[]; e=sqrt(epsilon02)*randn(K,1);
u0=[]; u0=mu_t0+sqrt(zeta2_t0)*randn(1);
u=[]; u=mu_t0*ones(K,1)+sqrt(zeta2_t0)*randn(K,1);

I=double(omega2~=Inf);
IJK=find(epsilon2~=Inf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through the Gibbs sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nn=1:NN
tic
    if mod(nn,10)==0 
        disp([num2str(nn),' of ',num2str(NN),' iterations done.']), 
    end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(beta|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=inv(Z_beta)*mu_beta;
    PSI=inv(Z_beta);
    for k=1:K
        if k==1
            sumrhow=[]; sumrhow=rho1*w0+rho2*wm1+rho3*wm2;
            sumrhot=[]; sumrhot=rho1*tr0+rho2*trm1+rho3*trm2;
        elseif k==2
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w0+rho3*wm1;
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr0+rho3*trm1;
        elseif k==3
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w(:,k-2)+rho3*w0;
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr0;
        else % k>3
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w(:,k-2)+rho3*w(:,k-3);
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3);
        end
        V=V+(w(:,k)-sumrhow)*(tr(k)-sumrhot)/sigma2;
     	PSI=PSI+(w(:,k)-sumrhow)*(w(:,k)'-sumrhow')/sigma2;
    end
    PSI=inv(PSI);
    beta=mvnrnd(PSI*V,PSI)';
    clear V PSI sumrho*
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(sigma_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for k=1:K
        if k==1
            sumrhow=[]; sumrhow=rho1*w0+rho2*wm1+rho3*wm2;
            sumrhot=[]; sumrhot=rho1*tr0+rho2*trm1+rho3*trm2;
        elseif k==2
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w0+rho3*wm1;
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr0+rho3*trm1;
        elseif k==3
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w(:,k-2)+rho3*w0;
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr0;
        else % k>3
            sumrhow=[]; sumrhow=rho1*w(:,k-1)+rho2*w(:,k-2)+rho3*w(:,k-3);
            sumrhot=[]; sumrhot=rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3);
        end
        DYKK=[];
        DYKK=tr(k)-sumrhot-(w(:,k)'-sumrhow')*beta;
        SUM_K=SUM_K+DYKK^2;           
    end
   	sigma2=1/randraw('gamma', [0,1/(chi_sigma2+1/2*SUM_K),...
     	(xi_sigma2+K/2)], [1,1]);
   	clear SUM_K DYKK sumrho*
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tau_2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUM_K=0;
    for k=1:K
        if k==1
            sumg=[]; sumg=u(k)-tr(k)-phi*(u0-tr0);
        else % k>1
            sumg=[]; sumg=u(k)-tr(k)-phi*(u(k-1)-tr(k-1));
        end
        SUM_K=SUM_K+sumg^2;           
    end
   	tau2=1/randraw('gamma', [0,1/(chi_tau2+1/2*SUM_K),...
     	(xi_tau2+K/2)], [1,1]);
   	clear SUM_K sumg
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho1|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=0;
    PSI=1/zeta2_rho;
    for k=1:K
        if k==1
            sumki=[]; sumki=tr0-w0'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho2*(trm1-wm1'*beta)+rho3*(trm2-wm2'*beta));
        elseif k==2
            sumki=[]; sumki=tr(k-1)-w(:,k-1)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho2*(tr0-w0'*beta)+rho3*(trm1-wm1'*beta));
        elseif k==3
            sumki=[]; sumki=tr(k-1)-w(:,k-1)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho2*(tr(k-2)-w(:,k-2)'*beta)+rho3*(tr0-w0'*beta));
        else % k>3
            sumki=[]; sumki=tr(k-1)-w(:,k-1)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho2*(tr(k-2)-w(:,k-2)'*beta)+rho3*(tr(k-3)-w(:,k-3)'*beta));
        end
        V=V+sumki*sumk/sigma2;
     	PSI=PSI+(sumki^2)/sigma2;
    end
    PSI=inv(PSI);
    rho1=mvnrnd(PSI*V,PSI)';
    clear V PSI sumk*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=0;
    PSI=1/zeta2_rho;
    for k=1:K
        if k==1
            sumki=[]; sumki=trm1-wm1'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr0-w0'*beta)+rho3*(trm2-wm2'*beta));
        elseif k==2
            sumki=[]; sumki=tr0-w0'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho3*(trm1-wm1'*beta));
        elseif k==3
            sumki=[]; sumki=tr(k-2)-w(:,k-2)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho3*(tr0-w0'*beta));
        else % k>3
            sumki=[]; sumki=tr(k-2)-w(:,k-2)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho3*(tr(k-3)-w(:,k-3)'*beta));
        end
        V=V+sumki*sumk/sigma2;
     	PSI=PSI+(sumki^2)/sigma2;
    end
    PSI=inv(PSI);
    rho2=mvnrnd(PSI*V,PSI)';
    clear V PSI sumk*
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(rho3|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=0;
    PSI=1/zeta2_rho;
    for k=1:K
        if k==1
            sumki=[]; sumki=trm2-wm2'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr0-w0'*beta)+rho2*(trm1-wm1'*beta));
        elseif k==2
            sumki=[]; sumki=trm1-wm1'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho2*(tr0-w0'*beta));
        elseif k==3
            sumki=[]; sumki=tr0-w0'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho2*(tr(k-2)-w(:,k-2)'*beta));
        else % k>3
            sumki=[]; sumki=tr(k-3)-w(:,k-3)'*beta;
            sumk=[]; sumk=tr(k)-w(:,k)'*beta-(rho1*(tr(k-1)-w(:,k-1)'*beta)+rho2*(tr(k-2)-w(:,k-2)'*beta));
        end
        V=V+sumki*sumk/sigma2;
     	PSI=PSI+(sumki^2)/sigma2;
    end
    PSI=inv(PSI);
    rho3=mvnrnd(PSI*V,PSI)';
    clear V PSI sumk*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(theta1|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(epsilon2~=Inf)~=0
        V=[]; PSI=[];
        V=0;
        PSI=1/zeta2_theta;
        for k=1:K
            if k==1
                V=V+e0*(y(k)-tr(k)-theta2*em1)/epsilon2(k);
                PSI=PSI+(e0^2)/epsilon2(k);
            elseif k==2
                V=V+e(k-1)*(y(k)-tr(k)-theta2*e0)/epsilon2(k);
                PSI=PSI+(e(k-1)^2)/epsilon2(k);
            else
                V=V+e(k-1)*(y(k)-tr(k)-theta2*e(k-2))/epsilon2(k);
            PSI=PSI+(e(k-1)^2)/epsilon2(k);
            end        
        end
        PSI=inv(PSI);
        theta1=mvnrnd(PSI*V,PSI)';
        clear V PSI sumk*
    else
        theta1=0;
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(theta2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(epsilon2~=Inf)~=0
        V=[]; PSI=[];
        V=0;
        PSI=1/zeta2_theta;
        for k=1:K
            if k==1
                V=V+em1*(y(k)-tr(k)-theta1*e0)/epsilon2(k);
                PSI=PSI+(em1^2)/epsilon2(k);
            elseif k==2
                V=V+e0*(y(k)-tr(k)-theta1*e(k-1))/epsilon2(k);
                PSI=PSI+(e0^2)/epsilon2(k);
            else
                V=V+e(k-2)*(y(k)-tr(k)-theta1*e(k-1))/epsilon2(k);
                PSI=PSI+(e(k-2)^2)/epsilon2(k);
            end        
        end
        PSI=inv(PSI);
        theta2=mvnrnd(PSI*V,PSI)';
        clear V PSI sumk*
    else
        theta2=0;
    end
    
  	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(phi|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=0;
    PSI=1/zeta2_phi;
    for k=1:K
        if k==1
            sumk=[]; sumk=u(k)-tr(k);
            sumkm1=[]; sumkm1=u0-tr0;
        else % k>2
            sumk=[]; sumk=u(k)-tr(k);
            sumkm1=[]; sumkm1=u(k-1)-tr(k-1);
        end
        V=V+sumkm1*sumk/tau2;
     	PSI=PSI+(sumkm1^2)/tau2;
    end
    PSI=inv(PSI);
    phi=mvnrnd(PSI*V,PSI)';
    clear V PSI sumk*
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(e0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(epsilon2~=Inf)~=0
        V=[]; PSI=[];
        V=theta1*(y(1)-tr(1)-theta2*em1)/epsilon2(1)+theta2*(y(2)-tr(2)-theta1*e(1))/epsilon2(2);
        PSI=1/(1/epsilon02+(theta1^2)/epsilon2(1)+(theta2^2)/epsilon2(2));
        e0=mvnrnd(PSI*V,PSI);
        clear V PSI
    else
        e0=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(em1|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(epsilon2~=Inf)~=0
        V=[]; PSI=[];
        V=theta2*(y(1)-tr(1)-theta1*e0)/epsilon2(1);
        PSI=1/(1/epsilon02+(theta2^2)/epsilon2(1));
        em1=mvnrnd(PSI*V,PSI);
        clear V PSI
    else
        em1=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(u0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=mu_t0/zeta2_t0+phi*(u(1)-(tr(1)-phi*tr0))/tau2;
    PSI=1/(1/zeta2_t0+(phi^2)/tau2);
    u0=mvnrnd(PSI*V,PSI);
    clear V PSI
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(uk|.), k=1,...,K-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=1:(K-1)
            V=[]; PSI=[];
            if k==1
                V=z(1)/omega2(1)+((1+phi^2)*tr(1)-phi*(tr0+tr(2))+phi*(u0+u(2)))/tau2;
                PSI=1/(1/omega2(1)+(1+phi^2)/tau2);         
            else
                V=z(k)/omega2(k)+((1+phi^2)*tr(k)-phi*(tr(k-1)+tr(k+1))+phi*(u(k-1)+u(k+1)))/tau2;
                PSI=1/(1/omega2(k)+(1+phi^2)/tau2);
            end
            u(k)=mvnrnd(PSI*V,PSI);
            clear V PSI
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(uK|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V=[]; PSI=[];
        V=z(K)/omega2(K)+(tr(K)+phi*(u(K-1)-tr(K-1)))/tau2;
        PSI=1/(1/omega2(K)+1/tau2);
        u(K)=mvnrnd(PSI*V,PSI);
        clear V PSI

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(trm2|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=mu_t0/zeta2_t0+rho3*(tr(1)-(w(:,1)'-(rho1*w0'+rho2*wm1'+rho3*wm2'))*beta-rho1*tr0-rho2*trm1)/sigma2;
    PSI=1/(1/zeta2_t0+(rho3^2)/sigma2);
    trm2=mvnrnd(PSI*V,PSI);
    clear V PSI
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(trm1|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=mu_t0/zeta2_t0+...
        rho2*(tr(1)-(w(:,1)'-(rho1*w0'+rho2*wm1'+rho3*wm2'))*beta-rho1*tr0-rho3*trm2)/sigma2+...
        rho3*(tr(2)-(w(:,2)'-(rho1*w(:,1)'+rho2*w0'+rho3*wm1'))*beta-rho1*tr(1)-rho2*tr0)/sigma2;
    PSI=1/(1/zeta2_t0+(rho2^2+rho3^2)/sigma2);
    trm1=mvnrnd(PSI*V,PSI);
    clear V PSI
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(tr0|.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=[]; PSI=[];
    V=mu_t0/zeta2_t0+phi*(tr(1)-(u(1)-phi*u0))/tau2+...
        rho1*(tr(1)-(w(:,1)'-(rho1*w0'+rho2*wm1'+rho3*wm2'))*beta-rho2*trm1-rho3*trm2)/sigma2+...
        rho2*(tr(2)-(w(:,2)'-(rho1*w(:,1)'+rho2*w0'+rho3*wm1'))*beta-rho1*tr(1)-rho3*trm1)/sigma2+...
        rho3*(tr(3)-(w(:,3)'-(rho1*w(:,2)'+rho2*w(:,1)'+rho3*w0'))*beta-rho1*tr(2)-rho2*tr(1))/sigma2;
    PSI=1/(1/zeta2_t0+(rho1^2+rho2^2+rho3^2)/sigma2+(phi^2)/tau2);
    tr0=mvnrnd(PSI*V,PSI);
    clear V PSI
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(trk|.), k=1,...,K-3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:K
            V=[]; PSI=[];
            if k==1
                V=x(k)/delta2(k)+(y(k)-theta1*e0-theta2*em1)/epsilon2(k)+...
                    (u(k)-phi*(u0-tr0))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w0'+rho2*wm1'+rho3*wm2'))*beta+rho1*tr0+rho2*trm1+rho3*trm2)/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w0'+rho3*wm1'))*beta-rho2*tr0-rho3*trm1)/sigma2+...
                    rho2*(tr(k+2)-(w(:,k+2)'-(rho1*w(:,k+1)'+rho2*w(:,k)'+rho3*w0'))*beta-rho1*tr(k+1)-rho3*tr0)/sigma2+...
                    rho3*(tr(k+3)-(w(:,k+3)'-(rho1*w(:,k+2)'+rho2*w(:,k+1)'+rho3*w(:,k)'))*beta-rho1*tr(k+2)-rho2*tr(k+1))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2+rho2^2+rho3^2)/sigma2+(1+phi^2)/tau2);
            elseif k==2
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e0)/epsilon2(k)+...
                    (u(k)-phi*(u(k-1)-tr(k-1)))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w0'+rho3*wm1'))*beta+rho1*tr(k-1)+rho2*tr0+rho3*trm1)/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w(:,k-1)'+rho3*w0'))*beta-rho2*tr(k-1)-rho3*tr0)/sigma2+...
                    rho2*(tr(k+2)-(w(:,k+2)'-(rho1*w(:,k+1)'+rho2*w(:,k)'+rho3*w(:,k-1)'))*beta-rho1*tr(k+1)-rho3*tr(k-1))/sigma2+...
                    rho3*(tr(k+3)-(w(:,k+3)'-(rho1*w(:,k+2)'+rho2*w(:,k+1)'+rho3*w(:,k)'))*beta-rho1*tr(k+2)-rho2*tr(k+1))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2+rho2^2+rho3^2)/sigma2+(1+phi^2)/tau2);
            elseif k==3
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e(k-2))/epsilon2(k)+...
                    (u(k)-phi*(u(k-1)-tr(k-1)))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w(:,k-2)'+rho3*w0'))*beta+rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr0)/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w(:,k-1)'+rho3*w(:,k-2)'))*beta-rho2*tr(k-1)-rho3*tr(k-2))/sigma2+...
                    rho2*(tr(k+2)-(w(:,k+2)'-(rho1*w(:,k+1)'+rho2*w(:,k)'+rho3*w(:,k-1)'))*beta-rho1*tr(k+1)-rho3*tr(k-1))/sigma2+...
                    rho3*(tr(k+3)-(w(:,k+3)'-(rho1*w(:,k+2)'+rho2*w(:,k+1)'+rho3*w(:,k)'))*beta-rho1*tr(k+2)-rho2*tr(k+1))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2+rho2^2+rho3^2)/sigma2+(1+phi^2)/tau2);
            elseif ((k>3)&&(k<(K-2)))
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e(k-2))/epsilon2(k)+...
                    (u(k)-phi*(u(k-1)-tr(k-1)))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w(:,k-2)'+rho3*w(:,k-3)'))*beta+rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3))/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w(:,k-1)'+rho3*w(:,k-2)'))*beta-rho2*tr(k-1)-rho3*tr(k-2))/sigma2+...
                    rho2*(tr(k+2)-(w(:,k+2)'-(rho1*w(:,k+1)'+rho2*w(:,k)'+rho3*w(:,k-1)'))*beta-rho1*tr(k+1)-rho3*tr(k-1))/sigma2+...
                    rho3*(tr(k+3)-(w(:,k+3)'-(rho1*w(:,k+2)'+rho2*w(:,k+1)'+rho3*w(:,k)'))*beta-rho1*tr(k+2)-rho2*tr(k+1))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2+rho2^2+rho3^2)/sigma2+(1+phi^2)/tau2);
            elseif k==(K-2)
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e(k-2))/epsilon2(k)+...
                    (u(k)-phi*(u(k-1)-tr(k-1)))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w(:,k-2)'+rho3*w(:,k-3)'))*beta+rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3))/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w(:,k-1)'+rho3*w(:,k-2)'))*beta-rho2*tr(k-1)-rho3*tr(k-2))/sigma2+...
                    rho2*(tr(k+2)-(w(:,k+2)'-(rho1*w(:,k+1)'+rho2*w(:,k)'+rho3*w(:,k-1)'))*beta-rho1*tr(k+1)-rho3*tr(k-1))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2+rho2^2)/sigma2+(1+phi^2)/tau2);
            elseif k==(K-1)
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e(k-2))/epsilon2(k)+...
                    (u(k)-phi*(u(k-1)-tr(k-1)))/tau2+phi*(tr(k+1)-(u(k+1)-phi*u(k)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w(:,k-2)'+rho3*w(:,k-3)'))*beta+rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3))/sigma2+...
                    rho1*(tr(k+1)-(w(:,k+1)'-(rho1*w(:,k)'+rho2*w(:,k-1)'+rho3*w(:,k-2)'))*beta-rho2*tr(k-1)-rho3*tr(k-2))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+(1+rho1^2)/sigma2+(1+phi^2)/tau2);
            else % k==(K)
                V=x(k)/delta2(k)+(y(k)-theta1*e(k-1)-theta2*e(k-2))/epsilon2(k)+(u(k)-phi*(u(k-1)-tr(k-1)))/tau2+...
                    ((w(:,k)'-(rho1*w(:,k-1)'+rho2*w(:,k-2)'+rho3*w(:,k-3)'))*beta+rho1*tr(k-1)+rho2*tr(k-2)+rho3*tr(k-3))/sigma2;
                PSI=1/(1/delta2(k)+1/epsilon2(k)+1/sigma2+1/tau2);
            end
            tr(k)=mvnrnd(PSI*V,PSI);
            clear V PSI
        
        % save out e
        if sum(epsilon2~=Inf)~=0
                if y(k)==0 % no data; just draw prior
                    dijk=abs(k-IJK); dijk=IJK(find(dijk==min(dijk),1,'first'));
                    e(k)=randn(1)*sqrt(epsilon2(dijk));
                else % data
                    if k==1
                        e(k)=y(k)-tr(k)-theta1*e0-theta2*em1;
                    elseif k==2
                        e(k)=y(k)-tr(k)-theta1*e(k-1)-theta2*e0;
                    else
                        e(k)=y(k)-tr(k)-theta1*e(k-1)-theta2*e(k-2);
                    end
                end
        else
            e=zeros(size(e));
        end
    end
        
    % define residuals
    d=x-tr;
    d(find(x==0))=nan;
    f=z-u;
    f(find(z==0))=nan;
    for k=1:K
        if k==1
            s(k)=tr(k)-w(:,k)'*beta-...
                rho1*(tr0-w0'*beta)-...
                rho2*(trm1-wm1'*beta)-...
                rho3*(trm2-wm2'*beta);
        elseif k==2
            s(k)=tr(k)-w(:,k)'*beta-...
                rho1*(tr(k-1)-w(:,k-1)'*beta)-...
                rho2*(tr0-w0'*beta)-...
                rho3*(trm1-wm1'*beta);
        elseif k==3
            s(k)=tr(k)-w(:,k)'*beta-...
                rho1*(tr(k-1)-w(:,k-1)'*beta)-...
                rho2*(tr(k-2)-w(:,k-2)'*beta)-...
                rho3*(tr0-w0'*beta);
        else
            s(k)=tr(k)-w(:,k)'*beta-...
                rho1*(tr(k-1)-w(:,k-1)'*beta)-...
                rho2*(tr(k-2)-w(:,k-2)'*beta)-...
                rho3*(tr(k-3)-w(:,k-3)'*beta);
        end
    end
    
    % define G
    for k=1:K
        if k==1
            g(k)=u(k)-tr(k)-phi*(u0-tr0);
        else
            g(k)=u(k)-tr(k)-phi*(u(k-1)-tr(k-1));
        end
    end
    
    % save out
    TR0(nn)=tr0;
    TRM1(nn)=trm1;
    TRM2(nn)=trm2;
    TR(nn,:)=tr;
    E0(nn)=e0;
    EM1(nn)=em1;
    E(nn,:)=e;
    U0(nn)=u0;
    U(nn,:)=u;
    BETA(nn,:)=beta;
    SIGMA2(nn)=sigma2;
    TAU2(nn)=tau2;
    RHO1(nn)=rho1;
    RHO2(nn)=rho2;
    RHO3(nn)=rho3;
    THETA1(nn)=theta1;
    THETA2(nn)=theta2;
    PHI(nn)=phi;
    S(nn,:)=s;
    F(nn,:)=f;
    D(nn,:)=d;
    G(nn,:)=g;

toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear burnin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR0(1:NN_burn)=[];
TRM1(1:NN_burn)=[];
TRM2(1:NN_burn)=[];
TR(1:NN_burn,:)=[];
E0(1:NN_burn)=[];
EM1(1:NN_burn)=[];
E(1:NN_burn,:)=[];
U0(1:NN_burn)=[];
U(1:NN_burn,:)=[];
BETA(1:NN_burn,:)=[];
SIGMA2(1:NN_burn)=[];
TAU2(1:NN_burn)=[];
RHO1(1:NN_burn)=[];
RHO2(1:NN_burn)=[];
RHO3(1:NN_burn)=[];
THETA1(1:NN_burn)=[];
THETA2(1:NN_burn)=[];
PHI(1:NN_burn)=[];
S(1:NN_burn,:)=[];
F(1:NN_burn,:)=[];
D(1:NN_burn,:)=[];
G(1:NN_burn,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thin chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR0=TR0(thin_period:thin_period:NN_post);
TRM1=TRM1(thin_period:thin_period:NN_post);
TRM2=TRM2(thin_period:thin_period:NN_post);
TR=TR(thin_period:thin_period:NN_post,:);
E0=E0(thin_period:thin_period:NN_post);
EM1=EM1(thin_period:thin_period:NN_post);
E=E(thin_period:thin_period:NN_post,:);
U0=U0(thin_period:thin_period:NN_post);
U=U(thin_period:thin_period:NN_post,:);
BETA=BETA(thin_period:thin_period:NN_post,:);
SIGMA2=SIGMA2(thin_period:thin_period:NN_post);
TAU2=TAU2(thin_period:thin_period:NN_post);
RHO1=RHO1(thin_period:thin_period:NN_post);
RHO2=RHO2(thin_period:thin_period:NN_post);
RHO3=RHO3(thin_period:thin_period:NN_post);
THETA1=THETA1(thin_period:thin_period:NN_post);
THETA2=THETA2(thin_period:thin_period:NN_post);
PHI=PHI(thin_period:thin_period:NN_post);
S=S(thin_period:thin_period:NN_post,:);
F=F(thin_period:thin_period:NN_post,:);
D=D(thin_period:thin_period:NN_post,:);
G=G(thin_period:thin_period:NN_post,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thin chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except TR* E* U* G BETA SIGMA2 TAU2 RHO* THETA* PHI S F D x y z delta2 epsilon2 omega2 t save_name w
TIME=datenum(1982,3,17+(1:numel(t)));
save(save_name,'TIME','TR*','E*','U*','BETA','SIGMA2','TAU2','RHO*','THETA*','PHI','S','F','D','G','x','y','z','delta2','epsilon2','omega2','t','w*')
return