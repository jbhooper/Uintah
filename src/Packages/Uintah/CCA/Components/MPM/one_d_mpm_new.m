function KK=one_d_mpm_new(problem_type,CFL,NN)

% One dimensional MPM

PPC=1;

E=1e8;
density = 1000.;
G=1e6;
K=1e6;

c = sqrt(E/density);

bar_length = 1.;
domain=1.;
area=1.;
dx=domain/(NN-1);
volp=dx/PPC;

dt=CFL*(1./40.)/c;

if problem_type==1  %impulsively loaded bar
    period = sqrt(16.*bar_length*bar_length*density/E);
    TipForce=10.;
    D=TipForce*bar_length/(area*E);
    M=4.*D/period;
end
if problem_type==2 %simple oscillator
    Mass=10000.;
    period=2.*3.14159/sqrt(E/Mass);
    v0=1.0;
    Amp=v0/(2.*3.14159/period);
end
if problem_type==3
    period=0.5;
end
if problem_type==4
    period=0.5;
end

tfinal=1.0*period;

ip=1;
xp(ip)=dx/(2.*PPC);

% create particles
while xp(ip)+dx/PPC < bar_length
    ip=ip+1;
    xp(ip)=xp(ip-1)+dx/PPC;
end

NP=ip  % Particle count

% initialize other particle variables
for ip=1:NP
    vol(ip)=volp;
    mp(ip)=volp*density;
    vp(ip)=0.;
    dp(ip)=0.;
    sigp(ip)=0.;
    Fp(ip)=1.;
    exFp(ip)=0.;
end

if problem_type==1
   exFp(NP)=TipForce;
end
if problem_type==2
   mp(NP)=Mass;
   vp(NP)=v0;
end

% create array of nodal locations, only used in plotting
for(ig=1:NN)
    xg(ig)=(ig-1)*dx;
    dug(ig)=0.;
end

% set up BCs
numBCs=1;
if problem_type==3
    numBCs=0;
end
if problem_type==4
    numBCs=0;
end

BCNode(1)=1;
BCNode(2)=NN;
BCValue(1)=0.;
BCValue(2)=1.;

t=0.;
tstep=0;

while t<tfinal
    tstep=tstep+1;
    t=t+dt;
    if problem_type==3
        exFp=mp.*t;
    end
    if problem_type==4
        exFp=mp;
    end

    % initialize arrays to be zero
    for ig=1:NN
        mg(ig)=0.;
        vg(ig)=0.;
        vg_nobc(ig)=0.;
        ag(ig)=0.;
        Feg(ig)=0.;
        Fig(ig)=0.;
    end

    % project particle data to grid
    for ip=1:NP
        [nodes,Ss]=findNodesAndWeights(xp(ip),dx);
        for ig=1:2
            mg(nodes(ig))=mg(nodes(ig))+mp(ip)*Ss(ig);
            vg(nodes(ig))=vg(nodes(ig))+mp(ip)*vp(ip)*Ss(ig);
            Feg(nodes(ig))=Feg(nodes(ig))+exFp(ip)*Ss(ig);
        end
    end

    % normalize by the mass
    vg=vg./mg;
    vg_nobc=vg;
    
    % set velocity BC
    for ibc=1:numBCs
      vg(BCNode(ibc))=BCValue(ibc);
    end
    
    %compute particle stress
    [sigp,vol,Fp]=computeStressFromVelocity(xp,dx,dt,vg,E,Fp,volp,NP);

    %compute internal force
    for ip=1:NP
        [nodes,Gs]=findNodesAndWeightGradients(xp(ip),dx);
        for ig=1:2
            Fig(nodes(ig))=Fig(nodes(ig))-(Gs(ig)/dx)*sigp(ip)*vol(ip);
        end
    end

    %compute acceleration
    ag=(Fig+Feg)./mg;
    vg_new=vg+ag.*dt;

    %set BCs again
    for ibc=1:numBCs
      vg_new(BCNode(ibc))=BCValue(ibc);
    end
    
    for ig=1:NN
      ag(ig)=(vg_new(ig)-vg_nobc(ig))/dt;
      dug(ig)=dug(ig)+vg(ig)*dt;
    end
    
    %project changes back to particles
    for ip=1:NP
        [nodes,Ss]=findNodesAndWeights(xp(ip),dx);
        dvp=0.;
        dxp=0.;
        for ig=1:2
            dvp=dvp+ag(nodes(ig))*dt*Ss(ig);
            dxp=dxp+vg(nodes(ig))*dt*Ss(ig);
        end
        vp(ip)=vp(ip)+dvp;
        xp(ip)=xp(ip)+dxp;
        dp(ip)=dp(ip)+dxp;
    end
    
    DX_tip(tstep)=dp(NP);
    T=t;%-dt;

    KE(tstep)=0.;SE(tstep)=0.;
    for ip=1:NP
        KE(tstep)=KE(tstep)+.5*mp(ip)*vp(ip)*vp(ip);
        SE(tstep)=SE(tstep)+.5*sigp(ip)*(Fp(ip)-1.)*vol(ip);
        TE(tstep)=KE(tstep)+SE(tstep);
    end

    if problem_type==1
        if(T<=period/2.)
            Exact_tip(tstep)=M*T;
        else
            Exact_tip(tstep)=4.*D-M*T;
        end
    end
    if problem_type==2
       Exact_tip(tstep)=Amp*sin(2.*3.14159*T/period);
    end
    if problem_type==3
       Exact_tip(tstep)=(1./6.)*T*T*T; 
    end
    if problem_type==4
       Exact_tip(tstep)=(1./2.)*T*T; 
    end
    
    TIME(tstep)=t;
end

close all;
subplot(2,1,1),plot(TIME,DX_tip,'bx');
hold on;
subplot(2,1,1),plot(TIME,Exact_tip,'r-');
%subplot(2,1,2),plot(TIME,KE,'g*');
%hold on;
%subplot(2,1,2),plot(TIME,SE,'r+');
subplot(2,1,2),plot(TIME,TE,'b-');

E_err=TE(1)-TE(tstep)

% compute error
err=abs(DX_tip(tstep)-Exact_tip(tstep))

length(TIME)

return;
