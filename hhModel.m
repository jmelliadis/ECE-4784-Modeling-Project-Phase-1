function hhModel(impulseTime)

simTime = 100000; %total simulation time 100 ms with 1/1000 ms steps
gK = 36;
gNa = 120;
gL = 0.3;
EK = -12;
ENa = 115;
EL = 10.6;
Vrest = -70;
Cm = 1.0;

V = [];
V(1) = 0; %assume resting potential is 0 and then shift total response array down by Vrest later
I = [50*ones(1,impulseTime*1000), zeros(1,100001-impulseTime*1000)]; %determine input current based off of impulseTime Parameter
Iion = [];
INa = [];
IK = [];
IL = [];

m = [];
n = [];
h = [];
alpham = 0.1*((25-V(1))/(exp((25-V(1))/10)-1));
betam = 4*exp(-V(1)/18);
alphan = 0.01*((10-V(1))/(exp((10-V(1))/10)-1));
betan = .125*exp(-V(1)/80);
alphah = .07*exp(-V(1)/20);
betah = 1/(exp((30-V(1))/10)+1);

m(1) = alpham/(alpham + betam); %initialize m
n(1) = alphan/(alphan + betan); %initialize n
h(1) = alphah/(alphah + betah); %initialize h

INa(1) = m(1)^3*gNa*h(1)*(V(1)-ENa);
IK(1) = n(1)^4*gK*(V(1)-EK);
IL(1) = gL*(V(1)-EL);

Iion(1) = I(1) - IK(1) - INa(1) - IL(1);

for tt=2:simTime+1, %loop through all values of time
    
    %calculate alpha and beta for m, n, h based off of Vm
    alpham = 0.1*((25-V(tt-1))/(exp((25-V(tt-1))/10)-1));
    betam = 4*exp(-V(tt-1)/18);
    alphan = .01*((10-V(tt-1))/(exp((10-V(tt-1))/10)-1));
    betan = .125*exp(-V(tt-1)/80);
    alphah = .07*exp(-V(tt-1)/20);
    betah = 1/(exp((30-V(tt-1))/10)+1);
    
    %calculate deltas for each diff. eq. (Vm, m, n, h)
    deltaV = Iion(tt-1)/Cm;
    deltam = alpham*(1-m(tt-1))-betam*m(tt-1);
    deltan = alphan*(1-n(tt-1))-betan*n(tt-1);
    deltah = alphah*(1-h(tt-1))-betah*h(tt-1);
    
    %increment each changing variable by the delta value
    V(tt) = V(tt-1) + .001*deltaV; %h=1/1000 (ms)
    m(tt) = m(tt-1) + .001*deltam; %h=1/1000 (ms)
    n(tt) = n(tt-1) + .001*deltan; %h=1/1000 (ms)
    h(tt) = h(tt-1) + .001*deltah; %h=1/1000 (ms)
    
    %calculate new current values for Na, K, L, and Iion
    INa(tt) = m(tt)^3*gNa*h(tt)*(V(tt)-ENa);
    IK(tt) = n(tt)^4*gK*(V(tt)-EK);
    IL(tt) = gL*(V(tt)-EL);
    Iion(tt) = I(tt) - IK(tt) - INa(tt) - IL(tt);
end

%plot Vm versus time
plot(0:1:simTime, V-70);

%plot gNa and gK versus time
figure
plot(0:1:simTime, m.^3*gNa.*h, 0:1:simTime, n.^4*gK)
%figure
%plot(0:1:simTime, m);
%figure
%plot(0:1:simTime, n);
%figure
%plot(0:1:simTime, h);
%figure
%plot(0:1:simTime, Iion);

end