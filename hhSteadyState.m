function hhSteadyState()

simTime = 100000;
gK = 36;
gNa = 120;
gL = 0.3;
EK = -12;
ENa = 115;
EL = 10.6;
Vrest = -70;
Cm = 1.0;

V = [];
V(1) = 0;
I = [5*ones(50), zeros(100000-50)];
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

m(1) = alpham/(alpham + betam);
n(1) = alphan/(alphan + betan);
h(1) = alphah/(alphah + betah);

INa(1) = m(1)^3*gNa*h(1)*(V(1)-ENa)
IK(1) = n(1)^4*gK*(V(1)-EK)
IL(1) = gL*(V(1)-EL)

Iion(1) = I(1) - IK(1) - INa(1) - IL(1)

for tt=2:simTime+1,
    alpham = 0.1*((25-V(tt-1))/(exp((25-V(tt-1))/10)-1));
    betam = 4*exp(-V(tt-1)/18);
    alphan = .01*((10-V(tt-1))/(exp((10-V(tt-1))/10)-1));
    betan = .125*exp(-V(tt-1)/80);
    alphah = .07*exp(-V(tt-1)/20);
    betah = 1/(exp((30-V(tt-1))/10)+1);
    
    deltaV = Iion(tt-1)/Cm;
    deltam = alpham*(1-m(tt-1))-betam*m(tt-1);
    deltan = alphan*(1-n(tt-1))-betan*n(tt-1);
    deltah = alphah*(1-h(tt-1))-betah*h(tt-1);
    
    V(tt) = V(tt-1) + .001*deltaV;
    m(tt) = m(tt-1) + .001*deltam;
    n(tt) = n(tt-1) + .001*deltan;
    h(tt) = h(tt-1) + .001*deltah;
    
    INa(tt) = m(tt)^3*gNa*h(tt)*(V(tt)-ENa);
    IK(tt) = n(tt)^4*gK*(V(tt)-EK);
    IL(tt) = gL*(V(tt)-EL);
    Iion(tt) = I(tt) - IK(tt) - INa(tt) - IL(tt);
end

plot(0:1:simTime, V-70);
figure
plot(0:1:simTime, m);
figure
plot(0:1:simTime, n);
figure
plot(0:1:simTime, h);
figure
plot(0:1:simTime, Iion);

end