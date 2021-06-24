gEph = 0.035; αG= 2.53; ainm = 2.46e-10
ainBohr = ainm/(5.29177e-11)
Nagib =sqrt(3)*ainBohr^2*αG^2*gEph^2/16
Eph = 0.7692307692307692
LD = 32.30769230769231
LG = 37.1213101534543
Dcut = LG;

function ImS(ω::Real)
  abs(-1 + ω)< 0.7692307692307692 ?  0 :  Nagib*abs(ω- sign(ω- 1)*Eph);
end

function ReS(ω)
0.0029194846882482015*(9.627235691523993 + (-0.7692307692307692 
+ ω)*log((1.7692307692307692 - ω)^2) + 
(0.7692307692307692 - ω)*log((32.30769230769231 - ω)^2) - 
0.7692307692307692*log((-0.23076923076923084 + ω)^2) - ω*log((-0.23076923076923084 + ω)^2) + 
1.5384615384615383*log((0.7692307692307692 + ω)^2) + 
2*ω*log((0.7692307692307692 + ω)^2) - 0.7692307692307692*log((37.1213101534543 + ω)^2) - 
ω*log((37.1213101534543 + ω)^2)) - log(abs(-32.30769230769231 + ω))/(
100000000000000000000*π) + log(abs(8.5 + ω))/(
100000000000000000000*π);
end

A(ω)= ω - ReS(ω)
B(ω) = ω - ReS(ω)
Ia(ω) = ImS(ω)
Ib(ω) = ImS(ω)

D1(ω, ω1) = 2*((A(ω1) - B(ω+ω1))^2 + (Ia(ω1) - Ib(ω+ω1))^2)*((A(ω1) - B(ω+ω1))^2 + (Ia(ω1) + Ib(ω+ω1))^2)
D2(ω, ω1) = 2*((A(ω1) + B(ω+ω1))^2 + (Ia(ω1) - Ib(ω+ω1))^2)*((A(ω1) + B(ω+ω1))^2 + (Ia(ω1) + Ib(ω+ω1))^2)

function pri(ω, ω1)
    A = JJDFTX.A(ω1)
    Ia = JJDFTX.Ia(ω1)
    B = JJDFTX.B(ω+ω1)
    Ib = JJDFTX.Ib(ω+ω1)
    D1 = JJDFTX.D1(ω, ω1)
    D2 = JJDFTX.D2(ω, ω1)

    Pribrojnik1 = 2*Ib/D1*(A^3 - 2*A^2*B - 2*B*Ia^2 + A*(B^2 + Ib^2 + Ia^2))*(2*atan(A/Ia) - atan((A + Dcut)/Ia) - atan((A - Dcut)/Ia));
    Pribrojnik2 = 2*Ia/D1*(B^3 - 2*B^2*A - 2*A*Ib^2 + B*(A^2 + Ia^2 + Ib^2))*(2*atan(A/Ib) - atan((A + Dcut)/Ib) - atan((A - Dcut)/Ib));
    Pribrojnik3 = 2*Ib/D2*(A^3 + 2*A^2*B + 2*B*Ia^2 + A*(B^2 + Ib^2 + Ia^2))*(2*atan(A/Ia) - atan((A + Dcut)/Ia) - atan((A - Dcut)/Ia));
    Pribrojnik4 = 2*Ia/D2*(B^3 + 2*B^2*A + 2*A*Ib^2 + B*(A^2 + Ia^2 + Ib^2))*(2*atan(A/Ib) - atan((A + Dcut)/Ib) - atan((A - Dcut)/Ib));
    Pribrojnik5 = (1/D1 + 1/D2)*Ia*Ib*(A^2 - B^2 + Ia^2 - Ib^2)*(2*log((A^2 + Ia^2)/(B^2 + Ib^2)) - log(((A - Dcut)^2 + Ia^2)/((B - Dcut)^2 + Ib^2)) - log(((A + Dcut)^2 + Ia^2)/((B + Dcut)^2 + Ib^2)));
    return Pribrojnik1+Pribrojnik2+Pribrojnik3+Pribrojnik4+Pribrojnik5
end

kernel(ω, ω1)= 4.0/π^2*pri(ω, ω1)/ω*(-heaviside(ω1 - 1) + heaviside(ω + ω1-1));
