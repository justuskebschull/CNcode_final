function [a,b]=getprobability(prob1,prob2,x1,y1)
syms P D1 D2 x y positive

eqn=y*(D1-(1-x)*P)/x+(1-y)*P==D2;
S=solve(eqn,P,'Real',true);
f=matlabFunction(S);
b=f(prob1,prob2,x1,y1);

syms P D1 D2 x y positive
eqn=y*P+(1-y)*(D1-x*P)/(1-x)==D2;
S=solve(eqn,P,'Real',true);
f=matlabFunction(S);
a=f(prob1,prob2,x1,y1);
