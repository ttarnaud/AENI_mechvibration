function OutP = LegendreArray(narray,x,OOI)
%narray =  degrees of interest
%x = values of interest
%OOI =  orders of interest (porbably [0,1])
OOI = OOI+1;
OutP = zeros(length(narray),length(x),length(OOI));
for i = 1:length(narray)
PI = legendre(narray(i),x);
OutP(i,:,OOI)=PI(OOI,:)';
end
end