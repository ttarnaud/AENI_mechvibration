function Out = GetGi(r,RBrain,RSkull,RScalp,Rair,sigmaBrain,sigmaSkull,sigmaScalp,sigmaAir,ix,posdip)
% calculate time independent components for 4 sphereical model, r=position
% poi can vary throughout whole model not necessarily on oute pshere
if ~isrow(sigmaBrain) || ~isrow(sigmaScalp) || ~isrow(sigmaSkull) || ~isrow(sigmaAir)
    error('sigmas wrong format')
elseif ~iscolumn(r)
    error('r wrong format')
end

R1 = RBrain;
R2 = RSkull;
R3 = RScalp;
R4 = Rair;
s1 = sigmaBrain;
s2 = sigmaSkull;
s3 = sigmaScalp;
s4 = sigmaAir;
b = posdip;
p=1;


A1 = b .^ (ix - 1) * p;

A2 = 0.2e1 .* (ix + 0.1e1 ./ 0.2e1) .* s2 .* R2 .^ 2 .* b .^ (ix - 0.1e1) .* (-((s2 + s3) .* ix + s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (0.2e1 .* ix - 0.1e1) .* R3 .^ (0.2e1 .* ix + 0.1e1) + R2 .^ (0.4e1 .* ix) .* ix .* (s2 - s3) .* (ix + 0.1e1) .* (s4 - s3)) .* p ./ (-((s4 + s3) .* ix + s4) .* (ix .* (s2 - s3) .* (s1 - s2) .* (ix + 0.1e1) .* R1 .^ (0.2e1 .* ix + 0.1e1) + ((s2 + s3) .* ix + s3) .* R2 .^ (0.2e1 .* ix + 0.1e1) .* ((s1 + s2) .* ix + s2)) .* R3 .^ (0.2e1 .* ix + 0.1e1) + (ix + 0.1e1) .* ix .* (R2 .^ (0.2e1 .* ix + 0.1e1) .* ((s2 + s3) .* ix + s2) .* (s1 - s2) .* R1 .^ (0.2e1 .* ix + 0.1e1) + (s2 - s3) .* ((s1 + s2) .* ix + s2) .* R2 .^ (0.4e1 .* ix + 0.2e1)) .* (s4 - s3));

A3 = 0.4e1 .* (ix + 0.1e1 ./ 0.2e1) .^ 2 .* s2 .* R3 .* R2 .* ((s4 + s3) .* ix + s4) .* s3 .* b .^ (ix - 0.1e1) .* p ./ ((-R2 .* ((s2 + s3) .* ix + s2) .* (s4 - s3) .* R3 .^ (-0.2e1 .* ix) + R3 .* (s2 - s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (-0.2e1 .* ix)) .* (ix + 0.1e1) .* ix .* (s1 - s2) .* R1 .^ (0.2e1 .* ix + 0.1e1) - ((s1 + s2) .* ix + s2) .* (R3 .^ (-0.2e1 .* ix) .* ix .* (s2 - s3) .* (ix + 0.1e1) .* (s4 - s3) .* R2 .^ (0.2e1 .* ix + 0.2e1) - ((s2 + s3) .* ix + s3) .* R3 .* R2 .* ((s4 + s3) .* ix + s4)));

A4 = R2 .* R3 .* s3 .* (2 .* ix + 1) .^ 3 .* s4 .* b .^ (ix - 1) .* p .* s2 ./ ((-R2 .* ((s2 + s3) .* ix + s2) .* (s4 - s3) .* R3 .^ (-2 .* ix) + R3 .* (s2 - s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (-2 .* ix)) .* (ix + 1) .* ix .* (s1 - s2) .* R1 .^ (2 .* ix + 1) - ((s1 + s2) .* ix + s2) .* (R3 .^ (-2 .* ix) .* ix .* (s2 - s3) .* (ix + 1) .* (s4 - s3) .* R2 .^ (2 .* ix + 2) - ((s2 + s3) .* ix + s3) .* R3 .* R2 .* ((s4 + s3) .* ix + s4)));

B1 = ((s1 - s2) .* (R3 .^ (2 .* ix + 1) .* ((s2 + s3) .* ix + s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (2 .* ix + 1) - ix .* R2 .^ (4 .* ix + 2) .* (s2 - s3) .* (ix + 1) .* (s4 - s3)) .* R1 .^ (-2 .* ix - 1) + ((s1 + s2) .* ix + s1) .* (-((s2 + s3) .* ix + s2) .* (s4 - s3) .* R2 .^ (2 .* ix + 1) + R3 .^ (2 .* ix + 1) .* (s2 - s3) .* ((s4 + s3) .* ix + s4))) .* (ix + 1) .* b .^ (ix - 1) .* p ./ ((ix + 1) .* ix .* (s1 - s2) .* (-((s2 + s3) .* ix + s2) .* (s4 - s3) .* R2 .^ (2 .* ix + 1) + R3 .^ (2 .* ix + 1) .* (s2 - s3) .* ((s4 + s3) .* ix + s4)) .* R1 .^ (2 .* ix + 1) + ((s1 + s2) .* ix + s2) .* (R3 .^ (2 .* ix + 1) .* ((s2 + s3) .* ix + s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (2 .* ix + 1) - ix .* R2 .^ (4 .* ix + 2) .* (s2 - s3) .* (ix + 1) .* (s4 - s3)));

B2 = 0.2e1 .* s2 .* (ix + 1) .* (-((s2 + s3) .* ix + s2) .* (s4 - s3) .* R2 .^ (2 .* ix + 1) + R3 .^ (2 .* ix + 1) .* (s2 - s3) .* ((s4 + s3) .* ix + s4)) .* (ix + 0.1e1 ./ 0.2e1) .* b .^ (ix - 1) .* p ./ ((ix + 1) .* ix .* (s1 - s2) .* (-((s2 + s3) .* ix + s2) .* (s4 - s3) .* R2 .^ (2 .* ix + 1) + R3 .^ (2 .* ix + 1) .* (s2 - s3) .* ((s4 + s3) .* ix + s4)) .* R1 .^ (2 .* ix + 1) + ((s1 + s2) .* ix + s2) .* (R3 .^ (2 .* ix + 1) .* ((s2 + s3) .* ix + s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (2 .* ix + 1) - ix .* R2 .^ (4 .* ix + 2) .* (s2 - s3) .* (ix + 1) .* (s4 - s3)));

B3 = -0.4e1 .* s2 .* R2 .* (ix + 1) .* (R3 .^ (-2 .* ix)) .* (ix + 0.1e1 ./ 0.2e1) .^ 2 .* s3 .* (s4 - s3) .* b .^ (ix - 1) .* p ./ ((-R2 .* ((s2 + s3) .* ix + s2) .* (s4 - s3) .* R3 .^ (-2 .* ix) + R3 .* (s2 - s3) .* ((s4 + s3) .* ix + s4) .* R2 .^ (-2 .* ix)) .* (ix + 1) .* ix .* (s1 - s2) .* R1 .^ (2 .* ix + 1) - ((s1 + s2) .* ix + s2) .* (R3 .^ (-2 .* ix) .* ix .* (s2 - s3) .* (ix + 1) .* (s4 - s3) .* R2 .^ (2 .* ix + 2) - ((s2 + s3) .* ix + s3) .* R3 .* R2 .* ((s4 + s3) .* ix + s4)));

G1i = 1./(4*pi*s1).*(A1./r.^(ix+1) + B1.*r.^ix);
G2i = 1./(4*pi*s2).*(A2./r.^(ix+1) + B2.*r.^ix);
G3i = 1./(4*pi*s3).*(A3./r.^(ix+1) + B3.*r.^ix);
G4i = double(s4~=0)*1./(4*pi*(s4+double(s4==0))).*A4./r.^(ix+1);

%% select correct Gi
Out = zeros(length(r),length(s1));
Gis = {G1i,G2i,G3i,G4i};
Rvals = double(b<r) + double(r>R1) + double(r>R2) + double(r>R3);
if any(Rvals==0)
    error('false POI pos')
end

%sigmas can be complex values
for i = 1:length(Rvals)
    
    Gi_temp =  Gis{Rvals(i)};
    rGi_temp = real(Gi_temp(i,:));imGi_temp = imag(Gi_temp(i,:));
    rGi_temp(isnan(rGi_temp)) = 0; rGi_temp(isinf(rGi_temp)) = 0;
    imGi_temp(isnan(imGi_temp)) = 0; imGi_temp(isinf(imGi_temp)) = 0;
    Gi_tempi = complex(rGi_temp,imGi_temp);
    Out(i,:) = Gi_tempi;
    
end

