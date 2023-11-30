function Out = GetGi(r,RBrain,RSkull,RScalp,Rair,sigmaBrain,sigmaSkull,sigmaScalp,sigmaAir,iSum,norm_PosCDInit)
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
b = norm_PosCDInit;

%% Maple code
G1i = 0.1e1 ./ pi ./ s1 .* (b .^ (iSum - 1) ./ r .^ (iSum + 1) + (iSum + 1) .* b .^ (iSum - 1) .* ((((s2 + s3) .* iSum + s3) .* ((s3 + s4) .* iSum + s4)...
    .* R3 .^ (2 .* iSum + 1) .* R2 .^ (2 .* iSum + 1) + R2 .^ (4 .* iSum + 2) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1)) .* (-s2 + s1) .* R1 .^ (-2 .* iSum - 1)...
    + (((s2 + s3) .* iSum + s2) .* (s3 - s4) .* R2 .^ (2 .* iSum + 1) + (-s3 + s2) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1)) .* (iSum .* (s1 + s2) + s1))...
    ./ (iSum .* (iSum + 1) .* (((s2 + s3) .* iSum + s2) .* (s3 - s4) .* R2 .^ (2 .* iSum + 1) + (-s3 + s2) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1)) .* (-s2 + s1)...
    .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2) .* (((s2 + s3) .* iSum + s3) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1) .* R2 .^ (2 .* iSum + 1)...
    + R2 .^ (4 .* iSum + 2) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1))) .* r .^ iSum) ./ 0.4e1;


G2i = 0.1e1 ./ pi ./ s2 .* (0.2e1 .* (iSum + 0.1e1 ./ 0.2e1) .* (R2 .^ 2) .* (((s2 + s3) .* iSum + s3) .* ((s3 + s4) .* iSum + s4) .* R2 .^ (2 .* iSum - 1)...
    .* R3 .^ (2 .* iSum + 1) + R2 .^ (4 .* iSum) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1)) .* s2 .* (b .^ (iSum - 1)) ./ (((-s3 + s2) .* (-s2 + s1)...
    .* (iSum + 1) .* iSum .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2) .* ((s2 + s3) .* iSum + s3) .* R2 .^ (2 .* iSum + 1)) .* ((s3 + s4) .* iSum + s4)...
    .* R3 .^ (2 .* iSum + 1) + iSum .* (iSum + 1) .* (s3 - s4) .* (((s2 + s3) .* iSum + s2) .* R2 .^ (2 .* iSum + 1) .* (-s2 + s1) .* R1 .^ (2 .* iSum + 1)...
    +(iSum .* (s1 + s2) + s2) .* (-s3 + s2) .* R2 .^ (4 .* iSum + 2))) ./ (r .^ (iSum + 1)) + 0.2e1 .* (iSum + 0.1e1 ./ 0.2e1) .* (iSum + 1) .* (b .^ (iSum - 1))...
    .* (((s2 + s3) .* iSum + s2) .* (s3 - s4) .* R2 .^ (2 .* iSum + 1) + (-s3 + s2) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1)) .* s2 ./ (iSum .* (iSum + 1) ...
    .* (((s2 + s3) .* iSum + s2) .* (s3 - s4) .* R2 .^ (2 .* iSum + 1) + (-s3 + s2) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1)) .* (-s2 + s1)...
    .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2) .* (((s2 + s3) .* iSum + s3) .* ((s3 + s4) .* iSum + s4) .* R3 .^ (2 .* iSum + 1) .* R2 .^ (2 .* iSum + 1)...
    + R2 .^ (4 .* iSum + 2) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1))) .* (r .^ iSum)) ./ 0.4e1;

G3i = 0.1e1 ./ pi ./ s3 .* (0.4e1 .* (iSum + 0.1e1 ./ 0.2e1) .^ 2 .* R2 .* s3 .* s2 .* ((s3 + s4) .* iSum + s4) .* (b .^ (iSum - 1)) .* R3 ./ (iSum .* (iSum + 1)...
    .* (((s2 + s3) .* iSum + s2) .* R2 .* (s3 - s4) .* R3 .^ (-2 .* iSum) + (-s3 + s2) .* R2 .^ (-2 .* iSum) .* ((s3 + s4) .* iSum + s4) .* R3) .* (-s2 + s1)...
    .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2) .* (R3 .^ (-2 .* iSum) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1) .* R2 .^ (2 .* iSum + 2)...
    + ((s2 + s3) .* iSum + s3) .* R2 .* ((s3 + s4) .* iSum + s4) .* R3)) ./ (r .^ (iSum + 1)) + (R3 .^ (-2 .* iSum) .* s3 .* (2 .* iSum + 1) .^ 2 .* b .^ (iSum - 1)...
    .* s2 .* (s3 - s4) .* (iSum + 1) .* R2 ./ (iSum .* (iSum + 1) .* (((s2 + s3) .* iSum + s2) .* R2 .* (s3 - s4) .* R3 .^ (-2 .* iSum) + (-s3 + s2)...
    .* R2 .^ (-2 .* iSum) .* ((s3 + s4) .* iSum + s4) .* R3) .* (-s2 + s1) .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2) .* (R3 .^ (-2 .* iSum)...
    .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1) .* R2 .^ (2 .* iSum + 2) + ((s2 + s3) .* iSum + s3) .* R2 .* ((s3 + s4) .* iSum + s4) .* R3)) .* r .^ iSum)) ./ 0.4e1;

G4i = 0.1e1 ./ pi .* R3 .* R2 .* s3 .* ((2 .* iSum + 1) .^ 3) .* (b .^ (iSum - 1)) .* s2 ./ (iSum .* (iSum + 1) .* (((s2 + s3) .* iSum + s2) .* R2 .* (s3 - s4)...
    .* R3 .^ (-2 .* iSum) + (-s3 + s2) .* R2 .^ (-2 .* iSum) .* ((s3 + s4) .* iSum + s4) .* R3) .* (-s2 + s1) .* R1 .^ (2 .* iSum + 1) + (iSum .* (s1 + s2) + s2)...
    .* (R3 .^ (-2 .* iSum) .* iSum .* (s3 - s4) .* (-s3 + s2) .* (iSum + 1) .* R2 .^ (2 .* iSum + 2) + ((s2 + s3) .* iSum + s3) .* R2 .* ((s3 + s4) .* iSum + s4)...
    .* R3)) ./ (r .^ (iSum + 1)) ./ 0.4e1;



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
            

end