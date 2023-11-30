function Out = calcfitime(CSource_Sel,CSink_Sel,iSum,nPOI,norm_PosCDInit)
PosCD = (CSource_Sel+CSink_Sel)./2; % Position centre of dipole
DM = 1e-6*(CSource_Sel-CSink_Sel); %dipole moment I = 1µA
if any(PosCD)
    nPosCD = PosCD/norm(PosCD); % normalized pos dipole
else
    nPosCD = CSource_Sel/norm(CSource_Sel); % normalized pos dipole
end
CosT = (nPosCD*nPOI')'; % cosine of angle between dipole and POI
CosT = min(CosT,ones(size(CosT)));CosT = max(CosT,-ones(size(CosT))); % ensure no numerical errors
nPosCD = repmat(nPosCD,size(nPOI,1),1);
ntPosCD = cross(cross(nPosCD,nPOI,2),nPosCD,2); % tangential of dipole [Nx3]
ntPosCD = ntPosCD./vecnorm(ntPosCD ,2,2);
ntPosCD(isnan(ntPosCD))=0;
Pi = legendre(iSum,CosT);
Out = (norm(PosCD(1,:))/norm_PosCDInit).^(iSum-1).*(iSum.*nPosCD.*Pi(1,:)'+(-1).*ntPosCD.*Pi(2,:)')*DM';
end