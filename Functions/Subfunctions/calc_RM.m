function A = calc_RM(r,t)
%A rotates r onto t
r = r./vecnorm(r,2,2); t = t./vecnorm(t,2,2);
v = (bsxfun(@cross,r',t'))';
s = vecnorm(v,2,2);
c = (r*t');
Vx = [zeros(1,1,size(r,1)),-permute(v(:,3),[3,2,1]),permute(v(:,2),[3,2,1]);permute(v(:,3),[3,2,1]),zeros(1,1,size(r,1)),-permute(v(:,1),[3,2,1]);...
    -permute(v(:,2),[3,2,1]),permute(v(:,1),[3,2,1]),zeros(1,1,size(r,1))];
Vx_cell = mat2cell(Vx,size(Vx,1),size(Vx,2),ones(size(Vx,3),1));
Vx2 = cell2mat(cellfun(@(x) x^2,Vx_cell,'UniformOutput',false));
A = bsxfun(@plus,eye(3),Vx+Vx2.*permute((1-c)./s.^2,[3,2,1]));
end