
POI = [0,1,0];
POI = POI/norm(POI);
dpOI = [1,1,0];
dpOI = dpOI/norm(dpOI);
DirusOI = [1,1,0];
DirusOI = DirusOI/norm(DirusOI);
ndpOI_POI = cross(POI,dpOI);
% if POI and dp on same vector
if ~any(ndpOI_POI)
    ndpOI_POI = cross(POIs(iPOI,:),CSource(idx_dpOI(idp),:));
    if ~any(ndpOI_POI)
        ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*dpOI';
        ndpOI_POI = ndpOI_POI';
        if ~any(ndpOI_POI)
            ndpOI_POI = [0,0,1;0,1,0;-1,0,0]*CSource(idx_dpOI(idp),:)';
            ndpOI_POI = ndpOI_POI';
        end
    end
end
ndpOI_POI = ndpOI_POI./norm(ndpOI_POI);

A = calc_RM(ndpOI_POI,[1,0,0]);
if any(isnan(A(:)))
    A = eye(3);
    if any(ndpOI_POI-ndpDB_POIDB)
        A = A.*(-2*double((ndpOI_POI-ndpDB_POIDB)~=0)+1);
    end
    
end

% rotate vibration direction
DirusOI_rot = (A*DirusOI')'
dpOInew = (A*dpOI')'


Aextra = calc_RM(dpOInew,[0,0,1])

DirusOI_rot = Aextra*DirusOI_rot'
dpOInew = Aextra*dpOInew'


function A = calc_RM(r,t)
        r = r./norm(r); t = t./norm(t);
        v = cross(r,t);
        s = norm(v);
        c = dot(r,t);
        Vx = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
        A = eye(3)+Vx+Vx^2.*(1-c)/s^2;
        
    end