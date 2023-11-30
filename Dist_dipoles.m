folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1224\Crop';
files = dir([folder_rawresult,'\R*.mat'])
T = struct2table(files); % convert the struct array to a table
sortedT = sortrows(T, 'date'); % sort the table by 'DOB'
files = table2struct(sortedT);
files = files([1:6,13:18])
%%
mindist_all = [];
for i=1:length(files)
    data_rr =  load(fullfile(folder_rawresult,files(i).name));
    Outall = data_rr.Outall;
    for j=1:size(Outall,2)
    CSource = Outall(1,j).Out.CSource;
    CSink = Outall(1,j).Out.CSink;
    dppos = (CSource+CSink)/2;
    DOI = dppos(1,:);
    mindist_all(i,j) = min(vecnorm(dppos(2:end,:)-DOI,2,2));
    end    
end
%%
diplim = 100
tic
for i=1:length(files)
    data_rr =  load(fullfile(folder_rawresult,files(i).name));
    Outall = data_rr.Outall;
    mindist = nan(1,size(Outall,2));
    mindist2 = nan(1,size(Outall,2));
    meandist = nan(1,size(Outall,2));
    maxdist = nan(1,size(Outall,2));

    for j=1:size(Outall,2)
        
        dist = [];
        CSource = Outall(1,j).Out.CSource;
        CSink = Outall(1,j).Out.CSink;
        dppos = (CSource+CSink)/2;
        if size(CSource,1)>diplim
            dist = [];
            kend = ceil(size(CSource,1)/diplim);
            mindist_intm=[];
            mindist2_intm=[];
            maxdist_intm=[];
            sumdist_intm = 0;
            for k=1:kend
                idx_start = (k-1)*diplim+1;
                if k==kend
                    idx_end = size(CSource,1);
                else
                    idx_end = k*diplim;
                end
                dist_intm = vecnorm(permute(dppos(idx_start:idx_end,:),[2,1])-permute(dppos,[2,3,1]),2,1);
                sumdist_intm = sumdist_intm+sum(dist_intm(:))/2;
                dist_intm = unique([mindist_intm;mindist2_intm;maxdist_intm;dist_intm(:)]); 
                dist_intm = dist_intm(dist_intm>0);              
                 
                
                mindist_intm = dist_intm(1);
                mindist2_intm = dist_intm(2);
                maxdist_intm = dist_intm(end);
            end
            mindist(1,j) = mindist_intm;
            mindist2(1,j) = mindist2_intm;
            meandist(1,j) = sumdist_intm/((size(CSource,1)-1)*size(CSource,1))*2;
            maxdist(1,j) = maxdist_intm;
        else          
            dist = vecnorm(permute(dppos,[2,1])-permute(dppos,[2,3,1]),2,1);
            sumdist_intm = sum(dist(:))/2;
            dist = unique(dist); dist = dist(dist>0);
            mindist(1,j) = dist(1);
            mindist2(1,j) = dist(2);
            meandist(1,j) = sumdist_intm/((size(CSource,1)-1)*size(CSource,1))*2;
            maxdist(1,j) = max(dist);
        end       
        
    end
    mindist_all(i,:) = mindist; 
    mindist2_all(i,:) = mindist2; 
    meandist_all(i,:) = meandist; 
    maxdist_all(i,:) = maxdist; 
end
