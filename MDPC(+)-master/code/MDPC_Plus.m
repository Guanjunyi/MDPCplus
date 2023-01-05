% Please kindly cite the paper Junyi Guan, Sheng li, Xiongxiong He, and Jiajia Chen 
%"Clustering by fast detection of main density peaks within a peak digraph" 
% Information Sciences,2023

% The code was written by Junyi Guan in 2022.
function [CL,NC,centers,runtime] = MDPC_Plus(data,k)
%% parameters
[n,dim]  = size(data); % n: the number of data points; d:dimensions
alpha = 0.5;
lambda = 2;
if nargin<2
    k = round(sqrt(n)); %k: the number of neighbors
end
k_b = min(2*floor(log(n)),k);
tic;

%% Fast KNN based on Kd-tree (when dimension is not large than 10)
if dim<=11;[knn,knn_dist] = knnsearch(data,data,'k',max(k,k_b));
else;DIST = squareform(pdist(data));
    [knn_dist,knn] = sort(DIST,2);
end

%% KNN-based Density
rho=sum(exp(-knn_dist(:,2:k)).^2,2)';

%% Allocation
[~,OrdRho]=sort(rho,'descend');

%% The identification of local density peaks and sub-clusters, and center-association degree
peaks = [];
for i=1:n
    point = OrdRho(i);
    all_neigh=knn(point,2:k);
    all_neigh_dist = knn_dist(point,2:k);
    big_neigh = all_neigh(rho(all_neigh)>rho(point));  
    
    bb = (rho(all_neigh)-rho(point))./rho(all_neigh);
    dd = all_neigh_dist./all_neigh_dist(end);
    bb = bb(bb>0);
    dd = dd(bb>0);
    bb=(bb-min(bb))./(max(bb)-min(bb));
    bb(isnan(bb))=0;
    dd=(dd-min(dd))./(max(dd)-min(dd));
    dd(isnan(dd))=0;
    bb_dd = alpha*bb + (1-alpha)*dd;
    
    if length(bb)>0
        best_big_index = find(bb_dd == min(bb_dd));
        best_big = big_neigh(best_big_index(1));
        neigh = best_big;
    else
         neigh = [];
    end
    if neigh
        phi(point) = (abs(rho(neigh)-rho(point))/rho(neigh))^lambda;
        pn(point)=neigh; %pn:parent node
    else
        phi(point) = 0;
        peaks = [point peaks];
    end
end

for i=1:n
    if pn(OrdRho(i))
        denisty_deviation_cost_to_peak(OrdRho(i)) = denisty_deviation_cost_to_peak(pn(OrdRho(i))) + phi(OrdRho(i));
    else
        denisty_deviation_cost_to_peak (OrdRho(i)) = 0;
    end
end

%% label initialization
for i=1:n
    sub_l(i)=-1; %% sub_l:sub-cluster label
end 

n_p = length(peaks);%% n_p:Number of peaks
sub_l(peaks) = (1:n_p);
for i=1:n
    if (sub_l(OrdRho(i))==-1)
        sub_l(OrdRho(i))=sub_l(pn(OrdRho(i)));
    end
end
%% edges matrix
rho_peaks= rho(peaks); [~,OrdRho_peaks]=sort(rho_peaks,'descend');
edges = Inf*ones(n_p,n_p);% edges: edge matrix between subtrees(graph theory)
for i=1:n
    BB = rho(i)*ones(1,k_b-1);
    CC = rho(knn(i,2:k_b));
    denisty_deviation_cost_of_link_set = (abs(BB-CC)./max(BB,CC)).^lambda;
    for j = 2:k_b
        jj = knn(i,j);
        AA = denisty_deviation_cost_to_peak(i)+denisty_deviation_cost_to_peak(jj); 
        if sub_l(i)~=sub_l(jj) & edges(sub_l(i),sub_l(jj))> AA
            if find(knn(jj,2:k)==i)
                denisty_deviation_cost_of_link = denisty_deviation_cost_of_link_set(j-1);
                edges(sub_l(i),sub_l(jj)) = AA+denisty_deviation_cost_of_link;
                edges(sub_l(jj),sub_l(i)) = AA+denisty_deviation_cost_of_link;
            end
        end
    end
end

%% PEAK-GRAPH
G=sparse(n_p,n_p);
for i=1:n_p
    for j = 1:n_p
        if edges(i,j) ~= Inf
            G(i,j) = edges(i,j);
            G(j,i) = edges(i,j);
        end
    end
end
[dim ~] = dijkstra(G,(1:n_p));
delta_peaks = Inf*ones(n_p,1);
pn_peaks  = -1*ones(n_p,1);
for i = 2:n_p
    ii = OrdRho_peaks(i);
    for j = 1:i-1
        jj = OrdRho_peaks(j);
        if delta_peaks(ii) > dim(ii,jj)
            delta_peaks(ii) = dim(ii,jj);
            pn_peaks(ii) = jj;
        end
    end
end

%% delta of peaks
delta = zeros(n,1);
delta(peaks) =  delta_peaks;
if n_p > 1
    Must_C = length(find(delta==Inf));
    if Must_C>1
        delta(delta==Inf) = max(delta(delta~=Inf))*1.2;
        delta_peaks = delta(peaks);
    else
        delta(delta==Inf) = max(delta(delta~=Inf)*1.2);
        delta_peaks = delta(peaks);
    end
end
time1 = toc;

%% center selection
figure;
plot(rho(peaks), delta(peaks),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
grid on;
axis([0 max(rho(peaks)) 0 max(delta(peaks))]);
title ('Decision Graph','FontSize',15.0);
xlabel ('\rho');
ylabel ('\delta');
rect = getrect;
rhomin = rect(1);
deltamin=rect(2);

%% center confirm
tic;
NC=0;
for i=1:n_p
    Cl_peaks(i)=-1; % clutser label of peaks
end
for i=1:n_p
    if rho_peaks(i)>rhomin & delta_peaks(i) > deltamin
        NC=NC+1;
        Cl_peaks(i)=NC;
        icl(NC)=i;
    end
end
for i=1:n_p
    if (Cl_peaks(OrdRho_peaks(i))==-1)
        Cl_peaks(OrdRho_peaks(i))=Cl_peaks(pn_peaks(OrdRho_peaks(i)));
    end
end

%% allocation
for i=1:n_p
    CL(sub_l== i) = Cl_peaks(i);
end
centers = peaks(icl);
time2 = toc;
runtime = time1+time2;


