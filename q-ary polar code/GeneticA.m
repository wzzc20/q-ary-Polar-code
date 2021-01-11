function [sampleM] = GeneticA(N, K, ebno1,ebno2,M,GA_itermax,alpha,beta)

addpath('Decoding_Index/')
addpath('GA/')

max_iter = 50;
max_err = 100;
max_runs = 1e3;

%generate M samples
    if(nchoosek(N,K)<=1)
        comb = combntns(1:N,K);
        combnum = size(comb,1);
        index = zeros(combnum,N);
        for i = 1:combnum
        index(i,comb(i,:)) = ones(1,K);
        end
        samplerows = randperm(combnum);
        samplerows = samplerows(1:M);
        sampleM = index(samplerows,:);
    else
        sampleM = zeros(M,N);
        for i = 1:M
            flag = 1;
            while(flag==1)
            sampleindex = randperm(N);
            sampleindex = sampleindex(1:K);
            index = zeros(1,N);
            index(1,sampleindex) = ones(1,K);
                for j = 1:i
                    flag = flag*(1-all(sampleM(i,:)==index));
                end
                flag = ~flag;
            end
            sampleM(i,:) = index;
        end
        %disp(sampleM)
    end
    % calculate the BER1*BER2
    %BER1 %BER2
    ber1 = zeros(1,M);
    ber2 = zeros(1,M);
    for i = 1:M
        [~, ber1(i)] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(i,:)==1), ebno1, N, K);
        [~, ber2(i)] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(i,:)==1), ebno2, N, K);
    end
    ber = ber1.*ber2;
    [~,order] = sort(ber);
    sampleM = sampleM(order,:);
    
    pdf = 1:M;
    pdf = exp(-alpha*pdf);
    pdf = pdf/sum(pdf);
    cdf = [0,cumsum(pdf)];
    %begin the Genetic algorithm
    for iter = 1:GA_itermax
        %generate the merged and mutated vector
        [~,p1] = histc(rand,cdf);
        [~,p2] = histc(rand,cdf);
        sample1 = sampleM(p1,:);
        sample2 = sampleM(p2,:);
        merge = 1-(1-sample1).*(1-sample2);
        merge_C = 1-merge;
        num = sum(merge_C);
        if(num*beta>=1)
        index_merge_C = find(merge_C==1);
        perm = randperm(num);
        perm = perm(1:floor(beta*num));
        index_merge_C_choose = index_merge_C(perm);
        merge(index_merge_C_choose)= 1;
        end
        
        numofmerge = sum(merge);        
        perm = randperm(numofmerge);
        perm = perm(1:K);   
        index_merge = find(merge==1);
        index_merge_choose = index_merge(perm);
        merge_selectk = 0*merge;
        merge_selectk(index_merge_choose) = 1;

        %test its performance
        [~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(merge_selectk==1), ebno1, N, K);
        [~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(merge_selectk==1), ebno2, N, K);
        ber_merge_selectk = ber1*ber2;
        ber_expand = [ber,ber_merge_selectk];
        sampleM_expand = [sampleM;merge_selectk];
        
        [ber_expand,order] = sort(ber_expand);
        sampleM_expand = sampleM_expand(order,:);
        ber = ber_expand(1:M);
        sampleM = sampleM_expand(1:M,:);
    end
    %disp(sampleM(1:5,:));
end
