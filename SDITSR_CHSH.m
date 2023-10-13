% This code compute a SDI lower bound on TSR given a temporal CHSH
% violation. I.e., Fig. 4 of https://arxiv.org/abs/2305.19548

clear all;
tic;

Bell_violation = 2*sqrt(2);

level = 1;
d = 2; % Uncommnet this for SDI
rank_constraint = 1; % 0: no rank constraint

nx = 2;
ny = 2;
na = 2;
nb = 2;
n_lambda = na^nx;

S = AMM_gen_xlevel_seq(ny,nb,level);

Sdag = [];
for i = 1:length(S)
    Sdag = [Sdag, proj_adjoint_poly(S(i))];
end


%%%%%% build \chi[I_{a|x}] (start)
chi_str_complex = AMM_string(S,Sdag,'complex');
chi_str_real = AMM_string(S,Sdag,'real');

uni_mono = unique(chi_str_real); % real or complex?

for x = 1:nx
    for a = 1:na
        chi_ax{a,x} = zeros(length(chi_str_complex));
    end
end

uni_mono(uni_mono==string('1'))=[];
uni_mono(uni_mono==string('0'))=[];


for idx = 1:length(uni_mono)
    
    tfMatrix = (chi_str_real==uni_mono(idx));
    
    for a = 1:na
        for x = 1:nx
            u{idx,a,x} = sdpvar(1,1,'hermitian','real');
            chi_ax{a,x} = chi_ax{a,x} + u{idx,a,x}.*tfMatrix;
        end
    end
    if mod(idx,50)==0
        disp(idx)
    end
end
%%%%%% build \chi[I_{a|x}] (end)


%%%%%% build \chi[G_\lambda] (start)
clear a
clear x
uni_mono_lambda = unique(chi_str_real); % real or complex?
for i_lambda = 1:n_lambda
    chi_lambda{i_lambda} = zeros(length(chi_str_complex));
end

%uni_mono_lambda(uni_mono_lambda==string('1'))=[];
% the uni_mono of string('1') should be kept when constructing chi_lambda
uni_mono_lambda(uni_mono_lambda==string('0'))=[];

for idx = 1:length(uni_mono_lambda)
    
    tfMatrix = (chi_str_real==uni_mono_lambda(idx));
    
    for i_lambda = 1:n_lambda
        
        v{idx,i_lambda} = sdpvar(1,1,'hermitian','real');
        chi_lambda{i_lambda} = chi_lambda{i_lambda} + v{idx,i_lambda}.*tfMatrix;
        
    end
    if mod(idx,50)==0
        disp(idx)
    end
end
%%%%%% build \chi[G_\lambda] (end)

D = detPax(nx,na);
for x = 1:nx
    for a = 1:na
        sum_lambda_D_chi_lambda{a,x} = 0;
        for i_lambda = 1:n_lambda
            sum_lambda_D_chi_lambda{a,x} = sum_lambda_D_chi_lambda{a,x} + D(i_lambda,a,x).*chi_lambda{i_lambda};
        end
    end
end
sum_lambda_chi_lambda = 0;
for i_lambda = 1:n_lambda
    sum_lambda_chi_lambda = sum_lambda_chi_lambda + chi_lambda{i_lambda};
end

% Uncomment the following for SDI
blk_gamma_SDP = [];
for x = 1:nx
    for a = 1:na
        blk_gamma_SDP = blkdiag(blk_gamma_SDP, chi_ax{a,x}); 
    end
end


for x = 1:nx
    pax{1,x} = sdpvar(1,1,'hermitian','real');
    pax{2,x} = 1 - pax{1,x};
end
for x = 1:nx
    for y = 1:ny
        pabxy{1,1,x,y} = sdpvar(1,1,'hermitian','real');
        pabxy{2,1,x,y} = sdpvar(1,1,'hermitian','real');
        
        pabxy{1,2,x,y} = pax{1,x} - pabxy{1,1,x,y};
        pabxy{2,2,x,y} = 1 - pabxy{1,1,x,y} - pabxy{2,1,x,y} - pabxy{1,2,x,y};
    end
end


for x = 1:nx
    for y = 1:ny
        E{x,y} = pabxy{1,1,x,y} + pabxy{2,2,x,y} - pabxy{1,2,x,y} - pabxy{2,1,x,y};
        
    end
end

SE = E{1,1}+E{1,2}+E{2,1}-E{2,2};


constr = [];

%%%%%% semi-DI constraint (start)
[S_num,~] = AMM_proj_seq_str2num(d,ny,nb,S,rank_constraint);
AMM_new = AMM_str2num_blk_temporal_ver2(d,nx,S_num,'real'); % randomnly generate AMM, in the block-diagonal form

G = AMM_new(:);
size_G = size(AMM_new,1)/(nx*na); % size of moment matrix, i.e., length of S
if size_G~=length(S)
    error('something may be wrong')
end

r = 1;
while r > 10^-7
    S_num = AMM_proj_seq_str2num(d,ny,nb,S,rank_constraint);
    AMM_new = AMM_str2num_blk_temporal_ver2(d,nx,S_num,'real');
    
    %  randomnly generate AMM (in the block-diagonal form) until the set
    %  of AMM (in the block-diagonal form) spans the space that AMM (...) lives in
    Gnew = AMM_new(:);
    for k=1:size(G,2)
        Gnew = Gnew - (Gnew'*G(:,k))*G(:,k); % Gram-Schmidt process
    end
    r = sqrt(Gnew'*Gnew);
    Gnew = Gnew./r;
    G = [G,Gnew];% set of orthonormal vectors (it forms an orthonormal basis when the loop is finished)
end

sum_mu_Gk = 0;
for kk=1:size(G,2)
    Gk{kk} = reshape(G(:,kk),sqrt(length(G(:,1))),[]);% reshape vectorized G into the the matrix form
    mu{kk} = sdpvar(1,1,'hermitian','real');
    sum_mu_Gk = sum_mu_Gk + mu{kk}.*Gk{kk};
end
constr = [constr, blk_gamma_SDP == sum_mu_Gk];

%%%%%% semi-DI constraint (end)


for y = 1:ny
    for b = 1:nb-1
        str_Eby = string(strcat('B_',num2str(b),'|',num2str(y)));
        for x = 1:nx
            for a = 1:na
                constr = [constr, pabxy{a,b,x,y} == u{uni_mono==str_Eby,a,x}];
            end
        end
    end
end

for x = 1:nx
    sum_a_gamma_SDP{x} = 0;
    for a = 1:na
        chi_ax{a,x}(1,1) = pax{a,x};
        sum_a_gamma_SDP{x} = sum_a_gamma_SDP{x} + chi_ax{a,x};
    end
end

for i_lambda = 1:n_lambda
    constr = [constr, chi_lambda{i_lambda} >= 0, chi_lambda{i_lambda}+chi_lambda{i_lambda}' >= 0];
end

for x = 1:nx
    for a = 1:na
        constr = [constr, sum_lambda_D_chi_lambda{a,x} >= chi_ax{a,x}];
        constr = [constr, chi_ax{a,x}>=0, (1/2).*(chi_ax{a,x}'+chi_ax{a,x})>=0];
    end
end

constr = [constr, SE == Bell_violation];

DITSR_var = sum_lambda_chi_lambda(1,1) - 1;
%opts = sdpsettings('verbose', 0);
sol=solvesdp(constr , DITSR_var);%, opts);
DITSR = double(DITSR_var)
sol

toc;

