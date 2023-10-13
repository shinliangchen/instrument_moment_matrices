%%% Computing the successful probabilities in the 2 --> 1 QRAC scenario
clear all;
tic;

m = input('enter 0 if no constraint on POVMs; enter 1 if rank of POVMs is 1; Ans:');

level = 1;
d = 2;

nx = 2;
ny = 2;
na = 2;
nb = 2;

S = AMM_gen_xlevel_seq(ny,nb,level);

Sdag = [];
for i = 1:length(S)
    Sdag = [Sdag, proj_adjoint_poly(S(i))];
end

gamma_str_complex = AMM_string(S,Sdag,'complex');
gamma_str_real = AMM_string(S,Sdag,'real');


uni_mono = unique(gamma_str_real); % real or complex?

for x = 1:nx
    for a = 1:na
        gamma_SDP{a,x} = zeros(length(gamma_str_complex));
        gamma_SDP{a,x} = gamma_SDP{a,x} + (gamma_str_complex == string('1'));
    end
end

uni_mono(uni_mono==string('1'))=[];
uni_mono(uni_mono==string('0'))=[];

for idx = 1:length(uni_mono)
    
    tfMatrix = (gamma_str_complex==uni_mono(idx));
    
    for a = 1:na
        for x = 1:nx
            if proj_adjoint_poly(uni_mono(idx))~=uni_mono(idx)
                tfMatrix2 = (gamma_str_complex==proj_adjoint_poly(uni_mono(idx)));
                u{idx,a,x} = sdpvar(1,1,'hermitian','real');
                gamma_SDP{a,x} = gamma_SDP{a,x} + u{idx,a,x}.*tfMatrix + u{idx,a,x}.*tfMatrix2;
            else
                u{idx,a,x} = sdpvar(1,1,'hermitian','real');
                gamma_SDP{a,x} = gamma_SDP{a,x} + u{idx,a,x}.*tfMatrix;
            end
        end
    end
    if mod(idx,50)==0
        disp(idx)
    end
end

blk_gamma_SDP = [];
for x = 1:nx
    for a = 1:na
        blk_gamma_SDP = blkdiag(blk_gamma_SDP, gamma_SDP{a,x});
    end
end

for a = 1:na
    for x = 1:nx
        for y = 1:ny
            p_B_baxy{1,a,x,y} = sdpvar(1,1,'hermitian','real');
            p_B_baxy{2,a,x,y} = 1 - p_B_baxy{1,a,x,y};
        end
    end
end

coun = 0;
SE = 0;
for x = 1:nx
    for a = 1:na
        for y = 1:ny
            for b = 1:nb
                xx = [a,x]; % x_1, x_2
                SE = SE + (1/8).*p_B_baxy{b,a,x,y}.*kronDel(b,xx(y));
                if b == xx(y)
                    coun = coun + 1;
                end
            end
        end
    end
end


constr = [];

%%%%%% semi-DI constraint (start)
[S_num,~] = AMM_proj_seq_str2num(d,ny,nb,S,m);
AMM_new = AMM_str2num_blk_temporal_normalized(d,nx,S_num,'real'); % randomnly generate AMM, in the block-diagonal form
%AMM_new = AMM_str2num_blk_temporal_normalized_ver2(d,nx,S_num,'real'); % randomnly generate AMM, in the block-diagonal form

G = AMM_new(:);
size_G = size(AMM_new,1)/(nx*na); % size of moment matrix, i.e., length of S
if size_G~=length(S)
    error('something may be wrong')
end

r = 1;
while r > 10^-7
    S_num = AMM_proj_seq_str2num(d,ny,nb,S,m);
    AMM_new = AMM_str2num_blk_temporal_normalized(d,nx,S_num,'real');
    %AMM_new = AMM_str2num_blk_temporal_normalized_ver2(d,nx,S_num,'real');
    
    %    randomnly generate AMM (in the block-diagonal form) until the set
    %    of AMM (in the block-diagonal form) spans the space that AMM (...) lives in
    Gnew = AMM_new(:);
    for k=1:size(G,2)
        Gnew = Gnew - (Gnew'*G(:,k))*G(:,k); %Gram-Schmidt process
    end
    r = sqrt(Gnew'*Gnew);
    Gnew = Gnew./r;
    %disp(r);
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
                constr = [constr, p_B_baxy{b,a,x,y} == u{uni_mono==str_Eby,a,x}];
            end
        end
    end
end

for x = 1:nx
    sum_a_gamma_SDP{x} = 0;
    for a = 1:na
        gamma_SDP{a,x}(1,1) = 1;%p_A_ax{a,x};
        sum_a_gamma_SDP{x} = sum_a_gamma_SDP{x} + gamma_SDP{a,x};
    end
end

for x = 1:nx
    for a = 1:na
        constr = [constr, gamma_SDP{a,x}>=0, (1/2).*(gamma_SDP{a,x}'+gamma_SDP{a,x})>=0];
    end
end


%opts = sdpsettings('verbose', 0);
sol=solvesdp(constr , -SE);%, opts);
double(SE)
sol

toc;
