function assemblage_moment_matrix_blk = AMM_str2num_blk_temporal_normalized(d,nx,S_num,real_or_complex)

% AMM_str2num_blk_temporal_normalized randomnly generate the temporal AMM (with assemblage being normalized
% ,with the entries obeying S_num)

% S_num: the sequence in the numerical form

na = 2;
%assemblage = randAssemblage_pure_temporal(nx,d);

for x = 1:nx
    for a = 1:na
        %assemblage(:,:,a,x) = (1/trace(assemblage(:,:,a,x))).*assemblage(:,:,a,x);
        
        psi = randPsi(d);
        assemblage(:,:,a,x) = psi*psi';
    end
end

CJ_MES = d.*MES(d)*MES(d)';
%rand_CJ = kron(eye(d),randU(d))*CJ_MES;
rand_CJ = kron(eye(d),eye(d))*CJ_MES;


for kk = 1:length(S_num)
    S_num_dag{kk} = S_num{kk}';
end

for x = 1:nx
    for a = 1:na
        
        for ii = 1:length(S_num_dag)
            for jj = 1:length(S_num)
                %assemblage_moment_matrix{a,x}(ii,jj) = trace(kron(assemblage(:,:,a,x),S_num_dag{ii}*S_num{jj})*rand_CJ);
                assemblage_moment_matrix{a,x}(ii,jj) = trace(assemblage(:,:,a,x)*S_num_dag{ii}*S_num{jj});
            end
        end
        
        if string(real_or_complex) == string('real')
            assemblage_moment_matrix{a,x} = (1/2).*(assemblage_moment_matrix{a,x} + assemblage_moment_matrix{a,x}');
        end
        
    end
end

assemblage_moment_matrix_blk = [];
for x = 1:nx
    for a = 1:na
        assemblage_moment_matrix_blk = blkdiag(assemblage_moment_matrix_blk, assemblage_moment_matrix{a,x});
    end
end

end