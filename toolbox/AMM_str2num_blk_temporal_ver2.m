function assemblage_moment_matrix_blk = AMM_str2num_blk_temporal_ver2(d,nx,S_num,real_or_complex)
% the difference is the way to create rand_CJ

% AMM_str2num_temporal randomnly generate the temporal AMM, with the entries obeying S_num

% S_num: the sequence in the numerical form

na = 2;
%assemblage = randAssemblage_pure_temporal(nx,d);
%CJ_MES = d.*MES(d)*MES(d)';
%rand_CJ = kron(eye(d),randU(d))*CJ_MES;
%rand_CJ = kron(eye(d),eye(d))*CJ_MES;

Iax = rand_CJmatrix_of_subUnitary(nx,d);
psi = randPsi(d);
rho_0 = psi*psi';

for kk = 1:length(S_num)
    S_num_dag{kk} = S_num{kk}';
end

for x = 1:nx
    for a = 1:na
        
        for ii = 1:length(S_num_dag)
            for jj = 1:length(S_num)
                assemblage_moment_matrix{a,x}(ii,jj) = trace(kron(eye(d),S_num_dag{ii}*S_num{jj})*Iax(:,:,a,x)*kron(rho_0,eye(d)));
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