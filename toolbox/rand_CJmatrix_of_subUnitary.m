function Iax = rand_CJmatrix_of_subUnitary(nx,d)

%%%%%%
%  randAssemblage_pure randomnly generates an pure temporal assemblage (with binary
%  outcome. By pure we mean its associated normalized states are all pure.
%  I.e., assemblage(:,:,a,x) = pax(a,x).*assemblage_norm(:,:,a,x) where
%  assemblage_norm(:,;,a,x) are all pure noramlized states.
%%%%%%

na = 2;

%psi = randPsi(d);
%rho_0 = psi*psi';

% for x = 1:nx
%     vec_Eax = randn(d,1)+1i*randn(d,1);
%     Eax{1,x} = vec_Eax*vec_Eax'./(vec_Eax'*vec_Eax);
% end

sub_MES = d.*MES(d)*MES(d)';

for x = 1:nx
    pax(1,x) = rand;
    %pax(1,x) = trace(Eax{1,x}*rho_0);
    pax(2,x) = 1 - pax(1,x);
    UU{x} = randU(d);
    for a = 1:na
        %Iax(:,:,a,x) = pax(a,x).*( kron(eye(d), UU{x}*sub_MES*UU{x}' ) ); % the 1st way of CJ matrix
        Iax(:,:,a,x) = pax(a,x).*transpose(( kron(eye(d),UU{x}) * sub_MES * kron(eye(d),UU{x})'  )); % the 3rd way of CJ matrix
    end
end

end