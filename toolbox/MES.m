function y = MES(d)
% MES generates the maximally entangled state with dimension d

if d<=1
    error('d must be an integer greater than 1')
end

y = 0;
for i = 1:d
    ket{i} = zeros(d,1);
    ket{i}(i) = 1;
    y = y + (1/sqrt(d)).*kron(ket{i},ket{i});
end

end

