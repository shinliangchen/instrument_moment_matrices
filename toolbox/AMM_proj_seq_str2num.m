function [S_num,Eby] = AMM_proj_seq_str2num(d,ny,nb,S,m)

% proj_seq_str2num converts the sequence from the string form into the
% numerical form
% m: 0 when no rank constraint is on POVM
%    1 when rank constraint is k=1

% level = 3;
% d = 2;
% 
% nx = 3;
% ny = 3;
% na = 2;
% nb = 2;
% 
% S = proj_gen_xlevel_seq(nx,ny,na,nb,level);
dd = 2;

%m = input('enter 0 if no constraint on POVMs; enter 1 if rank of POVMs is 1; Ans:');

if nb~=2
    error('the current code is only for nb=2')
end

if m == 0 % no rank constraint on POVM
    for y = 1:ny
        Ey = RandomPOVM(dd,nb);
        for b = 1:nb
            Eby{b,y} = Ey{b};
        end
    end
elseif m == 1 % rank constraint is k=1
    for y = 1:ny
        vec_Eby = randn(d,1)+1i*randn(d,1);
        Eby{1,y} = vec_Eby*vec_Eby'./(vec_Eby'*vec_Eby);
        Eby{2,y} = eye(dd) - Eby{1,y};
    end
else
    error('the type of generating POVM must be specified')
end
% for y = 1:ny
%     %%%%%% this part is for the QRAC scenario
%     %%%%%% or strictly speaking, rank constraint is k=1
%     vec_Eby = randn(d,1)+1i*randn(d,1);
%     Eby{1,y} = vec_Eby*vec_Eby'./(vec_Eby'*vec_Eby);
%     Eby{2,y} = eye(dd) - Eby{1,y};
% 
%     %%%%%%% this part is for temporal CHSH scenario
%     %%%%%%% or strictly speaking, no rank constraint on POVM
%     Ey = RandomPOVM(dd,nb);
%     for b = 1:nb
%         Eby{b,y} = Ey{b};
%     end
% end


for ii = 1:length(S)
    if S(ii) == string('Id')
        S_num{ii} = eye(d);
    else
        S_i_split = strsplit(S(ii),'*');
        S_num{ii} = eye(d);
        for jj = 1:length(S_i_split)
            A_or_B = strsplit(S_i_split(jj),'_');
            if A_or_B(1) == string('B')
                by = strsplit(A_or_B(2),'|');
                b = str2double(by(1));
                y = str2double(by(2));
                S_num{ii} = S_num{ii}*Eby{b,y};
            else
                error('something is wrong')
            end
        end
    end
end

if length(S_num)~=length(S)
    error('something is wrong 2')
end



