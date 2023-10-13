function seq_all = AMM_gen_xlevel_seq(ny,nb,level)

%if any(level == [1 2 3 4 5 6]) == 0
%    error('the level should be set as a positive integer lower or equal than 7')
%end



seq_S_B = [];
for y = 1:ny
    for b = 1:nb-1
        seq_S_B = [seq_S_B string(strcat('B_',num2str(b),'|',num2str(y)))];
    end
end

if level == 1
    seq_all = [string('Id') seq_S_B];
else
    seq_all = [string('Id') seq_S_B];
    seq_S_xB{1} = seq_S_B;
    %len_all = 0;
    %i_level = 2;
    for i_level = 2:level
        Sdag = seq_S_xB{i_level-1};
        S = seq_S_xB{1};
        
        idx = 1;
        for i = 1:length(Sdag)
            for j = 1:length(S)
                Sdag_i = strsplit(Sdag(i),'*');
                S_j = strsplit(S(j),'*');
                gamma = [Sdag_i S_j];
                
                gamma_B = gamma(~cellfun('isempty', strfind(gamma,'B')));
                % this is to only keep moments of B
                is_zero = false;
                
                if isempty(gamma_B) % length of gamma_B = 0
                    gamma_B_short = string('Id');
                elseif length(gamma_B) == 1
                    gamma_B_short = gamma_B;
                    % gamma_B_short = Bk, e.g., B3
                elseif length(gamma_B) >= 1
                    gamma_B_proj = gamma_B;
                    for i_B = 2:length(gamma_B)
                        by1 = strsplit(gamma_B_proj(i_B-1),'_');
                        by1 = strsplit(by1(2),'|');
                        b1 = str2num(char(by1(1)));
                        y1 = str2num(char(by1(2)));
                        by2 = strsplit(gamma_B_proj(i_B),'_');
                        by2 = strsplit(by2(2),'|');
                        b2 = str2num(char(by2(1)));
                        y2 = str2num(char(by2(2)));
                        
                        if and(y2==y1,b2~=b1)
                            is_zero = true;
                            break % stop the local loop
                        elseif gamma_B_proj(i_B) == gamma_B_proj(i_B-1)
                            gamma_B_proj(i_B) = string('Id');
                            % the two are combined to one if they are the same,
                            % e.g. B_1|1*B_1|1 = B_1|1*Id
                            gamma_B_proj = [gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'B')))...
                                gamma_B_proj(~cellfun('isempty', strfind(gamma_B_proj,'B')))];
                            % this is to move Id to the most-left part, e.g., B0IdIdB1=IdIdB0B1
                        end
                    end
                    gamma_B_short = gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'Id')));
                else
                    error('something is wrong')
                end
                
                
                if length(gamma_B_short)==i_level
                    
                    gamma_B_short(gamma_B_short == string('Id'))=[];
                    gamma_temp = gamma_B_short;
                    
                    if isempty(gamma_temp)==1
                        
                    elseif is_zero == 1
                        
                    else
                        seq_S_xB{i_level}(idx) = strjoin(gamma_temp,'*');
                        idx = idx + 1;
                    end
                    
                end
                
                
            end
        end
        seq_all = [seq_all, seq_S_xB{i_level}];
        %len_all = len_all + length(seq_S_xB{i_level});
    end
    
end

% seq_S_xB{1}
% seq_S_xB{2}
% seq_S_xB{3}
% seq_S_xB{4}

%len_all = len_all + 1 + length(seq_S_xB{1});
%len_all
%length(seq_S_xB{3})
%length(unique(seq_S_xB{3}))

%length(seq_S_xB{3})+length(seq_S_xB{2})+length(seq_S_xB{1})+1

end
