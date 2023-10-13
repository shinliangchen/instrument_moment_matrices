function string_moment_matrix = ...
    AMM_string(S,Sdag,real_or_complex)

for i = 1:length(Sdag)
    for j = 1:length(S)
        Sdag_i = strsplit(Sdag(i),'*');
        S_j = strsplit(S(j),'*');
        gamma = [Sdag_i S_j];
%        gamma_A = gamma(~cellfun('isempty', strfind(gamma,'A')));
        % this is to only keep moments of A
        gamma_B = gamma(~cellfun('isempty', strfind(gamma,'B')));
        % this is to only keep moments of B
        is_zero = false;
%         if isempty(gamma_A) % length of gamma_A = 0
%             gamma_A_short = string('Id');
%         elseif length(gamma_A) == 1
%             gamma_A_short = gamma_A;
%             % gamma_A_short = Ak, e.g., A3
%         elseif length(gamma_A) >= 1
%             gamma_A_proj = gamma_A;
%             for i_A = 2:length(gamma_A)
%                 ax1 = strsplit(gamma_A_proj(i_A-1),'_');
%                 ax1 = strsplit(ax1(2),'|');
%                 a1 = str2num(char(ax1(1)));
%                 x1 = str2num(char(ax1(2)));
%                 ax2 = strsplit(gamma_A_proj(i_A),'_');
%                 ax2 = strsplit(ax2(2),'|');
%                 a2 = str2num(char(ax2(1)));
%                 x2 = str2num(char(ax2(2)));
%                 
%                 if and(x2==x1,a2~=a1)
%                     is_zero = true;
%                     break % stop the local loop
%                 elseif gamma_A_proj(i_A) == gamma_A_proj(i_A-1)
%                     gamma_A_proj(i_A) = string('Id');
%                     % the two are combined to one if they are the same,
%                     % e.g. A_1|1*A_1|1 = A_1|1*Id
%                     gamma_A_proj = [gamma_A_proj(cellfun('isempty', strfind(gamma_A_proj,'A')))...
%                         gamma_A_proj(~cellfun('isempty', strfind(gamma_A_proj,'A')))];
%                     % this is to move Id to the most-left part, e.g., A0IdIdA1=IdIdA0A1
%                 end
%             end
%             gamma_A_short = gamma_A_proj(cellfun('isempty', strfind(gamma_A_proj,'Id')));
%         else
%             error('something is wrong')
%         end
        
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
        

        if all(gamma_B_short==string('Id'))
            string_moment_matrix(i,j) = string('1');
        elseif string(real_or_complex) == string('complex')
%            gamma_A_short(gamma_A_short == string('Id'))=[];
            gamma_B_short(gamma_B_short == string('Id'))=[];
%            gamma_temp = [gamma_A_short gamma_B_short];
            gamma_temp = gamma_B_short;
            
            if isempty(gamma_temp)==1
                string_moment_matrix(i,j) = string('1');
            elseif is_zero == 1
                string_moment_matrix(i,j) = string('0');
            else
                string_moment_matrix(i,j) = strjoin(gamma_temp,'*');
            end
            
        elseif string(real_or_complex) == string('real')
%            gamma_A_short(gamma_A_short == string('Id'))=[];
            gamma_B_short(gamma_B_short == string('Id'))=[];
%            gamma_temp = [gamma_A_short gamma_B_short];
            gamma_temp = gamma_B_short;
            gamma_temp_flip = flip(gamma_B_short);
            
            if isempty(gamma_temp) == 1
                string_moment_matrix(i,j) = string('1');
            elseif is_zero == 1
                string_moment_matrix(i,j) = string('0');
            elseif ~or(any(any(strjoin(gamma_temp,'*') == string_moment_matrix)) ,...
                    any(any(strjoin(gamma_temp_flip,'*') == string_moment_matrix)))
                string_moment_matrix(i,j) = strjoin(gamma_temp,'*');
            elseif any(any(strjoin(gamma_temp,'*') == string_moment_matrix))
                string_moment_matrix(i,j) = strjoin(gamma_temp,'*');
            elseif any(any(strjoin(gamma_temp_flip,'*') == string_moment_matrix))
                string_moment_matrix(i,j) = strjoin(gamma_temp_flip,'*');
            else
                error('something is wrong')
            end
            
        else
            error('something is wrong')
        end
        
        
    end
end

end