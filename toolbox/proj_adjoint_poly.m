function anti_comm_result_full = proj_adjoint_poly(string_of_poly)
%%%%%
%%% anti_comm_poly.m transform a polynomial composed of moments to its
%%% adjoint form. E.g. 1+A1A2-A1A2B1 = 1+A2A1-A2A1B1.
%%%
%%% Warning: This version is only for the coefficient of each term being
%%% either 1 or -1. E.g., -A1+A1B2-1+A3A2B1B4B5
%%%
%%% string_of_poly: the polynomial in the string form.
%%% E.g., string('1+A2-A1A2').
%%%
%%% list_of_string: Multiplication of observables. It can be either a
%%% a single string, e.g. string('B1') or be an array of string,
%%% e.g. [string('Id') string('A1') string('A2')]. Note that it is not
%%% necessary to put different parties' observables in an order in advance.
%%% For instance, either string('A1A2B3B2') or string('A1B3A2B2') is OK. 
%%%
%%% This version is only for the two-party scenario [A B]. 
%%%%%

poly_split = strsplit(string_of_poly,{'+','-'});

char_of_poly = char(string_of_poly);
plus_and_minus = char_of_poly(or(char_of_poly=='+',char_of_poly=='-'));

if poly_split(1)==string()
    poly_split(1) = [];
end

if length(poly_split) == length(plus_and_minus)+1
    % this means the leading coefficient is +, such as
    % 1+A1-A2 or A2-A1B2+B1
    
    plus_and_minus = ['+',plus_and_minus];
    
end

anti_comm_result_full = string();
for i = 1:length(poly_split)
    
gamma = poly_split(i);

if gamma == string('1')
    
    anti_comm_result(i) = string('1');
    
elseif gamma == string('-1')
    
    anti_comm_result(i) = string('-1');
    
else
    
    gamma_sep = strsplit(gamma,'*');
    
    if length(gamma_sep)==1
        anti_comm_result(i) = gamma_sep;
    else
        gamma_A = gamma_sep(~cellfun('isempty', strfind(gamma_sep,'A')));
        % this is to only keep projectors of A
        gamma_B = gamma_sep(~cellfun('isempty', strfind(gamma_sep,'B')));
        % this is to only keep projectors of B
        
        if ~isempty(gamma_A) % length of gamma_A =\= 0
            gamma_A_dag = flip(gamma_A);
            gamma_A_dag = strjoin(gamma_A_dag,'*');
        elseif isempty(gamma_A)
            gamma_A_dag = [];
        end
        
        if ~isempty(gamma_B) % length of gamma_B =\= 0
            gamma_B_dag = flip(gamma_B);
            gamma_B_dag = strjoin(gamma_B_dag,'*');
        elseif isempty(gamma_B)
            gamma_B_dag = [];
        end
        
        anti_comm_result(i) = strjoin(string([gamma_A_dag gamma_B_dag]),'*');
    end
    
end

anti_comm_result_full = strjoin([anti_comm_result_full...
        strjoin([string(plus_and_minus(i)) anti_comm_result(i)],'')],'');

end

anti_comm_result_full_char = char(anti_comm_result_full);
if anti_comm_result_full_char(1) == '+'
    
    anti_comm_result_full_char(1) = [];
    anti_comm_result_full = string(anti_comm_result_full_char);
    
end


end