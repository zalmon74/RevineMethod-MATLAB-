% coor_1 - ������ ��������� ���������� � ������� �������� ����� ��������
% coor_2 - ������ ��������� ���������� � ������� ���������� ����� ��������
% vec_sym - ������, ������� �������� �������, ������� ������� � �������
% func - �������

function coor = MethodRavine(coor_1, coor_2, eps, vec_sym, func)

% ���������� ��� ������������ ���������� �������
syms d
vec_sym_diff = vec_sym;
% ������� �������� ����������� � �������
for i = 1:length(vec_sym)
    vec_diff(i) = diff(func, vec_sym(i));
    % ���������� ���������� ������� � ������� �����������
    temp = symvar(vec_diff(i));
    temp(end+1:length(vec_sym)) = vec_sym(i);
    vec_sym_diff(i,:) = temp;
end

% �������� ����
for N = 1:9999999999
    
    if (N == 1)
        % ��������� ����� ����������
        new_coor_1 = CalcNewCoord(coor_1, vec_sym, vec_diff, vec_sym_diff, func);
        new_coor_2 = CalcNewCoord(coor_2, vec_sym, vec_diff, vec_sym_diff, func);
        % ��������� �������� ������� � ���������� �����������
        res_func_n_1 = CalcFuncCoor(new_coor_1, vec_sym, func);
        res_func_n_2 = CalcFuncCoor(new_coor_2, vec_sym, func);
    else
        % ����������� ���������� ��������� � coor_1, � � coor_2 ���������
        % ���� ���������� ��� � � coor_1 ������ �����������
        coor_1 = coor;
        coor_2 = ceil(coor);
        % ��������� ����� ����������
        new_coor_1 = CalcNewCoord(coor_1, vec_sym, vec_diff, vec_sym_diff, func);
        new_coor_2 = CalcNewCoord(coor_2, vec_sym, vec_diff, vec_sym_diff, func);
        % ��������� �������� ������� � ���������� �����������
        res_func_n_1 = CalcFuncCoor(new_coor_1, vec_sym, func);
        res_func_n_2 = CalcFuncCoor(new_coor_2, vec_sym, func);
    end

    % ���� �������� ���� ���������
    for i = 1:length(new_coor_1)
        % ���������� ����� ������� ����� ������ �������� � �� ������� � ��������� ������� ��� ����������� ����
        if (res_func_n_1 > res_func_n_2)
            func_new_coor(i) = new_coor_2(i)-d*(new_coor_1(i)-new_coor_2(i));
        else
            func_new_coor(i) = new_coor_1(i)-d*(new_coor_2(i)-new_coor_1(i));
        end
    end
    
    % ����������� ����� ��������� � �������� ������� ��� ���������� ����������
    new_func = InsertIntoMaineFunc(func_new_coor, vec_sym, func);
    % ������ ���������
    res_d = diff(new_func); % ����� ����������� �.�. �� ��������� ������� ����� ��, � ����� ������
    res_d = vpasolve(res_d==0);
    res_d = eval(res_d(1));
    % ����������� � ����� ������� ���������� �������� h
    new_coor = subs(func_new_coor, d, res_d);
    coor = eval(new_coor);
    % ������������ �������� � ������ ������������ � ���� �� ����� 0, �� ����������� 
    coor_grad = CalcCoorGrad(coor, vec_sym, vec_diff, vec_sym_diff);
    
    % ���� �������� ��������� ������ ������������ ��������, �� ���������
    flag_return = true;
    for i = 1:length(coor_grad)
        if (coor_grad(i) > eps)
            flag_return = false;
        end
    end
    if (flag_return == true)
        return;
    end
    
end

end

% ������� ���������� ��������� � �����������
% ��. ���������:
    % coor - ����������, � ������� ���������� ��������� �������� �������
    % vec_sym - ������, ������� �������� �������, ������� ������� � �������
    % vec_diff - ������� ����������� ������� (��� ��������)
    % vec_sym_diff - ���������� ������� � ������� �����������
% ���. ������:
    % vec_diff_val - �������� ��������� � ������������ �����������

function vec_diff_val_out = CalcCoorGrad(coor, vec_sym, vec_diff, vec_sym_diff)

% ���� �������� ���� ���������
for i = 1:length(vec_sym)
    vec_diff_val(i) = vec_diff(i);
    % ���� �������� ���������� �������� � ������� �����������
    for i_sym = 1:length(vec_sym_diff(i,:))
        % ���������� ������ �������� ������� � ����� ������� ��������
        ind_sym = find(vec_sym==vec_sym_diff(i, i_sym));
        % ����������� ���������� � ������� �����������
        vec_diff_val(i) = subs(vec_diff_val(i), vec_sym(ind_sym), coor(ind_sym));
        % ���������� �������� �� ������� � �������, ���� ���, ��
        % ����������� ������ ����
        temp = symvar(vec_diff_val(i));
        if (isempty(temp))
            break;
        end
    end
    % ��������������� ���������� �������� � double
    vec_diff_val_out(i) = eval(vec_diff_val(i));

end

end

% ������� ���������� �������� ������� � �����������
% ��. ���������:
    % coor - ����������, � ������� ���������� ��������� �������� �������
    % vec_sym - ������, ������� �������� �������, ������� ������� � �������
    % func - ������� � ������� ���������� ��������� ��������
% ���. ���������:
    % res - �������� �������
function res = CalcFuncCoor(coor, vec_sym, func)

res = func;
% ���� �������� ���� ��������
for i_sym = 1:length(vec_sym)
    res = subs(res, vec_sym(i_sym), coor(i_sym));
end
% ��������������� ���������� �������� � double
res = eval(res);

end

% ������� ������������ ����� ���. �������
% ��. ���������:
    % coor - ����������, ������� ���������� � ������������ ����� �������
    % grad_coor - �������� ���������, ������� � ������������ ����� �������
    % c_sym - ������, ��� ������ �������� ��������� ����� �������
% ���. ���������:
    % vec_new_func - ������, ������� ������ � ���� �������������� ����� �������
function vec_new_func = CreateNewFunc(coor, coor_grad, c_sym)

% ���� �������� ���� ���������
for i = 1:length(coor)

    % ������� ����� ������� new_func = coor(i)-grad_coor(i)*c_sym
    vec_new_func(i) = coor(i)-coor_grad(i)*c_sym;
    
end

end

% ������� ����������� ����� ������� � �������� �������
% ��. ���������:
    % vec_new_func - ������ � ���������, ������� ���������� ����������
    % vec_sym - ������, ������� �������� � ���� ������, ������� �������� � ���� �������� �������
    % func - �������� �������
% ���. ���������:
	% new_func - �������� �������, � ������� �������� ����� �������
function new_func = InsertIntoMaineFunc(vec_new_func, vec_sym, func)

new_func = func;

% ���� �������� ���� ��������
for i = 1:length(vec_sym)
    
    % ��������� ����� ������� ������ ������������� ������� � �������� �������
    new_func = subs(new_func, vec_sym(i), vec_new_func(i));
    
end

end

% ������� ���������� ����� ���������
% ��. ���������:
    % old_coor - ������ ����������
    % vec_sym - ������, ������� �������� �������, ������� ������� � �������
    % vec_diff - ������� ����������� ������� (��� ��������)
    % vec_sym_diff - ���������� ������� � ������� �����������
    % func - �������� �������
% ���. ���������:
    % coor - ����� ����������
function coor = CalcNewCoord(old_coor, vec_sym, vec_diff, vec_sym_diff, func)

syms h
% ��������� �������� ����������
coor_grad_1 = CalcCoorGrad(old_coor, vec_sym, vec_diff, vec_sym_diff);
% ��������� ����� �������
vec_new_func = CreateNewFunc(old_coor, coor_grad_1, h);
% ����������� � �������� ������� �������������� ����
new_func = InsertIntoMaineFunc(vec_new_func, vec_sym, func);
% ������ ���������
res_h = diff(new_func); % ����� ����������� �.�. �� ��������� ������� ����� ��, � ����� ������
res_h = vpasolve(res_h==0);
res_h = eval(res_h(1));
% ����������� � ����� ������� ���������� �������� h
new_coor_1 = subs(vec_new_func, h, res_h);
coor = eval(new_coor_1);

end
