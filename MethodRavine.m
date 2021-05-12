% coor_1 - Первые начальные координаты с которых начнется поиск минимума
% coor_2 - Вторые начальные координаты с которых начинается поиск минимума
% vec_sym - Вектор, который содержит символы, которые имеются в функции
% func - Функция

function coor = MethodRavine(coor_1, coor_2, eps, vec_sym, func)

% Переменная для формирования внутренний функций
syms d
vec_sym_diff = vec_sym;
% Находим часатные производные у функции
for i = 1:length(vec_sym)
    vec_diff(i) = diff(func, vec_sym(i));
    % Определяем оставшиеся символы в частных производных
    temp = symvar(vec_diff(i));
    temp(end+1:length(vec_sym)) = vec_sym(i);
    vec_sym_diff(i,:) = temp;
end

% Основной цикл
for N = 1:9999999999
    
    if (N == 1)
        % Вычисляем новые координаты
        new_coor_1 = CalcNewCoord(coor_1, vec_sym, vec_diff, vec_sym_diff, func);
        new_coor_2 = CalcNewCoord(coor_2, vec_sym, vec_diff, vec_sym_diff, func);
        % Вычисляем значение функции в полученных координатах
        res_func_n_1 = CalcFuncCoor(new_coor_1, vec_sym, func);
        res_func_n_2 = CalcFuncCoor(new_coor_2, vec_sym, func);
    else
        % Вычисленные координаты переходят в coor_1, а в coor_2 переходят
        % теже координаты что и в coor_1 только округленные
        coor_1 = coor;
        coor_2 = ceil(coor);
        % Вычисляем новые координаты
        new_coor_1 = CalcNewCoord(coor_1, vec_sym, vec_diff, vec_sym_diff, func);
        new_coor_2 = CalcNewCoord(coor_2, vec_sym, vec_diff, vec_sym_diff, func);
        % Вычисляем значение функции в полученных координатах
        res_func_n_1 = CalcFuncCoor(new_coor_1, vec_sym, func);
        res_func_n_2 = CalcFuncCoor(new_coor_2, vec_sym, func);
    end

    % Цикл перебора всех координат
    for i = 1:length(new_coor_1)
        % Определяем какая функция имеет меньше значение в ту сторону и формируем функцию для определения шага
        if (res_func_n_1 > res_func_n_2)
            func_new_coor(i) = new_coor_2(i)-d*(new_coor_1(i)-new_coor_2(i));
        else
            func_new_coor(i) = new_coor_1(i)-d*(new_coor_2(i)-new_coor_1(i));
        end
    end
    
    % Подставляем новые уравнения в основную функцию для нахождения экстремума
    new_func = InsertIntoMaineFunc(func_new_coor, vec_sym, func);
    % Решаем уравнение
    res_d = diff(new_func); % Берем производную т.к. по алгоритму сначала берем ее, а потом решаем
    res_d = vpasolve(res_d==0);
    res_d = eval(res_d(1));
    % Подставляем в новые функции полученное значение h
    new_coor = subs(func_new_coor, d, res_d);
    coor = eval(new_coor);
    % Рассчитываем градиент с новыми координатами и если он равен 0, то заканчиваем 
    coor_grad = CalcCoorGrad(coor, vec_sym, vec_diff, vec_sym_diff);
    
    % Если значение градиента меньше определенной точности, то завершаем
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

% Функция вычисления градиента в координатах
% ВХ. АРГУМЕНТЫ:
    % coor - координаты, в которых необходимо вычислить градиент функции
    % vec_sym - Вектор, который содержит символы, которые имеются в функции
    % vec_diff - Частные производные функции (сам градиент)
    % vec_sym_diff - оставшиеся символы в частных производных
% ВЫХ. ДАННЫЕ:
    % vec_diff_val - значение градиента в определенных координатах

function vec_diff_val_out = CalcCoorGrad(coor, vec_sym, vec_diff, vec_sym_diff)

% Цикл перебора всех симоволов
for i = 1:length(vec_sym)
    vec_diff_val(i) = vec_diff(i);
    % Цикл перебора оставшихся символов в частных произодвных
    for i_sym = 1:length(vec_sym_diff(i,:))
        % Определяем индекс текущего символа в общем векторе символов
        ind_sym = find(vec_sym==vec_sym_diff(i, i_sym));
        % Подставляем координаты в частные производные
        vec_diff_val(i) = subs(vec_diff_val(i), vec_sym(ind_sym), coor(ind_sym));
        % Определяем остались ли символы в функции, если нет, то
        % заканчиваем данный цикл
        temp = symvar(vec_diff_val(i));
        if (isempty(temp))
            break;
        end
    end
    % Преобразовываем символьное значение в double
    vec_diff_val_out(i) = eval(vec_diff_val(i));

end

end

% Функция вычисления значение функции в координатах
% ВХ. АРГУМЕНТЫ:
    % coor - координаты, в которых необходимо вычислить градиент функции
    % vec_sym - Вектор, который содержит символы, которые имеются в функции
    % func - функция в которой необходимо вычислить значение
% ВЫХ. АРГУМЕНТЫ:
    % res - значение функции
function res = CalcFuncCoor(coor, vec_sym, func)

res = func;
% Цикл перебора всех символов
for i_sym = 1:length(vec_sym)
    res = subs(res, vec_sym(i_sym), coor(i_sym));
end
% Преобразовываем символьное значение в double
res = eval(res);

end

% функция формирования новых мат. функций
% ВХ. АРГУМЕНТЫ:
    % coor - координаты, которые учавствуют в формировании новой функции
    % grad_coor - значения градиента, которые в формировании новой функции
    % c_sym - символ, при помощи которого создается новаы функция
% ВЫХ. АРГУМЕНТЫ:
    % vec_new_func - вектор, который хранит в себе сформированные новые функции
function vec_new_func = CreateNewFunc(coor, coor_grad, c_sym)

% Цикл перебора всех координат
for i = 1:length(coor)

    % Создаем новые функции new_func = coor(i)-grad_coor(i)*c_sym
    vec_new_func(i) = coor(i)-coor_grad(i)*c_sym;
    
end

end

% Функция подстановки новых функций в основную функцию
% ВХ. АРГУМЕНТЫ:
    % vec_new_func - вектор с функциями, которые необходимо подставить
    % vec_sym - вектор, который содержит в себе сиволы, которые содержит в себе основная функция
    % func - основная функция
% ВЫХ. АРГУМЕНТЫ:
	% new_func - основная функция, в которую вставили новые функции
function new_func = InsertIntoMaineFunc(vec_new_func, vec_sym, func)

new_func = func;

% Цикл перебора всех символов
for i = 1:length(vec_sym)
    
    % Вставляем новую функцию вместо определенного символа в основную функцию
    new_func = subs(new_func, vec_sym(i), vec_new_func(i));
    
end

end

% функция вычисления новых координат
% ВХ. АРГУМЕНТЫ:
    % old_coor - старые координаты
    % vec_sym - Вектор, который содержит символы, которые имеются в функции
    % vec_diff - Частные производные функции (сам градиент)
    % vec_sym_diff - оставшиеся символы в частных производных
    % func - основная функция
% ВЫХ. АРГУМЕНТЫ:
    % coor - новые координаты
function coor = CalcNewCoord(old_coor, vec_sym, vec_diff, vec_sym_diff, func)

syms h
% Вычисляем значение градиентов
coor_grad_1 = CalcCoorGrad(old_coor, vec_sym, vec_diff, vec_sym_diff);
% Формируем новые функции
vec_new_func = CreateNewFunc(old_coor, coor_grad_1, h);
% Подставляем в основную функцию сформированные выше
new_func = InsertIntoMaineFunc(vec_new_func, vec_sym, func);
% Решаем уравнение
res_h = diff(new_func); % Берем производную т.к. по алгоритму сначала берем ее, а потом решаем
res_h = vpasolve(res_h==0);
res_h = eval(res_h(1));
% Подставляем в новые функции полученное значение h
new_coor_1 = subs(vec_new_func, h, res_h);
coor = eval(new_coor_1);

end
