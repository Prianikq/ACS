[x, y, z] = ndgrid(-10:0.1:10, -10:0.1:10, -10:0.1:10);
u = 7*x.^2 - 14*x.*y + 9*y.^2 + 8*x.*z - 11*z.^2 + x + y + z - 10;

A = [ 
      7, -7, 4;
      -7, 9, 0;
      4, 0, -11 
      ];
b = [ 0.5, 0.5, 0.5 ];
X = [ 'x', 'y', 'z' ];
b = transpose(b);
a0 = -10;

solutions = linsolve(A, -b);
disp("Решение системы Ax+b=0: ");
for i = 1:rows(solutions)
  disp(cstrcat(X(i), " = ", num2str(solutions(i, 1))));
endfor

% Поскольку решение системы  Ax+b=0  существует,
% то данная система является совместной, 
% что приводит к 1-му случаю решения:

a0New = transpose(solutions) * b + a0;
disp(cstrcat("\nnew a0 = ", num2str(a0New)));

% Нахождение собственных значений

charPoly = poly(A);
valslin = roots(charPoly);

[vecs, vals] = eig(A, eigvalOption = "vector");

vars = sort(vals);
valslin = sort(valslin);

function res = printEigenvals(vals)
  for i = 1:rows(vals)
    disp(cstrcat("c. з. ", num2str(i), " =  ", num2str(vals(i, 1))));
  endfor
endfunction

disp("\nСобственные значения, найденные с помощью метода root: ");
printEigenvals(vals);

disp("\nСобственные значения, найденные с помощью метода eig: ");
printEigenvals(valslin);

precision = 4;

if mat2str(vals, precision) != mat2str(valslin, precision)
  disp("Собственные значения, найденные с помощью методов roots и eig не совпадают!");
else
  disp("Собственные значения, найденные с помощью методов roots и eig совпадают.");
endif

% Далее необходимо найти собственные векторы матрицы A,
% которые будут являться направляющими векторами для осей в новой системе координат:

disp("\nНайденные собственные векторы: ");
for i = 1:rows(vals)
    disp(cstrcat("для с.з. ", num2str(vals(i, 1)), ": vec = ", mat2str(transpose(vecs(:, i)), precision)));
endfor

%  Получено следующее каноническое уравнение в новой системе координат:
uCanonical = (vals(1)*x.^2 + vals(2)*y.^2 + vals(3)*z.^2 + a0New) / a0New;
isosurface(u, 0);
figure;
isosurface(uCanonical, 0);

