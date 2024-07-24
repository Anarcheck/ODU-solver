
Задаётся уравнение
y1''+p1(x)*y1'+q11(x)*y1+q12(x)y2=f1(x)

y2''+p2(x)*y2'+q21(x)*y1+q22(x)y2=f2(x)

с условиями y1(a)=y1a y1(b)=y1b y2(a)=y2a y2(b)=y2b на отрезке [a;b]. Отрезок [a;b] в начале делиться на k частей. 
Погрешность вычисления не должна превышать величину e.


Имя входного файла:"input.txt"
Формат данных во входном файле:
s(к-во уравнений)
a(левый конец отрезка) b(правый конец) eps(точность) k(начальное к-во точек) y1a y1b y2a y2b(Краевые условия)

Функции f1(x),f2(x),p1(x),p2(x),q11(x),q12(x),q21(x),q22(x) задаются вручную в файле func.c
В случае одного уравнения задаются только функции f1(x),p1(x),q11(x).

Для построения графика точное решение задается 
func1_acc(задаётся в return), либо func2_acc(задаются в res[0]-значение первой функции, res[1]-точное значение второй функции), в зависимости количества уравнений.

Результат выполнения программы - график, записывается в файл graph.jpg 
и файл points.txt, в котором в начале записываются значения x, затем наше найденное решение(1 или 2 строки, в зависимости от количества уравнений) 
и затем точное решение, заданное в функциях func1_acc либо func2_acc(Также в зависимости от количества уравнений).

Ключ для компиляции программы 

gcc -Wall main.c odu.c func.c -lm -o main 

Вызов программы main.exe
