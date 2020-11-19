# Ischemic-Stroke-Inflammation-Model-Solver
Этот проект призван описать инструменты работы с математической моделью био-хим явления - процессов острой фазы при ишемическом инсульте. 

Глобально репозиторий содержит реализацию двух программ: решение прямой задачи и решение обратной задачи в связи с первой.

## Прямая задача (*namespace StraightTask*)

Означенная заключается в *системе обыкновенных дифференциальных уравнений с отклоняющимися аргументами*, набором начальных условий, начальных множеств и начальных функций.  

(I) Projects/ST содержит код реализации решения прямой задачи на системе из ОДУ с запаздываниями, а именно метод последовательного интегрирования поверх численных методов решений задачи Коши.

Мат. модель представляет из себя систему ОДУ с запаздываниями(или без). Уравнения конкретной системы описаны в Projects/ST/ST/includes/*models or tests*/equations.h

Поверхность алгоритма решения прямой задачи описана в Projects/ST/ST/includes/base/Solver.h

Хедеры в ST/ST/include/*Folder* организованы и включаются друг в друга в заданном порядке, который работает.
В ST/ST/ST_main.cpp включаются только те хедеры, что лежат на вершине ST/ST/includes.

Ссылка на отчёты и наработки: 
https://www.overleaf.com/read/mtxqqwygypgp

## Обратная задача (*namespace ReverseTask*)
(II) Projects/RT содержит код реализации решения обратной покоэффициентной задачи поставленной на систему из ОДУ с запаздываниями, а именно стохастический метод BGA.

Input содержит необходимые входные данные для работы, такие как начальные условия, значения постоянных коэффициентов правой части уравнений по умолчанию...

Ссылка на отчёты и наработки: 
https://www.overleaf.com/read/fygdbgkzgxtz
