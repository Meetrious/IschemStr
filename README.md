# IschemStr

Проект призван описать инструменты для работы с математической моделью биохим. явления -- *процессов при острой фазе ишемического инсульта*. Подробнее о ней [здесь](https://www.overleaf.com/read/mtxqqwygypgp). Модель описана в терминах *задачи на обыкновенные дифференциальные уравнения (в дальнейшем ОДУ) первого порядка, разрешённые относительно производной*. 

Глобально репозиторий содержит реализацию двух программ: решение **прямой задачи** и решение **обратной задачи** в связи с первой.

## Прямая задача (*namespace StraightTask*)

### Постановка
Означенная состоит из *системы ОДУ с запаздывающими аргументами первого порядка*, набора начальных условий, и начальных множеств/функций для отклоняющихся аргументов. Путём использования *метода последовательного интегрирования* (в дальнейшем **метода шагов**) объявленная задача сводится к *задаче Коши*, приближенные решения которой могут быть получены с привлечением широко-известных численных методов, таких как например Рунге-Кутта, Адамса-Бэшфорда-Моултона(ABM), или Гира.

### Детали реализации
Директория [ST](projects/ST) содержит код реализации решения прямой задачи на системе из ОДУ с запаздываниями, а именно *метода шагов* поверх перечисленных выше методов решений задачи Коши.

// на очереди картинки с UML диаграммами, наглядно демонстрирующие структуру.

to be continued...

## Анализ чувствительности (*namespace SensAnalysisTask*)
Блок, чья задача численно высчитать относительную чувствительность решения от малейших изменений постоянных коэффициентов правой части уравнений системы.

### Постановка
Пусть есть некоторое *исходное* численное решение системы, которое нас устраивает (оно, как минимум, неперерывно на нужном нам конечном промежутке). И нас интересует вопрос: "Насколько сильно изментся его поведение (характер монотонности на промежутках, экстремумы, непрерывность), если внести малейшие изменения в коэффициенты правой части?" Изменения, конечно, можно вносить в несколько коэффициентов сразу, а после смотреть на изменение решений, что тоже даст нам некоторую информацию. Но мы же будем дёргать по одному коэффициенту за эксперимент, что, само собой, не всегда ответит нам на вопрос о влиянии совокупного изменения, но не всё же сразу.

### Детали реализации
Директория [SensAn](projects/SensAn) содержит код реализации.

#### Поверхностно на словах и пальцах
Выходные данные метода, который мы будем использовать - матрица *RS = (rs<sub>ij</sub>)* из *коэффициентов относительной чувствительности* i-того уравнения к изменению j-того коэффициента, где предварительно по (i,j) обеспечена некоторая индексация уравнений системы и коэффициентов правой части системы соответственно. Вносимые изменения в коэффициенты составляют ±5-20% от значений коэффициентов в *исходном* состоянии. *rs<sub>ij</sub>* главным образом определяется через разность решения i-того уравнения в *исходном* состоянии и изменённого решения i-того уравнения на некотором интересующем нас конечном промежутке через сеточную норму в l<sub>2</sub> (в средне-квадратичном смысле). Выбор знака вносимых изменений для j-того коэффициента зависит от того, с плюсом, или с минусом норма разностей решений была наибольшей.

#### Формально:
// Приведу картинку из LaTeX



## Обратная задача (*namespace ReverseTask*)

### Постановка
Означенная представляет из себя *задачу математического программирования* (a.k.a **задача оптимизации**). Её целевая функция - это функционал невязки между численным решением прямой задачи и экспериментальными данными. Её множество регулируемых параметров (относительно которых ищется минимум невязки) - это постоянные коэффициенты в уравнениях прямой задачи, а ограничения на их допустимое множество - "здравый смысл" и положительная определённость. Подробнее о ней [здесь](https://www.overleaf.com/read/fygdbgkzgxtz).


### Детали реализации 
Директория [RT](projects/RT) содержит код реализации решения поставленной покоэффициентной задачи, а именно вариации *стохастического метода* ***BGA***.

...

<!--
Input содержит необходимые входные данные для работы, такие как начальные условия, значения постоянных коэффициентов правой части уравнений по умолчанию...
-->

