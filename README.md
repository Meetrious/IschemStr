# IschemStr

Проект призван описать инструменты для работы с математической моделью биохим. явления -- *процессов при острой фазе ишемического инсульта*. Подробнее о ней [здесь](https://www.overleaf.com/read/mtxqqwygypgp). Модель описана в терминах *задачи на обыкновенные дифференциальные уравнения (в дальнейшем ОДУ) первого порядка, разрешённые относительно производной*. 

Глобально репозиторий содержит реализацию двух программ: решение **прямой задачи** и решение **обратной задачи** в связи с первой.

## Прямая задача (*namespace StraightTask*)

### Постановка
Означенная состоит из *системы ОДУ с запаздывающими аргументами первого порядка*, набора начальных условий, и начальных множеств/функций для отклоняющихся аргументов. Путём использования *метода последовательного интегрирования* (в дальнейшем **метода шагов**) объявленная задача сводится к *задаче Коши*, приближенные решения которой могут быть получены с привлечением широко-известных численных методов, таких как например Рунге-Кутта, Адамса-Бэшфорда-Моултона(ABM), или Гира.

### Детали реализации
Директория [ST](projects/ST/ST) содержит код реализации решения прямой задачи на системе из ОДУ с запаздываниями, а именно *метода шагов* поверх перечисленных выше методов решений задачи Коши.

// на очереди картинки с UML диаграммами, наглядно демонстрирующие структуру.

to be continued...

## Обратная задача (*namespace ReverseTask*)

### Постановка
Означенная представляет из себя *задачу математического программирования* (a.k.a **задача оптимизации**). Её целевая функция - это функционал невязки между численным решением прямой задачи и экспериментальными данными. Её множество регулируемых параметров (относительно которых ищется минимум невязки) - это постоянные коэффициенты в уравнениях прямой задачи, а ограничения на их допустимое множество - "здравый смысл" и положительная определённость. Подробнее о ней [здесь](https://www.overleaf.com/read/fygdbgkzgxtz).


### Детали реализации 
Директория ***Projects/RT*** содержит код реализации решения поставленной покоэффициентной задачи, а именно вариации *стохастического метода* ***BGA***.

...
<!--
Input содержит необходимые входные данные для работы, такие как начальные условия, значения постоянных коэффициентов правой части уравнений по умолчанию...
-->

