Программа находит совершенное паросочетание минимального веса.

В ней реализован функционал callback, который подсказывает солверу новое ограничение в случае, если решение получилось дробным. 

Формат ввода:

В первой строке заданы два целых числа n и m (2 ≤ n ≤ 500, 1 ≤ m ≤ 1000).

В последующих m строках заданы рёбра графа. Каждое ребро задаётся тройкой целых чисел ai, bi, wi (0 ≤ ai, bi ≤ n - 1, 1 ≤ wi ≤ 1000), где ai и bi — номера вершин, соединяемых ребром, а wi — его вес.

Формат вывода:

В первой строке минимально возможный вес совершенного паросочетания.

Во второй строке  n / 2 целых чисел — номера рёбер в совершенном паросочетании минимального веса. Рёбра пронумерованы от 0 до m-1 в порядке следования во входных данных.

