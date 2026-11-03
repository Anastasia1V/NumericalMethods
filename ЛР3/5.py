import math


OUTPUT_FILE = 'output5.txt'


def function(x):
    return (3 * x + 4) / (2 * x + 7)


def find_number_of_intervals(left_border, right_border, h):
    return round((right_border - left_border) / h)


def integrate_left(left_border, right_border, h):
    number = find_number_of_intervals(left_border, right_border, h)
    result = 0
    for i in range(number):
        x = left_border + i * h
        result += function(x)
    return h * result


def integrate_right(left_border, right_border, h):
    number = find_number_of_intervals(left_border, right_border, h)
    result = 0
    for i in range(1, number + 1):
        x = left_border + i * h
        result += function(x)
    return h * result


def integrate_middle(left_border, right_border, h):
    number = find_number_of_intervals(left_border, right_border, h)
    result = 0
    for i in range(number):
        x = left_border + (i + 0.5) * h
        result += function(x)
    return h * result


def integrate_trapezoid(left_border, right_border, h):
    number = find_number_of_intervals(left_border, right_border, h)
    result = (function(left_border) + function(right_border)) / 2
    for i in range(1, number):
        x = left_border + i * h
        result += function(x)
    return h * result


def integrate_simpson(left_border, right_border, h):
    number = find_number_of_intervals(left_border, right_border, h)
    result = function(left_border) + function(right_border)
    for i in range(1, number):
        x = left_border + i * h
        if i % 2 == 0:
            result += 2 * function(x)
        else:
            result += 4 * function(x)
    return h * result / 3


def runge_romberg(fine_value, coarse_value, order, ratio):
    return fine_value + (fine_value - coarse_value) / (ratio ** order - 1)


def estimate_error(fine_value, coarse_value, order, ratio):
    return abs(fine_value - coarse_value) / (ratio ** order - 1)


def integral_primitive(x):
    return 1.5 * x - 3.25 * math.log(abs(2 * x + 7))


def exact_integral(left_border, right_border):
    return integral_primitive(right_border) - integral_primitive(left_border)


def main():
    try:
        with open('input5.txt', 'r', encoding='utf-8') as file:
            values = list()
            for number in file.read().split():
                values.append(float(number))
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input5.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input5.txt: {e}')
    left_border = values[0]
    right_border = values[1]
    h1 = values[2]
    h2 = values[3]
    if h1 <= 0 or h2 <= 0:
        raise ValueError('Шаги интегрирования должны быть положительными')
    if right_border <= left_border:
        raise ValueError('Правый край должен быть больше левого')
    number_h1 = find_number_of_intervals(left_border, right_border, h1)
    number_h2 = find_number_of_intervals(left_border, right_border, h2)
    if abs(number_h1 * h1 - (right_border - left_border)) > 1e-7:
        raise ValueError('Шаг h1 должен точно разбивать отрезок')
    if abs(number_h2 * h2 - (right_border - left_border)) > 1e-7:
        raise ValueError('Шаг h2 должен точно разбивать отрезок')
    ratio = h1 / h2
    if ratio <= 1 or abs(ratio - round(ratio)) > 1e-7:
        raise ValueError('Шаг h1 должен быть кратен шагу h2')
    if number_h1 % 2 != 0 or number_h2 % 2 != 0:
        raise ValueError('Число интервалов должно быть чётным')
    if left_border <= -3.5 <= right_border:
        raise ValueError('Подынтегральная функция имеет разрыв на отрезке')
    methods = [
        ('Левые прямоугольники', integrate_left, 1),
        ('Правые прямоугольники', integrate_right, 1),
        ('Средние прямоугольники', integrate_middle, 2),
        ('Трапеции', integrate_trapezoid, 2),
        ('Симпсона', integrate_simpson, 4)
    ]
    results = list()
    for name, method, order in methods:
        value_h1 = method(left_border, right_border, h1)
        value_h2 = method(left_border, right_border, h2)
        refined_value = runge_romberg(value_h2, value_h1, order, ratio)
        error_estimate = estimate_error(value_h2, value_h1, order, ratio)
        results.append((
            name, value_h1, value_h2, order,
            refined_value, error_estimate
        ))
    exact_value = exact_integral(left_border, right_border)
    absolute_errors = list()
    for result in results:
        absolute_errors.append(abs(result[4] - exact_value))
    best_index = 0
    for i in range(1, len(results)):
        if absolute_errors[i] < absolute_errors[best_index]:
            best_index = i
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as file:
        file.write('Численное интегрирование' + '\n')
        file.write('y = (3x + 4) / (2x + 7)' + '\n')
        file.write(
            f'Пределы интегрирования: [{left_border:.3g}, '
            f'{right_border:.3g}]' + '\n'
        )
        file.write(f'h1 = {h1:.3g}, h2 = {h2:.3g}' + '\n')
        file.write(f'Число интервалов при h1: {number_h1}' + '\n')
        file.write(f'Число интервалов при h2: {number_h2}' + '\n\n')
        file.write('Результаты вычисления' + '\n')
        file.write(
            'Метод | F(h1) | F(h2) | p | F по Рунге–Ромбергу | '
            'Оценка погрешности при h2 | Абсолютная погрешность' + '\n'
        )
        for i in range(len(results)):
            result = results[i]
            file.write(
                f'{result[0]} | {result[1]:.3g} | {result[2]:.3g} | '
                f'{result[3]} | {result[4]:.3g} | {result[5]:.3g} | '
                f'{absolute_errors[i]:.3g}' + '\n'
            )
        file.write('\n')
        file.write(
            'Первообразная: F(x) = 3/2 * x - 13/4 * ln|2x + 7|' + '\n'
        )
        file.write('Точное значение интеграла' + '\n')
        file.write(f'F = {exact_value:.3g}' + '\n\n')
        file.write(
            f'Наиболее точный результат после уточнения получен методом '
            f'{results[best_index][0]}.' + '\n'
        )
        file.write(
            f'Абсолютная погрешность: {absolute_errors[best_index]:.3g}'
            + '\n'
        )


if __name__ == '__main__':
    main()
