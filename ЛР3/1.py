import math
import re
import matplotlib.pyplot as plt


OUTPUT_FILE = 'output1.txt'


def read_math(str):
    str = str.replace(' ', '')
    str = str.replace('pi', 'math.pi')
    str = re.sub(r'(\d)(math\.pi)', r'\1*\2', str)
    return eval(str, {'math': math, '__builtins__': {}})


def function(x):
    return math.cos(x) / math.sin(x)


def lagrange(x_nodes, y_nodes, x):
    result = 0
    for i in range(len(x_nodes)):
        s = y_nodes[i]
        for j in range(len(x_nodes)):
            if i != j:
                s *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        result += s
    return result


def divided_differences(x_nodes, y_nodes):
    size = len(x_nodes)
    differences = list()
    for i in range(size):
        differences.append([0] * size)
        differences[i][0] = y_nodes[i]
    for order in range(1, size):
        for i in range(size - order):
            differences[i][order] = (
                differences[i + 1][order - 1] - differences[i][order - 1]
            ) / (x_nodes[i + order] - x_nodes[i])
    return differences


def newton(x_nodes, differences, x):
    result = differences[0][0]
    for order in range(1, len(x_nodes)):
        s = differences[0][order]
        for i in range(order):
            s *= x - x_nodes[i]
        result += s
    return result


def lagrange_polynomial_to_str(x_nodes, y_nodes):
    polynomial = ''
    for i in range(len(x_nodes)):
        omega_derivative = 1
        for j in range(len(x_nodes)):
            if i != j:
                omega_derivative *= x_nodes[i] - x_nodes[j]
        coefficient = y_nodes[i] / omega_derivative
        if abs(coefficient) < 1e-12:
            continue
        polynomial_part = f'{abs(coefficient):.3g}'
        for j in range(len(x_nodes)):
            if i != j:
                polynomial_part += f' * (x - {x_nodes[j]:.3g})'
        if len(polynomial) == 0:
            if coefficient < 0:
                polynomial = '-'
            polynomial += polynomial_part
        elif coefficient < 0:
            polynomial += ' - ' + polynomial_part
        else:
            polynomial += ' + ' + polynomial_part
    if len(polynomial) == 0:
        polynomial = '0'
    return polynomial


def newton_polynomial_to_str(x_nodes, differences):
    polynomial = f'{differences[0][0]:.3g}'
    for order in range(1, len(x_nodes)):
        coefficient = differences[0][order]
        if abs(coefficient) < 1e-12:
            continue
        polynomial_part = f'{abs(coefficient):.3g}'
        for i in range(order):
            polynomial_part += f' * (x - {x_nodes[i]:.3g})'
        if coefficient < 0:
            polynomial += ' - ' + polynomial_part
        else:
            polynomial += ' + ' + polynomial_part
    return polynomial


def divided_differences_to_str(x_nodes, y_nodes, differences):
    output = list()
    header = ['i', 'x_i', 'f(x_i)']
    for order in range(1, len(x_nodes)):
        header.append(f'{order}-й порядок')
    output.append('  '.join(header))
    for i in range(len(x_nodes)):
        line = [
            str(i), f'{x_nodes[i]:.3g}', f'{y_nodes[i]:.3g}'
        ]
        for order in range(1, len(x_nodes)):
            if i + order < len(x_nodes):
                line.append(f'{differences[i][order]:.3g}')
            else:
                line.append('-')
        output.append('  '.join(line))
    return '\n'.join(output)


def make_graph(result, x_star, file_name, title_suffix):
    x_nodes = result[0]
    y_nodes = result[1]
    differences = result[2]
    left_border = x_nodes[0]
    right_border = x_nodes[-1]
    margin = (right_border - left_border) * 0.2
    left_border -= margin
    right_border += margin
    points_count = 1000
    list_x = list()
    for i in range(points_count + 1):
        x = left_border + (right_border - left_border) * i / points_count
        list_x.append(x)
    list_y = list()
    for x in list_x:
        list_y.append(function(x))
    lagrange_y = list()
    newton_y = list()
    for x in list_x:
        lagrange_y.append(lagrange(x_nodes, y_nodes, x))
        newton_y.append(newton(x_nodes, differences, x))
    figure, axes = plt.subplots(figsize=(10, 7))
    axes.plot(list_x, list_y, color='black', label='y = ctg(x)')
    axes.plot(list_x, lagrange_y, color='tab:blue', linestyle='-',
              linewidth=2, label='Многочлен Лагранжа')
    axes.plot(list_x, newton_y, color='tab:green', linestyle='--',
              linewidth=1.5, label='Многочлен Ньютона')
    axes.scatter(x_nodes, y_nodes, color='black', marker='o',
                 zorder=4, label='Узлы интерполяции')
    axes.scatter([x_star], [function(x_star)], color='red', marker='o',
                 zorder=5, label='Точка x*')
    axes.set_title(f'Интерполяция y = ctg(x), {title_suffix})')
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(file_name)
    plt.close(figure)


def main():
    try:
        with open('input1.txt', 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input1.txt')
    try:
        first_x_nodes = list()
        for expression in lines[0].split():
            first_x_nodes.append(read_math(expression))
        second_x_nodes = list()
        for expression in lines[1].split():
            second_x_nodes.append(read_math(expression))
        x_star = read_math(lines[2].strip())
    except Exception as e:
        raise ValueError(f'Некорректный input1.txt: {e}')
    results = list()
    node_sets = [first_x_nodes, second_x_nodes]
    for x_nodes in node_sets:
        if len(x_nodes) != 4:
            raise ValueError('Для интерполяции должны быть заданы четыре узла')
        y_nodes = list()
        for x in x_nodes:
            y_nodes.append(function(x))
        differences = divided_differences(x_nodes, y_nodes)
        lagrange_value = lagrange(x_nodes, y_nodes, x_star)
        newton_value = newton(x_nodes, differences, x_star)
        exact_value = function(x_star)
        results.append((
            x_nodes, y_nodes, differences, lagrange_value, newton_value,
            exact_value
        ))
    make_graph(results[0], x_star, '1a.png', 'а')
    make_graph(results[1], x_star, '1b.png', 'б')
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as file:
        file.write('Интерполяционные многочлены Лагранжа и Ньютона' + '\n')
        file.write('Функция: y = ctg(x)' + '\n')
        file.write(f'Заданная точка x* = {x_star:.3g}' + '\n')
        file.write(f'Точное значение f(x*) = {function(x_star):.3g}' + '\n\n')
        for i in range(len(results)):
            x_nodes, y_nodes = results[i][0], results[i][1]
            differences = results[i][2]
            lagrange_value, newton_value = results[i][3], results[i][4]
            exact_value = results[i][5]
            file.write(f'Набор узлов {"а" if i == 0 else "б"}' + '\n')
            file.write('Исходные узлы интерполяции' + '\n')
            file.write('i  x_i  y_i' + '\n')
            for j in range(len(x_nodes)):
                file.write(
                    f'{j}  {x_nodes[j]:.3g}  {y_nodes[j]:.3g}' + '\n'
                )
            file.write('\n')
            file.write('Интерполяционный многочлен Лагранжа' + '\n')
            file.write('L3(x) = ' + lagrange_polynomial_to_str(
                x_nodes, y_nodes
            ) + '\n\n')
            file.write('Таблица разделённых разностей' + '\n')
            file.write(divided_differences_to_str(
                x_nodes, y_nodes, differences
            ) + '\n\n')
            file.write('Интерполяционный многочлен Ньютона' + '\n')
            file.write('P3(x) = ' + newton_polynomial_to_str(
                x_nodes, differences
            ) + '\n\n')
            file.write(
                f'Значение многочлена Лагранжа в точке x*: '
                f'{lagrange_value:.3g}' + '\n'
            )
            file.write(
                f'Значение многочлена Ньютона в точке x*: '
                f'{newton_value:.3g}' + '\n'
            )
            file.write(f'Точное значение функции f(x*): {exact_value:.3g}')
            file.write('\n')
            file.write(
                f'Абсолютная погрешность интерполяции Лагранжа: '
                f'{abs(lagrange_value - exact_value):.3g}' + '\n'
            )
            file.write(
                f'Абсолютная погрешность интерполяции Ньютона: '
                f'{abs(newton_value - exact_value):.3g}' + '\n'
            )
            file.write(
                f'Разность между результатами методов: '
                f'{lagrange_value - newton_value:.3g}' + '\n'
            )
            file.write(
                f'Абсолютная разность между результатами методов: '
                f'{abs(lagrange_value - newton_value):.3g}' + '\n'
            )


if __name__ == '__main__':
    main()
