import matplotlib.pyplot as plt


GRAPH_FILE = '3.png'
OUTPUT_FILE = 'output3.txt'


def solve_linear_system(a_matrix, b_vector):
    size = len(b_vector)
    matrix = list()
    for i in range(size):
        line = list(a_matrix[i])
        line.append(b_vector[i])
        matrix.append(line)
    for i in range(size):
        lead_line = i
        max_value = abs(matrix[i][i])
        for j in range(i + 1, size):
            value = abs(matrix[j][i])
            if value > max_value:
                max_value = value
                lead_line = j
        if max_value < 1e-7:
            raise ValueError('Нормальная система вырождена')
        if lead_line != i:
            matrix[i], matrix[lead_line] = matrix[lead_line], matrix[i]
        for j in range(i + 1, size):
            multiplier = matrix[j][i] / matrix[i][i]
            for k in range(i, size + 1):
                matrix[j][k] -= multiplier * matrix[i][k]
    x = [0] * size
    for i in range(size - 1, -1, -1):
        s = 0
        for j in range(i + 1, size):
            s += matrix[i][j] * x[j]
        x[i] = (matrix[i][size] - s) / matrix[i][i]
    return x


def normal_system(x_nodes, y_nodes, degree):
    size = degree + 1
    a_matrix = list()
    b_vector = list()
    for i in range(size):
        line = list()
        for j in range(size):
            s = 0
            for x in x_nodes:
                s += x ** (i + j)
            line.append(s)
        a_matrix.append(line)
        s = 0
        for j in range(len(x_nodes)):
            s += y_nodes[j] * x_nodes[j] ** i
        b_vector.append(s)
    return a_matrix, b_vector


def polynomial(coefficients, x):
    result = 0
    for i in range(len(coefficients) - 1, -1, -1):
        result = result * x + coefficients[i]
    return result


def squared_errors(x_nodes, y_nodes, coefficients):
    phi = 0
    for i in range(len(x_nodes)):
        difference = polynomial(coefficients, x_nodes[i]) - y_nodes[i]
        phi += difference ** 2
    return phi


def polynomial_to_str(coefficients):
    polynomial = ''
    for i in range(len(coefficients)):
        coefficient = coefficients[i]
        if abs(coefficient) < 1e-7:
            continue
        polynomial_part = f'{abs(coefficient):.3g}'
        if i == 1:
            polynomial_part += ' * x'
        elif i > 1:
            polynomial_part += f' * x^{i}'
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


def find_approximation(x_nodes, y_nodes, degree):
    a_matrix, b_vector = normal_system(x_nodes, y_nodes, degree)
    coefficients = solve_linear_system(a_matrix, b_vector)
    phi = squared_errors(x_nodes, y_nodes, coefficients)
    mse = phi / len(x_nodes)
    return coefficients, phi, mse


def make_graph(x_nodes, y_nodes, approximations):
    left_border = x_nodes[0]
    right_border = x_nodes[-1]
    points_count = 1000
    list_x = list()
    for i in range(points_count + 1):
        x = left_border + (right_border - left_border) * i / points_count
        list_x.append(x)
    figure, axes = plt.subplots(figsize=(10, 7))
    axes.scatter(
        x_nodes, y_nodes, color='black', marker='o', zorder=4,
        label='Исходные данные'
    )
    colors = ['black', 'b', 'g', 'c', 'm', 'y', 'k']
    for i in range(len(approximations)):
        coefficients = approximations[i][0]
        list_y = list()
        for x in list_x:
            list_y.append(polynomial(coefficients, x))
        axes.plot(
            list_x, list_y, color=colors[i % len(colors)],
            label=f'Многочлен степени {i + 1}'
        )
    axes.set_title('Аппроксимация методом наименьших квадратов')
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(GRAPH_FILE)
    plt.close(figure)


def main():
    try:
        with open('input3.txt', 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input3.txt')
    try:
        indices = list()
        for number in lines[0].split():
            indices.append(int(number))
        x_nodes = list()
        for number in lines[1].split():
            x_nodes.append(float(number))
        y_nodes = list()
        for number in lines[2].split():
            y_nodes.append(float(number))
    except Exception as e:
        raise ValueError(f'Некорректный input3.txt: {e}')
    if len(indices) != len(x_nodes) or len(x_nodes) != len(y_nodes):
        raise ValueError('Некорректные данные в input3.txt')
    approximations = list()
    for degree in range(1, 4):
        approximations.append(find_approximation(x_nodes, y_nodes, degree))
    make_graph(x_nodes, y_nodes, approximations)
    best_degree = 1
    for i in range(1, len(approximations)):
        if approximations[i][1] < approximations[best_degree - 1][1]:
            best_degree = i + 1
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as file:
        file.write('Метод наименьших квадратов' + '\n')
        file.write('Исходные данные' + '\n')
        file.write('i  x_i  y_i' + '\n')
        for i in range(len(x_nodes)):
            file.write(
                f'{indices[i]}  {x_nodes[i]:.3g}  {y_nodes[i]:.3g}' + '\n'
            )
        file.write('\n')
        for i in range(len(approximations)):
            coefficients, phi, mse = approximations[i]
            degree = i + 1
            file.write(f'Коэффициенты многочлена степени {degree}' + '\n')
            for j in range(len(coefficients)):
                file.write(f'a{j} = {coefficients[j]:.3g}' + '\n')
            file.write(
                f'P{degree}(x) = {polynomial_to_str(coefficients)}' + '\n'
            )
            file.write(f'Сумма квадратов ошибок Φ{degree} = {phi:.3g}' + '\n')
            file.write(f'MSE{degree} = {mse:.3g}' + '\n\n')
        file.write('Сравнение полученных ошибок' + '\n')
        for i in range(len(approximations)):
            file.write(
                f'Φ{i + 1} = {approximations[i][1]:.3g}, '
                f'MSE{i + 1} = {approximations[i][2]:.3g}' + '\n'
            )
        file.write(
            f'Наилучшее приближение на исходных точках даёт '
            f'многочлен степени {best_degree}.' + '\n'
        )


if __name__ == '__main__':
    main()
