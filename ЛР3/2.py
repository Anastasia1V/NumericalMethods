import matplotlib.pyplot as plt


GRAPH_FILE = '2.png'
OUTPUT_FILE = 'output2.txt'


def solve_tridiagonal(a, b, c, d):
    size = len(b)
    p = [0] * size
    q = [0] * size
    if size == 0:
        raise ValueError('Система уравнений пуста')
    if abs(b[0]) < 1e-7:
        raise ValueError('Нулевой знаменатель в методе прогонки')
    if size > 1:
        p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]
    for i in range(1, size):
        denominator = b[i] + a[i] * p[i - 1]
        if abs(denominator) < 1e-7:
            raise ValueError('Нулевой знаменатель в методе прогонки')
        if i < size - 1:
            p[i] = -c[i] / denominator
        q[i] = (d[i] - a[i] * q[i - 1]) / denominator
    x = [0] * size
    x[size - 1] = q[size - 1]
    for i in range(size - 2, -1, -1):
        x[i] = p[i] * x[i + 1] + q[i]
    return x


def find_spline_coefficients(x_nodes, y_nodes):
    number_of_intervals = len(x_nodes) - 1
    h = list()
    slopes = list()
    for i in range(number_of_intervals):
        h.append(x_nodes[i + 1] - x_nodes[i])
        slopes.append((y_nodes[i + 1] - y_nodes[i]) / h[i])
    size = number_of_intervals - 1
    a = [0] * size
    b = [0] * size
    c = [0] * size
    d = [0] * size
    for i in range(size):
        node_number = i + 2
        h_left = h[node_number - 2]
        h_right = h[node_number - 1]
        b[i] = 2 * (h_left + h_right)
        if i > 0:
            a[i] = h_left
        if i < size - 1:
            c[i] = h_right
        d[i] = 3 * (slopes[node_number - 1] - slopes[node_number - 2])
    c_middle = solve_tridiagonal(a, b, c, d)
    c_values = [0] + c_middle
    a_values = list()
    b_values = list()
    d_values = list()
    for i in range(number_of_intervals):
        a_values.append(y_nodes[i])
        if i < number_of_intervals - 1:
            b_values.append(
                slopes[i] - h[i] * (c_values[i + 1] + 2 * c_values[i]) / 3
            )
            d_values.append(
                (c_values[i + 1] - c_values[i]) / (3 * h[i])
            )
        else:
            b_values.append(slopes[i] - 2 * h[i] * c_values[i] / 3)
            d_values.append(-c_values[i] / (3 * h[i]))
    return a_values, b_values, c_values, d_values, h


def find_interval(x_nodes, x):
    for i in range(len(x_nodes) - 1):
        if x_nodes[i] <= x <= x_nodes[i + 1]:
            return i
    raise ValueError('Точка x* находится вне отрезка интерполяции')


def evaluate_spline(x_nodes, a, b, c, d, x):
    interval = find_interval(x_nodes, x)
    delta_x = x - x_nodes[interval]
    return (
        a[interval] + b[interval] * delta_x
        + c[interval] * delta_x ** 2
        + d[interval] * delta_x ** 3
    )


def make_graph(x_nodes, y_nodes, a, b, c, d, h, x_star, y_star):
    figure, axes = plt.subplots(figsize=(10, 7))
    colors = ['b', 'g', 'b', 'g', 'm', 'y', 'k']
    for i in range(len(a)):
        list_x = list()
        list_y = list()
        points_count = 1000
        for j in range(points_count + 1):
            x = x_nodes[i] + h[i] * j / points_count
            delta_x = x - x_nodes[i]
            y = (
                a[i] + b[i] * delta_x + c[i] * delta_x ** 2
                + d[i] * delta_x ** 3
            )
            list_x.append(x)
            list_y.append(y)
        axes.plot(
            list_x, list_y, color=colors[i % len(colors)],
            label=f'Сплайн на [{x_nodes[i]:.3g}, {x_nodes[i + 1]:.3g}]'
        )
    axes.scatter(
        x_nodes, y_nodes, color='black', marker='o', zorder=4,
        label='Исходные узлы'
    )
    axes.scatter(
        [x_star], [y_star], color='red', marker='o', s=100, zorder=5,
        label='Точка (x*, S(x*))'
    )
    axes.set_title('Кубический сплайн')
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(GRAPH_FILE)
    plt.close(figure)


def main():
    try:
        with open('input2.txt', 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input2.txt')
    try:
        x_star = float(lines[0].strip())
        indices = list()
        for number in lines[1].split():
            indices.append(int(number))
        x_nodes = list()
        for number in lines[2].split():
            x_nodes.append(float(number))
        y_nodes = list()
        for number in lines[3].split():
            y_nodes.append(float(number))
    except Exception as e:
        raise ValueError(f'Некорректный input2.txt: {e}')
    if len(x_nodes) != len(y_nodes) or len(indices) != len(x_nodes):
        raise ValueError('Некорректные данные в input2.txt')
    a, b, c, d, h = find_spline_coefficients(x_nodes, y_nodes)
    interval = find_interval(x_nodes, x_star)
    spline_value = evaluate_spline(x_nodes, a, b, c, d, x_star)
    node_values = list()
    for x in x_nodes:
        node_values.append(evaluate_spline(x_nodes, a, b, c, d, x))
    left_second_derivative = 2 * c[0]
    last_interval = len(a) - 1
    right_second_derivative = (
        2 * c[last_interval] + 6 * d[last_interval] * h[last_interval]
    )
    make_graph(
        x_nodes, y_nodes, a, b, c, d, h, x_star, spline_value
    )
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as file:
        file.write(
            'Кубический сплайн с нулевой кривизной на концах отрезка' + '\n')
        file.write('Исходные узлы интерполяции' + '\n')
        file.write('i  x_i  f_i' + '\n')
        for i in range(len(x_nodes)):
            file.write(
                f'{indices[i]}  {x_nodes[i]:.3g}  {y_nodes[i]:.3g}' + '\n'
            )
        file.write('\n')
        file.write('Коэффициенты и интервалы сплайна' + '\n')
        file.write('i  [x_(i-1), x_i]  a_i  b_i  c_i  d_i' + '\n')
        for i in range(len(a)):
            file.write(
                f'{i + 1}  [{x_nodes[i]:.3g}, {x_nodes[i + 1]:.3g}]  '
                f'{a[i]:.3g}  {b[i]:.3g}  {c[i]:.3g}  {d[i]:.3g}' + '\n'
            )
        file.write('\n')
        file.write('Многочлены сплайна на каждом интервале' + '\n')
        for i in range(len(a)):
            file.write(
                f'S_{i + 1}(x) = {a[i]:.3g} '
                f'+ ({b[i]:.3g}) * (x - {x_nodes[i]:.3g}) '
                f'+ ({c[i]:.3g}) * (x - {x_nodes[i]:.3g})^2 '
                f'+ ({d[i]:.3g}) * (x - {x_nodes[i]:.3g})^3, '
                f'x в [{x_nodes[i]:.3g}, {x_nodes[i + 1]:.3g}]' + '\n'
            )
        file.write('\n')
        file.write('Значение сплайна в заданной точке' + '\n')
        file.write(f'x* = {x_star:.3g}' + '\n')
        file.write(f'Номер интервала: {interval + 1}' + '\n')
        file.write(f'S(x*) = {spline_value:.3g}' + '\n\n')
        file.write('Проверка прохождения сплайна через исходные узлы' + '\n')
        file.write('i  x_i  f_i  S(x_i)  |S(x_i) - f_i|' + '\n')
        for i in range(len(x_nodes)):
            file.write(
                f'{indices[i]}  {x_nodes[i]:.3g}  {y_nodes[i]:.3g}  '
                f'{node_values[i]:.3g}  '
                f'{abs(node_values[i] - y_nodes[i]):.3g}' + '\n'
            )
        file.write('\n')
        file.write('Значения второй производной на концах отрезка' + '\n')
        file.write(f'S"(x0) = {left_second_derivative:.3g}' + '\n')
        file.write(f'S"(x4) = {right_second_derivative:.3g}' + '\n')


if __name__ == '__main__':
    main()
