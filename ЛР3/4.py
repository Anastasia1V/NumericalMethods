OUTPUT_FILE = 'output4.txt'


def find_node_index(x_nodes, x):
    for i in range(len(x_nodes)):
        if abs(x_nodes[i] - x) < 1e-7:
            return i
    raise ValueError('Точка x* должна совпадать с внутренним узлом')


def find_derivatives(x_nodes, y_nodes, x):
    index = find_node_index(x_nodes, x)
    left_index = index - 1
    right_index = index + 1
    x_left = x_nodes[left_index]
    x_middle = x_nodes[index]
    x_right = x_nodes[right_index]
    y_left = y_nodes[left_index]
    y_middle = y_nodes[index]
    y_right = y_nodes[right_index]
    left_derivative = (y_middle - y_left) / (x_middle - x_left)
    right_derivative = (y_right - y_middle) / (x_right - x_middle)
    first_derivative = (
        left_derivative + (right_derivative - left_derivative)
        / (x_right - x_left) * (2 * x - x_left - x_middle)
    )
    second_derivative = (
        2 * (right_derivative - left_derivative)
        / (x_right - x_left)
    )
    return (
        index, left_derivative, right_derivative,
        first_derivative, second_derivative
    )


def main():
    try:
        with open('input4.txt', 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input4.txt')
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
        raise ValueError(f'Некорректный input4.txt: {e}')
    if len(indices) != len(x_nodes) or len(x_nodes) != len(y_nodes):
        raise ValueError('Некорректные данные в input4.txt')
    (index, left_derivative, right_derivative,
     first_derivative, second_derivative) = find_derivatives(
        x_nodes, y_nodes, x_star
    )
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as file:
        file.write('Численное дифференцирование' + '\n')
        file.write('Исходная таблица значений функции' + '\n')
        file.write('i  x_i  y_i' + '\n')
        for i in range(len(x_nodes)):
            file.write(
                f'{indices[i]}  {x_nodes[i]:.3g}  {y_nodes[i]:.3g}' + '\n'
            )
        file.write('\n')
        file.write(f'Точка вычисления x* = {x_star:.3g}' + '\n')
        file.write(
            'Узлы таблицы, использованные при вычислении производных' + '\n')
        for i in range(index - 1, index + 2):
            file.write(
                f'x_{indices[i]} = {x_nodes[i]:.3g}, '
                f'y_{indices[i]} = {y_nodes[i]:.3g}' + '\n'
            )
        file.write('\n')
        file.write('Первая производная по формуле (3.18)' + '\n')
        file.write(
            f'Левосторонняя производная: {left_derivative:.3g}' + '\n'
        )
        file.write(
            f'Правосторонняя производная: {right_derivative:.3g}' + '\n'
        )
        file.write(
            f"Первая производная по формуле (3.20): y'({x_star:.3g}) = "
            f'{first_derivative:.3g}' + '\n\n'
        )
        file.write('Вторая производная по формуле (3.21)' + '\n')
        file.write(
            f"y''({x_star:.3g}) = {second_derivative:.3g}" + '\n'
        )


if __name__ == '__main__':
    main()
