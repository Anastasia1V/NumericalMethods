EPSILON = 0.01


def line_to_left_and_right(line):
    left, right = line.replace(' ', '').split('=')
    i = 0
    left_part = dict()
    while i < len(left):
        is_positive = 1
        if i < len(left) and left[i] == '+':
            i += 1
        if i < len(left) and left[i] == '-':
            is_positive = -1
            i += 1
        number_index = i
        while i < len(left) and left[i].isdigit():
            i += 1
        if i < len(left) and left[i] == '.':
            i += 1
            while i < len(left) and left[i].isdigit():
                i += 1
        str_number = left[number_index:i]
        if len(str_number) > 0:
            number = float(str_number)
        else:
            number = 1
        i += 1
        index = i
        while i < len(left) and left[i].isdigit():
            i += 1
        left_part[int(left[index:i])] = is_positive * number
    right_part = float(right.replace(',', '.'))
    return left_part, right_part


def lines_to_a_and_b(lines):
    max_x = 0
    parts = list()
    for line in lines:
        line = line.replace(' ', '')
        left_part, right_part = line_to_left_and_right(line)
        parts.append((left_part, right_part))
        for index_x in left_part:
            if index_x > max_x:
                max_x = index_x
    a_matrix = list()
    b = list()
    for left_part, right_part in parts:
        line = list()
        for i in range(1, max_x + 1):
            if i in left_part:
                line.append(left_part[i])
            else:
                line.append(0)
        a_matrix.append(line)
        b.append(right_part)
    return a_matrix, b


def correct_diagonal(a_matrix, b_vector):
    size = len(a_matrix)
    need_to_change = False
    for i in range(size):
        if a_matrix[i][i] == 0:
            row = -1
            for j in range(i + 1, size):
                if row == -1 and a_matrix[j][i] != 0:
                    row = j
                    break
            if row != -1:
                a_matrix[i], a_matrix[row] = (
                    a_matrix[row], a_matrix[i]
                )
                b_vector[i], b_vector[row] = (
                    b_vector[row], b_vector[i]
                )
            else:
                for j in range(i):
                    if row == -1 and a_matrix[j][i] != 0:
                        row = j
                        break
                if row != -1:
                    a_matrix[i], a_matrix[row] = (
                        a_matrix[row], a_matrix[i]
                    )
                    b_vector[i], b_vector[row] = (
                        b_vector[row], b_vector[i]
                    )
                    need_to_change = True
                else:
                    raise ValueError(f'Остался ноль в строке {i + 1}')
    return need_to_change


def find_a_and_b(a_matrix, b_vector):
    size = len(a_matrix)
    iterations = 0
    need_to_change = True
    while need_to_change:
        need_to_change = False
        iterations += 1
        if iterations > size:
            raise ValueError(
                f'Нули на главной диагонали спустя {size} проходов'
            )
        need_to_change = correct_diagonal(a_matrix, b_vector)
    a = list()
    for i in range(size):
        line = list()
        for j in range(size):
            if i == j:
                line.append(0)
            else:
                line.append(-a_matrix[i][j] / a_matrix[i][i])
        a.append(line)
    b = list()
    for i in range(size):
        b.append(b_vector[i] / a_matrix[i][i])
    return a, b


def vector_norm(vector):
    n = abs(vector[0])
    for i in range(1, len(vector)):
        value = abs(vector[i])
        if value > n:
            n = value
    return n


def matrix_norm(matrix):
    n = 0
    for i in range(len(matrix)):
        s = 0
        for j in range(len(matrix[i])):
            s += abs(matrix[i][j])
        if s > n:
            n = s
    return n


def upper_triangle_norm(matrix):
    n = 0
    for i in range(len(matrix)):
        s = 0
        for j in range(i + 1, len(matrix[i])):
            s += abs(matrix[i][j])
        if s > n:
            n = s
    return n


def vector_subtraction(vector_a, vector_b):
    vector = list()
    for i in range(len(vector_a)):
        vector.append(vector_a[i] - vector_b[i])
    return vector


def find_x_simply(a, b, epsilon):
    norm = matrix_norm(a)
    x = list(b)
    iteration = 0
    max_iteration = 1000
    not_find = True
    while iteration < max_iteration:
        iteration += 1
        next_x = list()
        for i in range(len(b)):
            s = 0
            for j in range(len(b)):
                s += a[i][j] * x[j]
            next_x.append(b[i] + s)
        delta = vector_norm(vector_subtraction(next_x, x))
        if norm < 1:
            check = norm / (1.0 - norm) * delta
        else:
            check = delta
        x = next_x
        if check <= epsilon:
            not_find = False
            break
    if not_find:
        raise RuntimeError(
            'Метод простых итераций не сошёлся за '
            + str(max_iteration) + ' итераций'
        )
    return x, iteration


def find_x_by_seidel(a, b, epsilon):
    norm = matrix_norm(a)
    upper_norm = upper_triangle_norm(a)
    x = list(b)
    iteration = 0
    max_iteration = 1000
    not_find = True
    while iteration < max_iteration:
        iteration += 1
        next_x = list()
        for i in range(len(b)):
            next_x.append(x[i])
        for i in range(len(b)):
            s = 0
            for j in range(i):
                s += a[i][j] * next_x[j]
            for j in range(i + 1, len(b)):
                s += a[i][j] * x[j]
            next_x[i] = b[i] + s
        delta = vector_norm(vector_subtraction(next_x, x))
        if norm < 1:
            check = upper_norm / (1.0 - norm) * delta
        else:
            check = delta
        x = next_x
        if check <= epsilon:
            not_find = False
            break
    if not_find:
        raise RuntimeError(
            'Метод Зейделя не сошёлся за '
            + str(max_iteration) + ' итераций'
        )
    return x, iteration


def find_x_by_gauss(a_matrix, b):
    matrix = list()
    for i in range(len(a_matrix)):
        line = list(a_matrix[i])
        line.append(b[i])
        matrix.append(line)
    for i in range(len(a_matrix)):
        if matrix[i][i] == 0:
            row = -1
            for j in range(i + 1, len(a_matrix)):
                if matrix[j][i] != 0:
                    row = j
                    break
            if row == -1:
                raise ValueError(
                    f'В столбце {i + 1} ниже диагонали только нули'
                )
            matrix[i], matrix[row] = matrix[row], matrix[i]
        for j in range(i + 1, len(a_matrix)):
            multiplier = matrix[j][i] / matrix[i][i]
            for k in range(i, len(a_matrix) + 1):
                matrix[j][k] -= multiplier * matrix[i][k]
    x = [0] * len(a_matrix)
    for i in range(len(a_matrix) - 1, -1, -1):
        s = 0
        for j in range(i + 1, len(a_matrix)):
            s += matrix[i][j] * x[j]
        x[i] = (matrix[i][len(a_matrix)] - s) / matrix[i][i]
    return x


def vector_to_str(vector):
    line = list()
    for i in range(len(vector)):
        line.append(f'{vector[i]:.3g}')
    str_vector = ' '.join(line)
    return str_vector


def main():
    try:
        with open('input3.txt', 'r', encoding='utf-8') as file:
            a_matrix, b_vector = lines_to_a_and_b(file.readlines())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input3.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input3.txt: {e}')
    if len(a_matrix) == 0:
        raise ValueError('Отсутствует матрица')
    a, b = find_a_and_b(a_matrix, b_vector)
    with open('output3.txt', 'w', encoding='utf-8') as file:
        file.write('Заданная точность вычислений ε' + '\n')
        file.write(f'{EPSILON}' + '\n')
        file.write('\n')
        x_simply, simply_count = find_x_simply(a, b, EPSILON)
        file.write('Метод простых итераций' + '\n')
        file.write(f'Количество итераций: {simply_count}' + '\n')
        file.write(vector_to_str(x_simply) + '\n')
        file.write('\n')
        x_seidel, seidel_count = find_x_by_seidel(a, b, EPSILON)
        file.write('Метод Зейделя' + '\n')
        file.write(f'Количество итераций: {seidel_count}' + '\n')
        file.write(vector_to_str(x_seidel) + '\n')
        file.write('\n')
        x_gauss = find_x_by_gauss(a_matrix, b_vector)
        file.write('Эталонное решение для проверки результата' + '\n')
        file.write(vector_to_str(x_gauss) + '\n')
        file.write('\n')
        if seidel_count < simply_count:
            file.write(
                'Метод Зейделя сошёлся быстрее на '
                + str(simply_count - seidel_count) + ' итераций' + '\n'
            )
        elif seidel_count > simply_count:
            file.write(
                'Метод простых итераций сошёлся быстрее на '
                + str(seidel_count - simply_count) + ' итераций' + '\n'
            )
        else:
            file.write('Методы равны по скорости сходимости' + '\n')


if __name__ == '__main__':
    main()
