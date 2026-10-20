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


def lu_decomposition(a_matrix):
    u_matrix = list()
    for line in a_matrix:
        u_matrix.append(list(line))
    l_matrix = list()
    for i in range(len(a_matrix)):
        line = list()
        for j in range(len(a_matrix)):
            line.append(0)
        l_matrix.append(line)
    p_matrix = list()
    for i in range(len(a_matrix)):
        line = list()
        for j in range(len(a_matrix)):
            if i == j:
                line.append(1)
            else:
                line.append(0)
        p_matrix.append(line)
    det_is_positive = 1
    for i in range(len(a_matrix)):
        lead_line = i
        max_x = abs(u_matrix[i][i])
        for line in range(i + 1, len(a_matrix)):
            value = abs(u_matrix[line][i])
            if value > max_x:
                max_x = value
                lead_line = line
        if max_x < 1e-7:
            raise ValueError('Матрица вырожденная')
        if lead_line != i:
            u_matrix[i], u_matrix[lead_line] = (
                u_matrix[lead_line], u_matrix[i]
            )
            p_matrix[i], p_matrix[lead_line] = (
                p_matrix[lead_line], p_matrix[i]
            )
            for j in range(i):
                el = l_matrix[i][j]
                l_matrix[i][j] = l_matrix[lead_line][j]
                l_matrix[lead_line][j] = el
            det_is_positive = det_is_positive * -1
        l_matrix[i][i] = 1
        for line in range(i + 1, len(a_matrix)):
            mult = (u_matrix[line][i] / u_matrix[i][i])
            l_matrix[line][i] = mult
            for j in range(i, len(a_matrix)):
                u_matrix[line][j] = (u_matrix[line][j] - mult * u_matrix[i][j])
    return l_matrix, u_matrix, p_matrix, det_is_positive


def find_z_from_lz_equal_b(l_matrix, right_part):
    z = list()
    for i in range(len(right_part)):
        s = 0
        for j in range(i):
            s += l_matrix[i][j] * z[j]
        z.append(right_part[i] - s)
    return z


def find_x_from_ux_equal_z(u_matrix, right_part):
    x = [0] * len(right_part)
    for i in range(len(right_part) - 1, -1, -1):
        s = 0
        for j in range(i + 1, len(right_part)):
            s += u_matrix[i][j] * x[j]
        x[i] = (right_part[i] - s) / u_matrix[i][i]
    return x


def matrix_mult_b(matrix, b):
    v = list()
    for i in range(len(matrix)):
        s = 0
        for j in range(len(b)):
            s += matrix[i][j] * b[j]
        v.append(s)
    return v


def matrices_mult(matrix_a, matrix_b):
    matrix = list()
    for i in range(len(matrix_a)):
        line = list()
        for j in range(len(matrix_b[0])):
            s = 0
            for k in range(len(matrix_b)):
                s += matrix_a[i][k] * matrix_b[k][j]
            line.append(s)
        matrix.append(line)
    return matrix


def transposition(matrix):
    t_matrix = list()
    for i in range(len(matrix[0])):
        line = list()
        for j in range(len(matrix)):
            line.append(matrix[j][i])
        t_matrix.append(line)
    return t_matrix


def matrix_to_str(matrix):
    str_matrix = ''
    for i in range(len(matrix)):
        line = list()
        for j in range(len(matrix[i])):
            line.append(f'{matrix[i][j]:.3g}')
        str_matrix += ' '.join(line) + '\n'
    return str_matrix


def vector_to_str(vector):
    str_vector = ''
    line = list()
    for i in range(len(vector)):
        line.append(f'{vector[i]:.3g}')
    str_vector = ' '.join(line)
    return str_vector


def are_matrices_equal(matrix_a, matrix_b):
    for i in range(len(matrix_a)):
        for j in range(len(matrix_a[0])):
            if abs(matrix_a[i][j] - matrix_b[i][j]) > 1e-7:
                return False
    return True


def find_vector_x(l_matrix, u_matrix, p_matrix, right_part):
    p_matrix_mult_b = matrix_mult_b(p_matrix, right_part)
    z = find_z_from_lz_equal_b(l_matrix, p_matrix_mult_b)
    x = find_x_from_ux_equal_z(u_matrix, z)
    return x


def inverse_matrix(matrix, l_matrix, u_matrix, p_matrix):
    columns = list()
    for j in range(len(matrix)):
        e_j = list()
        for i in range(len(matrix)):
            if i == j:
                e_j.append(1)
            else:
                e_j.append(0)
        columns.append(find_vector_x(l_matrix, u_matrix, p_matrix, e_j))
    return transposition(columns)


def calc_det(matrix, det_is_positive):
    det = det_is_positive
    for i in range(len(matrix)):
        det *= matrix[i][i]
    return det


def main():
    try:
        with open('input1.txt', 'r', encoding='utf-8') as file:
            a_matrix, b = lines_to_a_and_b(file.readlines())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input1.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input1.txt: {e}')
    if len(a_matrix) == 0:
        raise ValueError('Отсутствует матрица')
    l_matrix, u_matrix, p_matrix, det_is_positive = lu_decomposition(a_matrix)
    with open('output1.txt', 'w', encoding='utf-8') as file:
        file.write('Матрица L' + '\n')
        file.write(matrix_to_str(l_matrix) + '\n')
        file.write('Матрица U' + '\n')
        file.write(matrix_to_str(u_matrix) + '\n')
        file.write('Результат умножения L · U' + '\n')
        lu_mult = matrices_mult(l_matrix, u_matrix)
        file.write(matrix_to_str(lu_mult) + '\n')
        file.write('Решение системы Ax = b' + '\n')
        vector_x = find_vector_x(l_matrix, u_matrix, p_matrix, b)
        file.write(vector_to_str(vector_x) + '\n')
        file.write('\n')
        file.write('Обратная матрица A⁻¹' + '\n')
        inverted = inverse_matrix(a_matrix, l_matrix, u_matrix, p_matrix)
        file.write(matrix_to_str(inverted) + '\n')
        det = calc_det(u_matrix, det_is_positive)
        file.write('Определитель матрицы A' + '\n')
        file.write(f'{det:.3g}' + '\n')
        file.write('\n')
        file.write('Проверка A · A⁻¹ = E' + '\n')
        matrices_check = matrices_mult(a_matrix, inverted)
        e_matrix = list()
        for i in range(len(a_matrix)):
            line = list()
            for j in range(len(a_matrix)):
                if i == j:
                    line.append(1)
                else:
                    line.append(0)
            e_matrix.append(line)
        inverse_check = are_matrices_equal(matrices_check, e_matrix)
        if inverse_check:
            file.write('Да' + '\n')
        else:
            file.write('Нет' + '\n')
        file.write('\n')
        file.write('Проверка L · U = P · A' + '\n')
        pa_mult = matrices_mult(p_matrix, a_matrix)
        decomposition_check = are_matrices_equal(lu_mult, pa_mult)
        if decomposition_check:
            file.write('Да' + '\n')
        else:
            file.write('Нет' + '\n')


if __name__ == "__main__":
    main()
