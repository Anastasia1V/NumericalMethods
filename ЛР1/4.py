import math


EPSILON = 0.01


def lines_to_matrix(lines):
    matrix = list()
    for line in lines:
        numbers = line.split()
        line_numbers = list()
        for number in numbers:
            line_numbers.append(float(number))
        matrix.append(line_numbers)
    return matrix


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


def non_diagonal_sqrt(matrix):
    s = 0
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            s += matrix[i][j] * matrix[i][j]
    return math.sqrt(s)


def rotation_method(matrix, epsilon):
    curr_matrix = list()
    for line in matrix:
        curr_matrix.append(list(line))
    vectors_matrix = list()
    for i in range(len(matrix)):
        line = list()
        for j in range(len(matrix)):
            if i == j:
                line.append(1)
            else:
                line.append(0)
        vectors_matrix.append(line)
    iterations = 0
    while non_diagonal_sqrt(curr_matrix) > epsilon:
        row = 0
        column = 1
        max_number = abs(curr_matrix[0][1])
        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                number = abs(curr_matrix[i][j])
                if number > max_number:
                    max_number = number
                    row = i
                    column = j
        if abs(curr_matrix[row][row] - curr_matrix[column][column]) < 1e-7:
            angle = math.pi / 4
        else:
            angle = 0.5 * math.atan(
                2 * curr_matrix[row][column]
                / (curr_matrix[row][row] - curr_matrix[column][column])
            )
        cosa = math.cos(angle)
        sina = math.sin(angle)
        rotation = list()
        for i in range(len(matrix)):
            line = list()
            for j in range(len(matrix)):
                if i == j:
                    line.append(1)
                else:
                    line.append(0)
            rotation.append(line)
        rotation[row][row] = cosa
        rotation[column][column] = cosa
        rotation[row][column] = -sina
        rotation[column][row] = sina
        curr_matrix = matrices_mult(
            matrices_mult(transposition(rotation), curr_matrix), rotation,
        )
        vectors_matrix = matrices_mult(vectors_matrix, rotation)
        iterations += 1
        if iterations > 1000:
            raise RuntimeError('Метод вращений не сошёлся')
    lambdas = list()
    for i in range(len(matrix)):
        lambdas.append(curr_matrix[i][i])
    return lambdas, vectors_matrix


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


def main():
    try:
        with open('input4.txt', 'r', encoding='utf-8') as file:
            matrix = lines_to_matrix(file.readlines())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input4.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input4.txt: {e}')
    if len(matrix) == 0:
        raise ValueError('Отсутствует матрица')
    lambdas, vectors_matrix = rotation_method(matrix, EPSILON)
    diagonal = list()
    for i in range(len(matrix)):
        line = list()
        for j in range(len(matrix)):
            if i == j:
                line.append(lambdas[i])
            else:
                line.append(0)
        diagonal.append(line)
    with open('output4.txt', 'w', encoding='utf-8') as file:
        file.write('Заданная точность вычислений ε' + '\n')
        file.write(f'{EPSILON}' + '\n')
        file.write('\n')
        file.write('Найденные собственные значения' + '\n')
        file.write(vector_to_str(lambdas) + '\n')
        file.write('\n')
        file.write('Найденные собственные векторы' + '\n')
        for j in range(len(lambdas)):
            column = list()
            for i in range(len(matrix)):
                column.append(vectors_matrix[i][j])
            file.write(f'x{j + 1} = ' + vector_to_str(column) + '\n')
        file.write('\n')
        file.write('Матрица собственных векторов' + '\n')
        file.write(matrix_to_str(vectors_matrix) + '\n')
        left = matrices_mult(matrix, vectors_matrix)
        right = matrices_mult(vectors_matrix, diagonal)
        is_correct = True
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if abs(left[i][j] - right[i][j]) > EPSILON:
                    is_correct = False
        file.write('Проверка A · V = V · Λ' + '\n')
        if is_correct:
            file.write('Да' + '\n')
        else:
            file.write('Нет' + '\n')


if __name__ == '__main__':
    main()
