import math
import numpy


EPSILON = 0.01


def lines_to_matrix(lines):
    matrix = list()
    for line in lines:
        line_numbers = list()
        for el in line.split():
            line_numbers.append(float(el))
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


def find_qr(matrix):
    r_matrix = list()
    for line in matrix:
        r_matrix.append(list(line))
    q_matrix = list()
    size = len(matrix)
    for i in range(size):
        line = list()
        for j in range(size):
            if i == j:
                line.append(1)
            else:
                line.append(0)
        q_matrix.append(line)
    for column in range(size - 1):
        vector = list()
        for row in range(column, size):
            vector.append(r_matrix[row][column])
        norm = 0
        for el in vector:
            norm += el * el
        norm = math.sqrt(norm)
        if norm < 1e-7:
            continue
        vector[0] += math.copysign(norm, vector[0])
        v_norm = 0
        for el in vector:
            v_norm += el * el
        v_norm = math.sqrt(v_norm)
        for i in range(len(vector)):
            vector[i] = vector[i] / v_norm
        h_matrix = list()
        for i in range(size):
            line = list()
            for j in range(size):
                if i == j:
                    line.append(1)
                else:
                    line.append(0)
            h_matrix.append(line)
        for i in range(column, size):
            for j in range(column, size):
                h_matrix[i][j] -= (
                    2 * vector[i - column] * vector[j - column]
                )
        r_matrix = matrices_mult(h_matrix, r_matrix)
        q_matrix = matrices_mult(q_matrix, h_matrix)
    return q_matrix, r_matrix


def split_into_blocks(matrix, epsilon):
    blocks = list()
    i = 0
    size = len(matrix)
    while i < size:
        if i == size - 1:
            blocks.append((i, 1))
            i += 1
        else:
            subdiag_norm = math.sqrt(
                sum(matrix[k][i] ** 2 for k in range(i + 1, size))
            )
            if subdiag_norm <= epsilon:
                blocks.append((i, 1))
                i += 1
            else:
                blocks.append((i, 2))
                i += 2
    return blocks


def find_lambdas_for_block(matrix, i):
    a11 = matrix[i][i]
    a12 = matrix[i][i + 1]
    a21 = matrix[i + 1][i]
    a22 = matrix[i + 1][i + 1]
    diag_sum = a11 + a22
    det = a11 * a22 - a12 * a21
    discriminant = diag_sum ** 2 - 4 * det
    if discriminant >= 0:
        lambda1 = (diag_sum + math.sqrt(discriminant)) / 2
        lambda2 = (diag_sum - math.sqrt(discriminant)) / 2
        return [lambda1, lambda2]
    lambda1 = complex(diag_sum / 2, math.sqrt(-discriminant) / 2)
    lambda2 = complex(diag_sum / 2, -1 * math.sqrt(-discriminant) / 2)
    return [lambda1, lambda2]


def find_lambdas_by_qr(matrix, epsilon):
    current_matrix = [list(row) for row in matrix]
    previous_lambdas = dict()
    blocks = list()
    for i in range(1000):
        q, r = find_qr(current_matrix)
        current_matrix = matrices_mult(r, q)
        blocks = split_into_blocks(current_matrix, epsilon)
        correct = True
        current_lambdas = dict()
        for start, size in blocks:
            if size == 2:
                lambdas = find_lambdas_for_block(current_matrix, start)
                current_lambdas[start] = lambdas
                if start not in previous_lambdas:
                    correct = False
                else:
                    for curr, prev in zip(lambdas, previous_lambdas[start]):
                        if abs(curr - prev) > epsilon:
                            correct = False
                            break
        previous_lambdas = current_lambdas
        if correct:
            break
    result = list()
    for start, size in blocks:
        if size == 1:
            result.append(current_matrix[start][start])
        else:
            result += find_lambdas_for_block(current_matrix, start)
    return result


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
        with open('input5.txt', 'r', encoding='utf-8') as file:
            matrix = lines_to_matrix(file.readlines())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input5.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input5.txt: {e}')
    if len(matrix) == 0:
        raise ValueError('Отсутствует матрица')
    q_matrix, r_matrix = find_qr(matrix)
    lambdas = find_lambdas_by_qr(matrix, EPSILON)
    example = [[0, -1, 0], [1, 0, 0], [0, 0, 3]]
    lambdas_by_qr = find_lambdas_by_qr(example, EPSILON)
    lambdas_by_numpy = numpy.linalg.eigvals(numpy.array(example))
    with open('output5.txt', 'w', encoding='utf-8') as file:
        file.write('Исходная матрица A' + '\n')
        file.write(matrix_to_str(matrix) + '\n')
        file.write('Матрица Q' + '\n')
        file.write(matrix_to_str(q_matrix) + '\n')
        file.write('Матрица R' + '\n')
        file.write(matrix_to_str(r_matrix) + '\n')
        file.write('Результат умножения Q · R' + '\n')
        qr_mult = matrices_mult(q_matrix, r_matrix)
        file.write(matrix_to_str(qr_mult) + '\n')
        file.write('Найденные собственные значения' + '\n')
        file.write(vector_to_str(lambdas) + '\n\n')
        file.write('Собственный пример матрицы с комплексно-сопряжёнными '
                   'собственными значениями' + '\n')
        file.write(matrix_to_str(example) + '\n')
        file.write('Собственные значения через QR-разложение' + '\n')
        file.write(vector_to_str(lambdas_by_qr) + '\n\n')
        file.write('Собственные значения через NumPy' + '\n')
        file.write(vector_to_str(lambdas_by_numpy))


if __name__ == '__main__':
    main()
