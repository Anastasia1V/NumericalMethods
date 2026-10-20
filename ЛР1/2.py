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


def lines_to_abcd(lines):
    a = list()
    b = list()
    c = list()
    d = list()
    i = 0
    for line in lines:
        line = line.replace(' ', '').replace('−', '-')
        left_part, right_part = line_to_left_and_right(line)
        i += 1
        if i - 1 in left_part:
            a.append(left_part[i - 1])
        else:
            a.append(0)
        if i in left_part:
            b.append(left_part[i])
        else:
            b.append(0)
        if i + 1 in left_part:
            c.append(left_part[i + 1])
        else:
            c.append(0)
        d.append(right_part)
    return a, b, c, d


def find_p_q_x(a, b, c, d):
    p = [0] * len(b)
    q = [0] * len(b)
    p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]
    for i in range(1, len(b) - 1):
        divisor = b[i] + a[i] * p[i - 1]
        p[i] = -c[i] / divisor
        q[i] = (d[i] - a[i] * q[i - 1]) / divisor
    p[len(b) - 1] = 0
    divisor = b[len(b) - 1] + a[len(b) - 1] * p[len(b) - 2]
    q[len(b) - 1] = (d[len(b) - 1] - a[len(b) - 1] * q[len(b) - 2]) / divisor
    x = [0] * len(b)
    x[len(b) - 1] = q[len(b) - 1]
    for i in range(len(b) - 2, -1, -1):
        x[i] = p[i] * x[i + 1] + q[i]
    return p, q, x


def main():
    try:
        with open('input2.txt', 'r', encoding='utf-8') as file:
            a, b, c, d = lines_to_abcd(file.readlines())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input2.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input2.txt: {e}')
    if len(b) == 0:
        raise ValueError('Отсутствует система')
    output = list()
    p, q, x = find_p_q_x(a, b, c, d)
    output.append('Все прогоночные коэффициенты')
    for i in range(len(p)):
        output.append(f'i = {i + 1}, P = {p[i]:.3g}, Q = {q[i]:.3g}')
    output.append('')
    output.append('Итоговое решение системы уравнений')
    line = list()
    for i in range(len(x)):
        line.append(f'{x[i]:.3g}')
    output.append(' '.join(line))
    with open('output2.txt', 'w', encoding='utf-8') as file:
        file.write('\n'.join(output))


if __name__ == '__main__':
    main()
