import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


LEFT_BORDER_1 = 0.84
RIGHT_BORDER_1 = 1.0
LEFT_BORDER_2 = 1.6
RIGHT_BORDER_2 = 1.7
X_BY_GRAPH = 0.92
Y_BY_GRAPH = 1.65
MAX_ITERATIONS = 1000
GRAPH_LEFT_X1 = 0.80
GRAPH_RIGHT_X1 = 1.00
GRAPH_LEFT_X2 = 1.6
GRAPH_RIGHT_X2 = 1.7
GRAPH_X1_STEP = 0.001
GRAPH_X2_STEP = 0.001
GRAPH_FILE = '2.png'


def function(vector):
    x1 = vector[0]
    x2 = vector[1]
    return x1 - math.cos(x2) - 1, x2 - math.log(x1 + 1) - 1


def equivalent_function(vector):
    x1 = vector[0]
    x2 = vector[1]
    return 1 + math.cos(x2), 1 + math.log(x1 + 1)


def vector_norm(vector):
    n = abs(vector[0])
    if abs(vector[1]) > n:
        n = abs(vector[1])
    return n


def subtract_vectors(vector_a, vector_b):
    return vector_a[0] - vector_b[0], vector_a[1] - vector_b[1]


def find_q():
    phi1_min = 1 + math.cos(RIGHT_BORDER_2)
    phi1_max = 1 + math.cos(LEFT_BORDER_2)
    phi2_min = 1 + math.log(LEFT_BORDER_1 + 1)
    phi2_max = 1 + math.log(RIGHT_BORDER_1 + 1)
    if phi1_min < LEFT_BORDER_1 or phi1_max > RIGHT_BORDER_1:
        raise ValueError('Не выполняется условие 1) теоремы 2.4')
    if phi2_min < LEFT_BORDER_2 or phi2_max > RIGHT_BORDER_2:
        raise ValueError('Не выполняется условие 1) теоремы 2.4')
    derivative_12 = abs(math.sin(LEFT_BORDER_2))
    if abs(math.sin(RIGHT_BORDER_2)) > derivative_12:
        derivative_12 = abs(math.sin(RIGHT_BORDER_2))
    derivative_21 = 1 / (LEFT_BORDER_1 + 1)
    q = max(derivative_12, derivative_21)
    if q >= 1:
        raise ValueError('Не выполняется условие 2) теоремы 2.4')
    return q


def simple_iter_method(epsilon, q):
    x = (X_BY_GRAPH, Y_BY_GRAPH)
    list_x = list()
    list_x.append(x)
    for i in range(MAX_ITERATIONS):
        next_x = equivalent_function(x)
        list_x.append(next_x)
        delta = vector_norm(subtract_vectors(next_x, x))
        error = q / (1 - q) * delta
        x = next_x
        if error <= epsilon:
            return x, list_x
    raise RuntimeError('Метод простой итерации не сошёлся')


def newton_method(epsilon, v):
    x = v
    list_x = list()
    list_x.append(x)
    for i in range(MAX_ITERATIONS):
        x1 = x[0]
        x2 = x[1]
        f1, f2 = function(x)
        det = 1 + math.sin(x2) / (x1 + 1)
        delta_x1 = (-f1 + math.sin(x2) * f2) / det
        delta_x2 = (-f2 - f1 / (x1 + 1)) / det
        next_x = (x1 + delta_x1, x2 + delta_x2)
        list_x.append(next_x)
        x = next_x
        if max(abs(delta_x1), abs(delta_x2)) <= epsilon:
            return x, list_x
    raise RuntimeError('Метод Ньютона не сошёлся')


def make_graph():
    x1_count = round((GRAPH_RIGHT_X1 - GRAPH_LEFT_X1) / GRAPH_X1_STEP)
    x2_count = round((GRAPH_RIGHT_X2 - GRAPH_LEFT_X2) / GRAPH_X2_STEP)
    list_x1 = list()
    for i in range(x1_count + 1):
        list_x1.append(GRAPH_LEFT_X1 + i * GRAPH_X1_STEP)
    list_x2 = list()
    for i in range(x2_count + 1):
        list_x2.append(GRAPH_LEFT_X2 + i * GRAPH_X2_STEP)
    list_f1 = list()
    list_f2 = list()
    for x2 in list_x2:
        row_f1 = list()
        row_f2 = list()
        for x1 in list_x1:
            f1, f2 = function((x1, x2))
            row_f1.append(f1)
            row_f2.append(f2)
        list_f1.append(row_f1)
        list_f2.append(row_f2)
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_title('Графики функций')
    axes.contour(list_x1, list_x2, list_f1, levels=[0],
                 colors=['tab:blue'])
    axes.contour(list_x1, list_x2, list_f2, levels=[0],
                 colors=['tab:green'])
    axes.plot([], [], color='tab:blue', label='x1 - cos(x2) - 1 = 0')
    axes.plot([], [], color='tab:green',
              label='x2 - ln(x1 + 1) - 1 = 0')
    axes.set_xlabel('x1')
    axes.set_ylabel('x2')
    axes.set_xlim(GRAPH_LEFT_X1, GRAPH_RIGHT_X1)
    axes.set_ylim(GRAPH_LEFT_X2, GRAPH_RIGHT_X2)
    axes.xaxis.set_major_locator(MultipleLocator(0.01))
    axes.yaxis.set_major_locator(MultipleLocator(0.01))
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(GRAPH_FILE)
    plt.close(figure)


def main():
    try:
        with open('input2.txt', 'r', encoding='utf-8') as file:
            epsilon = float(file.read().strip())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input2.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input2.txt: {e}')
    if epsilon <= 0:
        raise ValueError('Точность должна быть положительной')
    q = find_q()
    make_graph()
    simple_x, simple_list_x = simple_iter_method(epsilon, q)
    newton_x, newton_list_x = newton_method(epsilon, (X_BY_GRAPH, Y_BY_GRAPH))
    with open('output2.txt', 'w', encoding='utf-8') as file:
        file.write('x1 - cos(x2) - 1 = 0' + '\n')
        file.write('x2 - ln(x1 + 1) - 1 = 0' + '\n')
        file.write(f'Точность вычислений = {epsilon:.3g}' + '\n\n')
        file.write('Графически определено' + '\n')
        file.write(f'x1 на [{LEFT_BORDER_1}, {RIGHT_BORDER_1}]' + '\n')
        file.write(f'x2 на [{LEFT_BORDER_2}, {RIGHT_BORDER_2}]' + '\n')
        file.write('Начальное приближение x0 = '
                   f'({X_BY_GRAPH}, {Y_BY_GRAPH})' + '\n')
        file.write('\n' + 'Метод простой итерации' + '\n')
        file.write('i - номер итерации' + '\n')
        file.write('x1^(k) - текущее приближение x1' + '\n')
        file.write('x2^(k) - текущее приближение x2' + '\n')
        file.write('q/(1-q)*||x^(k+1)-x^(k)|| - оценка погрешности' + '\n')
        file.write('i  x1^(k)  x2^(k)  '
                   'q/(1-q)*||x^(k+1)-x^(k)||' + '\n')
        for i in range(len(simple_list_x) - 1):
            v_k = simple_list_x[i]
            v_k1 = simple_list_x[i + 1]
            delta = vector_norm(subtract_vectors(v_k1, v_k))
            error = q / (1 - q) * delta
            file.write(f'{i}  {v_k[0]:.3g}  {v_k[1]:.3g}  '
                       f'{error:.3g}' + '\n')
        file.write(f'Корень методом простой итерации: '
                   f'({simple_x[0]:.3g}, {simple_x[1]:.3g})' + '\n')
        file.write('\n' + 'Метод Ньютона' + '\n')
        file.write('Начальное приближение x0 = '
                   f'({X_BY_GRAPH}, {Y_BY_GRAPH})' + '\n')
        file.write('i - номер итерации' + '\n')
        file.write('x1^(k) - текущее приближение x1' + '\n')
        file.write('x2^(k) - текущее приближение x2' + '\n')
        file.write('||x^(k+1)-x^(k)|| - критерий остановки' + '\n')
        file.write('i  x1^(k)  x2^(k)  ||x^(k+1)-x^(k)||' + '\n')
        for i in range(len(newton_list_x) - 1):
            v_k = newton_list_x[i]
            v_k1 = newton_list_x[i + 1]
            delta = vector_norm(subtract_vectors(v_k1, v_k))
            file.write(f'{i}  {v_k[0]:.3g}  {v_k[1]:.3g}  '
                       f'{delta:.3g}' + '\n')
        file.write(f'Корень методом Ньютона: '
                   f'({newton_x[0]:.3g}, {newton_x[1]:.3g})')


if __name__ == '__main__':
    main()
