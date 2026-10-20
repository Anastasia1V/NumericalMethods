import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


LEFT_BORDER = 0.8
RIGHT_BORDER = 0.9
LAMBDA = 0.3
MAX_ITERATIONS = 1000
GRAPH_LEFT_X = 0
GRAPH_RIGHT_X = 1
GRAPH_LEFT_Y = -1
GRAPH_RIGHT_Y = 1
GRAPH_STEP = 0.001
GRAPH_FILE = '1.png'


def function(x):
    return x ** 3 + x ** 2 - x - 0.5


def first_derivative(x):
    return 3 * x ** 2 + 2 * x - 1


def equivalent_function(x):
    return x - LAMBDA * function(x)


def find_q():
    borders_values = [equivalent_function(LEFT_BORDER),
                      equivalent_function(RIGHT_BORDER)]
    if min(borders_values) < LEFT_BORDER or max(borders_values) > RIGHT_BORDER:
        raise ValueError('Не выполняется условие 1) теоремы 2.3')
    q = max(
        abs(1 - LAMBDA * first_derivative(LEFT_BORDER)),
        abs(1 - LAMBDA * first_derivative(RIGHT_BORDER))
    )
    if q >= 1:
        raise ValueError('Не выполняется условие 2) теоремы 2.3')
    return q


def simple_iter_method(epsilon, q):
    x = (LEFT_BORDER + RIGHT_BORDER) / 2
    list_x = list()
    list_x.append(x)
    for i in range(MAX_ITERATIONS):
        next_x = equivalent_function(x)
        list_x.append(next_x)
        error = q / (1 - q) * abs(next_x - x)
        x = next_x
        if error <= epsilon:
            return x, list_x
    raise RuntimeError('Метод простой итерации не сошёлся')


def newton_method(epsilon):
    x = RIGHT_BORDER
    list_x = list()
    list_x.append(x)
    for i in range(MAX_ITERATIONS):
        next_x = x - function(x) / first_derivative(x)
        list_x.append(next_x)
        if abs(next_x - x) <= epsilon:
            return next_x, list_x
        x = next_x
    raise RuntimeError('Метод Ньютона не сошёлся')


def make_graph():
    points_count = round((GRAPH_RIGHT_X - GRAPH_LEFT_X) / GRAPH_STEP)
    list_x = list()
    for i in range(points_count + 1):
        list_x.append(GRAPH_LEFT_X + i * GRAPH_STEP)
    list_y = list()
    for x in list_x:
        list_y.append(function(x))
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_title('График функции')
    axes.plot(list_x, list_y, label='y = x**3 + x**2 - x - 0.5')
    axes.set_xlabel('x')
    axes.set_ylabel('y')
    axes.set_xlim(GRAPH_LEFT_X, GRAPH_RIGHT_X)
    axes.set_ylim(GRAPH_LEFT_Y, GRAPH_RIGHT_Y)
    axes.xaxis.set_major_locator(MultipleLocator(0.05))
    axes.yaxis.set_major_locator(MultipleLocator(0.05))
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(GRAPH_FILE)
    plt.close(figure)


def main():
    try:
        with open('input1.txt', 'r', encoding='utf-8') as file:
            epsilon = float(file.read().strip())
    except FileNotFoundError:
        raise FileNotFoundError('Отсутствует input1.txt')
    except Exception as e:
        raise ValueError(f'Некорректный input1.txt: {e}')
    if epsilon <= 0:
        raise ValueError('Точность должна быть положительной')
    q = find_q()
    make_graph()
    simple_x, simple_list_x = simple_iter_method(epsilon, q)
    newton_x, newton_list_x = newton_method(epsilon)
    with open('output1.txt', 'w', encoding='utf-8') as file:
        file.write('x^3 + x^2 - x - 0.5 = 0' + '\n')
        file.write(f'Точность вычислений = {epsilon:.3g}' + '\n\n')
        file.write('Графически определено' + '\n')
        file.write(f'x на [{LEFT_BORDER}, {RIGHT_BORDER}]' + '\n')
        file.write('\n' + 'Метод простой итерации' + '\n')
        file.write('Начальное приближение x0 = '
                   f'{((LEFT_BORDER + RIGHT_BORDER) / 2):.3g}' + '\n')
        file.write('i - номер итерации' + '\n')
        file.write('x^(k) - текущее приближение корня' + '\n')
        file.write('q/(1-q)*|x^(k+1)-x^(k)| - оценка погрешности' + '\n')
        file.write('i  x^(k)  q/(1-q)*|x^(k+1)-x^(k)|' + '\n')
        for i in range(len(simple_list_x) - 1):
            x_k = simple_list_x[i]
            x_k1 = simple_list_x[i + 1]
            error = q / (1 - q) * abs(x_k1 - x_k)
            file.write(f'{i}  {x_k:.3g}  {error:.3g}' + '\n')
        file.write(f'Корень методом простой итерации: {simple_x:.3g}' + '\n')
        file.write('\n' + 'Метод Ньютона' + '\n')
        file.write(f'Начальное приближение x0 = {RIGHT_BORDER:.3g}' + '\n')
        file.write('i - номер итерации' + '\n')
        file.write('x^(k) - текущее приближение корня' + '\n')
        file.write('|x^(k+1)-x^(k)| - оценка погрешности' + '\n')
        file.write('i  x^(k)  |x^(k+1)-x^(k)|' + '\n')
        for i in range(len(newton_list_x) - 1):
            x_k = newton_list_x[i]
            x_k1 = newton_list_x[i + 1]
            delta = abs(x_k1 - x_k)
            file.write(f'{i}  {x_k:.3g}  {delta:.3g}' + '\n')
        file.write(f'Корень методом Ньютона: {newton_x:.3g}')


if __name__ == '__main__':
    main()
