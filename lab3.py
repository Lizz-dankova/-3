import timeit
import matplotlib.pyplot as plt

class GaloisField:
    def __init__(self, m, p):
        self.m = m
        self.p = p
        self.irreducible_poly = self.get_irreducible_poly()

    def get_irreducible_poly(self):
        return 0b10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001

    def add(self, a, b):
        # Додавання за модулем p
        result = a ^ b
        return result % self.p

    def multiply(self, a, b):
        # Множення за модулем незвідного многочлена
        product = 0
        while a and b:
            if b & 1:
                product ^= a
            a <<= 1
            if a & (1 << self.m):
                a ^= self.irreducible_poly
            b >>= 1
        return product % self.p

    def multiply_with_reduction(self, a, b):
        # Множення зі зведенням за модулем незвідного многочлена
        product = self.multiply(a, b)
        for _ in range(self.m, 2 * self.m - 1):
            if product & (1 << self.m):
                product ^= self.irreducible_poly
            product <<= 1
        for _ in range(self.m - 1):
            if product & (1 << self.m):
                product ^= self.irreducible_poly
            product <<= 1
        return product % self.p

    def square(self, a):
        # Піднесення до квадрату за модулем незвідного многочлена
        a_sq = 0
        for i in range(self.m):
            if a & (1 << i):
                a_sq ^= (1 << (2 * i))
        a_sq ^= a
        for _ in range(self.m - 1):
            if a_sq & (1 << self.m):
                a_sq ^= self.irreducible_poly
            a_sq <<= 1
        return a_sq % self.p

    def inverse_with_extended_euclidean(self, a):
        # Обчислення оберненого елемента за допомогою розширеного алгоритму Евкліда
        t, new_t = 0, 1
        r, new_r = self.p, a

        while new_r != 0:
            quotient = r // new_r
            t, new_t = new_t, t - quotient * new_t
            r, new_r = new_r, r - quotient * new_r

        if r > 1:
            raise ValueError(f"{a} is not invertible")
        if t < 0:
            t += self.p

        return t % self.p

    def zero_element(self):
        # Константа 0 – нейтральний елемент по операції «+»
        return 0

    def one_element(self):
        # Константа 1 – нейтральний елемент по операції «*»
        return 1

    def trace(self, a):
        # Слід елемента поля
        result = a
        for _ in range(1, self.m):
            a = self.multiply(a, a)
            result ^= a
        return result % self.p

    def power(self, a, n):
        # Піднесення до степеня
        result = 1
        while n > 0:
            if n & 1:
                result = self.multiply(result, a)
            a = self.square(a)
            n >>= 1
        return result % self.p

    def to_binary_string(self, a):
        # Конвертація елемента поля в m-бітний рядок
        return bin(a)[2:].zfill(self.m)

    def from_binary_string(self, s):
        # Конвертація m-бітного рядка в елемент поля
        return int(s, 2)

# Приклад використання:
m = 163
p = 0b1101100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
field = GaloisField(m, p)

# Приклади операцій:
a = 0b1101
b = 0b1010
c = 0b0110
d = 0b1100

# Додавання
sum_result = field.add(a, b)
print(f"Додавання: {field.to_binary_string(a)} + {field.to_binary_string(b)} = {field.to_binary_string(sum_result)}")

# Множення
mul_result = field.multiply(a, b)
print(f"Множення: {field.to_binary_string(a)} * {field.to_binary_string(b)} = {field.to_binary_string(mul_result)}")

# Множення зі зведенням за модулем незвідного многочлена
mul_with_reduction_result = field.multiply_with_reduction(a, b)
print(f"Множення зі зведенням: {field.to_binary_string(a)} * {field.to_binary_string(b)} mod p = {field.to_binary_string(mul_with_reduction_result)}")

# Піднесення до квадрату
square_result = field.square(a)
print(f"Піднесення до квадрату: {field.to_binary_string(a)}^2 mod p = {field.to_binary_string(square_result)}")

# Обернений елемент за допомогою extended Euclidean algorithm
inverse_result = field.inverse_with_extended_euclidean(a)
print(f"Inverse({field.to_binary_string(a)}) = {field.to_binary_string(inverse_result)}")

# Слід елемента поля
trace_result = field.trace(a)
print(f"Слід елемента поля {field.to_binary_string(a)}: {field.to_binary_string(trace_result)}")

# Піднесення до степеня
n = 5
power_result = field.power(a, n)
print(f"Піднесення до степеня {n}: {field.to_binary_string(a)}^{n} = {field.to_binary_string(power_result)}")

# Перевірка тотожностей

# Тотожність (a + b) * c = b * c + c * a
left_side = field.multiply(field.add(a, b), c)
right_side = field.add(field.multiply(b, c), field.multiply(c, a))
print(f"Тотожність 1: ({field.to_binary_string(a)} + {field.to_binary_string(b)}) * {field.to_binary_string(c)} = {field.to_binary_string(left_side)}")
print(f"Перевірка 1: {left_side == right_side}")

# Тотожність d^(2^m-1) = 1 для деякого d
left_side = field.power(d, (1 << m) - 1)
right_side = field.power(d, (1 << m) - 1)
print(f"Тотожність 2: {field.to_binary_string(d)}^(2^m-1) = {field.to_binary_string(left_side)}")
print(f"Перевірка 2: {left_side == right_side}")

# Вимір часу для додавання, множення та піднесення до степеня
num_trials = 10000
add_time = timeit.timeit(lambda: field.add(0b1101, 0b1010), number=num_trials) / num_trials
multiply_time = timeit.timeit(lambda: field.multiply_with_reduction(0b1101, 0b1010), number=num_trials) / num_trials
power_time = timeit.timeit(lambda: field.power(0b1101, 5), number=num_trials) / num_trials

# Вивід результатів та побудова графіка
print(f"Середній час для додавання: {add_time:.6f} секунд")
print(f"Середній час для множення: {multiply_time:.6f} секунд")
print(f"Середній час для піднесення до степеня: {power_time:.6f} секунд")

operations = ['Додавання', 'Множення', 'Піднесення до степеня']
times = [add_time, multiply_time, power_time]

plt.bar(operations, times)
plt.ylabel('Середній час (секунди)')
plt.title('Середній час для операцій в полі Галуа')
plt.show()
