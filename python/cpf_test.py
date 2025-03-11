import numpy as np
import timeit
from mHT.CPF import (
    cpf_fast,
    cpf_accurate,
    cpf_fast_vector,
    cpf_accurate_vector,
    cpf_fast_vector_based_on_numpy,
)

# Example parameters for CPF functions
x = 1.0  # Dimensionless
y = 1.0  # Dimensionless
x_ = np.array([x])
y_ = np.array([y])
number_of_calls = 10000

# Function to measure execution time
def measure_time(func, *args, number=1000):
    return timeit.timeit(lambda: func(*args), number=number)

# Compute and display results
print("=" * 50)
print("CPF FUNCTION OUTPUTS")
print("=" * 50)

# Scalar functions
cpf_acc_res = cpf_accurate(x, y)
cpf_fast_res = cpf_fast(x, y)

print(f"cpf_accurate({x}, {y}): {cpf_acc_res}")
print(f"cpf_fast({x}, {y}): {cpf_fast_res}")

# Vectorized functions
cpf_acc_vec_res = cpf_accurate_vector(x_, y_)
cpf_fast_vec_res = cpf_fast_vector(x_, y_)
cpf_fast_vec_np_res = cpf_fast_vector_based_on_numpy(x_, y_)

print(f"cpf_accurate_vector({x_}, {y_}): {cpf_acc_vec_res}")
print(f"cpf_fast_vector({x_}, {y_}): {cpf_fast_vec_res}")
print(f"cpf_fast_vector_based_on_numpy({x_}, {y_}): {cpf_fast_vec_np_res}")

# Measure execution times
print("\n" + "=" * 50)
print(f"CPF FUNCTION EXECUTION TIMES ({number_of_calls} runs)")
print("=" * 50)

time_acc = measure_time(cpf_accurate, x, y, number=number_of_calls)
time_fast = measure_time(cpf_fast, x, y, number=number_of_calls)
time_acc_vec = measure_time(cpf_accurate_vector, x_, y_, number=number_of_calls)
time_fast_vec = measure_time(cpf_fast_vector, x_, y_, number=number_of_calls)
time_fast_vec_np = measure_time(cpf_fast_vector_based_on_numpy, x_, y_, number=number_of_calls)

print(f"Execution time - cpf_accurate: {time_acc:.6f} s")
print(f"Execution time - cpf_fast: {time_fast:.6f} s")
print(f"Execution time - cpf_accurate_vector: {time_acc_vec:.6f} s")
print(f"Execution time - cpf_fast_vector: {time_fast_vec:.6f} s")
print(f"Execution time - cpf_fast_vector_based_on_numpy: {time_fast_vec_np:.6f} s")

print("\n" + "=" * 50)
print("PERFORMANCE COMPARISON")
print("=" * 50)
print(f"cpf_fast vs cpf_accurate speedup: {time_acc / time_fast:.2f}x faster")
print(f"cpf_fast_vector vs cpf_accurate_vector speedup: {time_acc_vec / time_fast_vec:.2f}x faster")
print(f"cpf_fast_vector_based_on_numpy vs cpf_fast_vector speedup: {time_fast_vec / time_fast_vec_np:.2f}x faster")

print("=" * 50)
print("END")
print("=" * 50)
