#!/usr/bin/env python

import numpy as np

# def check_csv_vectors_equal(file1, file2):
#     # Load matrices from csv files
#     matrix1 = np.loadtxt(file1, delimiter=',')
#     matrix2 = np.loadtxt(file2, delimiter=',')

#     # Check if matrices are equal
#     if np.array_equal(matrix1, matrix2):
#         print("Matrices are equal")
#     else:
#         print("Matrices are not equal")
#         # Show indices of different entries
#         rows = np.where(matrix1 != matrix2)
#         cols = np.zeros_like(rows)
#         print("Different entries:")
#         for row, col in zip(rows, cols):
#             print("({0}, {1})".format(row, col))
#     return

def check_csv_vectors_equal(file1, file2, rtol=1e-05, atol=1e-08):
    vector1 = np.loadtxt(file1, delimiter=',')
    vector2 = np.loadtxt(file2, delimiter=',')

    if np.isclose(vector1, vector2, rtol=rtol, atol=atol).all():
        print("Vectors are equal")
        return True
    else:
        diff_indices = np.where(~np.isclose(vector1, vector2, rtol=rtol, atol=atol))[0]
        diff_values1 = vector1[diff_indices]
        diff_values2 = vector2[diff_indices]
        print("Different entries:")
        for i in range(len(diff_indices)):
            print(f"({diff_indices[i]}): {diff_values1[i]} != {diff_values2[i]}")
        return False

# Example usage
check_csv_vectors_equal("../build/system_rhs_serial.csv","../build/system_rhs_parallel.csv")
