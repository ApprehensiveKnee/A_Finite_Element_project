#!/usr/bin/env python

import numpy as np

def check_csv_matrices_equal(file1, file2):
    # Load matrices from csv files
    matrix1 = np.loadtxt(file1, delimiter=',')
    matrix2 = np.loadtxt(file2, delimiter=',')

    # Check if matrices are equal
    if np.array_equal(matrix1, matrix2):
        print("Matrices are equal")
    else:
        print("Matrices are not equal")
        # Show indices of different entries
        rows, cols = np.where(matrix1 != matrix2)
        print("Different entries:")
        for row, col in zip(rows, cols):
            print("({0}, {1})".format(row, col))
    return

# Example usage
check_csv_matrices_equal("../build/system_mat.csv","../mesh/matlab_mesh_generator/Playground/A.csv")
