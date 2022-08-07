#include "s21_matrix.h"

// Создание матриц (create_matrix)
int s21_create_matrix(int rows, int columns, matrix_t *result) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double*));
    for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
    }
    return 0;
}

// Очистка матриц (remove_matrix)
void s21_remove_matrix(matrix_t *A) {
    if (A->matrix) {
        for (int j = 0; j < A->rows; j++) {
            free(A->matrix[j]);
        }
        free(A->matrix);
    }
}

// Сравнение матриц (eq_matrix)
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int flag_status = SUCCESS;
    if (A->rows != B->rows || A->columns != B->columns || A->rows <= 0 ||
        A->columns <= 0 || B->rows <= 0 || B->columns <= 0) {
        flag_status = FAILURE;
    } else {
        for (int i = 0; (i < A->rows) && flag_status; i++)
            for (int j = 0; (j < A->columns) && flag_status; j++)
                flag_status = (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07 ? FAILURE : SUCCESS);
    }
    return flag_status;
}

// Сложение (sum_matrix)
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int flag_status = OK;
    if ((A->rows != B->rows) || (A->columns != B->columns)) {
        flag_status = CALCULATION_ERROR;
    } else if ((A->matrix != NULL) && (B->matrix != NULL) && (flag_status == OK) && (A->rows > 0)) {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
            }
        }
    } else {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Вычитание матриц (sub_matrix)
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int flag_status = OK;
    if ((A->rows != B->rows) || (A->columns != B->columns)) {
        flag_status = CALCULATION_ERROR;
    } else if ((A->matrix != NULL) && (B->matrix != NULL) && (flag_status == OK) && (A->rows > 0)) {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
            }
        }
    } else {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Умножение матрицы на число (mult_number)
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int flag_status = (((A->rows) > 0 && (A->columns) > 0) ? OK : CALCULATION_ERROR);
    if ((A->matrix != NULL) && (flag_status == OK)) {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = A->matrix[i][j] * number;
            }
        }
    } else if (flag_status != 2) {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Умножение двух матриц (mult_matrix)
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int flag_status = OK;
    if ((!(A->rows > 0) && !(A->columns > 0)) || (!(B->rows > 0) && !(B->columns > 0))) {
        flag_status = CALCULATION_ERROR;
    }
    if (((A->matrix != NULL) && (flag_status == OK)) && ((B->matrix != NULL) && (flag_status == OK))) {
        s21_create_matrix(A->rows, B->columns, result);
        for (int i = 0; (i < A->rows); i++) {
            for (int j = 0; j < B->columns; j++) {
                for (int k = 0; k < A->columns; k++)
                    result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
            }
        }
    } else if (flag_status != 2) {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Транспонирование матрицы (transpose)
int s21_transpose(matrix_t *A, matrix_t *result) {
    int flag_status = (((A->rows) > 0 && (A->columns) > 0) ? OK : CALCULATION_ERROR);
    if ((A->matrix != NULL) && (flag_status == OK)) {
        s21_create_matrix(A->columns, A->rows, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[j][i] = A->matrix[i][j];
            }
        }
    }  else if (flag_status != 2) {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Минор матрицы и матрица алгебраических дополнений (calc_complements)
int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int flag_status = (A->matrix || (A->columns > 0) || (A->rows > 0) ? OK : CALCULATION_ERROR);
    if ((A->rows == A->columns) && (flag_status == OK)) {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                matrix_t temp_mtr;
                s21_create_matrix(A->rows - 1, A->columns - 1, &temp_mtr);
                int count_i = 0, count_j = 0;
                for (int k = 0; k < A->rows; k++) {
                    if (k == i) continue;
                    for (int l = 0; l < A->columns; l++) {
                        if (l == j) continue;
                        temp_mtr.matrix[count_i][count_j] = A->matrix[k][l];
                        count_j++;
                    }
                    count_j = 0;
                    count_i++;
                }
                double determinant_value;
                s21_determinant(&temp_mtr, &determinant_value);
                s21_remove_matrix(&temp_mtr);
                determinant_value *= pow(-1, 2 + i + j);
                result->matrix[i][j] = determinant_value;
            }
        }
    } else if (flag_status != 2) {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Определитель матрицы
int s21_determinant(matrix_t *A, double *result) {
    int flag_status = (A->columns != A->rows ? CALCULATION_ERROR : OK);
    if (result && (flag_status == OK)) {
        *result = (A->rows == 1 ? A->matrix[0][0] : *result);
        *result = (A->rows == 2 ? A->matrix[0][0] * A->matrix[1][1] -
        A->matrix[0][1] * A->matrix[1][0] : *result);
        if (A->rows != 1 && A->rows != 2) {
            *result = 0;
            for (int i = 0; i < A->columns; i++) {
                matrix_t temp_mtr;
                s21_create_matrix(A->rows - 1, A->columns - 1, &temp_mtr);
                int count_i = 0, count_j = 0;
                for (int k = 0; k < A->rows; k++) {
                    if (k == 0) continue;
                    for (int l = 0; l < A->columns; l++) {
                        if (l == i) continue;
                        temp_mtr.matrix[count_i][count_j] = A->matrix[k][l];
                        count_j++;
                    }
                    count_j = 0;
                    count_i++;
                }
                double tmp_det;
                s21_determinant(&temp_mtr, &tmp_det);
                s21_remove_matrix(&temp_mtr);
                double tmp_res = pow(-1, 2 + i) * A->matrix[0][i] * tmp_det;
                *result += tmp_res;
            }
        }
    } else if (flag_status != 2) {
        flag_status = INCORRECT_MTRX;
    }
    return flag_status;
}

// Обратная матрица (inverse_matrix)
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int flag_status = OK;
    if (!A || !A->matrix || !result || A->columns <= 1 || A->rows <= 1) {
        flag_status = INCORRECT_MTRX;
    } else if (A->rows != A->columns) {
        flag_status = CALCULATION_ERROR;
    } else {
        double determ = 0;
        matrix_t tmp;
        matrix_t tmp1;
        s21_determinant(A, &determ);
        if (fabs(determ) > 1e-7) {
            s21_calc_complements(A, &tmp);
            s21_transpose(&tmp, &tmp1);
            s21_remove_matrix(&tmp);
            s21_mult_number(&tmp1, 1. / fabs(determ), result);
            s21_remove_matrix(&tmp1);
            for (int i = 0; i < result->rows; i++) {
                for (int j = 0; j < result->columns; j++) {
                    result->matrix[i][j] *= -1;
                }
            }
        } else {
            flag_status = CALCULATION_ERROR;
        }
    }
    return flag_status;
}

// void s21_print_one_matrix(matrix_t A) {
//     for (int i = 0; i < A.rows; i++) {
//         for (int j = 0; j < A.columns; j++) {
//             if (j == 0) {
//                 printf("|%.lf ", A.matrix[i][j]);
//             } else if (j == (A.columns) -1) {
//                 printf(" %.lf|", A.matrix[i][j]);
//             } else {
//                 printf(" %.lf ", A.matrix[i][j]);
//             }
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

// void s21_print_two_matrix(matrix_t A, matrix_t B) {
//     printf("\n---------------------------------\n");
//     printf("\n--- A ---\n");
//     s21_print_one_matrix(A);
//     printf("\n--- B ---\n");
//     s21_print_one_matrix(B);
//     printf("\n---------------------------------\n");
// }

// void s21_print_three_matrix(matrix_t A, matrix_t B, matrix_t result) {
//     printf("\n---------------------------------\n");
//     printf("\n--- A ---\n");
//     s21_print_one_matrix(A);
//     printf("\n--- B ---\n");
//     s21_print_one_matrix(B);
//     printf("\n- result -\n");
//     s21_print_one_matrix(result);
//     printf("\n---------------------------------\n");
// }

void s21_gen_matrix(matrix_t *A) {
    double num = 1;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++, num++) {
            A->matrix[i][j] = num;
        }
    }
}

void s21_anti_gen_matrix(matrix_t *A) {
    double num = 32;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++, num--) {
            A->matrix[i][j] = num;
        }
    }
}
