#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <random>

constexpr int MINIMAL = 1;
constexpr int MAXIMAL = 20;
constexpr double REAL_ONE = 1.0;
constexpr double REAL_ZERO = 0.0;

enum TypeOfElement {
    MATRIX,
    VECTOR
};

enum TypeOfFiling {
    KEYBOARD,
    FROMFILE,
    RANDOM,
    ANOTHERCHOICE
};

gsl_matrix* mull_coll_vector_to_row_vector(gsl_vector* vector, gsl_vector* vectorTranspone) {
    gsl_matrix_view gsl_matrix1 = gsl_matrix_view_vector(vectorTranspone, 1, vectorTranspone->size);
    gsl_matrix_view gsl_matrix2 = gsl_matrix_view_vector(vector, vector->size, 1);
    gsl_matrix* matrix = gsl_matrix_alloc(vector->size, vector->size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, REAL_ONE, &gsl_matrix2.matrix, &gsl_matrix1.matrix, REAL_ZERO, matrix);
    
    return matrix;
}

gsl_matrix* mull_transponVector_to_matrix(gsl_vector* vector, gsl_matrix* matrix) {
    gsl_matrix_view gsl_row_vector = gsl_matrix_view_vector(vector, 1, vector->size);
    gsl_matrix* result = gsl_matrix_alloc(1, vector->size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, REAL_ONE, &gsl_row_vector.matrix, matrix, REAL_ZERO, result);

    return result;
}

void randomizeTheInputInFile(int size,const char* fileName) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(MINIMAL, MAXIMAL);

    FILE* filePtr;
    filePtr = fopen(fileName, "w");
    if (filePtr == NULL)
    {
        printf("Error!");
        exit(1);
    }

    for (int i = 0; i < size; ++i) {
        fprintf(filePtr, "%d ", dist(rng));
    }
    fclose(filePtr);
}

void writeResultInFile(gsl_matrix* matrix) {
    FILE* file_ptr;
    file_ptr = fopen("RESULT.TXT", "w");
    if (file_ptr == NULL) {
        printf("Error!");
        exit(1);
    }

    for (int i = 0; i < matrix->size1; i++) {
        for (int j = 0; j < matrix->size2; j++) {
            fprintf(file_ptr, "%lf ", gsl_matrix_get(matrix, i, j));
        }
        fprintf(file_ptr, "\n");
    }
    fclose(file_ptr);
}

void printfMatrix(gsl_matrix* matrix, int size, const char* nameOfElement = "Matrix") {
    printf("%s: \n", nameOfElement);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%lf ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
}

void printfMatrix(gsl_matrix* matrix, int size1, int size2, const char* nameOfElement = "Matrix") {
    printf("%s: \n", nameOfElement);
    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            printf("%lf ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
}

void printVector(gsl_vector* vector, int size, const char* nameOfElement = "Vector") {
    printf("%s: \n", nameOfElement);
    for (int i = 0; i < size; ++i) {
        printf("%lf\n", gsl_vector_get(vector, i));
    }
    printf("\n");
}

void inputVectorFromFile(gsl_vector* vector, int size, const char* fileName) {
    FILE* filePtr;
    const int sizeOfFileName = strlen(fileName) + 5;
 
    char *ch = new char[sizeOfFileName];
    strcpy(ch, fileName);
    strcat(ch, ".txt");

    filePtr = fopen(ch, "r");
    if (filePtr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    int val;
    for (int i = 0; i < size; ++i) {
        fscanf(filePtr, "%d", &val);
        gsl_vector_set(vector, i, val);
    }
    free(ch);
    fclose(filePtr);
}

void inputMatrixFromFile(gsl_matrix* matrix, int size, const char* fileName) {
    FILE* filePtr;

    const int sizeOfFileName = strlen(fileName) + 5;
    char* ch = new char[sizeOfFileName];
    strcpy(ch, fileName);
    strcat(ch, ".txt");

    filePtr = fopen(ch, "r");
    if (filePtr == NULL)
    {
        printf("Error!");
        exit(1);
    }
    int val;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            fscanf(filePtr, "%d", &val);
            gsl_matrix_set(matrix, i, j, val);
        }
    }
    fclose(filePtr);
}

void inputMatrixFromKeyboard(gsl_matrix* matrix, int size, const char * nameOfElement = "Matrix") {
    int val;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%s[%d][%d]: ", nameOfElement, i, j);
            scanf("%d", &val);
            gsl_matrix_set(matrix, i, j, val);
        }
    }
}

void inputMatrixWithRandom(gsl_matrix* matrix, int size, const char* nameOfElement = "Matrix") {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(MINIMAL, MAXIMAL);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            gsl_matrix_set(matrix, i, j, dist(rng));
        }
    }
}

void inputVectorFromKeyboard(gsl_vector* vector, int size, const char* nameOfElement = "Vector") {
    int value;
    for (int i = 0; i < size; ++i) {
        printf("%s[%d]: ", nameOfElement, i);
        scanf("%d", &value);
        gsl_vector_set(vector, i, value);
    }
}

void inputVectorWithRandom(gsl_vector* vector, int size, const char* nameOfElement = "Vector") {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(MINIMAL, MAXIMAL);

    for (int i = 0; i < size; ++i) {
        gsl_vector_set(vector, i, dist(rng));
    }
}

void wayOfFilling(void* element, int size, TypeOfElement type, const char* nameOfElement = "Matrix") {
    gsl_matrix* tmp_mat = NULL;
    gsl_vector* tmp_vec = NULL;

    if (type == MATRIX)
        tmp_mat = (gsl_matrix*)element;
    else
        tmp_vec = (gsl_vector*)element;

    printf("\nSelect an option for filling %s:\n", nameOfElement);
    printf("1. Fill array with random numbers\n");
    printf("2. Fill array with user input\n");
    printf("3. Fill array with data from file\n");
    printf("4. Quit\n");

    int choice = 0;
    bool runWhile = false;
    do {
        scanf("%d", &choice);
        switch (choice) {
        case 1:
            if (type == MATRIX) {
                inputMatrixWithRandom(tmp_mat, size, nameOfElement);
                printfMatrix(tmp_mat, size, nameOfElement);
            }
            else 
            {
                inputVectorWithRandom(tmp_vec, size, nameOfElement);
                printVector(tmp_vec, size, nameOfElement);
            }

            runWhile = false;
            break;
        case 2:
            if (type == MATRIX)
                inputMatrixFromKeyboard(tmp_mat, size, nameOfElement);
            else
                inputVectorFromKeyboard(tmp_vec, size, nameOfElement);
            runWhile = false;
            break;
        case 3:
            if (type == MATRIX)
            {
                inputMatrixFromFile(tmp_mat, size, nameOfElement);
                printfMatrix(tmp_mat, size, nameOfElement);
            }
            else
            {
                inputVectorFromFile(tmp_vec, size, nameOfElement);
                printVector(tmp_vec, size, nameOfElement);
            }
            runWhile = false;
            break;
        case 4:
            exit(1);
            break;
        default:
            printf("Error! Enter number is invalid please make ur choice one more time!\n");
            runWhile = true;
            break;
        }
    } while (runWhile);
}

TypeOfFiling mainMenu() {
    printf("Choose how to fill in the data:\n");
    printf("1. All data random\n");
    printf("2. All data from keyboard\n");
    printf("3. All data from file\n");
    printf("4. Each data is different\n");
    int choice;
    scanf("%d", &choice);
    if (choice == 1)
        return RANDOM;
    if (choice == 2)
        return KEYBOARD;
    if (choice == 3)
        return FROMFILE;
    return ANOTHERCHOICE;
}

gsl_vector* matrix_mul_vector(gsl_matrix* matrix, gsl_vector* vector) {
    gsl_vector* resultVector = gsl_vector_alloc(vector->size);
    gsl_blas_dgemv(CblasNoTrans, REAL_ONE, matrix, vector, REAL_ZERO, resultVector);
    return resultVector;
}

int printResult() {
    printf("Do you want to see all results?: ");
    char ch;
    scanf(" %c", &ch);
    switch(ch){
    case 'y':
    case 'Y':
    case '1':
        return 1;
    }
    return 0;
}

int main(int argc, char* argv[])
{
    int N;
    printf("Enter size of matrix n: ");
    scanf("%d", &N);

    
    //Allocate memory for data structs
    gsl_matrix* A = gsl_matrix_alloc(N, N);
    gsl_matrix* A1 = gsl_matrix_alloc(N, N);
    gsl_matrix* A2 = gsl_matrix_alloc(N, N);
    gsl_matrix* B2 = gsl_matrix_alloc(N, N);
    gsl_matrix* С2 = gsl_matrix_alloc(N, N);

    //new one
    gsl_matrix* y1_transp_mul_y1_mul_y3power2 = gsl_matrix_alloc(N, N);
    gsl_matrix* x = gsl_matrix_alloc(1, N);
    gsl_matrix* Y3 = gsl_matrix_alloc(N, N);
    gsl_matrix* Y3_Power_2 = gsl_matrix_alloc(N, N);
    gsl_matrix* Y3_Power_3 = gsl_matrix_alloc(N, N);
    gsl_matrix* Y3_Power_2_temp = gsl_matrix_alloc(N, N);
    gsl_matrix* y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp = gsl_matrix_alloc(N, N);
    gsl_matrix* y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3 = gsl_matrix_alloc(N, N);
    gsl_matrix* y1Transp_mull_y1_mull_Y3 = gsl_matrix_alloc(N, N);
    gsl_matrix* y1Transp_mull_y1_mull_Y3_add_Y3power3 = gsl_matrix_alloc(N, N);


    gsl_vector* b = gsl_vector_alloc(N);
    gsl_vector* b1 = gsl_vector_alloc(N);
    gsl_vector* c1 = gsl_vector_alloc(N);




    system("cls");

    //Generate numbers in file
    /*int doubleSize = N * N;
    randomizeTheInputInFile(doubleSize, "A.txt");
    randomizeTheInputInFile(doubleSize, "A1.txt");
    randomizeTheInputInFile(doubleSize, "A2.txt");
    randomizeTheInputInFile(doubleSize, "B2.txt");
    randomizeTheInputInFile(doubleSize, "C2.txt");
    randomizeTheInputInFile(N, "b.txt");
    randomizeTheInputInFile(N, "b1.txt");
    randomizeTheInputInFile(N, "c1.txt");*/

    //Fill the data
    TypeOfFiling tf = mainMenu();
    switch (tf) {
    case RANDOM:

        inputMatrixWithRandom(A, N, "A");
        inputMatrixWithRandom(A, N, "A");
        inputMatrixWithRandom(A1, N, "A1");
        inputMatrixWithRandom(A2, N, "A2");
        inputMatrixWithRandom(B2, N, "B2");
        

        inputVectorWithRandom(b1, N, "b1");
        inputVectorWithRandom(c1, N, "c1");

        break;
    case KEYBOARD:

        inputMatrixFromKeyboard(A, N, "A");
        inputMatrixFromKeyboard(A, N, "A");
        inputMatrixFromKeyboard(A1, N, "A1");
        inputMatrixFromKeyboard(A2, N, "A2");
        inputMatrixFromKeyboard(B2, N, "B2");
       

        inputVectorFromKeyboard(b1, N, "b1");
        inputVectorFromKeyboard(c1, N, "c1");
        break;
    case FROMFILE:

        inputMatrixFromFile(A, N, "A");
        inputMatrixFromFile(A, N, "A");
        inputMatrixFromFile(A1, N, "A1");
        inputMatrixFromFile(A2, N, "A2");
        inputMatrixFromFile(B2, N, "B2");
        

        inputVectorFromFile(b1, N, "b1");
        inputVectorFromFile(c1, N, "c1");
        break;

    case ANOTHERCHOICE:

        wayOfFilling(A, N, MATRIX, "A");
        wayOfFilling(A1, N, MATRIX, "A1");
        wayOfFilling(A2, N, MATRIX, "A2");
        wayOfFilling(B2, N, MATRIX, "B2");
        
        wayOfFilling(b1, N, VECTOR, "b1");
        wayOfFilling(c1, N, VECTOR, "c1");

        break;
    }

    system("cls");

    //print input data
    printfMatrix(A, N, "A"); printf("\n");
    printfMatrix(A1, N, "A1"); printf("\n");
    printfMatrix(A2, N, "A2"); printf("\n");
    printfMatrix(B2, N, "B2"); printf("\n");

    printVector(b1, N, "b1"); printf("\n");
    printVector(c1, N, "c1"); printf("\n");


    //Start of calculation
    //Calculate b
    for (int i = 0; i < N; ++i) {
        gsl_vector_set(b, i, (i + 1) * 9);
    }
    

    //Calculate C2
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            gsl_matrix_set(С2, i, j, (1.0 / ((i + 1) + (j + 1))));
        }
    }
    
    gsl_vector* y1 = matrix_mul_vector(A, b);
    
    gsl_vector_sub(b1, c1);
    
    gsl_vector* y2 = matrix_mul_vector(A1, b1);
    
    //B2 + C2
    gsl_matrix_add(B2, С2);
    
    
    //A2 (B2 + C2)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, REAL_ONE, A2, B2, REAL_ZERO, Y3);
    
    //gsl_vector_mul(a, b, c);
    double y1Transpon_mul_y1;
    gsl_blas_ddot(y1, y1, &y1Transpon_mul_y1);
    
    gsl_matrix* y2_mul_y2Transpon = mull_coll_vector_to_row_vector(y2, y2);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, REAL_ONE, Y3, Y3, REAL_ZERO, Y3_Power_2);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, REAL_ONE, Y3_Power_2, Y3, REAL_ZERO, Y3_Power_3);
    
    gsl_matrix_memcpy(Y3_Power_2_temp, Y3_Power_2);
    gsl_matrix_scale(Y3_Power_2_temp, y1Transpon_mul_y1);
    
    gsl_matrix_memcpy(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp, Y3_Power_2_temp);
    gsl_matrix_add(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp, y2_mul_y2Transpon);
    
    gsl_matrix_memcpy(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3, y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp);
    gsl_matrix_add(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3, Y3);
    
    gsl_matrix_memcpy(y1Transp_mull_y1_mull_Y3, Y3);
    gsl_matrix_scale(y1Transp_mull_y1_mull_Y3, y1Transpon_mul_y1);
    
    gsl_matrix_memcpy(y1Transp_mull_y1_mull_Y3_add_Y3power3, y1Transp_mull_y1_mull_Y3);
    gsl_matrix_add(y1Transp_mull_y1_mull_Y3_add_Y3power3, Y3_Power_3);
    
    gsl_matrix* last_left_side_of_example = mull_transponVector_to_matrix(y1, y1Transp_mull_y1_mull_Y3_add_Y3power3);
    
    gsl_matrix* last_right_side_of_example = mull_transponVector_to_matrix(y2, y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3);
    
    //find x
    gsl_matrix_memcpy(x, last_left_side_of_example);
    gsl_matrix_add(x, last_right_side_of_example);
    
    //write answer in the file

    if (printResult())
    {
        printVector(b, N, "b"); printf("\n");
        printfMatrix(С2, N, "C2"); printf("\n");
        printVector(y1, N, "y1"); printf("\n");
        printVector(b1, N, "b1_minus_c1"); printf("\n");
        printVector(y2, N, "y2"); printf("\n");
        printfMatrix(B2, N, "B2 + C2"); printf("\n");
        printfMatrix(Y3, N, "Y3"); printf("\n");
        printf("y1' * y1\n%lf", y1Transpon_mul_y1); printf("\n");
        printfMatrix(y2_mul_y2Transpon, N, "y2 * y2'"); printf("\n");
        printfMatrix(Y3_Power_2, N, "Y3^2"); printf("\n");
        printfMatrix(Y3_Power_3, N, "Y3^3"); printf("\n");
        printfMatrix(Y3_Power_2_temp, N, "y1' * y1 * Y3^2"); printf("\n");
        printfMatrix(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp, N, "y1' * y1 * Y3^2 + y2 * y2'"); printf("\n");
        printfMatrix(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3, N, "y1' * y1 * Y3^2 + y2 * y2' + Y3"); printf("\n");
        printfMatrix(y1Transp_mull_y1_mull_Y3, N, "y1' * y1 * Y3"); printf("\n");
        printfMatrix(y1Transp_mull_y1_mull_Y3_add_Y3power3, N, "y1' * y1 * Y3 + Y3^3"); printf("\n");
        printfMatrix(last_left_side_of_example, 1, N, "y1'(y1' * y1 * Y3 + Y3^3)"); printf("\n");
        printfMatrix(last_right_side_of_example, 1, N, "y2'(y1' * y1 * Y3^2 + y2 * y2' + Y3)"); printf("\n");
    }
    printf("\n\nRESULT:\n");
    printfMatrix(x, 1, N, "x = y1'(y1' * y1 * Y3 + Y3^3) + y2'(y1' * y1 * Y3^2 + y2 * y2' + Y3)"); printf("\n");
    writeResultInFile(x);

    //free allocated memory
    gsl_matrix_free(A);
    gsl_matrix_free(A1);
    gsl_matrix_free(A2);
    gsl_matrix_free(B2);
    gsl_matrix_free(С2);
    gsl_vector_free(b);
    gsl_vector_free(b1);
    gsl_vector_free(c1);
    gsl_vector_free(y1);
    gsl_vector_free(y2);
    gsl_matrix_free(Y3);
    gsl_matrix_free(Y3_Power_3);
    gsl_matrix_free(Y3_Power_2);
    gsl_matrix_free(y2_mul_y2Transpon);
    gsl_matrix_free(last_left_side_of_example);
    gsl_matrix_free(last_right_side_of_example);
    gsl_matrix_free(y1_transp_mul_y1_mul_y3power2);
    gsl_matrix_free(x);
    gsl_matrix_free(Y3_Power_2_temp);
    gsl_matrix_free(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp);
    gsl_matrix_free(y1Transp_mull_y1_mull_Y3power2_add_y2_mull_y2Transp_add_Y3);
    gsl_matrix_free(y1Transp_mull_y1_mull_Y3);
    gsl_matrix_free(y1Transp_mull_y1_mull_Y3_add_Y3power3);
	return 0;
}

