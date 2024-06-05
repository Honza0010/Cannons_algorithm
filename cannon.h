
#ifndef __CANNON__
#define __CANNON__

#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
#include <utility>
#include <fstream>
#include <random>
#include <cmath>

//Creation of matrix which is contiguous in memory  /   vytvoøení 2D matice, která má uložené elementy v pamìti za sebou
template <typename T>
T** createMatrix(int rows, int cols, const T& val = T())
{
    if (rows == 0 || cols == 0)
        throw std::invalid_argument("number of rows is 0");

    T** rowPtrs = nullptr;
    T* pool = nullptr;
    try
    {
        rowPtrs = new T * [rows];  //Allocation of rows (pointers)  /   Ukazatele na řádky
        pool = new T[rows * cols];  // allocation of columns (one array and row pointers will be set up uniformly in it)    /   vytvoření dat v jednom poli, kdy ukazatele na řádky budou rozděleny podle velikosti sloupce
        T* startpool = pool;  //Start of the pool   /   ukazatel na začátek poolu

        //the row pointers are set up to point to the appropriate positions in the memory pool  /   řádkové ukazatele jsou nastaveny tak, aby ukazovaly na správné místo v paměti
        for (int i = 0; i < rows; ++i, pool += cols)
            rowPtrs[i] = pool;

        // Initialize pool  /   zaplnění pole hodnotami
        std::fill(startpool, startpool + rows * cols, val);
        return rowPtrs;
    }
    catch (std::bad_alloc& ex)
    {
        delete[] rowPtrs;
        throw ex;  // memory allocation error
    }
}

template <typename T>
void deleteMatrix(T** matrix)
{
    delete[] matrix[0];  // remove the pool    /   smazání dat
    delete[] matrix;     // remove the pointers    /   smazání řádkových ukazatelů
}


template<typename T>
void printMatrix(T** A, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            //std::cout << std::setw(5) << A[i][j] << " ";
            std::cout << A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}


template<typename T, typename = void>
struct has_add_assign : std::false_type {};

template<typename T>
struct has_add_assign<T, std::void_t<decltype(std::declval<T&>() += std::declval<const T&>())>> : std::true_type {};

template<typename T>
std::enable_if_t<has_add_assign<T>::value, void> addMatrix(T** src, int rows, int cols, T** other)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            src[i][j] += other[i][j];
        }
    }
}


template<typename T>
void multiplyMatrices(T** A, T** B, int rows, int cols, T** res)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int val = 0;
            for (int k = 0; k < rows; k++)
            {
                val += A[i][k] * B[k][j];
            }
            res[i][j] = val;
        }
    }
}

template<typename T>
void copyMatrix(T** A, T** B, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            A[i][j] = B[i][j];
        }
    }
}


template<typename T, int blockSize = 16>
void multiplyMatricesBlockwise(T** A, T** B, int rows, int cols, T** res)
{
    if (rows % blockSize != 0 || cols % blockSize != 0)
    {
        throw std::out_of_range("The matrix cannot be divided into submatrices of size (blockSize x blockSize)");
    }

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            res[i][j] = 0;
        }
    }
    //std::cout << blockSize << std::endl;
    for (int i = 0; i < rows; i += blockSize)
    {
        for (int j = 0; j < cols; j += blockSize)
        {
            //int val = 0;
            for (int k = 0; k < rows; k += blockSize)
            {
                for (int bi = 0; bi < blockSize; bi++)
                {
                    for (int bj = 0; bj < blockSize; bj++)
                    {
                        for (int bk = 0; bk < blockSize; bk++)
                        {
                            //val += A[i+bi][k+bk] * B[k+bk][j+bj];
                            res[i + bi][j + bj] += A[i + bi][k + bk] * B[k + bk][j + bj];
                        }
                    }
                }
            }
            //res[i + bi][j + bj] = val;
        }
    }
}

template<typename T, int blockSize = 16>
void multiplyMatricesBlockwiseOptimalized(T** A, T** B, int rows, int cols, T** res)
{
    if (rows % blockSize != 0 || cols % blockSize != 0)
    {
        throw std::out_of_range("The matrix cannot be divided into submatrices of size (blockSize x blockSize)");
    }

    // Initialize the result matrix to zero
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            res[i][j] = 0;
        }
    }

    for (int i = 0; i < rows; i += blockSize)
    {
        for (int j = 0; j < cols; j += blockSize)
        {
            for (int k = 0; k < rows; k += blockSize)
            {
                T* resRow = &res[i][j];
                T* ARow = &A[i][k];
                for (int i1 = 0; i1 < blockSize; ++i1, resRow += rows, ARow += rows)
                {
                    T* BRow = &B[k][j];
                    for (int k1 = 0; k1 < blockSize; ++k1, BRow += rows)
                    {
                        for (int j1 = 0; j1 < blockSize; ++j1)
                        {
                            resRow[j1] += ARow[k1] * BRow[j1];
                        }
                    }
                }
            }
        }
    }
}


template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
T** readMatrix(const std::string& filename, int& rows_, int& cols_)     //Reads matrix from the file. The first row of the file has to contain two integers for number of rows and columns
{                                                                       //Načítání matice ze souboru. První řádek musí obsahovat 2 čísla - počet řádků a sloupců
    int rows, cols;

    std::ifstream file(filename, std::ios::in);

    if (!file.is_open())
    {
        throw std::ios_base::failure("The could not be opened");
    }

    file >> rows >> cols;
    rows_ = rows;
    cols_ = cols;
    T** matrix = createMatrix(rows, cols, 0);

    T val;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            file >> val;
            matrix[i][j] = val;
            if ((i < rows - 1 || j < cols - 1) && file.eof())
            {
                throw std::out_of_range("Too few elements in the file");
                file.close();
            }
            if (file.fail())
            {
                throw std::invalid_argument("Invalid data type");
                file.close();
            }
        }
    }
    
    file.close();

    return matrix;
}

double roundToDecimalPlaces(double value, int decimalPlaces)    //Rounds number to the number with maximum of decimalPlaces digits after decimal point
{                                                               //Zaokrouhlí číslo, aby mělo maximálně 4 desetinná místa
    double factor = std::pow(10, decimalPlaces);
    return std::round(value * factor) / factor;
}


template<typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
void randFillMatrix(T** matrix, int rows, int cols)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    if constexpr (std::is_integral_v<T>) {
        // uniform_int_distribution for integer types   /   rovnoměrné rozdělení celých čísel
        std::uniform_int_distribution<T> dis(0, 100); // Range [0, 100]
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = dis(gen);
            }
        }
    }
    else if constexpr (std::is_floating_point_v<T>) {
        // uniform_real_distribution for floating-point types   /   Rovnoměrné rozdělení reálných čísel
        std::uniform_real_distribution<T> dis(0.0, 100.0); // Range [0.0, 100.0]
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = roundToDecimalPlaces(dis(gen), 4);
            }
        }
    }
}


#endif // !__CANNON__

