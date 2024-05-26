#pragma once
#ifndef __CANNON__
#define __CANNON__

#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
#include <utility>


//Creation of matrix which is contiguous in memory
template <typename T>
T** createMatrix(int rows, int cols, const T& val = T())
{
    if (rows == 0 || cols == 0)
        throw std::invalid_argument("number of rows is 0");

    T** ptr = nullptr;
    T* pool = nullptr;
    try
    {
        ptr = new T * [rows];  // allocate pointers (can throw here)
        pool = new T[rows * cols];  // allocate pool (can throw here)
        T* startpool = pool;  // Remember the start of the pool

        // now point the row pointers to the appropriate positions in
        // the memory pool
        for (int i = 0; i < rows; ++i, pool += cols)
            ptr[i] = pool;

        // Initialize pool
        std::fill(startpool, startpool + rows * cols, val);
        return ptr;
    }
    catch (std::bad_alloc& ex)
    {
        delete[] ptr; // either this is nullptr or it was allocated
        throw ex;  // memory allocation error
    }
}

template <typename T>
void deleteMatrix(T** arr)
{
    delete[] arr[0];  // remove the pool
    delete[] arr;     // remove the pointers
}


template<typename T>
void printMatrix(T** A, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            std::cout << std::setw(4) << A[i][j] << " ";
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




#endif // !__CANNON__
