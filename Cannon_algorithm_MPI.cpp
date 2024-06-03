#include <iostream>
#include <mpi.h>
#include <vector>

#include "cannon.h"
#include "TimerRT.h"

#include "ArrayMatrix.h"


void testSequenceMultiplication(int matrixSize, double& seqentialTime)
{
    ArrayMatrix<int> A_(matrixSize), B_(matrixSize), C_(matrixSize);
    TimerRT timer;

    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            A_.setElement(i, j, i * matrixSize + j);
            B_.setElement(i, j, i * matrixSize + j);
        }
    }
    timer.reset();
    timer.start();
    if (matrixSize % 16 == 0)
    {
        C_.multiplyMatricesBlockwise3(A_, B_);
    }
    else
    {
        C_.multiplyMatrices(A_, B_);
    }
    //C_.multiplyMatrices(A_, B_);
    timer.stop();

    seqentialTime = timer.getTime();
    std::cout << "Sequential time: " << seqentialTime << std::endl;

    if (matrixSize <= 16)
    {
        std::cout << "Sequential multiplication result: " << std::endl;
        std::cout << C_ << std::endl;
    }
}


int main(int argc, char* argv[])
{
    double parallelTime, sequentialTime;
    TimerRT timer;
    int rows = 256*16;
    int cols = rows;     //Rozměry matice

    int rank, P;                //Cislo procesu / pocet procesu

    const int dimension = 2;    //Dimenze kartézského systému souřadnic procesů
    int dims[2];                //Velikost po osách x a y v systému procesů. dims[0] == dims[1] je nutnost.
    const int periodic[] = { 1,1 };     //osy x a y jsou zabaleny, tedy pokud např. proces 2 je úplně nejpravější proces po ose x a proces 1 nejlevější proces po ose x, pak pravý soused procesu 2 je proces 1 a levý soused procesu 1 je proces 2
    int blockSize;              //Velikost jednoho bloku (submatice) (rows / sqrt(P))
    int procDim;               //sqrt(P), rozměry kartézského systému souřadnic pro procesy

    std::vector<int> commonData(4, 0); // {rows, cols, blockDim, blocksize} - master process (0) přepošle všechna tato data pomocí MPI_Bcast do všech ostatních procesů
    int left_rank, right_rank, up_rank, down_rank;  //left/right/up/down sousedi současného procesu
    MPI_Comm newCommunicator;     //Komunikátor pro nový MPI grid

    int** A = nullptr, ** B = nullptr, ** C = nullptr;
    int** subA = nullptr, ** subB = nullptr, ** subC = nullptr;


    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (rank == 0)
    {
        if ((double)std::sqrt(P) - std::floor(std::sqrtf(P)) != 0)
        {
            std::cerr << "Number of processes is not perfect square!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
        int sqrtP = (int)std::sqrt(P);

        //Control if the matrices are square    /   Kontrola, zda jsou matice čtvercové
        if (rows != cols)
        {
            std::cerr << "Matrix must be square!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        //Control if the matrices can be divided equally into processes /   Kontrola, zda lze rozdělit matice do stejně velkých bloků
        if (rows % sqrtP != 0 || cols % sqrtP != 0)
        {
            std::cerr << "Number of rows / columns not divisible by square root of #Processes!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        procDim = sqrtP;
        blockSize = cols / sqrtP;

        std::cout << "Matrix size: " << rows << " x " << cols << std::endl;

        //Allocation of contiguous matrices /   alokace matic, které mají data uložena v paměti za sebou
        try
        {
            A = createMatrix<int>(rows, cols, 0);
            B = createMatrix<int>(rows, cols, 0);
            C = createMatrix<int>(rows, cols, 0);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    A[i][j] = i * rows + j;
                    B[i][j] = i * rows + j;
                }
            }

            //printMatrix(A, rows, cols);
        }
        catch (std::exception& e)
        {
            std::cerr << "Allocation of the matrices failed!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        commonData[0] = rows;
        commonData[1] = cols;
        commonData[2] = procDim;
        commonData[3] = blockSize;
    }
    timer.reset();
    timer.start();
    MPI_Bcast(&commonData[0], 4, MPI_INT, 0, MPI_COMM_WORLD);  //Master process B_CAST rows,cols,procDim,blockSize     /   hlavní proces rozešle tato data do všech ostatních procesů
    rows = commonData[0];
    cols = commonData[1];
    procDim = commonData[2];
    blockSize = commonData[3];

    dims[0] = procDim;         //Initializing cartesian dimensions for processes    /   dimenze kartézského systému souřadnic procesů
    dims[1] = procDim;

    MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periodic, 1, &newCommunicator);  //2D cartesian grid of processes wrapped  /   Kartézská soustava souřadnic procesů zaobalená do smyček

    try
    {
        subA = createMatrix<int>(blockSize, blockSize, 0);
        subB = createMatrix<int>(blockSize, blockSize, 0);
        subC = createMatrix<int>(blockSize, blockSize, 0);
    }
    catch (std::exception& e)
    {
        std::cerr << "Allocation of the matrices failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    MPI_Comm_rank(newCommunicator, &rank);    //Process getting a rank in new topology  /   procesu je přiřazen rank v nové topologii
    int coords[2];                      //Coordinates of the process in the grid    / souřadnice procesu v dané mřížce
    MPI_Cart_coords(newCommunicator, rank, 2, coords);

    //std::cout << "Process " << rank << " coords: [" << coords[0] << ", " << coords[1] << "]" << std::endl;

    MPI_Cart_shift(newCommunicator, 0, 1, &up_rank, &down_rank);    //Getting neighbours of current process along x axis    /   získání vertikálních sousedů
    MPI_Cart_shift(newCommunicator, 1, 1, &left_rank, &right_rank);     //Getting neighbours of current process along y axis    /   získání horizontálních sousedů

    //std::cout << "Process " << rank << " Neighbours: " << left_rank << " " << up_rank << " " << right_rank << " " << down_rank << std::endl;

    int* sendcounts = new int[P];   //sendCounts for MPI_Scatterv function  /   v tomto případě jedno políčko obsahuje počet submatic, kolik je posláno danému procesu
    int* displs = new int[P];    //Displacements for MPI_Scatterv function   /   v tomto případě jedno políčko obsahuje startovní pozici submatice, kterou daný proces dostane

    if (rank == 0) {
        for (int i = 0; i < P; i++) {
            sendcounts[i] = 1;
        }

        for (int i = 0; i < procDim; i++) {
            for (int j = 0; j < procDim; j++) {
                displs[i * procDim + j] = j + i * procDim * blockSize;   //The starting position is the number of the starting element among all elements.   /   Startovní pozice je pořadí daného políčka mezi všemi políčky. Využívá se tento vzorec, protože v poli displs se startovní pozice počítá tak, kolikrát musí přeskočit blocksize daného datového typu.
            }
        }
    }

    int matSize[2] = { rows, cols };    //Size of matrices A,B,C    /   velikost matic A,B,C
    int submatSize[2] = { blockSize, blockSize };   //Size of submatrices subA,subB,subC    /   velikost submatic subA,subB,subC
    int submatStart[2] = { 0,0 };
    MPI_Datatype type, submat;  //We have to create a MPI_Datatype for submatrices which will be scattered in processes     /   vytvoříme si datový typ pro MPI, který bude představovat danou submatici poslanou procesu
    MPI_Type_create_subarray(2, matSize, submatSize, submatStart, MPI_ORDER_C, MPI_INT, &type);
    MPI_Type_create_resized(type, 0, blockSize * sizeof(int), &submat);
    MPI_Type_free(&type);       //Nonresized type has to be deallocated     /   Nevyužitý datový typ smažeme
    MPI_Type_commit(&submat);   //We have to commit the MPI_Datatype to make it available for all processes     /   Musíme daný datový typ poslat všem procesům


    MPI_Scatterv(rank == 0 ? &(A[0][0]) : nullptr, sendcounts, displs, submat, &(subA[0][0]), blockSize * blockSize, MPI_INT, 0, MPI_COMM_WORLD);   //Sending submatrices into processes    /   rozdělování bloků matice do procesů
    MPI_Scatterv(rank == 0 ? &(B[0][0]) : nullptr, sendcounts, displs, submat, &(subB[0][0]), blockSize * blockSize, MPI_INT, 0, MPI_COMM_WORLD);

    /* if (rank == 3 && blockSize <= 8)
     {
         std::cout << "Process " << rank << " received localA" << std::endl;
         printMatrix(subA, blockSize, blockSize);
     }*/

     //Initial shifting allong x by i left    /   Počáteční rozmístění submatic do procesů
    for (int i = 0; i < coords[0]; i++)
    {
        MPI_Sendrecv_replace(&(subA[0][0]), blockSize * blockSize, MPI_INT, left_rank, 1, right_rank, 1, newCommunicator, MPI_STATUS_IGNORE);
    }
    // Initial shifting allong y by i up
    for (int i = 0; i < coords[1]; i++)
    {
        MPI_Sendrecv_replace(&(subB[0][0]), blockSize * blockSize, MPI_INT, up_rank, 1, down_rank, 1, newCommunicator, MPI_STATUS_IGNORE);
    }

    int** helpRes = nullptr;    //Local results, one part which will be added to subC, because subC is formed by addition of multiples from different submatrices   /   Pomocná matice pro výpočet subC, protože subC je tvořeno součtem násobků různých bloků
    try
    {
        helpRes = createMatrix<int>(blockSize, blockSize, 0);
    }
    catch (std::exception& e)
    {
        std::cerr << "Allocation of the matrices failed!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }



    /* int** pomA = nullptr, ** pomB = nullptr;
     pomA = createMatrix<int>(blockSize, blockSize,  0);
     pomB = createMatrix<int>(blockSize, blockSize, 0);*/

     //The algorithm has sqrt(P) steps
    for (int steps = 0; steps < procDim; steps++)
    {
        if (blockSize % 16 == 0)
        {
            multiplyMatricesBlockwiseOptimalized(subA, subB, blockSize, blockSize, helpRes);
        }
        else
        {
            multiplyMatrices(subA, subB, blockSize, blockSize, helpRes);
        }
        //multiplyMatricesBlockwiseOptimalized<int, 4>(subA, subB, blockSize, blockSize, helpRes);

        addMatrix(subC, blockSize, blockSize, helpRes);     //Adding local results to subC  /   sčítání sčítanců pro subC

        MPI_Sendrecv_replace(&(subA[0][0]), blockSize * blockSize, MPI_INT, left_rank, 1, right_rank, 1, newCommunicator, MPI_STATUS_IGNORE);    //Shifting by 1 left along x axis  /   posunutí o jedna vlevo podél osy x
        MPI_Sendrecv_replace(&(subB[0][0]), blockSize * blockSize, MPI_INT, up_rank, 1, down_rank, 1, newCommunicator, MPI_STATUS_IGNORE);       //Shifting by 1 up along y axis /   posunutí o jedna nahoru podél osy y
        //if (coords[1] == 0)     
        //{
        //    //sendVal = rank;
        //    MPI_Send(&(subA[0][0]), blockSize * blockSize, MPI_INT, left_rank, 0, newCommunicator);
        //    MPI_Recv(&(subA[0][0]), blockSize * blockSize, MPI_INT, right_rank, 0, newCommunicator, MPI_STATUS_IGNORE);
        //}
        //else
        //{
        //    copyMatrix(pomA, subA, blockSize, blockSize);
        //    MPI_Recv(&(subA[0][0]), blockSize* blockSize, MPI_INT, right_rank, 0, newCommunicator, MPI_STATUS_IGNORE);
        //    MPI_Send(&(pomA[0][0]), blockSize* blockSize, MPI_INT, left_rank, 0, newCommunicator);
        //}
        //if (coords[0] == 0)
        //{
        //    //sendVal = rank;
        //    MPI_Send(&(subB[0][0]), blockSize * blockSize, MPI_INT, up_rank, 0, newCommunicator);
        //    MPI_Recv(&(subB[0][0]), blockSize * blockSize, MPI_INT, down_rank, 0, newCommunicator, MPI_STATUS_IGNORE);
        //}
        //else
        //{
        //    copyMatrix(pomB, subB, blockSize, blockSize);
        //    MPI_Recv(&(subB[0][0]), blockSize * blockSize, MPI_INT, down_rank, 0, newCommunicator, MPI_STATUS_IGNORE);
        //    MPI_Send(&(pomB[0][0]), blockSize * blockSize, MPI_INT, up_rank, 0, newCommunicator);
        //}
    }

    //deleteMatrix(pomA);
    //deleteMatrix(pomB);

    //Creating final C matrix, C = AB.  /   Vytváření výsledné matice C = AB
    MPI_Gatherv(&(subC[0][0]), rows * cols / P, MPI_INT, rank == 0 ? &(C[0][0]) : nullptr, sendcounts, displs, submat, 0, MPI_COMM_WORLD);


    timer.stop();
    if (rank == 0)
    {
        parallelTime = timer.getTime();
        std::cout << "Number of processes: " << P << " time: " << parallelTime << std::endl;
    }

    if (rank == 0 && rows <= 16)
    {
        std::cout << "Parallel multiplication result: " << std::endl;
        printMatrix(C, rows, cols);
    }

    if (rank == 0)
    {
        testSequenceMultiplication(rows, sequentialTime);
        double speedUp = sequentialTime / parallelTime;
        std::cout << "Sequential time / parallel time = " << sequentialTime << " / " << parallelTime << " = " << speedUp << std::endl;
        std::cout << "Efficiency = " << speedUp << " / " << P << " = " << speedUp / P << std::endl;
    }

    /* if (rank == 0)
     {
         multiplyMatricesBlockwise<int, 2>(A, B, rows, cols, C);
         printMatrix(C, rows, cols);
     }*/

    deleteMatrix(subA);
    deleteMatrix(subB);
    deleteMatrix(subC);
    if (rank == 0)
    {
        deleteMatrix(A);
        deleteMatrix(B);
        deleteMatrix(C);
    }

    delete[] sendcounts;
    delete[] displs;

    MPI_Finalize();

    return 0;
}

