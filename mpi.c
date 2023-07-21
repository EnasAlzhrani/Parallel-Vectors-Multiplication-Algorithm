#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

double Sequential_code(int n)
{
    int i, j, total = 0;
    double start, end, time;
    int* vector1 = (int*)malloc(n * sizeof(int));
    int* vector2 = (int*)malloc(n * sizeof(int));
    int* multiplied = (int*)malloc(n * sizeof(int));

    for (i = 0; i < n; i++) { //Randomizing the values of each element in the two vectors
        vector1[i] = rand() % 100;
        vector2[i] = rand() % 100;

    }

    start = clock(); //Start the clock

    for (j = 0; j < n; j++) { //Multiplying the two vectors and storing the results in a new vector
        multiplied[j] = (vector1[j] * vector2[j]);
        total = total + multiplied[j]; //Calculating the addition of the values
    }

    end = clock();
    time = end - start;


    //printf("The total of the multiplied vectors = %d", total); //Printing the total 

    //printf("\nThe time the algorithm took is: %.2lf s", time); //Printing the time it took for the algorithm calculations to execute

    return time;
}


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, np;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned long long  n;
    if (argc == 2)
    {
        n = atoi(argv[1]);
    }
    else
    {
        if(rank==0)
            printf("Give size of the vector in command line argument\neg:(mpiexec -np 3 mpi.exe 10) 10 is the size of the vector the size could be any number\n");
        MPI_Finalize();
        exit(0);
    }

    double sequntial_time = 0;
    if (rank == 0)
        sequntial_time = Sequential_code(n);

    double start, end, time; //Calculating the time it takes to run the algorithm


    unsigned long long* vector1 = NULL, * vector2 = NULL, * multiplied = NULL;
    int* send_counts = NULL, * displs = NULL;
    unsigned long long  rank_size;
    MPI_Status status;
    if (rank == 0)
    {
        vector1 = (unsigned long long*)malloc(n * sizeof(unsigned long long));
        vector2 = (unsigned long long*)malloc(n * sizeof(unsigned long long));
        send_counts = (int*)malloc(np * sizeof(int));
        displs = (int*)malloc(np * sizeof(int));

        for (int i = 0; i < n; i++) { //Randomizing the values of each element in the two vectors
            vector1[i] = 1;
            vector2[i] = 1;

        }
        
        if (n % np != 0)
        {
            rank_size = n / np;
            unsigned long long  temp_size = n / np;
            send_counts[0] = temp_size;
            for (int i = 1; i < np; i++)
            {
                if (i == np - 1)
                    temp_size += n % np;

                MPI_Send(&temp_size, 1, MPI_UNSIGNED_LONG_LONG, i, 1, MPI_COMM_WORLD);
                send_counts[i] = temp_size;

            }
        }
        else
        {
            rank_size = n / np;
            send_counts[0] = rank_size;
            for (int i = 1; i < np; i++)
            {

                MPI_Send(&rank_size, 1, MPI_UNSIGNED_LONG_LONG, i, 1, MPI_COMM_WORLD);
                send_counts[i] = rank_size;

            }
        }
        displs[0] = 0;
        for (int i = 1; i < np; i++)
            displs[i] = displs[i - 1] + send_counts[i - 1];

    }
    else
    {
        MPI_Recv(&rank_size, 1, MPI_UNSIGNED_LONG_LONG, 0, 1, MPI_COMM_WORLD, &status);
    }

    
    unsigned long long* rank_vector1 = NULL, * rank_vector2 = NULL;
    unsigned long long  rank_total = 0;

    rank_vector1 = (unsigned long long*)malloc(rank_size * sizeof(unsigned long long));
    rank_vector2 = (unsigned long long*)malloc(rank_size * sizeof(unsigned long long));
    multiplied = (unsigned long long*)malloc(rank_size * sizeof(unsigned long long));

    MPI_Scatterv(vector1, send_counts, displs, MPI_UNSIGNED_LONG_LONG, rank_vector1, rank_size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Scatterv(vector2, send_counts, displs, MPI_UNSIGNED_LONG_LONG, rank_vector2, rank_size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
    {
        free(vector1); 
        free(vector2);
    }
    
    start = MPI_Wtime(); //Start the clock

    for (int j = 0; j < rank_size; j++) { //Multiplying the two vectors and storing the results in a new vector
        multiplied[j] = (rank_vector1[j] * rank_vector2[j]);
        rank_total = rank_total + multiplied[j]; //Calculating the addition of the values
    }
    


    end = MPI_Wtime();
    time = end - start;

    unsigned long long global_total;
    
    MPI_Reduce(&rank_total, &global_total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
    {
        printf("Total number of process = %d\n", np);
        printf("Total vector size = %lld\n", n);
        printf("The total of the multiplied vectors = %lld", global_total); //Printing the total 

        printf("\nThe mpi time the algorithm took is: %.8lf s", time); //Printing the time it took for the algorithm calculations to execute

        printf("\nSpeedup of the mpi time of the algorithm is: %.2lf sec \n", sequntial_time / time);

    }
    


    MPI_Finalize();
    return 0;
}
