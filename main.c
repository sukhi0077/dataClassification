//
//  main.c
//  MyTest
//
//  Created by Sukhwinder Singh on 14/06/2021.
//
/* Compiler command:
 * mpicc -O5 main.c -o fox
 *
 * Run command:
 * mpirun -n -4 ./fox
 */

#include <stdio.h>
#include "mpi.h"
#include <math.h>///
#include <stdlib.h>

typedef struct {
    int num_process;
    int g_order;
    int current_row;
    int current_col;
    int current_rank;
    MPI_Comm g_comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;

}
GRID_STRUCT;

#define MAX 1000
typedef struct {
    int n_bar;
    #define Order(A)((A) -> n_bar)
    float entries[MAX];
    #define Entry(A, i, j)( * (((A) -> entries) + ((A) -> n_bar) * (i) + (j)))
}
TEMP_MATRIX;

TEMP_MATRIX * Allocate_matrix_space(int n_bar);
void Clear_matrix_space(TEMP_MATRIX ** in_matrix_A);
void Read_matrix(char * prompt, TEMP_MATRIX * in_matrix_A,GRID_STRUCT * grid, int n);
void Print_matrix(char * title, TEMP_MATRIX * in_matrix_A,GRID_STRUCT * grid, int n);
void Reset_matrix(TEMP_MATRIX * in_matrix_A);
void Temp_matrix_calculation(TEMP_MATRIX * in_matrix_A,TEMP_MATRIX * in_matrix_B, TEMP_MATRIX * out_matrix_C);
void construct_m_type(TEMP_MATRIX * in_matrix_A);
MPI_Datatype local_matrix_mpi_t;
TEMP_MATRIX * temp_mat;

int main(int argc, char * argv[]) {
    int num_process;
    int current_rank;
    GRID_STRUCT grid;
    TEMP_MATRIX * in_matrix_A;
    TEMP_MATRIX * in_matrix_B;
    TEMP_MATRIX * out_matrix_C;
    int n;
    int n_bar;

    void Setup_grid(GRID_STRUCT * grid);
    void Fox(int n, GRID_STRUCT * grid, TEMP_MATRIX * in_matrix_A,
        TEMP_MATRIX * in_matrix_B, TEMP_MATRIX * out_matrix_C);

    MPI_Init( & argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & current_rank);

    Setup_grid( & grid);
    if (current_rank == 0) {
        printf("Size?\n");
        scanf("%d", & n);
    }

    MPI_Bcast( & n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    n_bar = n / grid.g_order;

    in_matrix_A = Allocate_matrix_space(n_bar);
    Order(in_matrix_A) = n_bar;
    Read_matrix("Enter first", in_matrix_A, & grid, n);

    in_matrix_B = Allocate_matrix_space(n_bar);
    Order(in_matrix_B) = n_bar;
    Read_matrix("Enter second", in_matrix_B, & grid, n);

    construct_m_type(in_matrix_A);
    temp_mat = Allocate_matrix_space(n_bar);

    out_matrix_C = Allocate_matrix_space(n_bar);
    Order(out_matrix_C) = n_bar;
    Fox(n, & grid, in_matrix_A, in_matrix_B, out_matrix_C);
    Print_matrix("The final matrix is", out_matrix_C, & grid, n);
    Clear_matrix_space( & in_matrix_A);
    Clear_matrix_space( & in_matrix_B);
    Clear_matrix_space( & out_matrix_C);

    MPI_Finalize();
}

void Setup_grid(
    GRID_STRUCT * grid) {
    int old_rank;
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int free_coords[2];

    MPI_Comm_size(MPI_COMM_WORLD, & (grid -> num_process));
    MPI_Comm_rank(MPI_COMM_WORLD, & old_rank);

    grid -> g_order = (int) sqrt((double) grid -> num_process);
    dimensions[0] = dimensions[1] = grid -> g_order;

    wrap_around[0] = wrap_around[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions,
        wrap_around, 1, & (grid -> g_comm));
    MPI_Comm_rank(grid -> g_comm, & (grid -> current_rank));
    MPI_Cart_coords(grid -> g_comm, grid -> current_rank, 2,
        coordinates);
    grid -> current_row = coordinates[0];
    grid -> current_col = coordinates[1];

    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid -> g_comm, free_coords, &
        (grid -> row_comm));

    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid -> g_comm, free_coords, &
        (grid -> col_comm));
}

void Fox(
    int n,
    GRID_STRUCT * grid,
    TEMP_MATRIX * in_matrix_A,
    TEMP_MATRIX * in_matrix_B,
    TEMP_MATRIX * out_matrix_C) {

    TEMP_MATRIX * temp_A;

    int stage;
    int bcast_root;
    int n_bar;
    int source;
    int dest;
    MPI_Status status;

    n_bar = n / grid -> g_order;
    Reset_matrix(out_matrix_C);

    source = (grid -> current_row + 1) % grid -> g_order;
    dest = (grid -> current_row + grid -> g_order - 1) % grid -> g_order;

    temp_A = Allocate_matrix_space(n_bar);

    for (stage = 0; stage < grid -> g_order; stage++) {
        bcast_root = (grid -> current_row + stage) % grid -> g_order;
        if (bcast_root == grid -> current_col) {
            MPI_Bcast(in_matrix_A, 1, local_matrix_mpi_t,
                bcast_root, grid -> row_comm);
            Temp_matrix_calculation(in_matrix_A, in_matrix_B,
                out_matrix_C);
        } else {
            MPI_Bcast(temp_A, 1, local_matrix_mpi_t,
                bcast_root, grid -> row_comm);
            Temp_matrix_calculation(temp_A, in_matrix_B,
                out_matrix_C);
        }
        MPI_Sendrecv_replace(in_matrix_B, 1, local_matrix_mpi_t,
            dest, 0, source, 0, grid -> col_comm, & status);

    }

}

TEMP_MATRIX * Allocate_matrix_space(int local_order) {
    TEMP_MATRIX * temp;

    temp = (TEMP_MATRIX * ) malloc(sizeof(TEMP_MATRIX));
    return temp;
}

void Clear_matrix_space(
    TEMP_MATRIX ** local_A_ptr) {
    free( * local_A_ptr);
}

void Read_matrix(
    char * prompt,
    TEMP_MATRIX * in_matrix_A,
    GRID_STRUCT * grid,
    int n) {

    int mat_row, mat_col;
    int grid_row, grid_col;
    int dest;
    int coords[2];
    float * temp;
    MPI_Status status;

    if (grid -> current_rank == 0) {
        temp = (float * ) malloc(Order(in_matrix_A) * sizeof(float));
        printf("%s\n", prompt);
        fflush(stdout);
        for (mat_row = 0; mat_row < n; mat_row++) {
            grid_row = mat_row / Order(in_matrix_A);
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid -> g_order; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid -> g_comm, coords, & dest);
                if (dest == 0) {
                    for (mat_col = 0; mat_col < Order(in_matrix_A); mat_col++)
                        scanf("%f",
                            (in_matrix_A -> entries) + mat_row * Order(in_matrix_A) + mat_col);
                } else {
                    for (mat_col = 0; mat_col < Order(in_matrix_A); mat_col++)
                        scanf("%f", temp + mat_col);
                    MPI_Send(temp, Order(in_matrix_A), MPI_FLOAT, dest, 0,
                        grid -> g_comm);
                }
            }
        }
        free(temp);
    } else {
        for (mat_row = 0; mat_row < Order(in_matrix_A); mat_row++)
            MPI_Recv( & Entry(in_matrix_A, mat_row, 0), Order(in_matrix_A),
                MPI_FLOAT, 0, 0, grid -> g_comm, & status);
    }

}

void Print_matrix(
    char * title,
    TEMP_MATRIX * in_matrix_A,
    GRID_STRUCT * grid,
    int n) {
    int mat_row, mat_col;
    int grid_row, grid_col;
    int source;
    int coords[2];
    float * temp;
    MPI_Status status;

    if (grid -> current_rank == 0) {
        temp = (float * ) malloc(Order(in_matrix_A) * sizeof(float));
        printf("%s\n", title);
        for (mat_row = 0; mat_row < n; mat_row++) {
            grid_row = mat_row / Order(in_matrix_A);
            coords[0] = grid_row;
            for (grid_col = 0; grid_col < grid -> g_order; grid_col++) {
                coords[1] = grid_col;
                MPI_Cart_rank(grid -> g_comm, coords, & source);
                if (source == 0) {
                    for (mat_col = 0; mat_col < Order(in_matrix_A); mat_col++)
                        printf("%4.1f ", Entry(in_matrix_A, mat_row, mat_col));
                } else {
                    MPI_Recv(temp, Order(in_matrix_A), MPI_FLOAT, source, 0,
                        grid -> g_comm, & status);
                    for (mat_col = 0; mat_col < Order(in_matrix_A); mat_col++)
                        printf("%4.1f ", temp[mat_col]);
                }
            }
            printf("\n");
        }
        free(temp);
    } else {
        for (mat_row = 0; mat_row < Order(in_matrix_A); mat_row++)
            MPI_Send( & Entry(in_matrix_A, mat_row, 0), Order(in_matrix_A),
                MPI_FLOAT, 0, 0, grid -> g_comm);
    }

}
void Reset_matrix(
    TEMP_MATRIX * in_matrix_A) {

    int i, j;

    for (i = 0; i < Order(in_matrix_A); i++)
        for (j = 0; j < Order(in_matrix_A); j++)
            Entry(in_matrix_A, i, j) = 0.0;

}
void construct_m_type(
    TEMP_MATRIX * in_matrix_A) {
    MPI_Datatype temp_mpi_t;
    int block_lengths[2];
    MPI_Aint displacements[2];
    MPI_Datatype typelist[2];
    MPI_Aint start_address;
    MPI_Aint address;

    MPI_Type_contiguous(Order(in_matrix_A) * Order(in_matrix_A),
        MPI_FLOAT, & temp_mpi_t);
    block_lengths[0] = block_lengths[1] = 1;

    typelist[0] = MPI_INT;
    typelist[1] = temp_mpi_t;

    MPI_Address(in_matrix_A, & start_address);
    MPI_Address( & (in_matrix_A -> n_bar), & address);

    displacements[0] = address - start_address;

    MPI_Address(in_matrix_A -> entries, & address);
    displacements[1] = address - start_address;

    MPI_Type_struct(2, block_lengths, displacements,
        typelist, & local_matrix_mpi_t);
    MPI_Type_commit( & local_matrix_mpi_t);
}

void Temp_matrix_calculation(
    TEMP_MATRIX * in_matrix_A,
    TEMP_MATRIX * in_matrix_B,
    TEMP_MATRIX * out_matrix_C) {
    int i, j, k;

    for (i = 0; i < Order(in_matrix_A); i++)
        for (j = 0; j < Order(in_matrix_A); j++)
            for (k = 0; k < Order(in_matrix_B); k++)
                Entry(out_matrix_C, i, j) = Entry(out_matrix_C, i, j) +
                Entry(in_matrix_A, i, k) * Entry(in_matrix_B, k, j);

}
