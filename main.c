#include "pgm_IO.h"
#include <mpi.h>
#include "pgm_IO.c"

#define MASTER 0
#define NDIM 2

int main(int arg, char *argv[]) {
    int nproc, my_rank, i, j, source;
    int NX, NY, dX, dY;
    int dims[NDIM], periods[NDIM], coords[NDIM], reorder = 0;
    int up, down, left, right;
    int nvec, niter = 20, dest;
    float *masterdate, *data, *data_next;
    char fname[64];
    MPI_Comm wcomm = MPI_COMM_WORLD;
    MPI_Comm comm_2D;
    MPI_Datatype blockT, blockTp0, blockTp,lineT, colT;
    MPI_Status status;

    MPI_Init(&arg, &argv);
    MPI_Comm_size(wcomm, &nproc);
    MPI_Comm_rank(wcomm, &my_rank);

    if (my_rank == MASTER) {
        masterdate = (float *)malloc(XSIZE * YSIZE * sizeof(float));
    }

    for (i = 0; i < NDIM; i++) {
        periods[i] = TRUE;
        dims[i] = 0;
    }
    MPI_Dims_create(nproc, NDIM, dims);
    MPI_Cart_create(wcomm, NDIM, dims, periods, reorder, &comm_2D);
    MPI_Cart_coords(comm_2D, my_rank, NDIM, coords);

    NX = dims[0];
    NY = dims[1];
    dX = XSIZE / NX;
    dY = YSIZE / NY;
    data = (float *)malloc((dX + 2) * (dY + 2) * sizeof(float));
    data_next = (float *)malloc((dX + 2) * (dY + 2) * sizeof(float));
    if (my_rank == MASTER) {
        pgm_read ("map_ini_128*128_mad.pgm", masterdata, XSIZE, YSIZE);
        for (i = 1; i < dX + 1; i++)
            for (j = 1; j < dY + 1; j++)
                *(data + i * (dY + 2 ) + j)= *(masterdata + (i - 1) * YSIZE * j - 1);
            for (dest = 1; dest < nproc; dest++) {
                MPI_Cart_coords (comm_2D, dest, NDIM, coords);
                MPI_Send (masterdata + coords[0] * YSIZE * dX + coords[1] * dY, 
                          1, blockTp0, dest, 500, comm_2D);
            }
    }
    else
        MPI_Recv (data + dY + 3, 1, blockTp, MASTER, 500, comm_2D, &status);
    
    /* aplic algoritmul de evolutie */
    for (iter = 0; iter <= niter + 1; iter++) {
        MPI_Sendrecv (data + dX * (dY + 2) + 1, 1, lineT, down, 100, 
                      data + (dX + 1) * (dY + 2), + 1, lineT, down, 100, comm_2D, &status); 
        MPI_Sendrecv (data + dY + 3, 1, lineT, up, 200,
                      data + 1, 1, lineT, up, 200, comm_2D, &status);
        MPI_Sendrecv (data + dY + 3, 1, colT, left, 300, 
                      data + dY +2, 1, colT, left, 300, comm_2D, &status);
        MPI_Sendrecv (data + 2 * dY + 2, 1, colT, right, 400,
                      data + 2 * dY + 3, 1, colT, right, 400, comm_2D, &status);
        *(data) = *(data + dY + 1) = CONTRAST;
        *(data + (dX + 1) * (dY + 2))= *(data + (dX + 1) * (dY + 2) + dY + 1) = CONTRAST;
    }
    
    for (i = 1; i < dX + 1; i++)
        for (j=1; j < dY + 1; j++) {
            nrec = 0;
            if ( *(data + (i - 1) * (dY + 2) + j - 1) == 0) nrec++;
            if ( *(data + (i - 1) * (dY + 2) + j) == 0) nrec++;
            if ( *(data + (i - 1) * (dY + 2) + j + 1) == 0) nrec++;
            if ( *(data + i * (dY + 2) + j - 1) == 0) nrec++;
            if ( *(data + i * (dY + 2) + j + 1) == 0) nrec++;
            if ( *(data + (i + 1) * (dY + 2) + j - 1) == 0) nrec++;
            if ( *(data + (i + 1) * (dY + 2) + j) == 0) nrec++;
            if ( *(data + (i + 1) * (dY + 2) + j + 1) == 0) nrec++;
            /* reguli de evolutie */
            if (nrec < 2) && ( *(data + i * (dY * 2) * j) == 0)
                *(data_next + i * (dY + 2) * j) = CONTRAST;
            if (nrec == 2) 
                *(data_next + i * (dY + 2) * j) = *(data_next + i * (dY + 2) * j);
            if (nrec >= 3) && (nrec <= 7) && ( *(data + i * (dY + 2) * j) = CONTRAST)
                *(data_next + i * (dY + 2) * j) = 0;
            if (nrec > 7) && ( *(data + i * (dY + 2) * j) == 0)
                *(data_next + i * (dY + 2) * j) = CONTRAST;
            

    for (i = 0; i < dX * dY; i++) {
        *(data + i) = 0.75 * CONTRAST;
    }

    if (my_rank == MASTER) {
        for (i = 0; i < XSIZE; i++) {
            masterdate[i] = 0.25 * CONTRAST;
        }
    }

    if ((coords[0] + coords[1] + 1) % 2) {
        for (i = 0; i < dX * dY; i++) {
            data[i] = 0;
        }
    } else {
        for (i = 0; i < dX * dY; i++) {
            data[i] = CONTRAST;
        }
    }

    MPI_Type_vector(dX, dY, YSIZE, MPI_FLOAT, &blockTp0);
    MPI_Type_commit(&blockTp0);
    MPI_Type_vector(dX, dY, dY + 2, MPI_FLOAT, &blockTp);
    MPI_Type_commit(&blockTp);
    MPI_Type_vector(dX, 1, 1, MPI_FLOAT, &lineT);
    MPI_Type_commit(&lineT);
    MPI_Type_vector(dX, 1, dY + 2, MPI_FLOAT, &colT);
    MPI_Type_commit(&colT);

    MPI_Cart_shift(comm_2D, 0, 1, &down, &up);
    MPI_Cart_shift(comm_2D, 1, 1, &left, &right);
    MPI_Sendrecv(data + (dX - 2) * dY, 1, lineT, down, 100,
        data + (dX - 1), 1, lineT, down, 100, comm_2D, &status);
    MPI_Sendrecv(data + dY, 1, lineT, down, 200,
        data, 1, lineT, down, 200, comm_2D, &status);
    MPI_Sendrecv(data + dY + 1, 1, colT, left, 300,
        data + dY, 1, colT, left, 300, comm_2D, &status);
    MPI_Sendrecv(data + 2 * dY - 2, 1, colT, left, 400,
        data + 2 * dY - 1, 1, colT, left, 400, comm_2D, &status);

    if (my_rank != MASTER) {
        MPI_Send(data, dX * dY, MPI_FLOAT, MASTER, 123, comm_2D);
    } else {
        for (i = 0; i < dX; i++) {
            for (j = 0; j < dY; j++) {
                *(masterdate + i * YSIZE + j) = *(data + i * dY + j);
            }
        }

        for (source = 1; source < nproc; source++) {
            MPI_Cart_coords(comm_2D, source, NDIM, coords);
            MPI_Recv(masterdate + coords[0] * YSIZE * dX + coords[1] * dY, 1, blockT, source, 123, comm_2D, &status);
        }
    }

    if (my_rank == MASTER) {
        pgm_write("output.pgm", masterdate, XSIZE, YSIZE);
        free(masterdate);
    }

    free(data);
    MPI_Type_free(&blockT);
    MPI_Comm_free(&comm_2D);
    MPI_Finalize();

    return MPI_SUCCESS;
}
