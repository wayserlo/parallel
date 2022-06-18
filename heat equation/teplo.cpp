#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

double iter(double a, double b, double c, double k, double dt, double h){//одна итерация, возвращает u_{i+1}
        return k*dt*(c - 2*b + a)/h/h + b;
}
void exact_sol(double T, double u0, double l, double k){//возвращает точное решение
        printf("exact solution on N = :\n");
        double ex_sol, i = 0;
        int m = 0;
        for(i; i < 11; i++){
                ex_sol = 0;
                for(m = 0; m < 3; m++){
                        ex_sol += 4*u0/M_PI*exp(-k*M_PI*M_PI*(2*m+1)*(2*m+1)*T/l/l)/(2*m+1)*sin(M_PI*(2*m+1)*(i/10)/l);
                }
                printf("u(%lf, %lf) = %lf\n", i/10, T, ex_sol);
        }
}

int main(int argc, char*argv[]){
        double u0 = 1, l = 1, k = 1;
        int N = atoi(argv[1]) + 1;
        double h = l/(N-1);
        double dt = h*h*0.5;
        double T = 500*dt;

        int myrank, size;
        double begin, end;
        MPI_Status Status;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        if(myrank == 0){
                printf("Solution on N = %d, p = %d, T = %lf\n", N, size, T);
                exact_sol(T, u0, l, k);
                printf("\n\n\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        begin = MPI_Wtime();
        double *layer, *next_layer;
        int i = 0, j = 0;
        int length = N/size;
        for(i = 0; i < N%size; i++){
                if(myrank == i){
                        length += 1;
                }
        }
        //if(myrank == size - 1){
        //      length += N%size;
        //}
        layer = (double*)malloc((length + 2) * sizeof(double));
        next_layer = (double*)malloc((length + 2) * sizeof(double));

        for(i = 0; i < length; i++){
                layer[i] = u0;
                next_layer[i] = u0;
        }
        if(myrank==0){
                layer[1] = 0;
                next_layer[1] = 0;
        }
        if(myrank==size-1){
                layer[length] = 0;
                next_layer[length] = 0;
        }
        printf("t = %d\n", (int) (T/dt));
        for(j=0; j <= (int) (T/dt) + 1; j++){
                //printf("j = %d\n", j);
                if(myrank%2 == 0 && myrank < size - 1){//сначала четные, кроме последнего, получают правую точку и отправляют свою правую  
                        MPI_Recv(layer+length +1, 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
                        MPI_Send(layer+length, 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
                }
                if(myrank%2 == 0 && myrank > 0){//четны, кроме первого, получают левую точку и отправляют левую
                        MPI_Recv(layer, 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
                        MPI_Send(layer+1, 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
                }
                if(myrank%2 == 1){//нечетные сначала отправляют левую точку и потом получают левую
                        MPI_Send(layer+1, 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
                        MPI_Recv(layer, 1, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD, &Status);
                }
                if(myrank%2 == 1 && myrank < size -1){//нечетные, кроме последнего, сначала отправляют правую точку и потом получают правую      
                        MPI_Send(layer+length, 1, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD);
                        MPI_Recv(layer+length+1, 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
                }
                /*for (myrank =0; myrank < size; myrank++){
                        if(myrank%2 == 1 && myrank < size -1){
                                MPI_Send(layer+1, 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
                                MPI_Recv(layer+length+1, 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
                        }
                        if(myrank%2 == 0 && myrank > 0){
                                MPI_Recv(layer+length+1, 1, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD, &Status);
                                MPI_Send(layer+1, 1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD);
                        }*/
                int left = (myrank == 0) ? 2 : 1;
                int right = (myrank == size - 1) ? length : length + 1;
                int n;
                for (n = left; n < right; n++) {
                        next_layer[n] = iter(layer[n-1], layer[n], layer[n+1], k, dt, h);
                        //printf("next(n) = %lf\n", next_layer[n]);
                }
                for(n = left; n < right; n++){
                        layer[n] = next_layer[n];
                }
        }
        int n;
        double  x;
        x = myrank * (N / size) * h;
        int left = 1, right = length + 1;
        for (n = left; n < right; n++){
                if (fabs((x * 10) - round(x * 10)) < 5 * h) {
                        printf("u(%lf, %lf) = %lf\n", x, T, layer[n]);
                }
                x += h;
        }
        free(layer);
        free(next_layer);

        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        if(myrank==0){
                printf("time = %lf\n", end - begin);
        }
        MPI_Finalize();
        return 0;
}
