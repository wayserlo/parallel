#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

double iter(double x, double y){//одна итерация
        return 0.5*(4/(1+x*x) + 4/(1+y*y))*(y-x);
}
double integ(double x, double  y, double h){//интегрирование на [x;y]
        double sum = 0;
        int i = 0;
        while(x+h < y){
                sum += iter(x, x+h);
                x += h;
        }
        sum += iter(x, y);
        return sum;
}

int main(int argc, char*argv[]){
        int i;
        int myrank, size;
        double begin, end;
        double I_serial = 0, I = 0, N = 10, x=0;
        double h = 1/N;
        double mes[2];
        MPI_Status Status;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        if(myrank == 0){
                //последовательный интеграл
                begin = MPI_Wtime();
                I_serial = integ(0.0, 1.0, h);
                end = MPI_Wtime();
                printf("I_serial = %lf, time = %lf\n", I_serial, end - begin);

                //параллельный интеграл
                begin = MPI_Wtime();
                x = 0;
                for(i=1; i < size; i++){//отправляем концы отрезков
                        mes[0] = x;
                        mes[1] = x + N*h/size;
                        x += N*h/size;
                        MPI_Send(&mes[0], 2, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
                }//считаем свой остаток интеграла
                I = integ(x, 1.0, h);
                printf("I_%d = %lf\n", myrank, I);

                for(i=1; i < size; i++){//принимаем посчитанные интегралы и суммируем
                        MPI_Recv(&mes[0], 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
                        I += mes[0];
                }
                end = MPI_Wtime();
                printf("I_parallel = %lf, time = %lf\n", I, end - begin);
                printf("I_serial - I_parallel = %lf\n", I_serial - I);
        }
        if(myrank != 0){
                MPI_Recv(&mes[0], 2, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &Status);
                mes[0] = integ(mes[0], mes[1], h);
                printf("I_%d = %lf\n", myrank, mes[0]);
                MPI_Send(&mes[0], 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
        }
        MPI_Finalize();
        return 0;
}
