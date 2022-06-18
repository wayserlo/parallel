#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <semaphore.h>
#include <time.h>
#include <semaphore.h>

sem_t mutex;

typedef struct Thread_Information {
    pthread_t thread_id;
    double sum;
    double x1, x2, y1, y2;
    int num;
    int N;
    double *S;
} ThreadInfo, *pThreadInfo;

int check(double x, double y) {//проверка, попадает ли в область
    return (x <= M_PI) && (x >= 0) && (y <= sin(x)) && (y >= 0);
}

void* thread_function(void *arg) {
    pThreadInfo thread = (pThreadInfo)arg;
    unsigned int x_k = thread->num, y_k = thread->num;
    int N = thread->N;
    double  x1 = thread->x1 ,
            x2 = thread->x2 ,
            y1 = thread->y1 ,
            y2 = thread->y2 ;

    int i = 0;
    double sum = 0;
    for (i = 0; i < N; i++) {
         double
                a = (double)rand_r(&x_k) / RAND_MAX * (x2 - x1) + x1 ,
                b = (double)rand_r(&y_k) / RAND_MAX * (y2 - y1) + y1 ;
                //printf("%lf, %lf\n", a, b);
        if ( check(a, b) == 1)
            sum += a * b;
    }
    sem_wait(&mutex);
    *(thread->S) += sum;
    sem_post(&mutex);
    thread->sum = sum;
    return (void *)thread;
}
int main(int argc, char *argv[]) {
    double x1 = 0;
    double x2 = M_PI;
    double y1 = 0;
    double y2 = 1;
    int N = 1e6;
    double S;
    int i;

    int p = atoi(argv[1]);

    sem_init(&mutex, 0, 1);
    pThreadInfo *Thread_Info = (pThreadInfo *)malloc(sizeof(pThreadInfo) * p);//массив под указатели

    for (i = 0; i < p; ++i) {
        pThreadInfo thread = (pThreadInfo)malloc(sizeof(ThreadInfo));//на каждуб нить
        Thread_Info[i] = thread;
    }
    struct timespec begin, end;

    clock_gettime(CLOCK_REALTIME, &begin);
    srand(time(NULL));
    for (i = 0; i < p; i++) {
        pThreadInfo thread = Thread_Info[i];
        thread->N = N/p;
        if ( N%p > 0 && i < N%p )
                thread->N += 1;
        thread->x1 = x1;
        thread->x2 = x2;
        thread->y1 = y1;
        thread->y2 = y2;
        thread->num = i;
        thread->S = &S;
        int rc = pthread_create(&(Thread_Info[i]->thread_id), NULL, thread_function, Thread_Info[i]);
    }


    for (i = 0; i < p; ++i){
        pthread_join(Thread_Info[i]->thread_id, NULL);
    }

    clock_gettime(CLOCK_REALTIME, &end);
    double parallel_time = end.tv_sec - begin.tv_sec;
    parallel_time += (end.tv_nsec - begin.tv_nsec) / (double)1e9;

    sem_destroy(&mutex);
    for (i = 0; i < p; ++i) {
        free(Thread_Info[i]);
    }

    printf("Numerical ans: %f \n", S * (x2 - x1) * (y2 - y1) / N);
    printf("Analytical ans: %f\n", M_PI * M_PI / 8);
    printf("time = %f for p = %d\n", parallel_time, p);
    return 0;
}
