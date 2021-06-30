#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
void strategy1(unsigned long long* a, int N, int Threads){
    long j = 0;
    unsigned long long* sum = (unsigned long long *)malloc(sizeof(unsigned long long)*((N+1)/2));
    unsigned long long ans = 0;
    double start_time = omp_get_wtime();
    int k = 0;
    while(N > 1){
        #pragma omp parallel num_threads(Threads) shared(sum, N, a) private(k)
        {

            #pragma omp for schedule(static)
            for(j=0; j < (((N+1)/2) - 1); j++){
                k = 2*j;
                sum[j] = a[k];
                k++;
                sum[j] += a[k];
            }
        }

        j = ((N+1)/2) - 1;
        int k = 2*j;
        sum[j] = a[k];
        if(k + 1 < N)
            sum[j] += a[k+1];
        
        unsigned long long* temp = a;
        a = sum;
        sum = temp;
        N = (N + 1) / 2;
    }

    ans = a[0];

    double time = omp_get_wtime() - start_time;
    // printf("time : %g\n", time);
    printf("%llu\n", ans);
    free(sum);
    free(a);
}
void strategy0(unsigned long long* a, int N, int Threads){

    double start_time = omp_get_wtime();
    long i = 0;
    unsigned long long sum = 0;
    unsigned long long localsum = 0;
    #pragma omp parallel num_threads(Threads) shared(sum) private(localsum)
    {
        localsum = 0;

        #pragma omp for nowait
        for(i=0; i < N; i++){
            localsum += a[i];
        }

        #pragma omp critical
        {
            sum += localsum;
        }

    } 
    double time = omp_get_wtime() - start_time;
    // printf("time : %g\n", time);
    printf("%llu\n", sum);
    free(a);
}
int main(int argc, char* argv[]){ 
    int strategy = atoi(argv[1]);
    int N = atoi(argv[2]);
    int Threads = atoi(argv[3]); 
    unsigned long long* a;
   
    a = (unsigned long long *)malloc(sizeof(unsigned long long)*N);
    a[0] = 1;

    for(int i=1;i< N;i++) 
        a[i] = i+1;
    if(strategy == 0){
        strategy0(a, N, Threads);
    }
    else{
        strategy1(a, N, Threads);   
    }

}