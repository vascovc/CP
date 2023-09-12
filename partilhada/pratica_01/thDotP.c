#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define NUM_THREADS 2
#define num 2000

struct vector{
	float *op1;
	float *op2;
	int n;
	float result;
};

void* Multiply(void *arr){
	struct vector v = *(struct vector*) arr;

	float partial_sum = 0;
	
	for(int i=0; i<v.n; i++){
		partial_sum += v.op1[i] * v.op2[i];
	}

	(*(struct vector*) arr).result = partial_sum;
	//printf("%f\n",v.result);
    printf("O resultado do 'dot product' = %f\n", partial_sum);
	float* resulatado_thread = (float *)malloc(sizeof(float));
	*resulatado_thread = partial_sum;
	pthread_exit(resulatado_thread);
}

int main(int argc, char *argv[]){
	
	float v_1[num];
	float v_2[num];
	
	for (int i=0;i<num;i++){
		v_1[i] = (float)(rand()%100);;
		v_2[i] = (float)(rand()%100);;
	}
	
	//struct vector vec_1;
	//vec_1.op1 = &v_1[0];
	//vec_1.op2 = &v_2[0];
	//vec_1.n = num/2;
	// na teoria o que esta em cima vai ser a mesma coisa
	//vec_1.op1 = v_1;
	//vec_1.op2 = v_2;
	//vec_1.n = num/2;
	//struct vector vec_2;
	//vec_2.op1 = v_1 + num/2;
	//vec_2.op2 = v_2 + num/2;
	//vec_2.n = num/2;

	pthread_t threads[NUM_THREADS];
	struct vector vec_s[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++){
		struct vector vec;
		vec.op1 = v_1 + i*num/NUM_THREADS;
		vec.op2 = v_2 + i*num/NUM_THREADS;
		vec.n = num/NUM_THREADS;
		vec_s[i] = vec;
	}

	int rc;
	for (int i=0; i<NUM_THREADS;i++){
		rc = pthread_create(&threads[i],NULL,Multiply,(void *)&vec_s[i]);
		if (rc) {
			printf("erro de threads");
			exit(1);
		}
	}

	float soma_elemento_thread = 0;
	for (int i=0; i<NUM_THREADS;i++){
		printf("resultado da thread %d: %f\n",i,vec_s[i].result);
		soma_elemento_thread += vec_s[i].result;
	}
	printf("resultado da soma, guadardando na struct: %f\n",soma_elemento_thread);

	float* thread_result;
	float soma_final = 0.0;
	for(int i = 0; i < NUM_THREADS; i++){

        // Função que se encarrega de assegurar que as threads esperam que as outras terminem
        // Passamos o endereço do ponteiro para que seja la escrito o resutlado retornado da thread
        pthread_join(threads[i], &thread_result);

        // Fazer a conta final
        soma_final += *thread_result;
        // Libertamos a memoria
        free(thread_result);

    }
	// Verificar a soma dos vetores
    float resultado_final = 0.0;
    for(int i = 0; i < num; i++){

        resultado_final += v_1[i] * v_2[i];

    }
	printf("\n\n");
    printf("Resultado final esperado: %f", resultado_final);
    printf("\n\n");
    // Printar o resultado final obtido
    printf("Resultado final obtido: %f", soma_final);
    printf("\n");
}

