#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


void *print_letter(void *letter){
	char l = (char) letter;
	printf("%c",l);
	pthread_exit(NULL);
}

int main (int argc, char *argv[]){
	if (argc != 2){
		printf("wrong number of parameters");
		exit(1);
	}
	char *c = argv[1];
	int count;
	
	for (count=0;c[count] !='\0';++count);
	
	pthread_t threads[count];

	int rc;
	int i;

	for (i=0;i<count;i++)
		rc = pthread_create(&threads[i],NULL, print_letter, (void *)c[i]);
		if(rc) {
			printf("deu erro na creacao da thread");
			exit(1);
		}
	
	for(i=0; i<count;i++){
		pthread_join(threads[i],NULL);
	}
	return(0);
}

