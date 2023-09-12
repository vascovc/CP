
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>



#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

#define HI(num)	(((num) & 0x0000FF00) << 8)
#define LO(num)	((num) & 0x000000FF)


// Global Var
unsigned int* h_idata=NULL;
unsigned int* h_odata=NULL;
unsigned int h, w;

float *filter;
unsigned int fh, fw;
int NUM_THREADS = 4; 

int loadPGM(char* fname, unsigned int** parray, unsigned int *pwidth, unsigned int *pheight, unsigned int *pmaxgray)
{
    FILE *pgmFile;
    char version[3];
    int i, j;
    int lo, hi;
    pgmFile = fopen(fname, "rb");
    if (pgmFile == NULL) {
        perror("cannot open file to read");
        return -1;
    }
    if(fgets(version, sizeof(version), pgmFile)==NULL) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if (strcmp(version, "P5")) {
        fprintf(stderr, "File format is not P5 pgm\n");
        return -1;
    }
    if(fscanf(pgmFile, "%d", pwidth)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if(fscanf(pgmFile, "%d", pheight)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    if(fscanf(pgmFile, "%d", pmaxgray)!=1) {
        fprintf(stderr, "error reading file\n");
        return -1;
    }
    fgetc(pgmFile);
 
    *parray = malloc( *pwidth * *pheight * sizeof(int) );
    if (*pmaxgray > 255) {
        for (i = 0; i < *pheight; ++i) {
            for (j = 0; j < *pwidth; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                (*parray)[i * *pwidth + j] = (hi << 8) + lo;
            }
        }
    }
    else {
        for (i = 0; i < *pheight; ++i) {
            for (j = 0; j < *pwidth; ++j) {
                lo = fgetc(pgmFile);
                (*parray)[i * *pwidth + j] = lo;
            }
        }
    }
 
    fclose(pgmFile);
    return 0;
 
}
 

int savePGM(char* fname, unsigned int* array, unsigned int width, unsigned int height, unsigned int maxgray)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;
 
    pgmFile = fopen(fname, "wb");
    if (pgmFile == NULL) {
        perror("cannot open file to write");
        return -1;
    }
 
    fprintf(pgmFile, "P5\n");
    fprintf(pgmFile, "%d\n%d\n", width, height);
    fprintf(pgmFile, "%d\n", maxgray);
 
    if (maxgray > 255) {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                hi = HI(array[i*width+j]);
                lo = LO(array[i*width+j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }
 
        }
    }
    else {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                lo = LO(array[i*width+j]);
                fputc(lo, pgmFile);
            }
        }
    }
 
    fclose(pgmFile);

    return 0;
}


// loads filter coefficients from file fname, 
// allocates memory through parray and stores width and height of filter through pwidth and pheight
int loadFilter(char* fname, float** parray, unsigned int *pwidth, unsigned int *pheight)
{
    FILE* fp;

    if( (fp=fopen(fname, "r")) == NULL)
    {    
        fprintf(stderr,"Failed to open filter file %s\n",fname);
        return -1;
    }

    if(fscanf(fp,"%u %u",pwidth,pheight)!=2) {
        fprintf(stderr,"Failed to read header of filter file %s\n",fname);
        return -1;
    }

    *parray = (float *) malloc((*pwidth)*(*pheight)*sizeof(float));

    int i;
    for(i=0;i<(*pwidth)*(*pheight);i++)
    {
        if(fscanf(fp,"%f",(*parray+i))!=1) {
           fprintf(stderr,"Failed to read data of filter file %s\n",fname);
           return -1;
        }
    }

    fclose(fp);
    
    return 0;
}


// filter image
void filterImage(unsigned int *h_idata, unsigned int w, unsigned int h, 
                float* filter, unsigned int fw, unsigned int fh, 
                unsigned int* h_odata)
{
    int i,j,k,l;

    int fh_2 = fh/2;
    int fw_2 = fw/2;

    for(i=0; i<h; i++) //height image
    {
        for(j=0; j<w; j++) //width image
        {
            float sum = 0;
            for(k=-fh_2; k<=fh_2; k++) //filter height
            {
                for(l=-fw_2; l<=fw_2; l++) //filter width
                {
                    if( (i+k >= 0) && (i+k < h))
                        if( (j+l >=0) && (j+l < w)) {
                            sum += h_idata[(i+k)*w + j+l]*filter[(k+fh/2)*fw + l+fw/2];         
                        }

                }
            }
            h_odata[i*w+j] = min(max(sum,0),255);
        }
    }
}   

// filter image w
void filterImage_w(unsigned int *h_idata, unsigned int w, unsigned int h,
                float* filter, unsigned int fw, unsigned int fh,
                unsigned int* h_odata, unsigned int w_start, unsigned int w_end)
{
    int i,j,k,l;

    int fh_2 = fh/2;
    int fw_2 = fw/2;

    for(i=0; i<h; i++) //height image
    {
        for(j=w_start; j<w_end; j++) //width image
        {
            float sum = 0;
            for(k=-fh_2; k<=fh_2; k++) //filter height
            {
                for(l=-fw_2; l<=fw_2; l++) //filter width
                {
                    if( (i+k >= 0) && (i+k < h))
                        if( (j+l >=0) && (j+l < w)) {
                            sum += h_idata[(i+k)*w + j+l]*filter[(k+fh/2)*fw + l+fw/2];
                        }

                }
            }
            h_odata[i*w+j] = min(max(sum,0),255);
        }
    }
}
// filter image h 
void filterImage_h(unsigned int *h_idata, unsigned int w, unsigned int h, 
                float* filter, unsigned int fw, unsigned int fh, 
                unsigned int* h_odata, unsigned int h_start, unsigned int h_end)
{
    int i,j,k,l;

    int fh_2 = fh/2;
    int fw_2 = fw/2;

    for(i=h_start; i<h_end; i++) //height image
    {
        for(j=0; j<w; j++) //width image
        {
            float sum = 0;
            for(k=-fh_2; k<=fh_2; k++) //filter height
            {
                for(l=-fw_2; l<=fw_2; l++) //filter width
                {
                    if( (i+k >= 0) && (i+k < h))
                        if( (j+l >=0) && (j+l < w)) {
                            sum += h_idata[(i+k)*w + j+l]*filter[(k+fh/2)*fw + l+fw/2];         
                        }

                }
            }
            h_odata[i*w+j] = min(max(sum,0),255);
        }
    }
}

// filter image thread -- by heigth
void * filterImage_th_h(void *arg)
{
    long t = (long) arg;
    int h_start = ceil(t*h/NUM_THREADS);
    int h_end = ceil((t + 1)*h/NUM_THREADS);
    filterImage_h(h_idata, w, h, filter, fw, fh, h_odata,h_start,h_end); 
    return NULL;
}

// filter image thread -- by width
void * filterImage_th_w(void *arg)
{
    long t = (long) arg;
    int w_start = ceil(t*w/NUM_THREADS);
    int w_end = ceil((t+1)*w/NUM_THREADS);
    filterImage_w(h_idata,w,h,filter,fw,fh,h_odata,w_start,w_end);
    return NULL;
}

// print command line format
void usage(char *command) 
{
    printf("Usage: %s [-h] [-d device] [-i inputfile] [-o outputfile] [-f filterfile]\n",command);
}

// main
int main( int argc, char** argv) 
{
    unsigned int maxgray; 

    
    // default command line options
    char *fileIn="lena.pgm",*fileOut="lenaOut.pgm",*fileFilter="filter.txt";

    // parse command line arguments
    int opt;
    while( (opt = getopt(argc,argv,"i:o:f:h:t:")) !=-1)
    {
        switch(opt)
        {
            case 'i':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }

                fileIn = strdup(optarg);
                break;
            case 'o':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                fileOut = strdup(optarg);
                break;
            case 'f':
                if(strlen(optarg)==0)
                {
                    usage(argv[0]);
                    exit(1);
                }
                fileFilter = strdup(optarg);
                break;
            case 'h':
                usage(argv[0]);
                exit(0);
                break;
	    case 't':
		if (strlen(optarg)==0)
		{
		    usage(argv[0]);
		    exit(1);
		}
		NUM_THREADS = atoi(optarg);
		break;

        }
    }

    //load pgm
    if (loadPGM(fileIn, &h_idata, &w, &h,&maxgray)==-1) {
        printf("Failed to load image file: %s\n", fileIn);
        exit(1);
    }

    //load filter
    if(loadFilter(fileFilter, &filter, &fw, &fh)==-1)
    {
       printf("Failed to load filter file: %s\n",fileFilter);
       exit(1);
    }

    // allocate mem for the result
    h_odata   = (unsigned int*) malloc( h*w*sizeof(unsigned int));
 
    struct timeval start, end;
    
   
    double speedup_total = 0;
    int n_filters = 100;
    for (int k = 0;k < n_filters;k++)
    {
    printf("Serial\n");
    gettimeofday(&start, NULL);

    // filter image serial 
    filterImage(h_idata, w, h, filter, fw, fh, h_odata);

    gettimeofday(&end, NULL);
    printf( "Processing time: %f (ms)\n", (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0);

    double time_serial = (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    fprintf(stderr, "NTH %d\n",NUM_THREADS);
    gettimeofday(&start, NULL);

    // filter image in NUM_THREAD threads
    pthread_t thread[NUM_THREADS];

    for (long t = 0; t < NUM_THREADS; t++)
    {
        // Choose one
	    pthread_create(&thread[t], NULL, filterImage_th_w, (void*)t);		
        //pthread_create(&thread[t], NULL, filterImage_th_h, (void*)t);		
    }    
    
    for (int i = 0; i<NUM_THREADS; i++)
    {
	    pthread_join(thread[i],NULL);
    }

    gettimeofday(&end, NULL);
    printf( "Processing time: %f (ms)\n", (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0);

    double time_parallel = (end.tv_sec-start.tv_sec)*1000.0 + ((double)(end.tv_usec - start.tv_usec))/1000.0;
    double speedup = time_serial/time_parallel;

    printf("SpeedUp = %f \n",speedup);
    speedup_total = speedup_total + speedup;
    }
    printf("Average SpeedUp = %f \n",speedup_total/n_filters);
    // save output image
    if (savePGM(fileOut, h_odata, w, h, maxgray)==-1) {
        printf("Failed to save image file: %s\n", fileOut);
        exit(1);
    }

    // cleanup memory
    free(h_idata);
    free(h_odata);
    free(filter);
}
