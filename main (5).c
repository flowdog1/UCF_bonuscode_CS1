#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
  extraMemoryAllocated += sz;
  size_t* ret = malloc(sizeof(size_t) + sz);
  *ret = sz;
  printf("Extra memory allocated, size: %ld\n", sz);
  return &ret[1];
}

void DeAlloc(void* ptr)
{
  size_t* pSz = (size_t*)ptr - 1;
  extraMemoryAllocated -= *pSz;
  printf("Extra memory deallocated, size: %ld\n", *pSz);
  free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
  return ((size_t*)ptr)[-1];
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapify(int arr[], int n, int i) {
    int largest = i;
    int left = 2*i + 1;
    int right = 2*i + 2; 

    if (left < n && arr[left] > arr[largest])
        largest = left;

    if (right < n && arr[right] > arr[largest])
        largest = right;

    if (largest != i) {
        int swap = arr[i];
        arr[i] = arr[largest];
        arr[largest] = swap;

        heapify(arr, n, largest);
    }
}

void heapSort(int arr[], int n) {
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i=n-1; i>0; i--) {
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        heapify(arr, i, 0);
    }
}
//End of heap sort

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        // Sort first & second halves
        mergeSort(pData, l, m);
        mergeSort(pData, m + 1, r);

        // Merge sorted halves
        int i, j, k;
        int n1 = m - l + 1;
        int n2 = r - m;

        // Allocate temp arrays
        int *L = (int *)Alloc(sizeof(int)*n1), *R = (int *)Alloc(sizeof(int)*n2);

        // Copy data to temp arrays L[] and R[]
        for (i = 0; i < n1; i++)
            L[i] = pData[l + i];
        for (j = 0; j < n2; j++)
            R[j] = pData[m + 1+ j];

        // Merge back into pData[l..r]
        i = 0; 
        j = 0; 
        k = l; 
        while (i < n1 && j < n2) {
            if (L[i] <= R[j]) {
                pData[k] = L[i];
                i++;
            }
            else {
                pData[k] = R[j];
                j++;
            }
            k++;
        }

        // Copy remaining elements of L[]
        while (i < n1) {
            pData[k] = L[i];
            i++;
            k++;
        }

        // Copy remaining elements of R[]
        while (j < n2) {
            pData[k] = R[j];
            j++;
            k++;
        }

        // Free temp arrays
        DeAlloc(L);
        DeAlloc(R);
    }
}

// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n) {
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = pData[i];
        j = i - 1;

        while (j >= 0 && pData[j] > key) {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        pData[j + 1] = key;
    }
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n) {
    int i, j, tmp;
    for (i = 0; i < n-1; i++)    
        for (j = 0; j < n-i-1; j++)
            if (pData[j] > pData[j+1]) {
                tmp = pData[j];
                pData[j] = pData[j+1];
                pData[j+1] = tmp;
            }
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n) {
    int i, j, minIdx, tmp;
    for (i = 0; i < n-1; i++) {
        minIdx = i;
        for (j = i+1; j < n; j++)
          if (pData[j] < pData[minIdx])
            minIdx = j;
        if (minIdx != i) {
            tmp = pData[i];
            pData[i] = pData[minIdx];
            pData[minIdx] = tmp;
        }
    }
}

// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
  FILE* inFile = fopen(inputFileName,"r");
  int dataSz = 0;
  *ppData = NULL;

  if (inFile)
  {
    fscanf(inFile,"%d\n",&dataSz);
    *ppData = (int *)Alloc(sizeof(int) * dataSz);
  for (int i = 0; i < dataSz; ++i) {
          fscanf(inFile, "%d", (*ppData) + i);
      }
      fclose(inFile);
  }

  return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
  int i, sz = dataSz - 100;
  printf("\tData:\n\t");
  for (i=0;i<100;++i)
  {
    printf("%d ",pData[i]);
  }
  printf("\n\t");

  for (i=sz;i<dataSz;++i)
  {
    printf("%d ",pData[i]);
  }
  printf("\n\n");
}

int main(void)
{
  clock_t start, end;
  int i;
    double cpu_time_used;
  char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};

  for (i=0;i<3;++i)
  {
    int *pDataSrc, *pDataCopy;
    int dataSz = parseData(fileNames[i], &pDataSrc);

    if (dataSz <= 0)
      continue;

    pDataCopy = (int *)Alloc(sizeof(int)*dataSz);

    printf("---------------------------\n");
    printf("Dataset Size : %d\n",dataSz);
    printf("---------------------------\n");

    printf("Selection Sort:\n");
    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
    extraMemoryAllocated = 0;
    start = clock();
    selectionSort(pDataCopy, dataSz);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
    printArray(pDataCopy, dataSz);

    printf("Insertion Sort:\n");
    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
    extraMemoryAllocated = 0;
    start = clock();
    insertionSort(pDataCopy, dataSz);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
    printArray(pDataCopy, dataSz);

    printf("Bubble Sort:\n");
    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
    extraMemoryAllocated = 0;
    start = clock();
    bubbleSort(pDataCopy, dataSz);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
    printArray(pDataCopy, dataSz);

    printf("Merge Sort:\n");
    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
    extraMemoryAllocated = 0;
    start = clock();
    mergeSort(pDataCopy, 0, dataSz - 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
    printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
    memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
    extraMemoryAllocated = 0;
    start = clock();
    heapSort(pDataCopy, dataSz); //Original line: heapSort(pDataCopy, 0, dataSz - 1); - returns error for too many arguments. Should have two but it has 3. Assuming that it just got mistakenly copied from merge sort and set the standard call.
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
    printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
    printArray(pDataCopy, dataSz);

    DeAlloc(pDataCopy);
    DeAlloc(pDataSrc);
  }

}