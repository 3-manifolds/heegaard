#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <termios.h>

/* #define DEBUG_MALLOC	*/

extern struct termios normal_termio;

void   SysBeep(int seconds);
char   **NewHandle(size_t);
void   DisposeHandle(char **);
size_t GetHandleSize(char **);
void   *NewPtr(size_t);
void   DisposePrt(void *);
size_t GetPtrSize(void *);
void   ReadDateTime(unsigned long *);
char   *ReadString(char *buffer, int size);

typedef struct Handle{
  void *data;
  size_t size;
} Handle_t;

typedef struct Pointer{
  size_t size;
} Pointer_t;
  
void SysBeep(int seconds){
	printf("\a");
}

#ifndef DEBUG_MALLOC
char **NewHandle(size_t size){
  Handle_t *handle;
  
  handle = (Handle_t *)malloc(sizeof(Handle_t));
  if (handle == NULL)
    return NULL;
  handle->data = malloc(size);
  handle->size = size;
  return (char **)handle;
}

void DisposeHandle(char **h){
  Handle_t *handle = (Handle_t *)h;
  if (handle == NULL)
    return;
  if (handle->data)
    free(handle->data);
  free(handle);
}

size_t GetHandleSize(char **h){
  Handle_t *handle = (Handle_t *)h;
  return handle->size;
}

void *NewPtr(size_t size){
  Pointer_t *pointer; 
  pointer = (Pointer_t *)malloc(size + sizeof(Pointer_t));
  if (pointer == NULL)
    return NULL;
  pointer->size = size;
  return (void *)pointer + sizeof(Pointer_t);
}

void DisposePtr(void *p){
  Pointer_t *pointer;
  if (p == NULL)
    return;
  pointer = (Pointer_t *)(p - sizeof(Pointer_t));
  free(pointer);
}

size_t GetPtrSize(void *p){
  Pointer_t *pointer = (Pointer_t *)(p - sizeof(Pointer_t));
  return pointer->size;
}
#endif

#ifdef DEBUG_MALLOC

static const char   test_padding[5] = "[]{}";
unsigned long		NumErrors;

char **NewHandle(size_t size){
  Handle_t *handle;
  char		*test_bytes;
  int		i;
  
  handle = (Handle_t *)malloc(sizeof(Handle_t));
  if (handle == NULL)
    return NULL;
  handle->data = malloc(size + 4);
  if(handle->data == NULL)
  	{
  	free(handle);
  	return NULL;
  	}
  handle->size = size;  
  test_bytes = (char *)handle->data + size;
  for(i = 0; i < 4; i++)
  	test_bytes[i] = test_padding[i];  
  return (char **)handle;
}

void DisposeHandle(char **h){
  char	*p,
  		*test_bytes;
  
  int 	i;
  
  if (h == NULL) return;
  Handle_t *handle = (Handle_t *)h;
  if (handle->data)
  	{
  	test_bytes = (char *)handle->data + handle->size;
  	for(i = 0; i < 4; i++) if(test_bytes[i] != test_padding[i])
  		{
   		NumErrors ++; 			
  		printf("\n\nError in DisposeHandle().\n");	
  		p = (char *)handle->data;
  		while(p < test_bytes + 4) putchar(*p++);
  		break;
  		}	
    free(handle->data);
    }
  free(handle);
}

size_t GetHandleSize(char **h){
  Handle_t *handle = (Handle_t *)h;
  return handle->size;
}
	
void *NewPtr(size_t size){
  char		*test_bytes;
  int		i;
  
  Pointer_t *pointer;  
  pointer = (Pointer_t *)malloc(size + sizeof(Pointer_t) + 4);
  if (pointer == NULL)
    return NULL;
  pointer->size = size;  
  test_bytes = (char *)pointer + sizeof(Pointer_t) + size;
  for(i = 0; i < 4; i++)
  	test_bytes[i] = test_padding[i];    
  return (void *)pointer + sizeof(Pointer_t);
}

void DisposePtr(void *p){
  char	*q,
  		*test_bytes;
  
  int 	i;
  		
  if (p == NULL) return;
  Pointer_t *pointer;
  pointer = (Pointer_t *)(p - sizeof(Pointer_t));  
  test_bytes = (char *)p + pointer->size;
  for(i = 0; i < 4; i++) if(test_bytes[i] != test_padding[i])
	{
	NumErrors ++; 			
	printf("\n\nError in DisposePtr().\n");	
	q = (char *)p;
	while(q < test_bytes + 4) putchar(*q++);
	break;
	} 
  free(pointer);
}

size_t GetPtrSize(void *p){
  Pointer_t *pointer = (Pointer_t *)(p - sizeof(Pointer_t));
  return pointer->size;
}
#endif

void Mem_Error(void)
{
printf("\n\nUnable to allocate memory. Heegaard will now quit.\n");
exit(0);
}

void   ReadDateTime(unsigned long *seconds){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  *seconds = tv.tv_sec;
}

char *ReadString(char *buffer, int size){
  struct termios this_termio;
  char *result;			
  
  tcgetattr(0, &this_termio);
  tcsetattr(0, 0, &normal_termio);
  result = fgets(buffer, size, stdin); 
  while (*buffer){
    if (*buffer == '\n') 
    	{
    	*buffer='\0';
    	break;
    	}
    buffer++;
  }
  tcsetattr(0, 0, &this_termio);
  return result;
}

/* 
main(){
  void **p;
  unsigned long time;

  ReadDateTime(&time);
  printf("%u seconds\n\r", time);
  printf("handle struct size is %d\n\r", sizeof(Handle_t));
  SysBeep(5);
  p = NewHandle(100);
  printf("Allocated %d bytes\n\r", GetHandleSize(p));
  ReallocateHandle(p, 200);
  printf("reallocated %d bytes\n\r", GetHandleSize(p));
  DisposeHandle(p);
  DisposeHandle(p);
  p = NewPtr(100);
  printf("Allocated %d bytes\n\r", GetPtrSize(p));
  DisposePtr(p);
}
*/
