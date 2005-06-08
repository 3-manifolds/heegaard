#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <termios.h>

extern struct termios normal_termio;

void   SysBeep(int seconds);
void   ObscureCursor(void);
void   **NewHandle(size_t);
void   DisposeHandle(void **);
void   ReallocateHandle(void **, size_t);
size_t GetHandleSize(void **);
void   SetHandleSize(void **, size_t);
void   EmptyHandle(void **);
void   HLock(void **);
void   HUnlock(void **);
void   *NewPtr(size_t);
void   DisposePrt(void *);
size_t GetPtrSize(void *);
void   ReadDateTime(unsigned long *);
char   *ReadString(char *buffer, int size);

typedef struct Handle{
  void *data;
  size_t size;
  int lock;
} Handle_t;

typedef struct Pointer{
  size_t size;
} Pointer_t;
  
void SysBeep(int seconds){
  printf("\007");
}

void ObscureCursor(void){
  return;
}

void **NewHandle(size_t size){
  Handle_t *handle;
  
  handle = (Handle_t *)malloc(sizeof(Handle_t));
  if (handle == NULL)
    return NULL;
  handle->data = malloc(size);
  handle->size = size;
  handle->lock = 0;
  return (void **)handle;
}

void ReallocateHandle(void **h, size_t size){
  Handle_t *handle = (Handle_t *)h;
  
  handle->data = realloc(handle->data, size);
  handle->size = size;
}

void DisposeHandle(void **h){
  Handle_t *handle = (Handle_t *)h;
  if (handle == NULL)
    return;
  if (handle->lock){
    fprintf(stderr, "Attempt to free locked handle!\n\r");
    return;
  }
  if (handle->data)
    free(handle->data);
  free(handle);
}

size_t GetHandleSize(void **h){
  Handle_t *handle = (Handle_t *)h;
  return handle->size;
}

void SetHandleSize(void **h, size_t size){
  Handle_t *handle = (Handle_t *)h;
  handle->data = realloc(handle->data, size);
  handle->size = size;
}

void EmptyHandle(void **h){
  SetHandleSize(h, 0);
}

void HLock(void **h){
  Handle_t *handle = (Handle_t *)h;
  handle->lock = 1;
}

void HUnlock(void **h){
  Handle_t *handle = (Handle_t *)h;
  handle->lock = 0;
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
    if (*buffer == '\n') *buffer='\0';
    buffer++;
  }
  tcsetattr(0, 0, &this_termio);
}

/*
main(){
  void **p;
  unsigned long time;

  ReadDateTime(&time);
  printf("%u seconds\n\r", time);
  printf("handle struct size is %d\n\r", sizeof(Handle_t));
  SysBeep(5);
  ObscureCursor();
  p = NewHandle(100);
  printf("Allocated %d bytes\n\r", GetHandleSize(p));
  ReallocateHandle(p, 200);
  printf("reallocated %d bytes\n\r", GetHandleSize(p));
  HLock(p);
  DisposeHandle(p);
  HUnlock(p);
  DisposeHandle(p);
  p = NewPtr(100);
  printf("Allocated %d bytes\n\r", GetPtrSize(p));
  SetPtrSize(p,200);
  printf("Now allocated %d bytes\n\r", GetPtrSize(p));
  DisposePtr(p);
}
*/
