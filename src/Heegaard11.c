#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L    9 User_Says_Quit(void)
L   65 Too_Many_Components_ALert(void)
********************************************************************************************/

int User_Says_Quit(void)
{
	unsigned char 	*ptr = NULL;
	
	unsigned long	newlimit;
	
GET_RESPONSE1:
	if(Batch == FALSE) SysBeep(5);
	ptr = (unsigned char *) NewPtr(100);
	if(ptr == NULL) Mem_Error();
	printf("\n\nMemory is getting rather low! It may be wise to save and quit!");	
	printf("\n\nHeegaard has used %lu K bytes of memory for storing presentations.",
		BytesUsed/1024);
	printf("\nThis exceeds the current limit of %lu K bytes reserved for this purpose.",
		BytesAvailable/1024);
	printf("\n\nHIT 'c' TO CONTINUE RUNNING THIS EXAMPLE.");	
	printf("\nHIT 't' TO TERMINATE THIS RUN.");
	printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS NOW IN MEMORY.");
	printf("\nHit 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY.");
	GET_RESPONSE2:	
	switch(WaitkbHit())
		{
		case 'c':
			GET_RESPONSE3:
			printf("\n\nENTER A NEW UPPER LIMIT ON THE AMOUNT OF MEMORY TO BE USED (in K) AND HIT 'return'.	");
			newlimit = 0L;
			ReadString((char *)ptr, GetPtrSize(ptr));
			sscanf((char *) ptr,"%lu",&newlimit);
			BytesAvailable = 1024*newlimit;
			printf("\nBytesAvailable = %lu K.",BytesAvailable/1024);
			if(BytesAvailable <= BytesUsed)
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\nThe number of bytes available must exceed the number already used!!");
				goto GET_RESPONSE3;
				}
			DisposePtr((char *) ptr);	
			return(FALSE);
		case 't':
			DisposePtr((char *) ptr);
			return(TRUE);
		case 'v':
			DisposePtr((char *) ptr);
			Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,NULL);
			goto GET_RESPONSE1;	
		 case 'w':
		 	DisposePtr((char *) ptr);
            printf("\n\n     Sorting presentations. . . .");
            Sort_Presentations_In_Memory(1);
            goto GET_RESPONSE1;
		default:
			if(Batch == FALSE) SysBeep(5);
			goto GET_RESPONSE2;		
		}			
}

void Too_Many_Components_ALert(void)
{
	printf("\n\n	When Heegaard finds diagrams which are not connected, it checks if presentations"); 
	printf("\n	of each component are on file. If a component's presentation is not on file, Heegaard");
	printf("\n	tries to create a new summand for that component. However, Heegaard can carry at most");
	printf("\n	%d components, and already has %d components/summands.",MAXNUMCOMPONENTS,TotalComp );
	if(BreadthFirstSearch)
		{
		printf("\n		Rerunning using Depth-First Search may help because Depth-First Search runs new");
		printf("\n	presentations earlier than Breadth-First Search.");
		}
	
	printf("\n\n	Hit 'v' or 'w' and 'v' to print the current presentations in memory, and scroll back");
	printf("\n	to see some details.");
}

