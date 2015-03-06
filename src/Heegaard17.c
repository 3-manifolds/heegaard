#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   17 Just_Delete_Primitives(char)
L  392 Reduce_The_Initial_Presentation_To_Minimal_Length_SPC(void)
L  449 Get_User_Input_SPC(char c)
L  562 OneGenerator_SPC(unsigned char NumMissingGenerators)
L  578 OneRelator_SPC(unsigned char NumMissingGenerators)
L  595 Defining_Relator_SPC(void)
L  773 Verify_Length(unsigned char ***MyRelators,int MyNumRelators)
L  788 Report_SPC(int *Table)
********************************************************************************************/
	
unsigned char  *NMG = NULL;	
		
int Just_Delete_Primitives(char F1)
{

char			c;

unsigned char   i,
				j,
				k,
				TheTestRelator,
				NumMissingGenerators,
				NumUncheckedRels,
				*p,
				*q,
				**MyTemp;

int             GUIRV,
				ii,
				jj,
				kk,
				NumEmptyRelators,
				SNumGenerators,
				SNumRelators;
				
long			NumRelatorsChecked,
				NumPrimitivesFound;

	NumRelatorsChecked = NumPrimitivesFound = 0L;
	NumRelTooLong = 0;

	NMG = (unsigned char*) NewPtr((sizeof(char)*MAX_SAVED_PRES));
	if(NMG == NULL) Mem_Error();

	if(Batch != 4 && Batch != 14 && Batch != 15 && Batch != 53)
	Turn_Micro_Print_On();	
	
	if(Micro_Print)
		{
		printf("\n\n/*******************************************************************************");	
		printf("\n This routine recursively reduces the number of generators of a presentation by");
		printf("\n deleting primitive relators.");	
		printf(" The presentation's realizability is not checked.");		
		printf("\n You can 'protect' a relator R which should not be deleted by replacing R with");
		printf("\n R^n, n > 1, in the input.");
		printf("\n*******************************************************************************/\n");
		}

	SNumRelators = Delete_Dups();
	if(SNumRelators < NumRelators)
		{
		if(Batch == FALSE) printf("\n");
		if(NumRelators - SNumRelators == 1)
			printf("\n Deleted a relator which is a cyclic conjugate of another relator or its inverse.");
		else
			printf("\n Deleted %d relators which are cyclic conjugates of other relators or their inverses."
			,NumRelators - SNumRelators);
		NumRelators = SNumRelators;
		}
	Vertices = 2*NumGenerators;	
	SNumGenerators = NumGenerators;
	switch(Reduce_The_Initial_Presentation_To_Minimal_Length_SPC())
		{
		case 0:
			break;
		case 3:
			DisposePtr((unsigned char*) NMG);
			return(3);
		case TOO_LONG:		
			DisposePtr((unsigned char*) NMG);
			return(TOO_LONG);
		}
	NMG[0] = SNumGenerators - NumGenerators;
	if(NMG[0] && Batch == FALSE) printf(" FGR %d",NMG[0]);
				
	while(1)
		{	
		  if( (c = mykbhit()) )
			 {
			 switch(GUIRV = Get_User_Input_SPC(c))
				{
				case 0:
					break;
				case 1:	
					DisposePtr((unsigned char*) NMG);
					return(1);
				case INTERRUPT:
					DisposePtr((unsigned char*) NMG);
					return(INTERRUPT);	
				}
			 }
	 
		if(NumFilled >= MAX_SAVED_PRES - 3)
			{
			if(Batch == FALSE) printf("\n");
			printf("\n Heegaard has saved the maximum number of presentations currently allowed.");
			if(Batch == FALSE) printf("\n");
			printf("\n Heegaard checked %ld relators, and found %ld primitives.",NumRelatorsChecked,NumPrimitivesFound);
			break;
			}
	
		if(F1) ReadPres = 0;
		else ReadPres = NumFilled - 1;
	
		while((NumUncheckedRels = *SUR[ReadPres][0][0]) == 0)
			{
			if(ReadPres == 0 && *SUR[0][0][0] == 0)
				{
				if(Batch == FALSE)
				printf("\n\n Done! Every relator of every presentation has been checked for primitivity.");
				if(Batch == FALSE) printf("\n");
				if(NumRelatorsChecked == 1)
					{
					if(NumPrimitivesFound == 1)
						{
						printf("\n Heegaard checked 1 relator. It was primitive.");
						if(Batch == 15 && H_Results != NULL)
							fprintf(H_Results,"\n Relators checked 1, Primitives 1");
						}
					else
						{
						printf("\n Heegaard checked 1 relator. It was not primitive.");
						if(Batch == 15 && H_Results != NULL)
							fprintf(H_Results,"\n Relators checked 1, Primitives 0");
						}
					}
				else
					{
					if(F1)
						{
						printf("\n Heegaard checked %ld relator(s) and found %ld primitive(s).",
						NumRelatorsChecked,NumPrimitivesFound);
						if(Batch == 15 && H_Results != NULL)
							fprintf(H_Results,"\n\n%s \nRelators checked %ld, Primitives %ld",
							PresName,NumRelatorsChecked,NumPrimitivesFound);					
						}
					else
						printf("\n Heegaard checked %ld relator(s) of %u presentation(s) and found %ld primitive(s).",
						NumRelatorsChecked,NumFilled,NumPrimitivesFound);				
					}
				goto END;
				}
			ReadPres = FR[ReadPres];
			}
		NumUncheckedRels = *SUR[ReadPres][0][0];
		
		k = abs(rand()) % NumUncheckedRels;
		k++;
		p = *SUR[ReadPres][0];
		for(i = 1, j = 0; ; ) 
			{
			j+= p[i];
			if(j == k) break;
			i++;
			}
		TheTestRelator = i;	
		p[i] = 0;
		p[0]--;
		OnStack--;
	
		if(Micro_Print) 
			{
			if(NumFilled > 1)
				{
				printf("\n\nPresentation %d is:",ReadPres + 1);
				Print_Relators(SUR[ReadPres],NR[ReadPres]);
				}
			printf("\n Checking relator %d of %d unchecked relators of %d relators of Presentation %d for primitivity.\n"
			,TheTestRelator, NumUncheckedRels, NR[ReadPres], ReadPres + 1);
			if(TheTestRelator > 1) printf("\n Swapped Relators[%d] and Relators[1].",TheTestRelator);
			}
		
		if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
		Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][TheTestRelator]));
		if(Relators[1] == NULL) Mem_Error();
		p = *Relators[1];	
		q = *SUR[ReadPres][TheTestRelator];
		while( (*p++ = *q++) ) ;
	
		if(Micro_Print)
			{
			printf("\n Relators[1] is:");
			Print_Relators(Relators,1);
			}
	
		NumGenerators 			= NG[ReadPres];
		NumMissingGenerators 	= NMG[ReadPres];
		Vertices 				= 2*NumGenerators;
		NumRelators 			= 1;
		NumRelatorsChecked ++;
		jj = CheckPrimitivity( );
		NumRelators 			= NR[ReadPres];	
	
		if(Micro_Print && jj != 1) 
		printf("\n\n Relator %d of Presentation %d is not primitive.",TheTestRelator, ReadPres + 1);

		if(jj == 1)
			{
			NumPrimitivesFound ++;
			if(Get_Relators_From_SUR(ReadPres)) 
				{
				DisposePtr((unsigned char*) NMG);
				return TOO_LONG;
				}

			if(TheTestRelator > 1)
				{
				MyTemp       			 = Relators[1];
				Relators[1] 			 = Relators[TheTestRelator];
				Relators[TheTestRelator] = MyTemp;
				}
		
			kk = Find_Primitives(1);
		
			if(kk != 1) continue;
		
			/************************************************************************
				Call Defining_Relator_SPC(), but first deal with the special cases 
						NumGenerators = 1 and NumRelators = 1.
			*************************************************************************/
		
			if(NumGenerators == 1)
				{
				NMG[NumFilled] = NumMissingGenerators;
				OneGenerator_SPC(NumMissingGenerators);
				if(Batch == FALSE) printf("\n");
				if(F1)
					{
					printf("\n Heegaard checked %ld relator(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumPrimitivesFound);
					if(Batch == 15 && H_Results != NULL)
						fprintf(H_Results,"\n\n%s \nRelators checked %ld, Primitives %ld",
						PresName,NumRelatorsChecked,NumPrimitivesFound);				
					}			
				else
					printf("\n Heegaard checked %ld relator(s) of %u presentation(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumFilled,NumPrimitivesFound);
				goto END;
				}
			if(NumRelators == 1)
				{
				NumMissingGenerators += NumGenerators - 1;
				NMG[NumFilled] = NumMissingGenerators;
				OneRelator_SPC(NumMissingGenerators);
				if(Batch == FALSE) printf("\n");
				if(F1)
					{
					printf("\n Heegaard checked %ld relator(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumPrimitivesFound);
					if(Batch == 15 && H_Results != NULL)
						fprintf(H_Results,"\n\n%s \nRelators checked %ld, Primitives %ld",
						PresName,NumRelatorsChecked,NumPrimitivesFound);				
					}			
				else
					printf("\n Heegaard checked %ld relator(s) of %u presentation(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumFilled,NumPrimitivesFound);
				goto END;
				}
		
			if(Defining_Relator_SPC() == TOO_LONG) continue;
			NumEmptyRelators = Freely_Reduce();
		
			if(Micro_Print && NumEmptyRelators)
				{
				if(NumEmptyRelators == 1)
					printf("\n\n One relator reduced to an empty word.");
				else
					printf("\n\n %d relators reduced to empty words.", NumEmptyRelators);
				}

			if(NumEmptyRelators == TOO_LONG) continue;
			Length 			= OrigLength; 
			SNumGenerators 	= NumGenerators; 	
			Rewrite_Input();
			NumMissingGenerators += SNumGenerators - NumGenerators;
			if(NumGenerators == 0) 										/*	The presentation is empty! 	*/
				{
				NMG[NumFilled] = NumMissingGenerators;
				OneRelator_SPC(NumMissingGenerators);
				if(Batch == FALSE) printf("\n");
				if(F1)
					{
					printf("\n Heegaard checked %ld relator(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumPrimitivesFound);	
					if(Batch == 15 && H_Results != NULL)
						fprintf(H_Results,"\n\n%s \nRelators checked %ld, Primitives %ld",
						PresName,NumRelatorsChecked,NumPrimitivesFound);				
					}		
				else
					printf("\n Heegaard checked %ld relator(s) of %u presentation(s) and found %ld primitive(s).",
					NumRelatorsChecked,NumFilled,NumPrimitivesFound);
				goto END;
				}
			
			Length 			= OrigLength;
			SNumGenerators  = NumGenerators;
			ii = Find_Flow_A(NORMAL,FALSE);
			if(ii == TOO_LONG) continue;
			if(ii == 1) 
				{
				SNumGenerators = NumGenerators;	
				Rewrite_Input();        	
				NumMissingGenerators += SNumGenerators - NumGenerators;			
				}
			
			if(Micro_Print)
				{
				if(Length == OrigLength)
					printf("\n\nThis presentation has minimal length.");
				else
					{
					printf("\n\n%lu automorphism(s) reduced the length from %ld to %lu.\n",
						Automorphisms,OrigLength,Length);
					if(SNumGenerators > NumGenerators)
						printf("\nAnd %d generator(s) no longer appear in the relators.\n",
						SNumGenerators - NumGenerators);
					} 
				 }
			 
			Freely_Reduce();
			Length 			= OrigLength;	 
			SNumRelators 	= Delete_Dups();
			if(SNumRelators < NumRelators)
				{
				NumRelators = SNumRelators;
				Vertices = 2*NumGenerators;
				for(i = 1,Length = 0L; i <= NumRelators; i++) Length += GetHandleSize((char **) Relators[i]) - 1;
				ii = Find_Flow_A(NORMAL,FALSE);
				if(ii == TOO_LONG) continue;
				if(ii == 1) 
					{
					SNumGenerators = NumGenerators;	
					Rewrite_Input();
					NumMissingGenerators += SNumGenerators - NumGenerators;				
					}
				}
					 
			if(Batch != 4 && On_File() == NumFilled) 
				{
				NMG[NumFilled] = NumMissingGenerators;
				if(Save_Pres(ReadPres,0,Length,1,6,1,1,0))
					{
					DisposePtr((unsigned char*) NMG);
					return(TOO_LONG);
					}
				if(NMG[NumFilled - 1] && Batch == FALSE) printf(" FGR %d",NMG[NumFilled - 1]);
				}
			}
		}
	END:
		if(Batch == 4 || Batch == 53)
			{
			DisposePtr((unsigned char*) NMG);
			if(NumPrimitivesFound > 0) 
				return(1);
			else
				return(0);
			}
		if(Batch) printf("\n");	
		Band_Sums = NumDiagrams = 0;
		for(ii = 0; ii < NumFilled; ii++) if(NMG[ii])
			{
			if(Batch == FALSE) printf("\n");
			printf("\nNote that a Presentation with FGR n means the original group G is a free product of a");
			printf(" free group of rank n, and a group G' on the generators explicitly appearing in the relators.\n");
			break;
			}
		if(NumRelTooLong == 1)
			printf("\nDeleting primitives produced one relator which was too long. Scroll back for details.");				
		if(NumRelTooLong > 1)
			printf("\nDeleting primitives produced %u relators which were too long. Scroll back for details.",NumRelTooLong);     
		NumRelTooLong = 0; 
		Sort_Presentations_In_Memory(2);
		DisposePtr((unsigned char*) NMG);
		return(1);
}

int Reduce_The_Initial_Presentation_To_Minimal_Length_SPC(void)
{
    int   	i,
    		SNumGenerators;
    
    long    Scratch;

	for(i = 1, Length = 0L; i <= NumRelators; i++) Length += GetHandleSize((char **) Relators[i]);
	Length -= NumRelators;
	
    Scratch = Length;
    SNumGenerators = NumGenerators;
    
    switch(Find_Flow_A(NORMAL,FALSE))
        {
        case 1:
            Rewrite_Input();
            if(NumGenerators < SNumGenerators) 
            	{
            	printf("\n\nAfter automorphisms some generators no longer appear in the relators.");
            	printf("\nHence the original group G is a free product of a free group of rank >= %d,",
            	SNumGenerators - NumGenerators);
            	printf(" and a group G' on the generators explicitly appearing in the relators.\n");
            	}
            break;
        case TOO_LONG:
            printf("\n\n     This presentation may be too long for Heegaard to handle. Sorry!");
            return(TOO_LONG);    
        } 
    
    if(Micro_Print)
    	{
		if(Length == Scratch)
			printf("\n\nThis presentation has minimal length.");
		else
			{
			printf("\n\n%lu automorphism(s) reduced the length from %ld to %lu.\n",
				Automorphisms,Scratch,Length);
			if(SNumGenerators > NumGenerators)
				printf("\nand reduced the number of generators from %d to %d.\n",
				SNumGenerators,NumGenerators);
    		} 
    	 }
    	         
    Canonical_Rewrite(Relators,FALSE,FALSE);
    if(Batch == 4 && Length < Scratch)
    	{
    	printf("\nAutomorphisms reduced the length of this presentation from %ld to %ld.",Scratch, Length);
    	printf("\nHence this presentation is not a Heegaard Splitting Rep!");
		if(H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep!",PresName);    	
    	return(3);
    	}
    if(Save_Pres(ReadPres,0,Length,1,2,1,1,0)) return(TOO_LONG);
    if(Micro_Print) Report(0,0,0,0,0,0,0,0,1,NULL); 
    return(0);        
}

int Get_User_Input_SPC(char c)
{
int		SMicro_Print,
		SMicro_Print_F;
		
	switch(c)
		{
		#ifdef DEBUGGING
		case 'd':
			Debug();
			break;
		case 'r':
			Print_Relators(Relators,NumRelators);
			break;
		#endif    
		case 's':
			if(Batch)
				printf("\n  Status: NumFilled %u, NumPresExamined %u",NumFilled,NumPresExamined);
			else	
				printf("\n  Status: NumFilled %u",NumFilled);
			break;                        
		case ' ':
			if(Batch == FALSE)
				{
				SMicro_Print = Micro_Print;
				SMicro_Print_F = Micro_Print_F;
				Micro_Print = FALSE;
				Micro_Print_F = FALSE;
				LIST_OPTIONS2:
				printf("\n\nHIT 't' TO TERMINATE THIS RUN.");
				printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS.");
				if(NumFilled > 1)
					{
					printf("\nHIT 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY BY SUMMAND NUMBER,");
					printf("\n        NUMGENERATORS, NUMRELATORS, LENGTH AND 'LEXICOGRAPHIC' ORDER.");
					}
				printf("\nHIT 'x' TO TOGGLE MICRO_PRINTING ON AND OFF.");
				printf("\nHIT 'r' TO RESUME RUNNING THIS EXAMPLE.");
				GET_RESPONSE1:		
				switch(WaitkbHit())
					{ 
					case 'r':
						break;    
					case 't':
						return(1);
					case 'v':
						Report(Band_Sums,NumDiagrams,OnStack,0,0,0,1,0,1,NULL);
						SaveMinima = TRUE;
						Input = RESET;
						goto LIST_OPTIONS2;
					case 'w':
						printf("\n\n     Sorting presentations. . . .");
						Sort_Presentations_In_Memory(2);
						SaveMinima = TRUE;
						Input = RESET;
						goto LIST_OPTIONS2;    
					case 'x':
						printf("\n    HIT 's' TO MICRO_PRINT TO THE SCREEN.");
						if(SMicro_Print)
							printf("\n    HIT 'o' TO TURN MICRO_PRINTING COMPLETELY OFF.");        
						GET_RESPONSE2:				
						switch(WaitkbHit())
							{
							case 's':
								SMicro_Print = TRUE;
								SMicro_Print_F = FALSE;
								break;                            
							case 'o':
								if(!SMicro_Print)
									{
									SysBeep(5);
									goto GET_RESPONSE1;                                
									}
								SMicro_Print = FALSE;
								SMicro_Print_F = FALSE;
								break;                        
							default:
								SysBeep(5);
								goto GET_RESPONSE2;
							}    
						break;                
					default:
						SysBeep(5);
						goto GET_RESPONSE1;    
					}
				printf("\n\n     Resumed running. . . .\n");
				NoReport = TRUE;
				Micro_Print = SMicro_Print;
				Micro_Print_F = SMicro_Print_F;
				break;
				}
			if(Batch)
				{
				printf("\nStatus: NumFilled %u, NumPresExamined %u",NumFilled,NumPresExamined);
				printf("\nHIT 't' TO TERMINATE THIS RUN.");
				printf("\nHIT 'r' TO RESUME RUNNING.");
				GET_RESPONSE3:
				switch(WaitkbHit())
					{
					case 't':
						return(INTERRUPT);                           
					case 'r':
						break;                        
					default:
						SysBeep(5);
						goto GET_RESPONSE3;
					}
				break;
				}
		}
	return(0);	
}

int OneGenerator_SPC(unsigned char NumMissingGenerators)
{
	/****************************************************************************************
		This routine is called when only one generator appears in Relators[] and 
		Relators[1] is a defining relator for that generator. This implies the current
		presentation is an empty presentation of the trivial group.
	****************************************************************************************/
	
	NumGenerators = NumRelators = 0;
	Length = 0L;
	if(Save_Pres(ReadPres,0,Length,1,6,1,2,0)) return(0);
	if(NumMissingGenerators && Batch == FALSE) printf(" FGR %d",NumMissingGenerators);
	if(Micro_Print) printf("\n\n Presentation %d is empty, and 'presents' the trivial group.",NumFilled);
	return(1);
}

int OneRelator_SPC(unsigned char NumMissingGenerators)
{
	/****************************************************************************************
		This routine is called when only Relator[1] appears in Relators[]. Since Relators[1]
		is a defining relator, this implies the current presentation is a presentation of a
		free group on (NumGenerators - 1) generators.
	****************************************************************************************/
	
	NumGenerators = 0;
	NumRelators = 0;
	Length = 0L;
	if(Save_Pres(ReadPres,0,Length,1,6,1,2,0)) return(0);
	if(NumMissingGenerators && Batch == FALSE) printf(" FGR %d",NumMissingGenerators);
	if(Micro_Print) printf("\n\n Presentation %d is empty, and 'presents' the trivial group.",NumFilled);
	return(1);
}

int Defining_Relator_SPC()
{        
    /******************************************************************************************
    	Defining_Relator_SPC() is called when Relators[1] is a defining relator, NumGenerators
    	> 1 and NumRelators > 1. It makes the appropriate substitutions in the remaining 
    	relators, reduces the number of generators and the number of relators and calls ???
    ******************************************************************************************/
    
    unsigned char	s,
					t,
					x,
					y,
					z,
					*p,
					*q,
					*r,           
					**Temp;
                            
	unsigned int    C[125],
					j,
					k;
                            
    unsigned long	M;
      
      
    if(Micro_Print)
    	{
    	printf("\n\n Going into Defining_Relator_SPC() the presentation is:");
    	Print_Relators(Relators,NumRelators);
    	}  
    	
	/**************************************************************************************            
		Look for and select a random generator which appears only once in Relators[1].
	**************************************************************************************/
	
	for(x = 'A',y = 'a'; x < 'A' + NumGenerators; x++,y++)    C[x] = C[y] = 0;
	p = *Relators[1];
	while( (z = *p++) ) C[z]++;
	for(x = 'A',y = 'a',j = 0; x < 'A' + NumGenerators; x++,y++) if(C[x] + C[y] == 1) j++;
	if(j == 0) return(TOO_LONG);
	k = abs(rand()) % j;
	k++;
	for(x = 'A',y = 'a',j = 0; x < 'A' + NumGenerators; x++,y++)
		if(C[x] + C[y] == 1 && ++j == k) break;
		
	/**********************************************************************************             
		Relators[1] is a defining relator, and 'x' together with 'y = x^{-1}' appear
							only once in Relators[1].
	**********************************************************************************/
	
								
	/**********************************************************************************     
			Create two strings which will be substituted for the defined
						 generator 'x' and its inverse 'y'.                                            
	**********************************************************************************/

	if(Micro_Print)
		{
		printf("\n\nRelators[1] is a defining relator which was used to eliminate generator %c.",x);
		printf("\nSwapped Relator %d with Relator %d and reduced the number of relators.",1,NumRelators);
		}
		
	p = *Relators[1];
	while( (z = *p++) ) if(z == x || z == y) break;
	if(Temp6 != NULL) DisposeHandle((char **) Temp6);
	Temp6 = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[1]) -1);
	if(Temp6 == NULL) Mem_Error();
	q = *Temp6;
	while(*p) *q++ = *p++;
	p = *Relators[1];
	while( (z != *p) ) *q++ = *p++;
	*q = EOS;
	p = *Temp6;
	if(Temp7 != NULL) DisposeHandle((char **) Temp7);
	Temp7 = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[1]) -1);
	if(Temp7 == NULL) Mem_Error();
	q = *Temp7;
	while( (*q++ = *p++) ) ;
	if(z == x) Inverse(*Temp6);
	else Inverse(*Temp7);
		
	/**********************************************************************************
		 Exchange Relator[NumRelators] and Relator[1].
	 *********************************************************************************/
	 
		Temp = Relators[1];
		Relators[1] = Relators[NumRelators];
		Relators[NumRelators] = Temp;
		
	/**********************************************************************************
					Make the appropriate substitutions in the relators.         
	**********************************************************************************/
	
	for(j = 1; j < NumRelators; j++)
		{
		k = 0;    
		p = *Relators[j];
		while( (z = *p++) ) if(z == x || z == y) k++;
		if(k)
			{
			M = GetHandleSize((char **) Relators[j]) + k*GetHandleSize((char **) Relators[NumRelators]);
			M -= 3*k;
			if(M > MAXLENGTH)
				{
				if(Batch == FALSE)
				printf("\nA relator derived from Pres %d is too long after deleting a primitive.",ReadPres + 1);
				NumRelTooLong ++;
				return(TOO_LONG);
				}
			if(Temp4 != NULL) DisposeHandle((char **) Temp4);
			Temp4 = (unsigned char **) NewHandle(M);	
			if(Temp4 == NULL) Mem_Error(); 
			q = *Temp4;   
			p = *Relators[j];
			while( (z = *p++) )
				{
				*q++ = z;
				if(z == x)
					{
					q--;
					r = *Temp6;
					while( (*q++ = *r++) ) ;
					q--;
					}
				if(z == y)
					{
					q--;
					r = *Temp7;
					while( (*q++ = *r++) ) ;
					q--;
					}
				}
			*q = EOS;    
			Temp = Relators[j];
			Relators[j] = Temp4;
			Temp4 = Temp;                    
			}        
		}
		
	/**********************************************************************************
				If x is not the largest generator, replace the largest generator
				throughout with x and its inverse with y.                                
	**********************************************************************************/
			 
	if(x < NumGenerators + 64)
		{    
		s = NumGenerators + 64;
		t = NumGenerators + 96;
		for(j = 1; j < NumRelators; j++)
			{
			p = *Relators[j];
			while( (z = *p) )
				{
				if(z == s) *p = x;
				if(z == t) *p = y;
				p++;
				}
			}
		if(Micro_Print)
			printf("\n\nRewrote the relators by replacing %c with %c.",s,x);
		}    
		
	/**********************************************************************************
						Update NumGenerators and NumRelators.
	**********************************************************************************/
	
	NumGenerators --;
	NumRelators --;
	
	if(Micro_Print)
		{
		printf("\n\n In Defining_Relator_SPC(). The presentation is:");
		Print_Relators(Relators, NumRelators);
		}
	
	return(0);             
}

int Verify_Length(unsigned char ***MyRelators,int MyNumRelators)
{
	unsigned char 	*p;
	
    int    			i;
        
    long			MyLength;
    
    for(i = 1, MyLength = 0L; i <= MyNumRelators; i++)    
        for(p = *MyRelators[i]; *p; p++) MyLength ++;    
	if(MyLength == Length) return(FALSE);
	printf("\n\n MyLength = %ld, Length = %ld.",MyLength, Length);
	return(TRUE);
}

int Report_SPC(int *Table)
{
    /******************************************************************************************
        Report_SPC() is an output routine called on termination of the sorting
        routine which sorts the presentations presently in memory after SPC() runs.
    ******************************************************************************************/  
                       
    unsigned int    i,
                    j,
                    MyMinNumGenerators,
                    MyMinNumRelators,
                    n;                                
    
    unsigned long    Length;
    
    if(Batch == 0 || B14B15PrintPres) for(n = 0; n < NumFilled; n++)
        {
        i = Table[n];
        NumRelators = NR[i];
        Length = SURL[i];
        printf("\n\nPresentation %d  of Summand %u:  Gen  %d  Length  %lu  From Pres %u  NumHits %d  ",
        i+1,ComponentNum[i],NG[i],Length,FR[i]+1,SURNumX[i]);
        printf("FP  FGR %d", NMG[i]);
                
        printf("\n");
        for(j = 1; j <= NumRelators; j++) printf("\n    %s",*SUR[i][j]);
        } 
        
    if(Batch == 14 && !B14B15PrintPres)
    	{
    	MyMinNumGenerators = NG[Table[NumFilled -1]];
    	MyMinNumRelators = NR[Table[NumFilled -1]];
    	if(MyMinNumGenerators > 1 && MyMinNumRelators > 0) return(0);
    	MyMinNumGenerators ++;
    	MyMinNumRelators ++;
    	for(n = 0; n < NumFilled; n++)
    		{
			i = Table[n];
			if((NG[i] > MyMinNumGenerators) && (NR[i] > MyMinNumRelators)) continue;
			NumRelators = NR[i];
			Length = SURL[i];
			printf("\n\nPresentation %d  of Summand %u:  Gen  %d  Length  %lu  From Pres %u  NumHits %d  ",
			i+1,ComponentNum[i],NG[i],Length,FR[i]+1,SURNumX[i]);
			printf("FP  FGR %d", NMG[i]);
				
			printf("\n");
			for(j = 1; j <= NumRelators; j++) printf("\n    %s",*SUR[i][j]);
 			}   		
		}                           
    return(0);    
}
