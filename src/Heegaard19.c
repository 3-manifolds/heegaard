#include "Heegaard.h"
#include "Heegaard_Dec.h"

unsigned char	NumAEXPS,
				NumBEXPS;
				
int				Alphas[4],
				Betas[4],
				EXPA1_SF[3],
				EXPA2_SF[6],
				EXPB1_SF[3],
				EXPB2_SF[6],
				NEXA1_SF[3],
				NEXA2_SF[6],
				NEXB1_SF[3],
				NEXB2_SF[6];
				

/****************************** function prototypes *****************************************
L   33 Init_Genus_Two_Seifert_Fibered(int*,int,int)
L  120 Genus_Two_Seifert_Fibered(int OrbitNum)
L  486 Get_Genus_2_SF_EXPS1(void)
L  614 Get_Genus_2_SF_EXPS2(void)
L  734 Do_Aut_SF(int Num)
L  923 Get_SF_Alphas1(int Num)
L 1173 Get_SF_Invariants(int OrbitNum)
L 1343 Get_SF_Alphas2(int Num)
L 1815 SF_Sort_And_Print(int H1, int n, int A1, int A2, int A3, int B1, int B2, int B3, 
		int NumSolns, int* SolV)
L 2055 Test_Transverse(void)
********************************************************************************************/				

int Init_Genus_Two_Seifert_Fibered(int* MyTable,int MyStart,int MyCompNum)
{
	unsigned char	*p,
					*q;
					
	int				i,
					m,
					n,
					MultipleSolns,
					NumSFChecked,
					NumSFFound;

	/****** Check which genus two presentations are Seifert Fibered. ******/	

	for(n = MyStart,NumSFChecked = NumSFFound = MultipleSolns = 0; n >= 0; n--) 
		{
		ReadPres 		= MyTable[n];
		if(CS[ReadPres] > 0) continue;
		if(ComponentNum[ReadPres] < MyCompNum) continue;
		if(ComponentNum[ReadPres] > MyCompNum) return(n);
		NumGenerators 	= NG[ReadPres];
		NumRelators 	= NR[ReadPres];
		if(NumGenerators != 2) continue;
		if(NumRelators > 2) continue;
		Vertices		= 2*NumGenerators;
		Length 			= SURL[ReadPres];
	    for(i = 1; i <= NumRelators; i++)
			{
			if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
			Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][i]));
			if(Relators[i] == NULL) Mem_Error();
			q = *Relators[i];    
			p = *SUR[ReadPres][i];
			while( (*q++ = *p++) ) ;
			}
		NumSFChecked ++;
		if(NumSFChecked > 100) return(n);
		SFSolV[1] = ReadPres;
		SFSolV[2] = MyCompNum;
		SFSolV[3] = NumRelators;
		m = Genus_Two_Seifert_Fibered(ReadPres + 1);
		if(B10B11Finite && FoundFiniteSF == FALSE && (Batch == 10 || Batch == 11)) return(n);
		if(m == 13 || m == 14) NumSFFound ++;
		if(m == 14) MultipleSolns ++;	
		if(NumSFFound > MultipleSolns) 
			{
			if(SFSols[MyCompNum]) DisposePtr((int*) SFSols[MyCompNum]);
			SFSols[MyCompNum] = (int*) NewPtr(sizeof(int)*20);
			if(SFSols[MyCompNum] == NULL) Mem_Error();
			for(i = 0; i < 20; i++) SFSols[MyCompNum][i] = SFSolV[i];
			return(n);
			}
		if(MultipleSolns == 1)
			{
			if(SFSols[MyCompNum]) DisposePtr((int*) SFSols[MyCompNum]);
			SFSols[MyCompNum] = (int*) NewPtr(sizeof(int)*20);
			if(SFSols[MyCompNum] == NULL) Mem_Error();
			for(i = 0; i < 20; i++) SFSols[MyCompNum][i] = SFSolV[i];			
			}	
		}
			
	if(NumSFFound == 1)
		{
		if(NumSFChecked == 1)
			printf("\n\n Heegaard checked one presentation, which was SF.");			
		else	
			printf("\n\n Heegaard found one SF presentation in the %d presentations checked.",NumSFChecked);
		}
	if(NumSFFound > 1)
		{
		printf("\n\n Heegaard found %d SF presentations in the %d presentations checked.",NumSFFound,NumSFChecked);
		}
	if(MultipleSolns && NumSFFound == MultipleSolns)
		{
		printf("\n\n		A Potential Seifert Fibration Ambiguity Exists!");
		printf("\n These ambiguities arise in the following way: Suppose M = ±SF(0;e;B1/A1,B2/A2,B3/A3).");
		printf("\n Then, given A1,A2,A3,B1, and B2, which Heegaard computes, there must exist integers"); 
		printf("\n B3 and e with gcd(A3,B3) = 1 such that at least one of:");
		printf("\n 	|H1(M)| = B1*A2*A3 + A1*B2*A3 + A1*A2*B3 - e*A1*A2*A3,");
		printf("\n 	|H1(M)| = (A1-B1)*A2*A3 + A1*(A2-B2)*A3 + A1*A2*B3 - e*A1*A2*A3");
		printf("\n is satisfied. Generally, there is only one solution with gcd(A3,B3) = 1.");
		printf("\n However, there may by two solutions when 2*A3(A1*B2+A2*B1) = 0 mod A1*A2.");
		}
	
	return(n);		
}
		
int Genus_Two_Seifert_Fibered(int OrbitNum)
{
char			CheckedRelatorTwo,
				SF;

unsigned char	*p,
				*q,
				**Copy_Of_Rel_3[3],
				**Temp,
				x;
							
int				i,
				NumSolns,
				NumTries;
				
unsigned int	*ptr1,
				*ptr2;				
				
long			H1;

unsigned long	SLength;			
					
	SF = FALSE;
				
	if(NumRelators > 2) 	return(1);
	if(NumGenerators != 2) 	return(2);
	
	/** SFSolV[] has the form: [Flag,ReadPres,MyCompNum,NumRelators,e,B1,A1,B2,A2,B3,A3,e',B1',A1',B2',A2',B3',A3'] **/
	
	SFSolV[0] = 0;

	/**********************************************************************************
				Save a copy of the current relators in Copy_Of_Rel_3[]. 
	**********************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		Copy_Of_Rel_3[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(Copy_Of_Rel_3[i] == NULL) Mem_Error();
		p = *Copy_Of_Rel_3[i];
		q = *Relators[i];
		while( (*p++ = *q++) ) ;                    
		}
		
	/*********************************************************************************
				If there are two relators, check for S^1 X S^2.
	*********************************************************************************/
	
	if(NumRelators == 2)
		{
		ptr1 = (unsigned int *)NewPtr(sizeof(int)*100);
		if(ptr1 == NULL) Mem_Error();
		ptr2 = (unsigned int *)NewPtr(sizeof(int)*100);
		if(ptr2 == NULL) Mem_Error();

		ptr1['A'] = ptr1['B'] = ptr1['a'] = ptr1['b'] = 0;
		p = *Relators[1];
		while( (x = *p++) ) ptr1[x] ++;
		ptr2['A'] = ptr2['B'] = ptr2['a'] = ptr2['b'] = 0;
		p = *Relators[2];
		while( (x = *p++) ) ptr2[x] ++;
			
		H1 = (ptr1['A'] - ptr1['a'])*(ptr2['B'] - ptr2['b']) - (ptr2['A'] - ptr2['a'])*(ptr1['B'] - ptr1['b']);
	  	
		if((ptr1['A'] == 0 || ptr1['a'] == 0) && (ptr1['B'] == 0 || ptr1['b'] == 0))
			{
			Vertices = 2*NumGenerators;
			SLength = Length;
			Length = GetHandleSize((char **)Relators[1]) - 1;
			Find_Flow_A(NORMAL,1);			
			Length = SLength;
			if(Length1 == 1) 
				{
				if(H1 == 0)
					{
					FoundSF = TRUE;
					printf("\n\nA Relator is primitive, and the manifold of Orbit %d is S^1 X S^2.",OrbitNum);
					if(B10B11Recognized) 
						{
						SFSolV[0] = 1;
						SFSolV[4] = 1;
						SFSolV[5] = 0;						
						}
					DisposePtr((char *) ptr1);
					DisposePtr((char *) ptr2);  		
					return(13);
					}
				if(Lens_Space()) 
					{
					DisposePtr((char *) ptr1);
					DisposePtr((char *) ptr2);						
					return(1);
					}
				FoundSF = TRUE;
				FoundFiniteSF = TRUE;
				printf("\n\nA Relator is primitive, so the manifold of Orbit %d is a Lens space: L(%lu,%lu)",OrbitNum,P,Q);
				if(B10B11Finite || B10B11Recognized) 
					{
					SFSolV[0] = 2;
					SFSolV[4] = Q;
					SFSolV[5] = P;
					}
				DisposePtr((char *) ptr1);
				DisposePtr((char *) ptr2);					
				return(13);
				}
			}
		if((ptr2['A'] == 0 || ptr2['a'] == 0) && (ptr2['B'] == 0 || ptr2['b'] == 0))
			{
			Vertices = 2*NumGenerators;
			if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
			Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Rel_3[2]));
			if(Relators[1] == NULL) Mem_Error();
			p = *Relators[1];
			q = *Copy_Of_Rel_3[2];
			while( (*p++ = *q++) ) ;
			if(Relators[2] != NULL) DisposeHandle((char **) Relators[2]);
			Relators[2] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Rel_3[1]));
			if(Relators[2] == NULL) Mem_Error();
			p = *Relators[2];
			q = *Copy_Of_Rel_3[1];
			while( (*p++ = *q++) ) ;
			SLength = Length;
			Length = GetHandleSize((char **)Relators[1]) - 1;
			Find_Flow_A(NORMAL,1);			
			Length = SLength;
			if(Length1 == 1) 
				{
				if(H1 == 0)
					{
					FoundSF = TRUE;
					printf("\n\nA Relator is primitive, and the manifold of Orbit %d is S^1 X S^2.",OrbitNum);
					if(B10B11Recognized) 
						{
						SFSolV[0] = 3;
						SFSolV[4] = 1;
						SFSolV[5] = 0;
						}					
					DisposePtr((char *) ptr1);
					DisposePtr((char *) ptr2);  		
					return(13);
					}
				if(Lens_Space()) 
					{
					DisposePtr((char *) ptr1);
					DisposePtr((char *) ptr2);					
					return(1);
					}
				FoundSF = TRUE;
				FoundFiniteSF = TRUE;
				printf("\n\nA Relator is primitive, so the manifold of Orbit %d is a Lens space: L(%lu,%lu)",OrbitNum,P,Q);
				if(B10B11Finite || B10B11Recognized)				
					{
					SFSolV[0] = 4;
					SFSolV[4] = Q;
					SFSolV[5] = P;
					}				 
				DisposePtr((char *) ptr1);
				DisposePtr((char *) ptr2);					
				return(13);					
				}
			}
							
		if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
		Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Rel_3[1]));
		p = *Relators[1];
		q = *Copy_Of_Rel_3[1];
		while( (*p++ = *q++) ) ;
		if(Relators[2] != NULL) DisposeHandle((char **) Relators[2]);
		Relators[2] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Rel_3[2]));
		p = *Relators[2];
		q = *Copy_Of_Rel_3[2];
		while( (*p++ = *q++) )
		    ;
    	DisposePtr((char *) ptr1);
    	DisposePtr((char *) ptr2);
    	}
		
	CheckedRelatorTwo = FALSE;
	CHECKRELATORTWO:
						
		Vertices = 2*NumGenerators;
		SLength = Length;
		Length = GetHandleSize((char **)Relators[1]) - 1;
		i = Find_Flow_A(NORMAL,1);	
		Length = SLength;

		NumTries = 0;
	TRYAGAIN:
		if(NumTries++ > 1) goto END;
	
		if(Get_Genus_2_SF_EXPS1()) return(5);

		if(NumRelators == 2) i = Get_Genus_2_SF_EXPS2();
		if(i) return(6);
	/*
	Print_Relators(Relators,NumRelators);
		printf("\n");
		for(i = 0; i < 3; i++) printf(" EXPA1_SF[%d] = %3d",i,EXPA1_SF[i]);
		printf("\n");
		for(i = 0; i < 3; i++) printf(" NEXA1_SF[%d] = %3d",i,NEXA1_SF[i]);
		printf("\n");
		for(i = 0; i < 3; i++) printf(" EXPB1_SF[%d] = %3d",i,EXPB1_SF[i]);
		printf("\n");
		for(i = 0; i < 3; i++) printf(" NEXB1_SF[%d] = %3d",i,NEXB1_SF[i]);
		printf("\n");
		printf(" NumAEXPS = %d, NumBEXPS  = %d",NumAEXPS, NumBEXPS);	
	
	if(NumRelators > 1)
		{	
		printf("\n\n");
		for(i = 0; i < 6; i++) printf(" EXPA2_SF[%d] = %3d",i,EXPA2_SF[i]);
		printf("\n");
		for(i = 0; i < 6; i++) printf(" NEXA2_SF[%d] = %3d",i,NEXA2_SF[i]);
		printf("\n");
		for(i = 0; i < 6; i++) printf(" EXPB2_SF[%d] = %3d",i,EXPB2_SF[i]);
		printf("\n");
		for(i = 0; i < 6; i++) printf(" NEXB2_SF[%d] = %3d",i,NEXB2_SF[i]);
		}																		*/																				

		if(NumAEXPS == 1 && NumBEXPS == 1) 
			{
			if(NEXA1_SF[0] == 1 && NEXB1_SF[0] == 1)
				{
				SF = TRUE;
				Get_SF_Alphas1(1);
				NumSolns = 0;
				if(Get_SF_Invariants(OrbitNum) == 5) NumSolns = 2;
				goto END;
				}
			}
		if(NumAEXPS == 2 && NumBEXPS == 1)
			{
			if(abs(EXPA1_SF[0] - EXPA1_SF[1]) == 1)
				{
				SF = TRUE;
				Get_SF_Alphas1(2);
				NumSolns = 0;
				if(Get_SF_Invariants(OrbitNum) == 5) NumSolns = 2;
				goto END;
				}
			if((NEXB1_SF[0] == 2) && (abs(EXPB1_SF[0]) == 1) && (EXPA1_SF[0]*EXPA1_SF[1] < 0))
				{
				if(Do_Aut_SF(1) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;
				}
			if((NEXB1_SF[0] == 2) && (EXPA1_SF[0]*EXPA1_SF[1] == -1))
				{
				SF = TRUE;			/* R = AB^saB^S: SF over the Mobius band. */
				FoundSF = TRUE;
				printf("\n\nA relator has the form AB^SaB^S ==> SF over the Mobius band with one exceptional fiber of index %d.",abs(EXPB1_SF[0])); 		
				if(B10B11Recognized)
					{
					if(NumRelators == 1) SFSolV[0] = 5;
					else SFSolV[0] = 6;
					}	
				}		
			}
		if(NumAEXPS == 1 && NumBEXPS == 2)
			{
			if(abs(EXPB1_SF[0] - EXPB1_SF[1]) == 1)
				{
				SF = TRUE;
				Get_SF_Alphas1(3);
				NumSolns = 0;
				if(Get_SF_Invariants(OrbitNum) == 5) NumSolns = 2;
				goto END;
				}
			if((NEXA1_SF[0] == 2) && (abs(EXPA1_SF[0]) == 1) && (EXPB1_SF[0]*EXPB1_SF[1] < 0))
				{
				if(Do_Aut_SF(2) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;
				}
			if((NEXA1_SF[0] == 2) && (EXPB1_SF[0]*EXPB1_SF[1] == -1))
				{
				SF = TRUE;			/* R = A^PBA^Pb: SF over the Mobius band. */
				FoundSF = TRUE;
				printf("\n\nA relator has the form A^PBA^Pb ==> SF over the Mobius band with one exceptional fiber of index %d.",abs(EXPA1_SF[0])); 		
				if(B10B11Recognized)
					{
					if(NumRelators == 1) SFSolV[0] = 5;
					else SFSolV[0] = 6;
					}	
				}		
			}
		if(NumAEXPS == 2 && NumBEXPS == 2)
			{
			if(EXPB1_SF[0]*EXPB1_SF[1] == 2 && EXPA1_SF[0]*EXPA1_SF[1] < 0 && (abs(EXPA1_SF[0]) == 1 || abs(EXPA1_SF[1]) == 1))
				{
				if(Do_Aut_SF(3) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;
				}
			if(EXPA1_SF[0]*EXPA1_SF[1] == 2 && EXPB1_SF[0]*EXPB1_SF[1] < 0 && (abs(EXPB1_SF[0]) == 1 || abs(EXPB1_SF[1]) == 1))
				{
				if(Do_Aut_SF(4) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;
				}		
			}	
		if(NumAEXPS == 3 && NumBEXPS == 1 && abs(EXPB1_SF[0]) == 1)
			{
			if(abs(EXPA1_SF[0] - EXPA1_SF[1]) == 1 && EXPA1_SF[0]*EXPA1_SF[2] < 0)
				{
				if(Do_Aut_SF(5) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}
			if(abs(EXPA1_SF[0] - EXPA1_SF[2]) == 1 && EXPA1_SF[0]*EXPA1_SF[1] < 0)
				{
				if(Do_Aut_SF(6) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}
			if(abs(EXPA1_SF[1] - EXPA1_SF[2]) == 1 && EXPA1_SF[0]*EXPA1_SF[1] < 0)
				{
				if(Do_Aut_SF(7) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}								
			}
		if(NumAEXPS == 1 && NumBEXPS == 3 && abs(EXPA1_SF[0]) == 1)
			{
			if(abs(EXPB1_SF[0] - EXPB1_SF[1]) == 1 && EXPB1_SF[0]*EXPB1_SF[2] < 0)
				{
				if(Do_Aut_SF(8) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}
			if(abs(EXPB1_SF[0] - EXPB1_SF[2]) == 1 && EXPB1_SF[0]*EXPB1_SF[1] < 0)
				{
				if(Do_Aut_SF(9) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}
			if(abs(EXPB1_SF[1] - EXPB1_SF[2]) == 1 && EXPB1_SF[0]*EXPB1_SF[1] < 0)
				{
				if(Do_Aut_SF(10) == TOO_LONG) return(TOO_LONG);
				goto TRYAGAIN;			
				}								
			}

END:
		
		if(NumRelators == 2 && SF == FALSE && CheckedRelatorTwo == FALSE)
			{
			CheckedRelatorTwo = TRUE;
			for(i = 1; i <= NumRelators; i++)
				{
				if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
				Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Rel_3[i]));
				if(Relators[i] == NULL) Mem_Error();
				p = *Relators[i];	
				q = *Copy_Of_Rel_3[i];
				while( (*p++ = *q++) ) ;                    
				}
			Temp 		= Relators[1];
			Relators[1] = Relators[2];
			Relators[2] = Temp;
			goto CHECKRELATORTWO;		
			}	

		for(i = 1; i <= NumRelators; i++) if(Copy_Of_Rel_3[i]) 
			{
			if(Copy_Of_Rel_3[i] != NULL) DisposeHandle((char **) Copy_Of_Rel_3[i]);
			Copy_Of_Rel_3[i] = NULL;
			}

		if(NumSolns == 2) return(14);
		if(SF == TRUE) return(13);
		return(0);
}


int Get_Genus_2_SF_EXPS1()
{
    /*********************************************************************************************
        This routine determines the exponents with which each generator appears in Relators[1].
    These exponents are left in the arrays EXPA1_SF[] and EXPB1_SF[]. The routine also counts
    the number of times each exponent of each generator appears in Relators[1], and saves these
    values in the arrays NEXA1_SF[] and NEXB1_SF[]. 
    	If Relators[1] has minimal length and Relators[1] is an SF relator, then no generator
    can appear in Relators[1] with more than 3 exponents. (Note exponents are signed here, so 
    e and -e are different exponents.)
    *********************************************************************************************/    
    						
register unsigned char  i,
						*p,
						x,
						y;
						
register int   			ex,
						sex,
						*q,
						*r;
    
    if(GetHandleSize((char**) Relators[1]) < 3) return(1);  /* Relators[1] is either empty, or a primitive! */
                           
    for(i = 0; i < 3; i++) NEXA1_SF[i] = NEXB1_SF[i] = EXPA1_SF[i] = EXPB1_SF[i] = 0;
	
	p = *Relators[1];
	x = *p++;
	if(!x) return(2);			/* Relators[1] is empty!  */
	ex = 1;
	while(*p == x)
		{
		ex++;
		p++;
		}
	sex = ex;
	ex = 0;
	if(*p) x = *p;
	else return(3); 			/* Relators[1] is a proper power!  */
	while( (y = *p++) )
		{
		if(x == y)
			ex++;
		else
			{
			if(x == 'A' || x == 'a') 
				{
				q = EXPA1_SF;
				r = NEXA1_SF;
				}
			else
				{
				q = EXPB1_SF;
				r = NEXB1_SF;
				}
			if(x > 96) ex = -ex;		
			for(x = 0; x < 3 && ex != q[x] && q[x]; x++) ;
			if(x == 3) return(4);							/* Not Genus 2 SF! */
			q[x] = ex;
			r[x]++;
			ex = 1;    
			}
		x = y;
		}
	y = **Relators[1];        
	if(x == y)
		{
		if(x == 'A' || x == 'a') 
			{
			q = EXPA1_SF;
			r = NEXA1_SF;
			}
		else
			{
			q = EXPB1_SF;
			r = NEXB1_SF;
			}
		ex += sex;	
		if(x > 96) ex = -ex;		
		for(x = 0; x < 3 && ex != q[x] && q[x]; x++) ;
		if(x == 3) return(5);							/* Not Genus 2 SF! */
		q[x] = ex;
		r[x]++;
		}
	else
		{
		if(x == 'A' || x == 'a') 
			{
			q = EXPA1_SF;
			r = NEXA1_SF;
			}
		else
			{
			q = EXPB1_SF;
			r = NEXB1_SF;
			}
		if(x > 96) ex = -ex;		
		for(x = 0; x < 3 && ex != q[x] && q[x]; x++) ;
		if(x == 3) return(6);							/* Not Genus 2 SF! */
		q[x] = ex;
		r[x]++;
		
		ex = sex;
		x = y;
		
		if(x == 'A' || x == 'a') 
			{
			q = EXPA1_SF;
			r = NEXA1_SF;
			}
		else
			{
			q = EXPB1_SF;
			r = NEXB1_SF;
			}
		if(x > 96) ex = -ex;		
		for(x = 0; x < 3 && ex != q[x] && q[x]; x++) ;
		if(x == 3) return(7);							/* Not Genus 2 SF! */
		q[x] = ex;
		r[x]++;
		}
    
    for(i = NumAEXPS = 0; i < 3; i++) if(EXPA1_SF[i]) NumAEXPS ++;
    for(i = NumBEXPS = 0; i < 3; i++) if(EXPB1_SF[i]) NumBEXPS ++;
    
    return(0);            
}

int Get_Genus_2_SF_EXPS2()
{
    /******************************************************************************************
        This routine does two things. First, it determines the exponents with which each 
    generator appears in Relators[2]. It leaves these exponents in the six columns of the 
    arrays EXPA2_SF[] and EXPB2_SF[].
        Second, the routine counts the number of appearances of each exponent and leaves
    these counts in the arrays NEXA2_SF[] and NEXB2_SF[].
    ******************************************************************************************/    
    						
register unsigned char  i,
						*p,
						x,
						y;
						
register int   			ex,
						sex,
						*q,
						*r;
                               
    for(i = 0; i < 6; i++) NEXA2_SF[i] = NEXB2_SF[i] = EXPA2_SF[i] = EXPB2_SF[i] = 0;
	
	p = *Relators[2];
	x = *p++;
	if(!x) return(2);			/* Relators[2] is empty!  */
	ex = 1;
	while(*p == x)
		{
		ex++;
		p++;
		}
	sex = ex;
	ex = 0;
	if(*p) x = *p;
	while( (y = *p++) )
		{
		if(x == y)
			ex++;
		else
			{
			if(x == 'A' || x == 'a') 
				{
				q = EXPA2_SF;
				r = NEXA2_SF;
				}
			else
				{
				q = EXPB2_SF;
				r = NEXB2_SF;
				}
			if(x > 96) ex = -ex;		
			for(x = 0; x < 6 && ex != q[x] && q[x]; x++) ;
			if(x == 6) return(4);							/* Too Many Exponents! */
			q[x] = ex;
			r[x]++;
			ex = 1;    
			}
		x = y;
		}
	y = **Relators[2];        
	if(x == y)
		{
		if(x == 'A' || x == 'a') 
			{
			q = EXPA2_SF;
			r = NEXA2_SF;
			}
		else
			{
			q = EXPB2_SF;
			r = NEXB2_SF;
			}
		ex += sex;	
		if(x > 96) ex = -ex;		
		for(x = 0; x < 6 && ex != q[x] && q[x]; x++) ;
		if(x == 6) return(5);							/* Too Many Exponents! */
		q[x] = ex;
		r[x]++;
		}
	else
		{
		if(x == 'A' || x == 'a') 
			{
			q = EXPA2_SF;
			r = NEXA2_SF;
			}
		else
			{
			q = EXPB2_SF;
			r = NEXB2_SF;
			}
		if(x > 96) ex = -ex;		
		for(x = 0; x < 6 && ex != q[x] && q[x]; x++) ;
		if(x == 6) return(6);							/* Too Many Exponents! */
		q[x] = ex;
		r[x]++;
		
		ex = sex;
		x = y;
		
		if(x == 'A' || x == 'a') 
			{
			q = EXPA2_SF;
			r = NEXA2_SF;
			}
		else
			{
			q = EXPB2_SF;
			r = NEXB2_SF;
			}
		if(x > 96) ex = -ex;		
		for(x = 0; x < 6 && ex != q[x] && q[x]; x++) ;
		if(x == 6) return(7);							/* Too Many Exponents! */
		q[x] = ex;
		r[x]++;
		}
    
return(0);            
}

int Do_Aut_SF(int Num)
{
char			x,
				y;
			
unsigned char	*p,
				**Temp;			
			
int				i,
				ii,
				j,
				k;	
		
unsigned long	SLength;				

	if(Num == 1)
		{
		if(EXPB1_SF[0] > 0) 
			x = 'B';
		else 
			x = 'b';
		if(abs(EXPA1_SF[0]) < abs(EXPA1_SF[1]))
			i = 0;
		else
			i = 1;
		if(EXPA1_SF[i] > 0)
			y = 'A';
		else
			y = 'a';
		k = abs(EXPA1_SF[i]);
		}
	if(Num == 2)
		{
		if(EXPA1_SF[0] > 0) 
			x = 'A';
		else 
			x = 'a';
		if(abs(EXPB1_SF[0]) < abs(EXPB1_SF[1]))
			i = 0;
		else
			i = 1;
		if(EXPB1_SF[i] > 0)
			y = 'B';
		else
			y = 'b';
		k = abs(EXPB1_SF[i]);
		}	
	if(Num == 3)
		{
		if(EXPB1_SF[0] > 0) 
			x = 'B';
		else 
			x = 'b';
		if(abs(EXPA1_SF[0]) > abs(EXPA1_SF[1]))
			i = 0;
		else
			i = 1;
		if(EXPA1_SF[i] > 0)
			y = 'A';
		else
			y = 'a';
		k = abs(EXPA1_SF[i]);
		}
	if(Num == 4)
		{
		if(EXPA1_SF[0] > 0) 
			x = 'A';
		else 
			x = 'a';
		if(abs(EXPB1_SF[0]) > abs(EXPB1_SF[1]))
			i = 0;
		else
			i = 1;
		if(EXPB1_SF[i] > 0)
			y = 'B';
		else
			y = 'b';
		k = abs(EXPB1_SF[i]);
		}
	if(Num == 5)
		{
		if(EXPB1_SF[0] > 0) 
			x = 'B';
		else 
			x = 'b';
		if(EXPA1_SF[2] > 0)
			y = 'A';
		else
			y = 'a';
		k = abs(EXPA1_SF[2]);
		}
	if(Num == 6)
		{
		if(EXPB1_SF[0] > 0) 
			x = 'B';
		else 
			x = 'b';
		if(EXPA1_SF[1] > 0)
			y = 'A';
		else
			y = 'a';
		k = abs(EXPA1_SF[1]);
		}					
	if(Num == 7)
		{
		if(EXPB1_SF[0] > 0) 
			x = 'B';
		else 
			x = 'b';
		if(EXPA1_SF[0] > 0)
			y = 'A';
		else
			y = 'a';
		k = abs(EXPA1_SF[0]);
		}		
	if(Num == 8)
		{
		if(EXPA1_SF[0] > 0) 
			x = 'A';
		else 
			x = 'a';
		if(EXPB1_SF[2] > 0)
			y = 'B';
		else
			y = 'b';
		k = abs(EXPB1_SF[2]);
		}
	if(Num == 9)
		{
		if(EXPA1_SF[0] > 0) 
			x = 'A';
		else 
			x = 'a';
		if(EXPB1_SF[1] > 0)
			y = 'B';
		else
			y = 'b';
		k = abs(EXPB1_SF[1]);
		}
	if(Num == 10)
		{
		if(EXPA1_SF[0] > 0) 
			x = 'A';
		else 
			x = 'a';
		if(EXPB1_SF[0] > 0)
			y = 'B';
		else
			y = 'b';
		k = abs(EXPB1_SF[0]);
		}										
	
	/****************************************************************************************************
		In order to get an SF relator in which one of the generators appears only with exponents ±2, we
		create an additional relator of the form xy^kxy^k, which automorphisms will reduce to x^2.
		This will force the SF relator to have standard form in which one of the generators appears
		only with exponents ±2.
	*****************************************************************************************************/
				
	ii = NumRelators + 1;
	if(Relators[ii] != NULL) DisposeHandle((char **) Relators[ii]);
	Relators[ii] = (unsigned char **) NewHandle((2*k + 3)*sizeof(char));
	if(Relators[ii] == NULL) Mem_Error();			
	p = *Relators[ii];
	*p++ = x;
	for(j = 0; j < k; j++) *p++ = y;
	*p++ = x;
	for(j = 0; j < k; j++) *p++ = y;
	*p = EOS;
	
	Temp 			= Relators[ii];
	Relators[ii] 	= Relators[ii-1];
	Relators[ii-1] 	= Temp;			
	
	NumGenerators 	= 2;
	NumRelators 	= ii;
	Vertices 		= 4;
	SLength 		= Length;	
	for(i = 1, Length = 0L; i < NumRelators; i++) Length += GetHandleSize((char **) Relators[i]) - 1;	
	Find_Flow_A(NORMAL,2);
	Temp 			= Relators[ii];
	Relators[ii] 	= Relators[ii-1];
	Relators[ii-1] 	= Temp;	
	NumRelators 	= ii - 1;
	Length 			= SLength;
	
return(0);
}

int Get_SF_Alphas1(int Num)
{
char		F1,
			FA,
			FB;

int			Alp3A,
			Alp3B,
			i,
			j,
			MA,
			MB,
			NA,
			NB;

	if(Num == 1)
		{
		Alphas[1] = abs(EXPA1_SF[0]);
		Alphas[2] = abs(EXPB1_SF[0]);
		if(NumRelators == 1) return(0);
		if(NumRelators == 2)
			{			
			/******************************************
			 	Get tentative Betas[1] and Betas[2]. 
			******************************************/
			
			for(i = 0, Betas[1] = 0; i < 6; i++) if(EXPA2_SF[i] && abs(EXPA2_SF[i]) != Alphas[1])
				{
				Betas[1] = EXPA2_SF[i] % Alphas[1];
				if(Betas[1] < 0) Betas[1] += Alphas[1];
				break;
				}
			for(i = 0, Betas[2] = 0; i < 6; i++) if(EXPB2_SF[i] && abs(EXPB2_SF[i]) != Alphas[2])
				{
				Betas[2] = EXPB2_SF[i] % Alphas[2];
				if(Betas[2] < 0) Betas[2] += Alphas[2];
				break;
				}
			for(i = 0, F1  = TRUE; i < 6; i++)
				{
				if(abs(EXPA2_SF[i]) == Alphas[1])
					{
					F1 = FALSE;
					break;
					}
				if(abs(EXPB2_SF[i]) == Alphas[2])
					{
					F1 = FALSE;
					break;
					}	
				}
			if(F1 == TRUE)
				{
				for(i = 0, Alphas[3] = 0; i < 6; i++) Alphas[3] += NEXA2_SF[i];
				
				/*********************************************************************************
					If both generators appear in the SF relator Relators[1] with exponents of 
					the same sign, Betas[2] needs to be replaced by -Betas[2] mod Alphas[2].
				*********************************************************************************/
				
				if(EXPA1_SF[0]*EXPB1_SF[0] > 0) Betas[2] = Alphas[2] - Betas[2];
				return(0);
				}
				
			/*** Check if an 'A' exponent has absolute value > 2 in R1 or R2. ***/
			
			if(Alphas[1] > 2)
				FA = TRUE;
			else for(i = 0, FA  = FALSE; i < 6; i++)
				{
				if(abs(EXPA2_SF[i]) > 2)
					{
					FA = TRUE;
					break;
					}
				}
				
			if(FA)
				{
				for(i = MA = 0; i < 6; i++)
					{
					j = abs(EXPA2_SF[i]);
					if(j && j != Alphas[1] && j > MA) MA = j;
					}
				if(MA == 0)
					{
					Alphas[3] = 0;
					return(0);					
					}	
				for(i = NA = 0; i < 6; i++)
					{
					j = abs(EXPA2_SF[i]);
					if(j && j != Alphas[1] && j != MA)
						{
						NA = j;
						break;
						}
					}
				if(NA == 0)
					{
					for(i = Alp3A = 0; i < 6; i++)
						{
						j = EXPA2_SF[i];
						if(j ==  MA) Alp3A += NEXA2_SF[i];
						if(j == -MA) Alp3A -= NEXA2_SF[i];
						}		
					}
				if(NA && (MA < Alphas[1]))
					{
					for(i = Alp3A = 0; i < 6; i++)
						{
						j = EXPA2_SF[i];
						if(j ==  MA) Alp3A += NEXA2_SF[i];
						if(j == -MA) Alp3A -= NEXA2_SF[i];
						if(j ==  NA) Alp3A -= NEXA2_SF[i];
						if(j == -NA) Alp3A += NEXA2_SF[i];
						}							
					}
				if(NA && (MA > Alphas[1]))
					{
					for(i = Alp3A = 0; i < 6; i++)
						{
						j = EXPA2_SF[i];
						if(j ==  MA) Alp3A += NEXA2_SF[i];
						if(j == -MA) Alp3A -= NEXA2_SF[i];
						if(j ==  NA) Alp3A += NEXA2_SF[i];
						if(j == -NA) Alp3A -= NEXA2_SF[i];
						}							
					}
				Alphas[3] = abs(Alp3A);					
				}	
				
			/*** Check if a 'B' exponent has absolute value > 2 in R1 or R2. ***/
			
			if(Alphas[2] > 2)
				FB = TRUE;
			else for(i = 0, FB  = FALSE; i < 6; i++)
				{
				if(abs(EXPB2_SF[i]) > 2)
					{
					FB = TRUE;
					break;
					}
				}
				
			if(FB)
				{
				for(i = MB = 0; i < 6; i++)
					{
					j = abs(EXPB2_SF[i]);
					if(j && j != Alphas[2] && j > MB) MB = j;
					}
				if(MB == 0)
					{
					Alphas[3] = 0;
					return(0);					
					}	
				for(i = NB = 0; i < 6; i++)
					{
					j = abs(EXPB2_SF[i]);
					if(j && j != Alphas[2] && j != MB)
						{
						NB = j;
						break;
						}
					}
				if(NB == 0)
					{
					for(i = Alp3B = 0; i < 6; i++)
						{
						j = EXPB2_SF[i];
						if(j ==  MB) Alp3B += NEXB2_SF[i];
						if(j == -MB) Alp3B -= NEXB2_SF[i];
						}		
					}
				if(NB && (MB < Alphas[2]))
					{
					for(i = Alp3B = 0; i < 6; i++)
						{
						j = EXPB2_SF[i];
						if(j ==  MB) Alp3B += NEXB2_SF[i];
						if(j == -MB) Alp3B -= NEXB2_SF[i];
						if(j ==  NB) Alp3B -= NEXB2_SF[i];
						if(j == -NB) Alp3B += NEXB2_SF[i];
						}							
					}
				if(NB && (MB > Alphas[2]))
					{
					for(i = Alp3B = 0; i < 6; i++)
						{
						j = EXPB2_SF[i];
						if(j ==  MB) Alp3B += NEXB2_SF[i];
						if(j == -MB) Alp3B -= NEXB2_SF[i];
						if(j ==  NB) Alp3B += NEXB2_SF[i];
						if(j == -NB) Alp3B -= NEXB2_SF[i];
						}							
					}
				Alphas[3] = abs(Alp3B);
				}
			
			if(Alphas[1] > 2 && Alphas[2] > 2)
				{
				Betas[1] = MA % Alphas[1];
				Betas[2] = MB % Alphas[2];
							
				/*********************************************************************************
					If both generators appear in the SF relator Relators[1] with exponents of 
					identical (resp different) signs and ALp3A and Alp3B also have identical 
					(resp different) signs, Betas[2] needs to be replaced by -Betas[2] mod 
					Alphas[2]. 		
				*********************************************************************************/
				
				if(EXPA1_SF[0]*EXPB1_SF[0] > 0 && Alp3A*Alp3B > 0) Betas[2] = Alphas[2] - Betas[2];
				if(EXPA1_SF[0]*EXPB1_SF[0] < 0 && Alp3A*Alp3B < 0) Betas[2] = Alphas[2] - Betas[2];
				
				}
						
			if(FA || FB) return(0);
					
			/*** Neither generator appears in R1 or R2 with exponents greater than 2! ***/
							
			Get_SF_Alphas2(1);					
			}
		}
		
	if(Num == 2)
		{
		for(i = 0, Alphas[1] = 0; i < 3; i++) Alphas[1] += NEXA1_SF[i]*(abs(EXPA1_SF[i]));
		Alphas[2] = abs(EXPB1_SF[0]);
		if(NumRelators == 1) return(0);
		if(NumRelators == 2)
			{
			for(i = 0, Betas[1] = 0; i < 3; i++) Betas[1] += NEXA1_SF[i];		
			Get_SF_Alphas2(2);
			}					
		}
	if(Num == 3)
		{
		for(i = 0, Alphas[2] = 0; i < 3; i++) Alphas[2] += NEXB1_SF[i]*(abs(EXPB1_SF[i]));
		Alphas[1] = abs(EXPA1_SF[0]);
		if(NumRelators == 1) return(0);
		if(NumRelators == 2)
			{
			for(i = 0, Betas[2] = 0; i < 3; i++) Betas[2] += NEXB1_SF[i];
			Get_SF_Alphas2(3);
			}					
		}				
	return(0);
}

int Get_SF_Invariants(int OrbitNum)
{
int			A1,
			A2,
			B1,
			B2,
			H1,
			i,
			j,
			k,
			n,
			NumSolns,
			PSum1,
			PSum2,
			PSum3,
			SolV[7],
			SumAR1,
			SumBR1,
			SumAR2,
			SumBR2;
	
	printf("\n\nThe Manifolds of Orbit %4d are:", OrbitNum);
	
	if(NumRelators == 1)
		{
		A1 = Alphas[1];
		A2 = Alphas[2];
		if(A1 > A2)
			{
			i = A1;
			A1 = A2;
			A2 = i;
			}
		FoundSF = TRUE;	
		printf(" SF(0;m/%d,n/%d), 0 < m < %d, 0 < n < %d, gcd(m,%d) = gcd(n,%d) = 1.",
		A1,A2,A1,A2,A1,A2);
		if(B10B11Recognized)
			{
			SFSolV[0] = 7;			
			SFSolV[5] = A1;
			SFSolV[7] = A2;
			}
		return(1);
		}
			
	if(Alphas[1] == 0 || Alphas[2] == 0) return(TOO_LONG);

	if(Alphas[3] == 0)
		{
		A1 = Alphas[1];
		A2 = Alphas[2];
		B1 = Betas[1];
		B2 = Betas[2];
		if(B1) B1 = Recip_Mod_P(A1,B1);
		if(B2) B2 = Recip_Mod_P(A2,B2);
		if(A1 > A2 || (A1 == A2 && B1 > B2))
			{
			i = A1;
			A1 = A2;
			A2 = i;
			j = B1;
			B1 = B2;
			B2 = j;
			}
		if(B1 == 0 && B2 == 0)
			{
			FoundSF = TRUE;
			printf(" Connected sums of lens spaces: L(%d,Q1) # L(%d,Q2).",A1,A2);
			if(B10B11Recognized)
				{
				SFSolV[0] = 8;			
				SFSolV[5] = A1;
				SFSolV[7] = A2;
				}			
			return(2);
			}
		if(B1 == 0)
			{
			FoundSF = TRUE;
			printf(" Connected sums of lens spaces: L(%d,Q) # L(%d,%d).",A1,A2,B2);
			if(B10B11Recognized)
				{
				SFSolV[0] = 9;			
				SFSolV[5] = A1;
				SFSolV[6] = B2;
				SFSolV[7] = A2;
				}			
			return(2);
			}
		if(B2 == 0)
			{
			FoundSF = TRUE;
			printf(" Connected sums of lens spaces: L(%d,%d) # L(%d,Q).",A1,B1,A2);
			if(B10B11Recognized)
				{
				SFSolV[0] = 10;
				SFSolV[4] = B1;			
				SFSolV[5] = A1;
				SFSolV[6] = B2;
				SFSolV[7] = A2;
				}					
			return(2);
			}
		FoundSF = TRUE;	
		printf(" Connected sums of lens spaces: L(%d,%d) # L(%d,%d).",A1,B1,A2,B2);
		if(B10B11Recognized)
			{
			SFSolV[0] = 11;
			SFSolV[4] = B1;			
			SFSolV[5] = A1;
			SFSolV[6] = B2;
			SFSolV[7] = A2;
			}				
		return(2);
		}
		
	for(i = 0, SumAR1 = 0; i < 3; i++)	SumAR1 += NEXA1_SF[i]*EXPA1_SF[i];
	for(i = 0, SumBR1 = 0; i < 3; i++)	SumBR1 += NEXB1_SF[i]*EXPB1_SF[i];
	for(i = 0, SumAR2 = 0; i < 6; i++)	SumAR2 += NEXA2_SF[i]*EXPA2_SF[i];
	for(i = 0, SumBR2 = 0; i < 6; i++)	SumBR2 += NEXB2_SF[i]*EXPB2_SF[i];
		
	H1 = abs(SumAR1*SumBR2-SumAR2*SumBR1);
	
	/*** If Alphas[3] = 1, we have a lens space. ***/
	
	NumSolns = 0;	
	PSum1 = Betas[1]*Alphas[2]*Alphas[3]+Alphas[1]*Betas[2]*Alphas[3];
	PSum2 = H1-PSum1;
	j = Alphas[1]*Alphas[2];
	k = PSum2 % j;
	if(k == 0)
		{
		PSum3 = PSum2/j;
		Betas[3] = PSum3 % Alphas[3];
		if(Betas[3] < 0) Betas[3] = Betas[3] + Alphas[3];
		PSum3 -= Betas[3];	
		n = PSum3/Alphas[3];
		n = -n;	
		if(GCD(Alphas[3],Betas[3]) == 1) 
			{
			SF_Sort_And_Print(H1,n,Alphas[1],Alphas[2],Alphas[3],Betas[1],Betas[2],Betas[3],NumSolns,SolV);		
			NumSolns ++;
			}
		}
		
	if(Alphas[1] > 2 || Alphas[2] > 2)
		{
		PSum1 = (Alphas[1]-Betas[1])*Alphas[2]*Alphas[3]+Alphas[1]*(Alphas[2]-Betas[2])*Alphas[3];
		PSum2 = H1-PSum1;
		j = Alphas[1]*Alphas[2];
		k = PSum2 % j;
		if(k == 0)
			{
			PSum3 = PSum2/j;
			Betas[3] = PSum3 % Alphas[3];
			if(Betas[3] < 0) Betas[3] = Betas[3] + Alphas[3];
			PSum3 -= Betas[3];
			n = PSum3/Alphas[3];
			n = -n;	
			if(GCD(Alphas[3],Betas[3]) == 1)
				{
				if(SF_Sort_And_Print(H1,n,Alphas[1],Alphas[2],Alphas[3],Alphas[1]-Betas[1],Alphas[2]-Betas[2],Betas[3],NumSolns,SolV))
				NumSolns ++;
				}
			}
		}
if(NumSolns == 2) return(5);			
return(0);		
}

int Get_SF_Alphas2(int Num)
{
    /******************************************************************************************
   		The routine Get_SF_Alphas1() can easily determine the values of Alphas 1, 2 and 3
   	when each generator appears in the relators with an exponent of absolute value greater 
   	than 2.
   		The routine Get_SF_Alphas2() is used in the remaining cases when a generator appears
   	appears only with exponents of absolute value less than 3.
    ******************************************************************************************/    
    
    char		*p,
				*q,
    			**Rel = NULL,
    			Sign;
    						
    unsigned char  	x,
					y;

                            
    int			i,
	   			AExpP,
    			AExpQ,
    			AExpR,
    			BExpS,
    			BExpT,
    			BExpU,
    			ex,
    			j,
    			Mex,
				sex;
		
    p = (char *)*Relators[2];
	x = *p++;
	if(!x) return(2);			/* Relators[2] is empty!  */
	
	Rel = (char **) NewHandle(GetHandleSize((char **) Relators[2]));
	if(Rel == NULL) Mem_Error();
	q = *Rel;
	
	if(Num == 2)
		{
		for(i = 0, Mex = 0; i < 3; i++) if(abs(EXPA1_SF[i]) > Mex) Mex = abs(EXPA1_SF[i]);
		Mex--;
		for(i = 0, BExpS = Alphas[2], BExpT = BExpU = 0; i < 6; i++)
			{
			j = abs(EXPB2_SF[i]);
			if(j == 0 || j == BExpS) continue;
			if(j == BExpT) continue;
			if(BExpT == 0) 
				{
				BExpT = j;
				continue;
				}
			if(j == BExpU) continue;
			if(BExpU == 0) BExpU = j;
			}
		}
	if(Num == 3)
		{
		for(i = 0, Mex = 0; i < 3; i++) if(abs(EXPB1_SF[i]) > Mex) Mex = abs(EXPB1_SF[i]);
		Mex--;
		for(i = 0, AExpP = Alphas[1], AExpQ = AExpR = 0; i < 6; i++)
			{
			j = abs(EXPA2_SF[i]);
			if(j == 0 || j == AExpP) continue;
			if(j == AExpQ) continue;
			if(AExpQ == 0) 
				{
				AExpQ = j;
				continue;
				}
			if(j == AExpR) continue;
			if(AExpR == 0) AExpR = j;
			}
		}
				
	ex = 1;
	while(*p == x)
		{
		ex++;
		p++;
		}
	sex = ex;
	ex = 0;
	if(*p) x = *p;
	while( (y = *p++) )
		{
		if(x == y)
			ex++;
		else
			{
			if(Num == 1)
				{
				if(abs(ex) == 1) *q++ = x;
				if(ex == 2)
					{
					if(x == 'A') *q++ = 'C';
					if(x == 'a') *q++ = 'c';					
					if(x == 'B') *q++ = 'D';
					if(x == 'b') *q++ = 'd';
					}
				}
			if(Num == 2)
				{
				if(x == 'A')
					{
					if(ex == Mex) *q++ = 'c';
					else *q++ = 'C';
					}
				if(x == 'a')
					{
					if(ex == Mex) *q++ = 'C';
					else *q++ = 'c';
					}
				if(x == 'B')
					{
					if(ex == BExpS) *q++ = 'S';
					if(ex == BExpT) *q++ = 'T';
					if(ex == BExpU) *q++ = 'U';
					}
				if(x == 'b')
					{
					if(ex == BExpS) *q++ = 's';
					if(ex == BExpT) *q++ = 't';
					if(ex == BExpU) *q++ = 'u';
					}								
				}
			if(Num == 3)
				{
				if(x == 'A')
					{
					if(ex == AExpP) *q++ = 'P';
					if(ex == AExpQ) *q++ = 'Q';
					if(ex == AExpR) *q++ = 'R';
					}
				if(x == 'a')
					{
					if(ex == AExpP) *q++ = 'p';
					if(ex == AExpQ) *q++ = 'q';
					if(ex == AExpR) *q++ = 'r';
					}								
				if(x == 'B')
					{
					if(ex == Mex) *q++ = 'd';
					else *q++ = 'D';
					}
				if(x == 'b')
					{
					if(ex == Mex) *q++ = 'D';
					else *q++ = 'd';
					}	
				}						

			ex = 1;    
			}
		x = y;
		}
	y = **Relators[2];        
	if(x == y)
		{
		ex += sex;	
		if(Num == 1)
			{
			if(abs(ex) == 1) *q++ = x;
			if(ex == 2)
				{
				if(x == 'A') *q++ = 'C';
				if(x == 'a') *q++ = 'c';					
				if(x == 'B') *q++ = 'D';
				if(x == 'b') *q++ = 'd';
				}
			}
		if(Num == 2)
			{
			if(x == 'A')
				{
				if(ex == Mex) *q++ = 'c';
				else *q++ = 'C';
				}
			if(x == 'a')
				{
				if(ex == Mex) *q++ = 'C';
				else *q++ = 'c';
				}
			if(x == 'B')
				{
				if(ex == BExpS) *q++ = 'S';
				if(ex == BExpT) *q++ = 'T';
				if(ex == BExpU) *q++ = 'U';
				}
			if(x == 'b')
				{
				if(ex == BExpS) *q++ = 's';
				if(ex == BExpT) *q++ = 't';
				if(ex == BExpU) *q++ = 'u';
				}									
			}
		if(Num == 3)
			{
			if(x == 'A')
				{
				if(ex == AExpP) *q++ = 'P';
				if(ex == AExpQ) *q++ = 'Q';
				if(ex == AExpR) *q++ = 'R';
				}
			if(x == 'a')
				{
				if(ex == AExpP) *q++ = 'p';
				if(ex == AExpQ) *q++ = 'q';
				if(ex == AExpR) *q++ = 'r';
				}												
			if(x == 'B')
				{
				if(ex == Mex) *q++ = 'd';
				else *q++ = 'D';
				}
			if(x == 'b')
				{
				if(ex == Mex) *q++ = 'D';
				else *q++ = 'd';
				}	
			}												
		}
	else
		{
		if(Num == 1)
			{
			if(abs(ex) == 1) *q++ = x;
			if(ex == 2)
				{
				if(x == 'A') *q++ = 'C';
				if(x == 'a') *q++ = 'c';					
				if(x == 'B') *q++ = 'D';
				if(x == 'b') *q++ = 'd';
				}
			}
		if(Num == 2)
			{
			if(x == 'A')
				{
				if(ex == Mex) *q++ = 'c';
				else *q++ = 'C';
				}
			if(x == 'a')
				{
				if(ex == Mex) *q++ = 'C';
				else *q++ = 'c';
				}
			if(x == 'B')
				{
				if(ex == BExpS) *q++ = 'S';
				if(ex == BExpT) *q++ = 'T';
				if(ex == BExpU) *q++ = 'U';
				}
			if(x == 'b')
				{
				if(ex == BExpS) *q++ = 's';
				if(ex == BExpT) *q++ = 't';
				if(ex == BExpU) *q++ = 'u';
				}												
			}
		if(Num == 3)
			{
			if(x == 'A')
				{
				if(ex == AExpP) *q++ = 'P';
				if(ex == AExpQ) *q++ = 'Q';
				if(ex == AExpR) *q++ = 'R';
				}
			if(x == 'a')
				{
				if(ex == AExpP) *q++ = 'p';
				if(ex == AExpQ) *q++ = 'q';
				if(ex == AExpR) *q++ = 'r';
				}							
			if(x == 'B')
				{
				if(ex == Mex) *q++ = 'd';
				else *q++ = 'D';
				}
			if(x == 'b')
				{
				if(ex == Mex) *q++ = 'D';
				else *q++ = 'd';
				}	
			}					
		
		ex = sex;
		x = y;
		
		if(Num == 1)
			{
			if(abs(ex) == 1) *q++ = x;
			if(ex == 2)
				{
				if(x == 'A') *q++ = 'C';
				if(x == 'a') *q++ = 'c';					
				if(x == 'B') *q++ = 'D';
				if(x == 'b') *q++ = 'd';
				}
			}
		if(Num == 2)
			{
			if(x == 'A')
				{
				if(ex == Mex) *q++ = 'c';
				else *q++ = 'C';
				}
			if(x == 'a')
				{
				if(ex == Mex) *q++ = 'C';
				else *q++ = 'c';
				}
			if(x == 'B')
				{
				if(ex == BExpS) *q++ = 'S';
				if(ex == BExpT) *q++ = 'T';
				if(ex == BExpU) *q++ = 'U';
				}
			if(x == 'b')
				{
				if(ex == BExpS) *q++ = 's';
				if(ex == BExpT) *q++ = 't';
				if(ex == BExpU) *q++ = 'u';
				}									
			}
		if(Num == 3)
			{
			if(x == 'A')
				{
				if(ex == AExpP) *q++ = 'P';
				if(ex == AExpQ) *q++ = 'Q';
				if(ex == AExpR) *q++ = 'R';
				}
			if(x == 'a')
				{
				if(ex == AExpP) *q++ = 'p';
				if(ex == AExpQ) *q++ = 'q';
				if(ex == AExpR) *q++ = 'r';
				}						
			if(x == 'B')
				{
				if(ex == Mex) *q++ = 'd';
				else *q++ = 'D';
				}
			if(x == 'b')
				{
				if(ex == Mex) *q++ = 'D';
				else *q++ = 'd';
				}	
			}								
		}
		
    *q = EOS;
    
    if(Num == 1)
    	{
    	for(p = *Rel,Sign = 1,ex = 0; *p; p++)
    		{
    		if(*p == 'A' || *p == 'a')
    			{
    			if(Sign > 0) ex++;
    			if(Sign < 0) ex--;
    			}
    		if(*p == 'C' || *p == 'c' || *p == 'D' || *p == 'd') Sign = -Sign;	
    		}	
    	Alphas[3] = abs(ex);
    	}
    
    if(Num == 2)
    	{
    	for(p = *Rel, q = p + 1, ex = 0; *p && *q; p++, q++) 
    	if(*p == 'T' || *p == 'U' || *p == 't' || *p == 'u')
    		{
    		if(*q == 'C') ex++;
    		if(*q == 'c') ex--;
    		}
    	q = *Rel;
    	if(*p == 'T' || *p == 'U' || *p == 't' || *p == 'u')
    		{
    		if(*q == 'C') ex++;
    		if(*q == 'c') ex--;
    		}
    	Alphas[3] = abs(ex);	
    	}
    if(Num == 3)
    	{
    	for(p = *Rel, q = p + 1, ex = 0; *p && *q; p++, q++) 
    	if(*p == 'Q' || *p == 'R' || *p == 'q' || *p == 'r')
    		{
    		if(*q == 'D') ex++;
    		if(*q == 'd') ex--;
    		}
    	q = *Rel;
    	if(*p == 'Q' || *p == 'R' || *p == 'q' || *p == 'r')
    		{
    		if(*q == 'D') ex++;
    		if(*q == 'd') ex--;
    		}
    	Alphas[3] = abs(ex);	
    	}    	
 	
 	if(Num == 2 && Alphas[2] > 2) 		/*** Set Betas[2]  ***/
 		{
 		for(p = *Rel, q = p + 1; *p && *q; p++, q++)
 			{
 			if(*p == 'C')
 				{
 				if(*q == 'S' || *q == 's') continue;
 				if(*q == EOS) q = *Rel;
 				if(*q == 'T') j =  BExpT;
 				if(*q == 'U') j =  BExpU; 
 				if(*q == 't') j = -BExpT;
 				if(*q == 'u') j = -BExpU; 
 				}
 			if(*p == 'c')
 				{
 				if(*q == 'S' || *q == 's') continue;
 				if(*q == EOS) q = *Rel;
 				if(*q == 'T') j = -BExpT;
 				if(*q == 'U') j = -BExpU;
 				if(*q == 't') j =  BExpT;
 				if(*q == 'u') j =  BExpU;
 				}	
 			}
 		
 		j = j % Alphas[2];
 		if(j < 0) j += Alphas[2];
 		if(EXPA1_SF[0]*EXPB1_SF[0] > 0) j = Alphas[2] - j;
 		Betas[2] = j;
 		}
 	if(Alphas[2] == 2) Betas[2] = 1;	
 	
 	if(Num == 3 && Alphas[1] > 2) 		/*** Set Betas[1]  ***/
 		{
 		for(p = *Rel, q = p + 1; *p && *q; p++, q++)
 			{
 			if(*p == 'D')
 				{
 				if(*q == 'P' || *q == 'p') continue;
 				if(*q == EOS) q = *Rel;
 				if(*q == 'Q') j =  AExpQ;
 				if(*q == 'R') j =  AExpR; 
 				if(*q == 'q') j = -AExpQ;
 				if(*q == 'r') j = -AExpR; 
 				}
 			if(*p == 'd')
 				{
 				if(*q == 'P' || *q == 'p') continue;
 				if(*q == EOS) q = *Rel;
 				if(*q == 'Q') j = -AExpQ;
 				if(*q == 'R') j = -AExpR;
 				if(*q == 'q') j =  AExpQ;
 				if(*q == 'r') j =  AExpR;
 				}	
 			}
 		
 		j = j % Alphas[1];
 		if(j < 0) j += Alphas[1];
 		if(EXPA1_SF[0]*EXPB1_SF[0] > 0) j = Alphas[1] - j;
 		Betas[1] = j;
 		}
 	if(Alphas[1] == 2) Betas[1] = 1;			
 
    if(Rel != NULL)
    	{
    	DisposeHandle((char **) Rel);
    	Rel = NULL;
    	}
    return(0);            
}

int SF_Sort_And_Print(int H1, int n, int A1, int A2, int A3, int B1, int B2, int B3, int NumSolns, int* SolV)
{
char	Finite;

int		i,
		Q,
		R;

	if(A3 == 1)		/*** We have a lens space. ***/
		{
		i = n*A1 - B1;
		Q = GCD(A1,abs(i));
		if(i < 0)
			Q = A2*Recip_P-B2*Recip_Q;
		else
			Q = A2*Recip_P+B2*Recip_Q;			
		if(Q < 0) Q = -Q;
		if(H1 == 0)
			{
			FoundSF = TRUE;
			printf(" S^1 X S^2");
			if(B10B11Recognized) 
				{
				SFSolV[0] = 12;
				SFSolV[4] = 1;			
				SFSolV[5] = 0;
				}				
			return(0);
			}
		Q = Q % H1;	
		R = H1 - Q;
		GCD(H1,Q);
		if(labs(Recip_Q) < Q) Q = labs(Recip_Q);
		GCD(H1,R);
		if(labs(Recip_Q) < R) R = labs(Recip_Q);
		if(R < Q) Q = R;
		}

	if(A1 > A2)
		{
		i = A1;
		A1 = A2;
		A2 = i;
		i = B1;
		B1 = B2;
		B2 = i;				
		}
	if(A2 > A3)
		{
		i = A2;
		A2 = A3;
		A3 = i;
		i = B2;
		B2 = B3;
		B3 = i;				
		}
	if(A1 > A2)
		{
		i = A1;
		A1 = A2;
		A2 = i;
		i = B1;
		B1 = B2;
		B2 = i;				
		}
	if(A1 == A2 && B1 > B2)
		{
		i = B1;
		B1 = B2;
		B2 = i;
		}
	if(A2 == A3 && B2 > B3)
		{
		i =  B2;
		B2 = B3;
		B3 = i;
		}
	if(A1 == A2 && B1 > B2)
		{
		i = B1;
		B1 = B2;
		B2 = i;
		}
	if(Alphas[3] == 1)  /*** We have a lens space. ***/
		{
		if(H1 == 0)
			{
			FoundSF = TRUE;
			printf(" S^1 X S^2");
			if(B10B11Recognized) 
				{
				SFSolV[0] = 13;
				SFSolV[4] = 1;			
				SFSolV[5] = 0;
				}				
			}
		else
			{
			FoundSF = TRUE;
			FoundFiniteSF = TRUE;
			printf(" Lens spaces: SF(0;%d;%d/%d,%d/%d) = L(%d,%d)",n,B2,A2,B3,A3,H1,Q);
			if(B10B11Finite || B10B11Recognized) 
				{
				SFSolV[0] = 14;
				SFSolV[6] = B2;
				SFSolV[7] = A2;
				SFSolV[8] = B3;
				SFSolV[9] = A3;
				SFSolV[10] = n;
				SFSolV[18] = H1;
				SFSolV[19] = Q;
				}				
			}
		}
	else
		{	
		if(NumSolns == 0)
			{
			SolV[0] = n;
			SolV[1] = B1;
			SolV[2] = A1;
			SolV[3] = B2;
			SolV[4] = A2;
			SolV[5] = B3;
			SolV[6] = A3;
			FoundSF = TRUE;
			printf(" SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",n,B1,A1,B2,A2,B3,A3);
			printf(" SF(0;%d;%d/%d,%d/%d,%d/%d)",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3);
			if(B10B11Recognized)
				{
				SFSolV[0] = 15;
				SFSolV[4] = B1;
				SFSolV[5] = A1;
				SFSolV[6] = B2;
				SFSolV[7] = A2;
				SFSolV[8] = B3;
				SFSolV[9] = A3;
				SFSolV[10] = n;
				}					
			Finite = FALSE;	
			if(A1 == 2)
				{
				switch(A2)
					{
					case 2:
						Finite = TRUE;
						printf(" Finite!");
						break;
					case 3:
						switch(A3)
							{
							case 3:
							case 4:
							case 5:
								Finite = TRUE;
								printf(" Finite!");
								break;
							default:
								break;	
							}
					default:
						break;
					}			
				}
			if(Finite) FoundFiniteSF = TRUE;	
			if(Finite && B10B11Finite) 
				{
				SFSolV[0] = 16;
				SFSolV[4] = B1;
				SFSolV[5] = A1;
				SFSolV[6] = B2;
				SFSolV[7] = A2;
				SFSolV[8] = B3;
				SFSolV[9] = A3;
				SFSolV[10] = n;
				}								
			}						
		else
			{
			if(SolV[0] != n || SolV[1] != B1 || SolV[2] != A1 || SolV[3] != B2 || SolV[4] != A2 || SolV[5] != B3 || SolV[6] != A3 )
				{
				FoundSF = TRUE;
				printf("\n                     or perhaps:");			
				printf(" SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",n,B1,A1,B2,A2,B3,A3);
				printf(" SF(0;%d;%d/%d,%d/%d,%d/%d)",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3);
				if(B10B11Recognized)
					{
					SFSolV[0] = 17;
					SFSolV[11] = B1;
					SFSolV[12] = A1;
					SFSolV[13] = B2;
					SFSolV[14] = A2;
					SFSolV[15] = B3;
					SFSolV[16] = A3;
					SFSolV[17] = n;
					}								
				Finite = FALSE;	
				if(A1 == 2)
					{
					switch(A2)
						{
						case 2:
							Finite = TRUE;
							printf(" Finite!");
							break;
						case 3:
							switch(A3)
								{
								case 3:
								case 4:
								case 5:
									Finite = TRUE;
									printf(" Finite!");
									break;
								default:
									break;	
								}
						default:
							break;
						}			
					}
				if(Finite) FoundFiniteSF = TRUE;	
				if(Finite && B10B11Finite)
					{
					SFSolV[0] = 18;
					SFSolV[11] = B1;
					SFSolV[12] = A1;
					SFSolV[13] = B2;
					SFSolV[14] = A2;
					SFSolV[15] = B3;
					SFSolV[16] = A3;
					SFSolV[17] = n;
					}							
				return(1);
				}
			}
		}
return(0);
}

void Test_Transverse()
{
unsigned char	*ptr = NULL;

unsigned int	Whitehead_Graph();

	NumGenerators = 2;
	NumRelators   = 2;
	Vertices      = 4;
	
	if(Find_Flow_A(NORMAL,FALSE)) return;
	
	Whitehead_Graph();
	
	ptr = (unsigned char*) NewPtr(500);
	if(ptr == NULL) Mem_Error();
	Transverse(ptr);
	
	Print_Relators(Relators,NumRelators);
	printf("\n %s",ptr);
	DisposePtr((char*) ptr);	
}	
