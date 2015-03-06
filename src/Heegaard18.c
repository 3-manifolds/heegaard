#include "Heegaard.h"
#include "Heegaard_Dec.h"
#include <ctype.h>
#include <string.h>

/*******************************************************************************************
	The routines in this file allow Heegaard to find the canonical representative
	presentation in the orbit of a given presentation.
*******************************************************************************************/

/****************************** function prototypes *****************************************
L   38 Find_Canonical_Orbit_Reps(int* MyTable, int F1)
L  891 In_File2(int Test, unsigned char ***MyRelators)
L  986 Save_Pres2(void)
L 1019 Delete_Old_PresentationsSLP(void)
L 1043 Delete_Old_PresentationsSMGP(int MyNumSavedPres,unsigned int* SUR_Num)
L 1064 qksort2(int first, int last, int NumRelators, unsigned int* SUR_Num)
L 1096 qkst_compare2(int i,int j, int NumRelators, unsigned int* SUR_Num)
L 1131 qkst_swap2(i,j)
L 1142 ID_PMQPM(int MyNumSavedPres, char* PMQPML, unsigned int* SUR_Num)
L 1207 ID_A_PMQPM(unsigned int i)
L 1265 Rewrite_Orbit_Reps(int MyNumRelators,int NumOrbits,unsigned int* OrbitNum2SLRNum)
L 1334 MergeHegSpl(unsigned int,unsigned int)
L 1379 Display_HS_Diagrams(int NumHSReps,int* HSRepL)
L 1435 Check_HS_Uniqueness(int NumHSReps,int* HSRepL)
L 1510 Check_HS_Uniqueness_Sub1(int MyHSNum,int MyPresNum)
L 1694 Check_HS_Simple_Circuits(int NumHSReps,int* HSRepL)
L 1791 Find_Simple_Circuits(void)
L 2259 CHSP_Check_Simple_Circuits(unsigned int,int*,int,unsigned char**,unsigned char**)
L 2349 Check_HS_Reps(int NumHSReps,int* HSRepL)
L 2424 Get_Next_Presentation_From_File()
L 2564 Is_IP_In_HS_Reps(int NumHSReps,int* HSRepL)
********************************************************************************************/

int     *Table2 = NULL,
		*Table3 = NULL;

int Find_Canonical_Orbit_Reps(int* MyTable, int F1)
{
char			*PMQPML = NULL,
				ReTry;

unsigned char	*p,
				*q,
				*NewRep = NULL,
				**Temp;				
				
int				HitSum,
				*HitSumL = NULL,
				*HSRepL  = NULL,				
				LoopSum,
				*LoopSumL = NULL,
				i,
				j,
				k,
				l,
				m,
				MissingPres,
				MissingCanonicalRep,
				MultipleSolns,				
				MyMinNumGenerators,
				MyMinNumRelators,
				MyNumSavedPres,
				*NonSFL = NULL,
				NumHSReps,
				NumOrbits,
				NumSFChecked,
				NumSFFound,
				NumSplittings,
				n,
				NumInOrbit;

unsigned int	FOHSNum,
				FOLength,
				*HSL = NULL,
				*HSN = NULL,
				*HegSplNumCount = NULL,
				*HegSplNumFirstOrbit = NULL,
				MyNode,
				MyOrbitSize,
				MyRep,
				*OrbitLength = NULL,
				*OrbitNum2SLRNum = NULL,
				*OrbitSize = NULL,
				SMergers,
				SNumFilled,
				SSNumFilled;
				
	if(TotalComp > 1)
		{
		printf("\n\n Heegaard won't find Canonical Rep Presentations for connected sums!");
		return(TOO_LONG);
		}							
				
	MyMinNumGenerators = NG[MyTable[NumFilled -1]];
	if(MyMinNumGenerators == 0)
		{
		printf("\n\nThe initial presentation was: %s",PresName);
		printf("\n\n Heegaard won't find Canonical Rep Presentations for the trivial group!");
		return(TOO_LONG);
		}
	if(MyMinNumGenerators == 1)
		{
		printf("\n\nThe initial presentation was: %s",PresName);
		printf("\n\n Heegaard won't find Canonical Rep Presentations for one generator groups!");
		return(TOO_LONG);
		}
	MyMinNumRelators = NR[MyTable[NumFilled - 1]];
	if(MyMinNumRelators == 0)
		{
		printf("\n\nThe initial presentation was: %s",PresName);
		printf("\n\n Heegaard won't find Canonical Rep Presentations for Free groups!");
		return(TOO_LONG);
		}
	
	MissingPres = FALSE;
	MissingCanonicalRep = FALSE;
	SNumFilled = NumFilled;

	/******************************************************************************************
		Copy the Presentations of the current component on MyMinNumGenerators into SMGP[ ].
	******************************************************************************************/
	ReTry = FALSE;
_RETRY:
	for(n = NumFilled - 1, MyNumSavedPres = 0; n >= 0; n--)
		{
		if(MyNumSavedPres >= MAX_MIN_GEN_PRES)
			{
			printf("\n\n Find_Canonical_Orbit_Reps has saved %d presentations, which is the current maximal number allowed.",
			MAX_MIN_GEN_PRES);
			printf("\n Will proceed to find Canonical Orbit Representatives of this truncated set of presentations.");
			break;
			}
		ReadPres = MyTable[n];
		NumGenerators = NG[ReadPres];
		if(NumGenerators > MyMinNumGenerators) break;
		if(ComponentNum[ReadPres] > 1) break;    /* Skip presentations with ComponentNums > 1. */
		NumRelators = NR[ReadPres];
		if(F1 != 2 && !ReTry && ID_A_PMQPM(ReadPres) == FALSE) continue;
		if(NumRelators == 0)
			{
			printf("\n\nThe initial presentation was: %s",PresName);
			printf("\n\n Heegaard won't find Canonical Rep Presentations for Free groups!");
			return(TOO_LONG);
			}	
	
		for(i = 1; i <= NumRelators; i++)
			{
			if(SMGP[MyNumSavedPres][i] != NULL) DisposeHandle((char **) SMGP[MyNumSavedPres][i]);
			SMGP[MyNumSavedPres][i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][i]));
			if(SMGP[MyNumSavedPres][i] == NULL) Mem_Error();
			q = *SMGP[MyNumSavedPres][i];
			p = *SUR[ReadPres][i];
			while( (*q++ = *p++) ) ;
			}
		
		SUR_Num[MyNumSavedPres] = ReadPres;
		MyNumSavedPres++;
		}
		
	/*********************************************************************************************** 
		Normally, only presentations on MyMinNumGenerators and MyMinNumRelators marked PM or QPM by 
	Heegaard are used here. If there are no presentations marked PM or QPM, we abandon that 
	requirement and go to RETRY.
	***********************************************************************************************/
		
	if(MyNumSavedPres == 0 && ++ReTry == 1) goto _RETRY;
	if(MyNumSavedPres == 0)
		{
		printf("\n\nThe initial presentation was: %s",PresName);
		printf("\n\nNo presentations of component 1 meet the requirements of the Canonical Rep routine.");
		printf("\nSorry!");
		return(TOO_LONG);
		}
	
	OrbitSize = (unsigned int*) NewPtr((sizeof(int)*(MyNumSavedPres+1)));
	if(OrbitSize == NULL) Mem_Error();
			
	OrbitNum2SLRNum = (unsigned int*) NewPtr((sizeof(int)*(MyNumSavedPres+1)));
	if(OrbitNum2SLRNum == NULL) Mem_Error();
	
	NewRep = (unsigned char*) NewPtr((sizeof(char)*(MyNumSavedPres+1)));
	if(NewRep == NULL) Mem_Error();
		
	PMQPML = (char*) NewPtr((sizeof(char)*(MyNumSavedPres+1)));
	if(PMQPML == NULL) Mem_Error();			

	if(F1 != 2) ID_PMQPM(MyNumSavedPres, PMQPML, SUR_Num);
			
	for(i = 0; i < MyNumSavedPres; i++) 
		{
		OrbitSize[i] = 0;
		NewRep[i] = FALSE;
		BeenChecked[i] = FALSE;
		}

	for(NumOrbits = n = 0; n < MyNumSavedPres; n++) if(BeenChecked[n] == FALSE)
		{
		ReadPres 		= SUR_Num[n];
		NumGenerators 	= NG[ReadPres];
		NumRelators 	= NR[ReadPres];
		Vertices 		= 2*NumGenerators;
		Length 			= SURL[ReadPres];
		SLength 		= Length;
		NumOrbits ++;
	
		/********************************************************************
				Set up the root presentation of the orbit under level-
				transformations.
		********************************************************************/
			
		SLP[0] = (unsigned char ***) NewPtr(sizeof(long)*(NumRelators + 1));
		if(SLP[0] == NULL) Mem_Error();
	
		for(i = 1; i <= NumRelators; i++)
			{
			SLP[0][i] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[n][i]));            
			if(SLP[0][i] == NULL) Mem_Error();
			q = *SLP[0][i];	
			p = *SMGP[n][i];
			while( (*q++ = *p++) ) ;                                    
			}
		NumFilled = 0;    
		ReadPres = NumFilled;    
		NumFilled ++;

		do
			{
			j = Find_Level_Transformations(FALSE,2);
			if(j == NOT_CONNECTED)
				{
				printf("\n\nThe Whitehead Graph of Presentation %u is not connected!",SUR_Num[n]);
				printf("\nHeegaard will only find orbits of presentations with connected Whitehead Graphs.");
				if(OrbitSize) 		DisposePtr((unsigned int*) OrbitSize);
				if(OrbitNum2SLRNum)	DisposePtr((unsigned int*) OrbitNum2SLRNum);
				if(NewRep)			DisposePtr((unsigned char*) NewRep);
				if(PMQPML)			DisposePtr((char*) PMQPML); 
				Delete_Old_PresentationsSLP();
				NumFilled = SNumFilled;
				return(TOO_LONG);
				}
			if(j == FULL_HOUSE)
				{
				printf("\n\nThe orbit of Presentation %u contains at least %u members.",
				SUR_Num[n], MAX_SAVED_PRES - 3);
				printf("\nStopping because the identify of the Orbit's Canonical Rep can't be guaranteed!");
				Compute_Stabilizers = FALSE;
				if(OrbitSize)		DisposePtr((unsigned int*) OrbitSize);
				if(OrbitNum2SLRNum)	DisposePtr((unsigned int*) OrbitNum2SLRNum);
				if(NewRep)			DisposePtr((unsigned char*) NewRep);
				if(PMQPML)			DisposePtr((char*) PMQPML);
				Delete_Old_PresentationsSLP();
				NumFilled = SNumFilled;
				return(TOO_LONG);
				}
			if(j == 5)
				{
				printf("\n\nWe have run out of memory set aside for orbits under level-presentations!!!!");
				Compute_Stabilizers = FALSE;
				if(OrbitSize)		DisposePtr((unsigned int*) OrbitSize);
				if(OrbitNum2SLRNum)	DisposePtr((unsigned int*) OrbitNum2SLRNum);
				if(NewRep)			DisposePtr((unsigned char*) NewRep);
				if(PMQPML)			DisposePtr((char*) PMQPML);
				Delete_Old_PresentationsSLP();
				NumFilled = SNumFilled;
				return(TOO_LONG);
				}
			if(j == TOO_LONG) 
				{
				printf("\n\nNumFilled = %u, Find_Level_Transformations() returned TOO_LONG.", NumFilled);
				if(OrbitSize)		DisposePtr((unsigned int*) OrbitSize);
				if(OrbitNum2SLRNum)	DisposePtr((unsigned int*) OrbitNum2SLRNum);
				if(NewRep)			DisposePtr((unsigned char*) NewRep);
				if(PMQPML)			DisposePtr((char*) PMQPML);
				Delete_Old_PresentationsSLP();
				NumFilled = SNumFilled;
				return(TOO_LONG);
				}  
			ReadPres ++;
			}
		while(ReadPres < NumFilled);
	
    	OrbitSize[n] = NumFilled;
	
		/*****************************************************************
			Find the node of the canonical representative in the orbit.
		*****************************************************************/	
	
		for(MyNode = 0; Left[MyNode] < INFINITE; MyNode = Left[MyNode]) ;
	
		BeenChecked[n] = NumOrbits + MAX_MIN_GEN_PRES;
	
		/******************************************************************
			Check if any current unchecked presentations in SMGP[] lie in 
			the orbit of SMGP[n]. 
		******************************************************************/
	
		for(i = 0; i < MyNumSavedPres; i++)
			{
			ReadPres = SUR_Num[i];
			if((BeenChecked[i] == FALSE) && (NR[ReadPres] == NumRelators) && (SURL[ReadPres] == SLength)) 
			if(In_File2(TRUE, SMGP[i]) < INFINITE) BeenChecked[i] = NumOrbits;
			}
		
		/*******************************************************************
			If MyNode > 0, swap presentation SMGP[n] with SLP[MyNode].
							Set NewRep[n] = TRUE.
		*******************************************************************/
	
		if(MyNode)
			{
			for(i = 1; i <= NumRelators; i++)
				{
				Temp 			= SMGP[n][i];
				SMGP[n][i] 		= SLP[MyNode][i];
				SLP[MyNode][i] 	= Temp;
				}
			NewRep[n] = TRUE;
			MissingCanonicalRep = TRUE;
			}
	
		Delete_Old_PresentationsSLP();
		switch(mykbhit()) 
			{
			case ' ':
				printf("\n	Status: In 'Find_Canonical_Orbit_Reps()'. Processed %d of %d presentations.",
					n,MyNumSavedPres); 
				printf("\n  Hit 'c' to continue. Hit 'q' to abort.");
				GET_RESPONSE2:			
				switch(WaitkbHit())
					{
					case 'c':
						break;
					case 'q':
						{
						if(OrbitSize) 		DisposePtr((unsigned int*) OrbitSize);
						if(OrbitNum2SLRNum) DisposePtr((unsigned int*) OrbitNum2SLRNum);
						if(NewRep) 			DisposePtr((unsigned char*) NewRep);
						if(PMQPML) 			DisposePtr((char*) PMQPML);
						Delete_Old_PresentationsSLP();
						NumFilled = SNumFilled;
						return(INTERRUPT);		
						}
					default: goto GET_RESPONSE2;
					}
				break;						
			case 's':
				printf("\n	Status: In 'Find_Canonical_Orbit_Reps()'. Processed %d of %d presentations.",
				n,MyNumSavedPres); 
				break;
			default:
				break;
			}
		}

	/*******************************************************************
							Report Results.
	*******************************************************************/
	
	if(Batch == FALSE)
		{	
		if(MyNumSavedPres == 1)
			{
			printf("\n\n Heegaard found one Presentation of Component 1 on %d Generators:",
			MyMinNumGenerators);
			if(F1 != 2)
			printf("\n Note xx' indicates Presentation xx is pseudo-minimal or quasi-pseudo-minimal.\n\n");
			}
		else
			{
			printf("\n\n Heegaard found the following %d Presentations of Component 1 on %d Generators:",
			MyNumSavedPres, MyMinNumGenerators);
			if(F1 != 2)
			printf("\n Note xx' indicates Presentation xx is pseudo-minimal or quasi-pseudo-minimal.\n\n");
			}
		printf("{");
		for(i = 0; i < MyNumSavedPres - 1; i++) 
			{
			if(F1 != 2 && PMQPML[i])
				printf("%d',",SUR_Num[i] + 1);
			else		
				printf("%d,",SUR_Num[i] + 1);
			}
		if(F1 != 2 && PMQPML[i])
			printf("%d'}",SUR_Num[i] + 1);
		else		
			printf("%d}",SUR_Num[i] + 1); 
		}
		
	Table2 = (int*) NewPtr((sizeof(int)*NumOrbits));
	if(Table2 == NULL) Mem_Error();
	Table3 = (int*) NewPtr((sizeof(int)*(NumOrbits + 1)));
	if(Table3 == NULL) Mem_Error();
	HitSumL = (int*) NewPtr((sizeof(int)*NumOrbits));
	if(HitSumL == NULL) Mem_Error();
	LoopSumL = (int*) NewPtr((sizeof(int)*NumOrbits));
	if(LoopSumL == NULL) Mem_Error();	

	for(i = j = 0; i < MyNumSavedPres; i++) if(BeenChecked[i] > MAX_MIN_GEN_PRES)
		{
		Table2[j] = BeenChecked[i] - MAX_MIN_GEN_PRES;
		Table3[Table2[j]] = i;
		j++;
		}
	 	
	qksort2(0, NumOrbits, MyMinNumRelators, SUR_Num);

	if(Batch == FALSE)
		{
		if(NumOrbits == 1)
			printf("\n\n These all belong to one orbit.");
		else	 
			printf("\n\n These fall into the following %d orbits under level-transformations:\n",NumOrbits);
		}
		
	HSN 					= (unsigned int*) NewPtr((sizeof(int)*NumOrbits));
	if(HSN 					== NULL) Mem_Error();
	HegSplNumCount 			= (unsigned int*) NewPtr((sizeof(int)*SNumFilled));
	if(HegSplNumCount 		== NULL) Mem_Error();	
	HegSplNumFirstOrbit 	= (unsigned int*) NewPtr((sizeof(int)*SNumFilled));
	if(HegSplNumFirstOrbit 	== NULL) Mem_Error();
	OrbitLength 			= (unsigned int*) NewPtr((sizeof(int)*NumOrbits));
	if(OrbitLength 			== NULL) Mem_Error();	
	for(i = 0; i < SNumFilled; i++) HegSplNumCount[i] = 0;
	SMergers = Mergers;
	
	for(k = NumOrbits - 1; k >= 0; k--)
		{
		j = Table2[k];
		for(i = 0; i < MyNumSavedPres; i++) if(BeenChecked[i] == j + MAX_MIN_GEN_PRES)
			{
			MyOrbitSize = OrbitSize[i];
			OrbitLength[k] = SURL[SUR_Num[i]];
			break;
			}	
		for(i = NumInOrbit = 0; i < MyNumSavedPres; i++) 
			if((BeenChecked[i] == j) || (BeenChecked[i] == j + MAX_MIN_GEN_PRES)) NumInOrbit++;
		if(MyOrbitSize > NumInOrbit) MissingPres = TRUE;	
		for(i = HitSum = LoopSum = 0; i < MyNumSavedPres; i++) if((BeenChecked[i] == j) || (BeenChecked[i] == j + MAX_MIN_GEN_PRES)) 
			{
			HitSum  += SURNumX[SUR_Num[i]];
			LoopSum += NumLoops[SUR_Num[i]];
			}
		HitSumL[k]  = HitSum;
		LoopSumL[k] = LoopSum;		

		for(i = 0; i < MyNumSavedPres; i++) if((BeenChecked[i] == j) || (BeenChecked[i] == j + MAX_MIN_GEN_PRES))
			{
			l = i;
			break;
			}
		for(i = l + 1; i < MyNumSavedPres; i++) if((BeenChecked[i] == j) || (BeenChecked[i] == j + MAX_MIN_GEN_PRES))
			{
			MergeHegSpl(SUR_Num[l],SUR_Num[i]);
			l = i;
			}
		HSN[k] = HegSplNum[SUR_Num[l]];
		HegSplNumCount[HSN[k]] ++; 		/** Counts the number of orbits with the same HSN. **/
		HegSplNumFirstOrbit[HSN[k]] = k;

	if(Batch == FALSE)		
		printf("\nOrbit %4d, Size %6u, Length %6u, OrbHits %6d, OrbLoops %6d, HegSplNum %u, {",
		k + 1, MyOrbitSize, OrbitLength[k],HitSum,LoopSum,HSN[k]);
		
		for(i = 0; i < MyNumSavedPres; i++) if((BeenChecked[i] == j) || (BeenChecked[i] == j + MAX_MIN_GEN_PRES)) 
			{
			NumInOrbit--;
			if(Batch == FALSE)
				{
				if(NumInOrbit)
					{
					if(F1 != 2 && PMQPML[i])
						printf("%d',",SUR_Num[i] + 1);
					else		
						printf("%d,",SUR_Num[i] + 1);
					}
				else
					{
					if(F1 != 2 && PMQPML[i])
						printf("%d'}",SUR_Num[i] + 1);
					else		
						printf("%d}",SUR_Num[i] + 1); 
					}
				}		
			}
		}
	if(Batch == FALSE) printf("\n");
		
	for(k = NumOrbits - 1; k >= 0; k--)
		{
		if(NumOrbits == 1)
			j = BeenChecked[0] - MAX_MIN_GEN_PRES;
		else
			j = Table2[k]; 
		for(i = 0; i < MyNumSavedPres; i++) if(BeenChecked[i] == j + MAX_MIN_GEN_PRES)
			{
			ReadPres 		= SUR_Num[i];
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			MyOrbitSize     = OrbitSize[i];
			if(NewRep[i]) 
				MyRep = INFINITE;
			else 
				MyRep = ReadPres + 1;
			if(Batch == FALSE)
				{
				if(NumRelators == 1)
					{
					if(MyRep == INFINITE)
						printf("\nOrbit %d: Size %u, Canonical Rep Pres ??, Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
						k + 1, MyOrbitSize, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
					else
						{
						if(F1 != 2 && PMQPML[i])
							printf("\nOrbit %d: Size %u, Canonical Rep Pres %u', Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
							k + 1, MyOrbitSize, MyRep, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
						else	
							printf("\nOrbit %d: Size %u, Canonical Rep Pres %u, Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
							k + 1, MyOrbitSize, MyRep, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
						}
					}
				else
					{
					if(MyRep == INFINITE)
						printf("\nOrbit %d: Size %u, Canonical Rep Pres ??, Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
						k + 1, MyOrbitSize, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
					else
						{
						if(F1 != 2 && PMQPML[i])
							printf("\nOrbit %d: Size %u, Canonical Rep Pres %u', Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
							k + 1, MyOrbitSize, MyRep, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
						else	
							printf("\nOrbit %d: Size %u, Canonical Rep Pres %u, Gen %d, Rel %d, Length %lu, OrbHits %d, HegSplNum %u",
							k + 1, MyOrbitSize, MyRep, NumGenerators, NumRelators, Length, HitSumL[k],HSN[k]);
						}
					}
				}	
			OrbitNum2SLRNum[k+1] = i;	
			if(Batch == FALSE) Print_Relators(SMGP[i], NumRelators);
			}
		}

	/* Count the number of Heegaard Splittings. */
	
	for(i = NumSplittings = 0; i < SNumFilled; i++) if(HegSplNumCount[i]) NumSplittings ++;
	HSL = (unsigned int*) NewPtr((sizeof(int)*(NumSplittings)));
	if(HSL == NULL) Mem_Error();
	for(i = NumSplittings = 0; i < SNumFilled; i++) if(HegSplNumCount[i]) 
		HSL[NumSplittings++] = HegSplNumFirstOrbit[i];
	HSRepL = (int*) NewPtr((sizeof(int)*(NumOrbits + NumSplittings)));
	if(HSRepL == NULL) Mem_Error();		

	/* Sort the Heegaard Splittings by FirstOrbitNumbers. */
	
	m = 1;
	do
		{
		for(i = j = 0; i < NumSplittings - m; i++) if(HSL[i] > HSL[i+1])
			{
			k = HSL[i];
			HSL[i] = HSL[i+1];
			HSL[i+1] = k;
			j++;
			}
		m++;	
		}		
	while(j);		
	
	/* For each Heegaard splitting, print the orbit reps of that splitting with minimal length. */
	if(Batch == FALSE)
		{	
		printf("\n******* Below are the minimal genus Heegaard splittings and the minimal length *******");
		printf("\n******* orbit representatives which Heegaard found for each splitting.         *******\n");
		}		
	for(i = NumHSReps = 0; i < NumSplittings; i++)
		{
		FOLength = OrbitLength[HSL[i]];
		FOHSNum  = HSN[HSL[i]];
		for(k = 0, m = 0; k < NumOrbits; k++) if(OrbitLength[k] == FOLength && HSN[k] == FOHSNum)
			{
			if(NumOrbits == 1)
				j = BeenChecked[0] - MAX_MIN_GEN_PRES;
			else
				j = Table2[k]; 
			for(l = 0; l < MyNumSavedPres; l++) if(BeenChecked[l] == j + MAX_MIN_GEN_PRES)
				{
				ReadPres 		= SUR_Num[l];
				NumGenerators 	= NG[ReadPres];
				NumRelators 	= NR[ReadPres];
				Length 			= SURL[ReadPres];
				m++;
				if(m == 1) HSRepL[NumHSReps ++] = -(i+1);
				if(NewRep[l]) HSRepL[NumHSReps ++] = l + INFINITE;
				else HSRepL[NumHSReps ++] = ReadPres;
				printf("\n%s HS %u, P %d, L %lu, Gen %d, Rel %d ",PresName,i+1,m,Length,
					NumGenerators,NumRelators);
				Print_Relators(SMGP[l],NumRelators);	
				if((Batch == 10 || Batch == 11) && H_Results != NULL && B10B11HSReps == TRUE)
					{
					fprintf(H_Results,"\n%s HS %u, P %d, L %lu, Gen %d, Rel %d ",PresName,i+1,m,Length,
						NumGenerators,NumRelators);
					Print_Relators2(SMGP[l],NumRelators);
					}								
				}		
			}
		}
	if(Batch == FALSE || Batch == 53)
		{
		if(Batch == FALSE)
			{
			printf("\n******* Above are the minimal genus Heegaard splittings and the minimal length *******");
			printf("\n******* orbit representatives which Heegaard found for each splitting.         *******\n\n");	
			}
		if(Batch == 53) printf("\n\n");	
		printf("(  HS,  HSNum, FirstOrbit, Length, Orbits, FirstOrbitHits, Total HS Hits, Total HS Loops)\n");
		for(l = 0; l < NumSplittings; l++)
		for(i = 0; i < SNumFilled; i++) if(HegSplNumCount[i] && HegSplNumFirstOrbit[i] == HSL[l]) 
			{
			printf("(%4u, %6d, %10u, %6u, %6u, %14d,",l+1,i,HSL[l]+1,OrbitLength[HSL[l]],HegSplNumCount[i],HitSumL[HSL[l]]);
			for(k = j = m = 0; k < NumOrbits; k++) if(HSN[k] == i) 
				{
				j += HitSumL[k];
				m += LoopSumL[k];
				}
			printf(" %13d, %14d)\n",j,m);	
			}
		if(Batch == FALSE)
			{	
			printf("\n******* The table above gives more detailed info about the splittings Heegaard found. *******\n");
			SSNumFilled = SNumFilled;
			if(PRIM[SNumFilled - 1] == 70 || PRIM[SNumFilled - 1] == 170) SSNumFilled --;
			if(SSNumFilled == 0) SSNumFilled = 1;
			if(BreadthFirstSearch)
				{
				printf("\n Checking this list of candidate HS Reps by copying and pasting it into \042Input_Presentations\042\n");
				printf(" and rerunning individual presentations may reveal that \042Splittings\042 with relatively few hits\n");
				printf(" or HSNum > %u merge with others. In addition, further processing may reveal new or additional\n",FR[SSNumFilled - 1] + 1);
				printf(" minimal length representative presentations for a splitting.");
				}
			else
				{
				printf("\n Checking this list of candidate HS Reps by copying and pasting it into \042Input_Presentations\042\n");
				printf(" and rerunning individual presentations may reveal that \042Splittings\042 with relatively few hits\n");		
				printf(" merge with others. In addition, further processing may reveal new or additional minimal\n");
				printf(" length representative presentations for a splitting.");
				}
			}
		}
	
	if(Batch != 53) Check_HS_Uniqueness(NumHSReps,HSRepL);
	
	printf("\n\nThe initial presentation was: %s",PresName);

	if(Batch == FALSE)
		{			
		if(NumOrbits == 1)
			{
			if(MyNumSavedPres > 1)
				{
				printf("\n\nNote %d presentations on %d generators form one orbit under level-transformations.\n",
				MyNumSavedPres, MyMinNumGenerators);
				printf("These all lie on one Heegaard surface.");
				}
			}
		else
			{
			printf("\n\nNote %d presentations on %d generators form %d orbits under level-transformations.\n",
			MyNumSavedPres, MyMinNumGenerators, NumOrbits);
			if(NumSplittings == 1)
				printf("These all lie on one Heegaard surface.");
			else
				printf("These form %u equivalence classes under bandsums; so may represent %u Heegaard splittings."
				,NumSplittings, NumSplittings);
			}	
	
		i = Mergers - SMergers;	
		if(i == 1)
			printf("\n\n One additional merger was performed.");			
		if(i > 1)
			printf("\n\n %d additional mergers were performed.",i);	
	
		if(MyNumSavedPres >= MAX_MIN_GEN_PRES)
			printf("\n\n Note Heegaard only partitioned the first %d presentations into orbits.",MyNumSavedPres);
		
		printf("\n\n Note xx' indicates Presentation xx is pseudo-minimal or quasi-pseudo-minimal.");	
		if(MissingPres)
			{
			printf("\n\n Note: Initially Heegaard does not look for every presentation in an orbit under");
			printf("\n level-transformations. Hence an orbit may contain presentations not listed.");
			}
		if(MissingCanonicalRep)
			printf("\n In particular, the Canonical Rep Presentation may be missing from an orbit list.");

		if(F1 == 2)
			{
			printf("\n\n Note: The routine Just_Delete_Primitives() deletes primitives from the initial presentation,");
			printf("\n but does not check its realizability. If the initial presentation is not realizable, none of");
			printf("\n Heegaard's behavior and results beyond this point are guaranteed. However, if the initial");
			printf("\n presentation is realizable, all should be well.");
			printf("\n\n HIT 'p' TO PROCEED ANYWAY. HIT ANY OTHER KEY TO ABORT.");

			if(Batch == FALSE)
				{
				switch(WaitkbHit())
					{
					case 'p':
						F1 = 1;
						break;
					default:
						break;        
					}
				}		
			}
		}	
	
	if(F1 == 1)
		{
		/****** Check which genus two presentations are Seifert Fibered. And locate meridian reps. ******/
	
		NonSFL = (int*) NewPtr((sizeof(int)*NumOrbits));
		if(NonSFL == NULL) Mem_Error();

		MultipleSolns = 0;	

		for(k = 0, NumSFChecked = NumSFFound = 0; k < NumOrbits ; k++)
			{
			if(NumOrbits == 1)
				j = BeenChecked[0] - MAX_MIN_GEN_PRES;
			else
				j = Table2[k]; 
			for(i = 0; i < MyNumSavedPres; i++) if(BeenChecked[i] == j + MAX_MIN_GEN_PRES)
				{
				ReadPres 		= SUR_Num[i];
				NumGenerators 	= NG[ReadPres];
				NumRelators 	= NR[ReadPres];
				Vertices		= 2*NumGenerators;
				Length 			= SURL[ReadPres];
				MyOrbitSize     = OrbitSize[i];
				if(NumGenerators > 2) continue;
				if(NumRelators > 2) continue;
				if(NewRep[i] || PMQPML[i] || (k == 0))
					{
					if(NumGenerators == 2 && NumRelators == 1)
						{
						if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
						Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[i][1]));
						if(Relators[1] == NULL) Mem_Error();
						p = *Relators[1];
						q = *SMGP[i][1];
						while( (*p++ = *q++) ) ;
						Length = GetHandleSize((char **) Relators[1]) - 1;
						Genus_Two_Meridian_Reps(k+1,0);
						}
					
				   for(l = 1; l <= NumRelators; l++)
						{
						if(Relators[l] != NULL) DisposeHandle((char **) Relators[l]);
						Relators[l] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[i][l]));
						if(Relators[l] == NULL) Mem_Error();
						p = *Relators[l];
						q = *SMGP[i][l];
						while( (*p++ = *q++) ) ;
						}					
					NumSFChecked ++;
					m = Genus_Two_Seifert_Fibered(k + 1);
					if(m == 13 || m == 14)
						NumSFFound ++;
					else 
						NonSFL[NumSFChecked - NumSFFound - 1] = k + 1;
					if(m == 14) MultipleSolns ++;	
					if(NumSFFound > MultipleSolns) goto FOUND_GOOD_SOLN;
					}
				}
			}	
		}

FOUND_GOOD_SOLN:	
	if(F1 == 1)
		{
		if(NumSFFound == 1)
			{
			if(NumSFChecked == 1)
				printf("\n\n Heegaard checked one orbit rep presentation, which was SF.");			
			else	
				printf("\n\n Heegaard found one SF presentation in the %d orbit rep presentations checked.",NumSFChecked);
			}
		if(NumSFFound > 1)
			{
			printf("\n\n Heegaard found %d SF presentations in the %d orbit rep presentations checked.",NumSFFound,NumSFChecked);
			}		
		if(NumSFChecked > NumSFFound)
			{
			printf("\n\n Orbits of rep presentations checked but not SF: {");
			for(i = 0; i < NumSFChecked - NumSFFound -1; i++) printf("%d,",NonSFL[i]);
			printf("%d}",NonSFL[i]);
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
		if(Batch == FALSE && MyMinNumGenerators == 2 && MyMinNumRelators <= 2)
			{
			printf("\n\n Note 1) Heegaard only checks for Seifert fibrations of manifolds in Orbit 1 and manifolds");
			printf("\n in other orbits for which the 'Canonical Rep Presentation' is marked 'pseudo-minimal' or");
			printf("\n 'quasi-pseudo-minimal' or the 'Canonical Rep Presentation' is missing from the orbit list.");
			printf("\n Note 2) Heegaard recognizes a Seifert manifold M from special features of M's presentation. Since");
			printf("\n only some presentations of M may exhibit these features, Heegaard looks at multiple presentations.");
			}
		}
	
	if(Batch != 53)
		{
		Check_HS_Simple_Circuits(NumHSReps,HSRepL);
		if(CheckHSReps) Check_HS_Reps(NumHSReps,HSRepL);
		}
	else Is_IP_In_HS_Reps(NumHSReps,HSRepL);

	if(Batch == FALSE)
		{
		printf("\n\nDISPLAY DIAGRAMS OF HEEGAARD SPLITTING REPS ? HIT 'y' OR 'n'.");
		GET_RESPONSE5:
		switch(WaitkbHit())
			{
			case 'y':
				printf("\n");
				Display_HS_Diagrams(NumHSReps,HSRepL);
				break;
			case 'n':
				printf("\n");
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE5;
			}
		}
		
	if(Batch == FALSE)
		{	
		printf("\n\nREWRITE SELECTED ORBIT REPS USING EXPONENTS ? HIT 'y' OR 'n'.");
		GET_RESPONSE6:
		switch(WaitkbHit())
			{
			case 'y':
				printf("\n");
				Rewrite_Orbit_Reps(MyMinNumRelators,NumOrbits,OrbitNum2SLRNum);
				break;
			case 'n':
				printf("\n");
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE6;
			}	
		}
		
	/*******************************************************************
				Delete presentations saved in SMGP[ ].
	*******************************************************************/
	
	Delete_Old_PresentationsSMGP(MyNumSavedPres, SUR_Num);
	
	if(Table2) 				DisposePtr((int*) Table2);
	if(Table3) 				DisposePtr((int*) Table3);
	if(HSL) 				DisposePtr((unsigned int*) HSL);
	if(HSN)					DisposePtr((unsigned int*) HSN);
	if(HSRepL)				DisposePtr((int*) HSRepL);
	if(HegSplNumCount) 		DisposePtr((unsigned int*) HegSplNumCount);
	if(HegSplNumFirstOrbit)	DisposePtr((unsigned int*) HegSplNumFirstOrbit);
	if(HitSumL)				DisposePtr((int*) HitSumL);
	if(LoopSumL)			DisposePtr((int*) LoopSumL);
	if(OrbitLength)			DisposePtr((unsigned int*) OrbitLength);
	if(OrbitSize)			DisposePtr((unsigned int*) OrbitSize);
	if(OrbitNum2SLRNum)		DisposePtr((unsigned int*) OrbitNum2SLRNum);
	if(NewRep)				DisposePtr((unsigned char*) NewRep);
	if(PMQPML)				DisposePtr((char*) PMQPML);
	if(F1 == 1) 			DisposePtr((int*) NonSFL);

	NumFilled = SNumFilled;
	
return(0);
}

unsigned int In_File2(int Test, unsigned char ***MyRelators)
{    
    register unsigned char  *p,
    						*q;
    
    unsigned char           *r;
                            
    int                     i,
    						Result;
    
    unsigned int            Node;
    
    Size                    HSP,
    						HSQ;     
    
    if(Test)
    	{
    	for(i = 1; i <= NumRelators; i++) LR[i] = GetHandleSize((char **) MyRelators[i]) - 1;
    	}
    else
    	Canonical_Rewrite(Relators,FALSE,FALSE);

    Node = 0;
    while(1)
        {
        for(i = 1,Result = 0; i <= NumRelators; i++)
            {
         	HSP = LR[i] + 1;
            HSQ = GetHandleSize((char **) SLP[Node][i]);
            if(HSP > HSQ)
                {
                Result = 1;
                break;
                }
            if(HSP < HSQ)
                {
                Result = -1;
                break;
                }
            }
        if(Result == 0)    for(i = 1; i <= NumRelators; i++)
            {
            r = *MyRelators[i] + LR[i];
            *r = 125;
            for(p = *MyRelators[i],q = *SLP[Node][i]; *p == *q; p++,q++) ;
      	    *r = EOS;
            if(*p < *q)
                {
                Result = 1;
                break;
                }
            if(*p > *q)
                {
                Result = -1;
                break;
                }
            }        
        switch(Result)
            {
            case 1:
                if(Left[Node] == INFINITE)
                    {
                    if(Test) return(INFINITE);
                    if(Compute_Stabilizers)
                         printf("  %d -> %d",ReadPres + 1,NumFilled + 1);
                    if(Save_Pres2()) return(TOO_LONG);
                    Left[Node] = NumFilled - 1;
                    Left[NumFilled - 1] = Right[NumFilled - 1] = INFINITE;
                    return(NumFilled - 1);
                    }
                else
                    Node = Left[Node];
                break;        
            case 0:
                if(Compute_Stabilizers)
                     printf("  %d -> %d",ReadPres + 1,Node + 1);
                return(Node);
            case -1:
                if(Right[Node] == INFINITE)
                    {
                    if(Test) return(INFINITE);
                    if(Compute_Stabilizers)
                         printf("  %d -> %d",ReadPres + 1,NumFilled + 1);
                    if(Save_Pres2()) return(TOO_LONG);
                    Right[Node] = NumFilled - 1;
                    Left[NumFilled - 1] = Right[NumFilled - 1] = INFINITE;
                    return(NumFilled - 1);
                    }
                else
                    Node = Right[Node];
                break;            
            }
        }
}

int Save_Pres2(void)
{
    /******************************************************************************************
        Save_Pres2() is called from InFile2() when Heegaard has determined that a 
        presentation should be saved. It saves a copy of the presentation in the array
        SLP[][].
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q;
                    
    register int            i;
                        
	SLP[NumFilled] = (unsigned char ***) NewPtr(sizeof(long)*(NumRelators + 1));
	if(SLP[NumFilled] == NULL) Mem_Error();
			
    for(i = 1; i <= NumRelators; i++)
        {
        SLP[NumFilled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));            
        if(SLP[NumFilled][i] == NULL) Mem_Error();
        q = *SLP[NumFilled][i];    
        p = *Relators[i];
        while( (*q++ = *p++) ) ;                                    
        }

    NumFilled ++;
        
    if(Micro_Print)
        printf("\n\nSaved the current presentation as: Presentation %u.\n",NumFilled);

    return(NO_ERROR);        
}	

void Delete_Old_PresentationsSLP(void)
{
    unsigned int    i,
                    j;
                            
    for(i = 0; i < NumFilled; i++)
        {
        for(j = 1; j <= NumRelators; j++) if(SLP[i][j] != NULL) 
        	{
        	DisposeHandle((char **) SLP[i][j]);
        	SLP[i][j] = NULL;
        	}
        if(SLP[i]) 
        	{
        	DisposePtr((char ***) SLP[i]);
        	SLP[i] = NULL;
        	}
        }
    NumFilled         = 0;
    BytesAvailable    = BYTES_AVAILABLE;
    BytesUsed         = 0L;
    UserSaidQuit      = FALSE;                    
}

void Delete_Old_PresentationsSMGP(int MyNumSavedPres,unsigned int* SUR_Num)
{
    unsigned int    i,
                    j;
                            
    for(i = 0; i < MyNumSavedPres; i++)
        {
        ReadPres = SUR_Num[i];
        NumRelators = NR[ReadPres];
        for(j = 1; j <= NumRelators; j++) if(SMGP[i][j] != NULL) 
        	{
        	DisposeHandle((char **) SMGP[i][j]);
        	SMGP[i][j] = NULL;
        	}
        }
    NumFilled         = 0;
    BytesAvailable    = BYTES_AVAILABLE;
    BytesUsed         = 0L;
    UserSaidQuit      = FALSE;                    
}

void qksort2(int first, int last, int NumRelators, unsigned int* SUR_Num)
{
	int 		i;		/*  "static" to save stack space  */
	int 		j;
	int         qkst_compare2();
    void        qkst_swap2();
 
	while (last - first > 1) 
		{
		i = first;
		j = last;
		for (;;) 
			{
			while (++i < last  && qkst_compare2(i, first, NumRelators, SUR_Num) < 0) 	;
			while (--j > first && qkst_compare2(j, first, NumRelators, SUR_Num) > 0)	;
			if (i >= j)	break;
			qkst_swap2(i, j);
			}
		qkst_swap2(first, j);
		if (j - first < last - (j + 1)) 
			{
			qksort2(first, j, NumRelators, SUR_Num);
			first = j + 1;					/*  qsort2(j + 1, last, NumRelators, SUR_Num);  */
			}
		else 
			{
			qksort2(j + 1, last, NumRelators, SUR_Num);
			last = j;						/*  qsort2(first, j, NumRelators, SUR_Num);  */
			}
		}
}

int qkst_compare2(int i,int j, int NumRelators, unsigned int* SUR_Num)
{
    register unsigned char  	*p,
                            	*q;
                            
    int                        	Pi,
                            	Pj,
                            	k;
              
    Pi = Table3[Table2[i]];
    Pj = Table3[Table2[j]];
    if(Pi == Pj) return(0);

    if(SURL[SUR_Num[Pi]] < SURL[SUR_Num[Pj]]) return(-1);
    if(SURL[SUR_Num[Pi]] > SURL[SUR_Num[Pj]]) return(1);
    for(k = 1; k <= NumRelators; k++)
        {
        if(GetHandleSize((char **) SMGP[Pi][k]) > GetHandleSize((char **) SMGP[Pj][k])) return(-1);
        if(GetHandleSize((char **) SMGP[Pi][k]) < GetHandleSize((char **) SMGP[Pj][k])) return(1);
        }														
    for(k = 1; k <= NumRelators; k++)
        {
        p = *SMGP[Pi][k];
        q = *SMGP[Pj][k];
        while(*p && *p == *q)
            {
            p++;
            q++;
            }
        if(*p < *q) return(-1);
        if(*p > *q) return(1);
        }									
    return(0);
}	

void qkst_swap2(i,j)
int       	i,
            j;
{
    int            Temp;
    
    Temp          = Table2[i];
    Table2[i]     = Table2[j];
    Table2[j]     = Temp;
}

void ID_PMQPM(int MyNumSavedPres, char* PMQPML, unsigned int* SUR_Num)
{
unsigned int	i,
				j,
				n;

for(n = 0; n < MyNumSavedPres; n++) PMQPML[n] = FALSE;
for(n = 0; n < MyNumSavedPres; n++)
	{
	i= SUR_Num[n];
	if(NR[i] == 1) PMQPML[n] = TRUE;
	else
		switch(UDV[i])
			{
			case SPLIT:		
			case GENERIC_LENS_SPACE:        	
			case THREE_SPHERE:	
			case NOT_CONNECTED:		
			case S1_X_S2:	
			case S1_X_D2:	
			case S1_X_X2:
			case MISSING_GEN_DONE2:		
			case MISSING_GEN_DONE1:		
			case KNOWN_LENS_SPACE:
				break;	
			case SEP_PAIRS:
				if(PRIM[i] >= 100) PMQPML[n] = TRUE;	                    
				break;	
			case ANNULUS_EXISTS:	
			case V2_ANNULUS_EXISTS:                	
			case DELETED_RELATOR:	
			case NON_UNIQUE_4:	
			case NON_UNIQUE_3:	
			case NON_UNIQUE_2:	
			case NON_UNIQUE_1:	
			case DUPLICATE:
				break;        	
			default:
				{
				j = PRIM[i];
				switch(j)
					{
					case 8:
					case 108:
						PMQPML[n] = TRUE;
						break;
					case 70:
					case 75:
						if(QPM[i]) PMQPML[n] = TRUE;
						break;    
					case 170:
					case 175:
						PMQPML[n] = TRUE;
						break;    
					default:
						if(j >= 100) PMQPML[n] = TRUE;
						if(QPM[i]) PMQPML[n] = TRUE;    
						break;
					}
				break;
				}                                                                                
			}
	}
}

int ID_A_PMQPM(unsigned int i)
{
unsigned int	j;

	if(NR[i] == 1) return(TRUE);
	else switch(UDV[i])
		{
		case SPLIT:		
		case GENERIC_LENS_SPACE:        	
		case THREE_SPHERE:	
		case NOT_CONNECTED:		
		case S1_X_S2:	
		case S1_X_D2:	
		case S1_X_X2:
		case MISSING_GEN_DONE2:		
		case MISSING_GEN_DONE1:		
		case KNOWN_LENS_SPACE:
			break;	
		case SEP_PAIRS:
			if(PRIM[i] >= 100) return(TRUE);	                    
			break;	
		case ANNULUS_EXISTS:	
		case V2_ANNULUS_EXISTS:                	
		case DELETED_RELATOR:	
		case NON_UNIQUE_4:	
		case NON_UNIQUE_3:	
		case NON_UNIQUE_2:	
		case NON_UNIQUE_1:	
		case DUPLICATE:
			break;        	
		default:
			{
			j = PRIM[i];
			switch(j)
				{
				case 8:
				case 108:
					return(TRUE);
					break;
				case 70:
				case 75:
					if(QPM[i]) return(TRUE);
					break;    
				case 170:
				case 175:
					return(TRUE);
					break;    
				default:
					if(j >= 100) return(TRUE);
					if(QPM[i]) return(TRUE);    
					break;
				}
			break;
			}                                                                                
		}
	return(FALSE);	
}

int Rewrite_Orbit_Reps(int MyNumRelators,int NumOrbits,unsigned int* OrbitNum2SLRNum)
{	
	char			*ptr = NULL;
	
	unsigned char 	***MyRelators = NULL,
					*p,
					x,
					y;

	int				MyOrbit;

	unsigned int	ex,
					i,
					j;
	
	ptr = (char*) NewPtr(sizeof(char)*100);	
    if(ptr == NULL) Mem_Error();

START:             
    printf("\n\nEnter an orbit rep from 1 to %d to rewrite and hit 'return' or enter '0' and hit 'return' to exit. ", 
    	NumOrbits);
GET_RESPONSE1:        
    ReadString((char *)ptr, GetPtrSize(ptr));
    sscanf((char *) ptr,"%d",&MyOrbit);
    if(MyOrbit == 0) 
    	{
    	DisposePtr((char *) ptr);
    	return(0);
    	}
    if(MyOrbit < 1 || MyOrbit > NumOrbits) goto GET_RESPONSE1;
    
    j = OrbitNum2SLRNum[MyOrbit];
    MyRelators = SMGP[j];
    printf("\nRep of Orbit %d:", MyOrbit);
	for(i = 1; i <= MyNumRelators; i++)
		{
		printf("\n");
		p = *MyRelators[i];
		x = *p++;
		if(!x) continue;			/* Relators[i] is empty!  */
		ex = 1;
		while(*p == x)
			{
			ex++;
			p++;
			}
		printf("%c",x);	
		if(ex > 1) printf("^%d",ex);
		ex = 0;
		if(*p) x = *p;
		else continue; 			/* Relators[i] is a proper power!  */
		while( (y = *p++) )
			{
			if(x == y)
				ex++;
			else
				{
				printf("%c",x);
				if(ex > 1) printf("^%d",ex);	
				ex = 1;    
				}
			x = y;
			}
		printf("%c",x);
		if(ex > 1) printf("^%d",ex);
		}
	goto START;	
}

void MergeHegSpl(unsigned int i,unsigned int j)
{ 
unsigned int	ii,
				jj,
				k;
				
	if(NG[i] == NG[j] && ComponentNum[i] == ComponentNum[j])
		{
		if(HegSplNum[i] < HegSplNum[j])
			{
			ii = HegSplNum[i];
			jj = HegSplNum[j];
			k = j;
			while(1)
				{
				if(HegSplNum[k] == ii) break;
				HegSplNum[k] = ii;
				k = HegSplNxt[k];
				}
			ii = HegSplNxt[i];
			jj = HegSplNxt[j];
			HegSplNxt[j] = ii;
			HegSplNxt[i] = jj;	
			Mergers ++;
			}
		if(HegSplNum[i] > HegSplNum[j])
			{
			ii = HegSplNum[i];
			jj = HegSplNum[j];
			k = i;
			while(1)
				{
				if(HegSplNum[k] == jj) break;
				HegSplNum[k] = jj;
				k = HegSplNxt[k];
				}
			ii = HegSplNxt[i];
			jj = HegSplNxt[j];
			HegSplNxt[j] = ii;
			HegSplNxt[i] = jj;
			Mergers ++;
			}			
		}        
}

int Display_HS_Diagrams(int NumHSReps,int* HSRepL)
{
unsigned char	*p,
				*q;
				
int				i,
				j,
				k,
				l,
				m;
	
	printf("\n");
	for(i = 0; i < NumHSReps; i++)
		{
		if(HSRepL[i] < 0) 
			{
			j = -HSRepL[i];
			k = 1;
			continue;
			}
		if(HSRepL[i] >= INFINITE) 
			{
			/* The canonical rep presentation is new and needs to be checked for annuli and uniqueness. */
			printf("\nPres %d of HS %d is new.",k,j);
			l = HSRepL[i] - INFINITE;
			ReadPres 		= SUR_Num[l];
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			Vertices 		= 2*NumGenerators;
			WhichInput		= MAX_SAVED_PRES - 1;
			UDV[WhichInput] = 0;
			
			for(m = 1; m <= NumRelators; m++)
				{
				if(Relators[m] != NULL) DisposeHandle((char **) Relators[m]);
				Relators[m] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[l][m]));
				if(Relators[m] == NULL) Mem_Error();
				p = *SMGP[l][m];
				q = *Relators[m];
				while( (*q++ = *p++) ) ;
				}
			if(Display_A_Diagram(TRUE,k,j) == 2) return(2);	
			k++;		
			}
		else 
			{
			WhichInput = HSRepL[i];
			if(Display_A_Diagram(TRUE,k,j) == 2) return(2);
			k++;
			}
		}
		
	return(0);	
}

void Check_HS_Uniqueness(int NumHSReps,int* HSRepL)
{
unsigned char	*p,
				*q;
				
int				i,
				j,
				k,
				l,
				m,
				NumBad;
				
	printf("\n\n	Making sure each Heegaard Splitting Representative Presentation has been checked...\n");
	
	for(i = NumBad = 0; i < NumHSReps; i++)
		{
		if(HSRepL[i] < 0) 
			{
			j = -HSRepL[i];
			k = 1;
			continue;
			}
		if(HSRepL[i] >= INFINITE) 
			{
			/* The canonical rep presentation is new and needs to be checked for annuli and uniqueness. */
			l 				= HSRepL[i] - INFINITE;
			ReadPres 		= SUR_Num[l];
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			Vertices 		= 2*NumGenerators;
			
			for(m = 1; m <= NumRelators; m++)
				{
				if(Relators[m] != NULL) DisposeHandle((char **) Relators[m]);
				Relators[m] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[l][m]));
				if(Relators[m] == NULL) Mem_Error();
				p = *SMGP[l][m];
				q = *Relators[m];
				while( (*q++ = *p++) ) ;
				}
				
			if(Check_HS_Uniqueness_Sub1(j,k)) NumBad ++;
			k++;		
			}
		else 
			{
			ReadPres 		= HSRepL[i];
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			Vertices 		= 2*NumGenerators;
			
			for(m = 1; m <= NumRelators; m++)
				{
				if(Relators[m] != NULL) DisposeHandle((char **) Relators[m]);
				Relators[m] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][m]));
				if(Relators[m] == NULL) Mem_Error();
				p = *SUR[ReadPres][m];
				q = *Relators[m];
				while( (*q++ = *p++) ) ;
				}
			
			if(Check_HS_Uniqueness_Sub1(j,k)) NumBad ++;
			k++;
			}
		}
		
	if(NumBad == 0) 
		{
		if(Batch == FALSE) printf("\n\n");
		printf("	Each Heegaard Splitting Representative Presentation has a unique realization.");
		}	
}

int Check_HS_Uniqueness_Sub1(int MyHSNum,int MyPresNum)
{
int				i,
				j,
				k,
				LTRV;
				
unsigned int	DMRV;
			
long			HSS;

unsigned int	Diagram_Main();			
	
	if(Length == 0 || NumGenerators == 0 || NumRelators == 0) 
		{
		printf("\nPres %d of Heegaard-Splitting %d is empty!",MyPresNum,MyHSNum);
		return(1);
		}
		
	if(NumGenerators == 1)
		{
		if(NumRelators == 1 && GetHandleSize((char **) Relators[1]) > 4)
			{
			printf("\nPres %d of Heegaard-Splitting %d does not have a unique realization.",MyPresNum,MyHSNum);
			return(1);
			}
		if(NumRelators > 1)
			{
			HSS = GetHandleSize((char **) Relators[1]);
			for(i = 2; i <= NumRelators; i++) if(HSS != GetHandleSize((char **) Relators[i]))
				{
				printf("\nPres %d of Heegaard-Splitting %d is not realizable.",MyPresNum,MyHSNum);
				return(1);
				}
			if(HSS > 4)
				{
				printf("\nPres %d of Heegaard-Splitting %d does not have a unique realization.",MyPresNum,MyHSNum);
				return(1);
				}				
			}	
		}	

	Num_Level_Slides = 0;
	
RETRY:

	Fill_A(NumRelators);
    
    /****************************************************************************************** 
            Set VWG[i] equal to the valence of vertex i in the "reduced" Whitehead graph. 
            Set NumEdges equal to the number of edges in the "reduced" Whitehead graph.
            Use A[][] to setup AJ1[][].
    ******************************************************************************************/
    
    for(i = NumEdges = 0;i < Vertices; i ++)
        {
        A[i][i] = 0;
        for(j = k = 0; j < Vertices; j ++) if(A[i][j])
            {
            AJ1[i][k] = j;
            k++;
            }
        VWG[i] = k;
        NumEdges += k;
        AJ1[i][k] = VERTICES;
        }                            
    NumEdges /= 2;
	
	for(i = 0; i < Vertices; i++) ZZ[i] = 0;
 	if(Connected_(0,0) == FALSE)
 		{
		printf("\nThe Whitehead Graph of Pres %d of Heegaard-Splitting %d is not connected.",MyPresNum,MyHSNum);
		return(1);		
		}
		
	SepPairs = Sep_Pairs(0,0,1);
	
	if(SepPairs)
		{
		printf("\nThe Whitehead Graph of Pres %d of Heegaard-Splitting %d has a separating pair of vertices.",MyPresNum,MyHSNum);
		Num_Level_Slides = 0;
		Num_Saved_LPres  = 0;
		TestRealizability4 = TRUE;
		LTRV = Level_Transformations(0,0,1);
		TestRealizability4 = FALSE;
		switch(LTRV)
			{
			case 0: return(1);
			case 2: goto RETRY;
			case 3:
				{
				printf("\nPres %d of Heegaard-Splitting %d is not realizable.",MyPresNum,MyHSNum);
                return(1);				
				}
			case 13:
				{
				if(Num_Level_Slides == 0) printf("\nHeegaard found an annulus in Pres %d of Heegaard-Splitting %d.",MyPresNum,MyHSNum);
				if(Num_Level_Slides == 1) 
					printf("\nAfter %lu Sep-Vert-Slide, an annulus exists in a presentation obtained from Pres %d of Heegaard-Splitting %d."
					,Num_Level_Slides,MyPresNum,MyHSNum);
				if(Num_Level_Slides > 1)
					printf("\nAfter %lu Sep-Vert-Slides, an annulus exists in a presentation obtained from Pres %d of Heegaard-Splitting %d."
					,Num_Level_Slides,MyPresNum,MyHSNum);	
				printf("\nRerun Pres %d of Heegaard-Splitting %d for details.",MyPresNum,MyHSNum);
				return(1);
				}				
			case 6:
				{
				printf("\nThe search for level-transformations of Pres %d of Heegaard-Splitting %d was interrupted.",MyPresNum,MyHSNum);
				return(1);
				}
			case TOO_LONG:
				{
				printf("\nAn error occurred when checking Pres %d of Heegaard-Splitting %d. Perhaps the presentation is too long.",MyPresNum,MyHSNum);
                return(1);
                }		
			}			
		}
	else
		{
		if(Planar(TRUE,FALSE))
			{
			printf("\nThe Whitehead Graph of Pres %d of Heegaard-Splitting %d is not planar.",MyPresNum,MyHSNum);
			return(1);			
			}
		TestRealizability4 = TRUE;
		DMRV = Diagram_Main();
		TestRealizability4 = FALSE;
		switch(DMRV)
			{
			case NO_ERROR:
				{
				if(Num_Level_Slides == 1)
					printf("\nAfter %lu Sep-Vert-Slide, Pres %d of Heegaard-Splitting %d is uniquely realizable.",
					Num_Level_Slides,MyPresNum,MyHSNum);				
				if(Num_Level_Slides > 1)
					printf("\nAfter %lu Sep-Vert-Slides, Pres %d of Heegaard-Splitting %d is uniquely realizable.",
					Num_Level_Slides,MyPresNum,MyHSNum);
				return(0);
				}
			case NON_UNIQUE_1:
				{
				printf("\nThe diagram of Pres %d of Heegaard-Splitting %d is not unique because there is a generator which",MyPresNum,MyHSNum);
                printf("\nappears with only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
				return(1);
				}
			case NON_UNIQUE_2:
				{
				printf("\nThe diagram of Pres %d of Heegaard-Splitting %d is not unique because there is a generator which",MyPresNum,MyHSNum);
                printf("\nappears only with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
				return(1);
				}
			case NON_UNIQUE_3:
				{
				printf("\nThe diagram of Pres %d of Heegaard-Splitting %d is not unique because there is a generator which",MyPresNum,MyHSNum);
                printf("\nappears only with exponent 5.");
                return(1);
				}
			case NON_UNIQUE_4:
				{
				printf("\nThe diagram of Pres %d of Heegaard-Splitting %d is not unique because there is a generator which",MyPresNum,MyHSNum);
                printf("\nappears with only one exponent and that exponent is greater than 6.");
                return(1);
				}
			case V2_ANNULUS_EXISTS:
				{
				printf("\nThe diagram of Pres %d of Heegaard-Splitting %d is not unique because it has a valence two annulus.",MyPresNum,MyHSNum);
                return(1);				
				}
			case FATAL_ERROR:
				{
				printf("\nPres %d of Heegaard-Splitting %d is not realizable.",MyPresNum,MyHSNum);
                return(1);				
				}	
			case TOO_LONG:
				{
				printf("\nAn error occurred when checking Pres %d of Heegaard-Splitting %d. Perhaps the presentation is too long.",MyPresNum,MyHSNum);
                return(1);				
				}	
			}
		}
	return(0);			
}

int Check_HS_Simple_Circuits(int NumHSReps,int* HSRepL)
{
	unsigned char	Flag,
					*p,
					*q;
				
	int				i,
					j,
					k,
					l,
					m;
	
	printf("\n");
	Flag = FALSE;
	for(i = 0; i < NumHSReps; i++)
		{
		if(HSRepL[i] < 0) 
			{
			j = -HSRepL[i];
			k = 1;
			continue;
			}
		if(HSRepL[i] >= INFINITE)
			{			
			/* The canonical rep presentation is new and needs to be checked for annuli and uniqueness. */
			l = HSRepL[i] - INFINITE;
			ReadPres 		= SUR_Num[l];		
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			Vertices 		= 2*NumGenerators;
			WhichInput		= MAX_SAVED_PRES - 1;
			UDV[WhichInput] = 0;						
			for(m = 1; m <= NumRelators; m++)
				{
				if(Relators[m] != NULL) DisposeHandle((char **) Relators[m]);
				Relators[m] = (unsigned char **) NewHandle(GetHandleSize((char **) SMGP[l][m]));
				if(Relators[m] == NULL) Mem_Error();
				p = *SMGP[l][m];
				q = *Relators[m];
				while( (*q++ = *p++) ) ;
				}
			if(Find_Simple_Circuits()) 
				{
				if(Flag)
					printf(", Pres %d of HS %d",k,j);
				else
					{
					printf("\nNote: There are primitives, proper-powers, or curves not of full rank disjoint from the relators of:");
					printf(" Pres %d of HS %d",k,j);
					Flag = TRUE;
					}
				}
			k++;		
			}
		else 
			{
			ReadPres 		= HSRepL[i];
			NumGenerators 	= NG[ReadPres];
			NumRelators 	= NR[ReadPres];
			Length 			= SURL[ReadPres];
			Vertices 		= 2*NumGenerators;
			WhichInput 		= HSRepL[i];
			
			for(m = 1; m <= NumRelators; m++)
				{
				if(Relators[m] != NULL) DisposeHandle((char **) Relators[m]);
				Relators[m] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][m]));
				if(Relators[m] == NULL) Mem_Error();
				p = *SUR[ReadPres][m];
				q = *Relators[m];
				while( (*q++ = *p++) ) ;
				}
				
			if(Find_Simple_Circuits())
				{
				if(Flag)
					printf(", Pres %d of HS %d",k,j);
				else
					{
					printf("\nNote: There are primitives, proper-powers, or curves not of full rank disjoint from the relators of:");
					printf(" Pres %d of HS %d",k,j);
					Flag = TRUE;
					}
				}
			k++;
			}
		}
		
	if(Flag == FALSE) 
		printf("\nNote: Heegaard found no simple circuits, disjoint from the relators, with less than full rank.");	
	else
		printf("\nCheck paths and simple circuits of diagrams of the HS presentations for details.");		
		
	return(0);		
}

int Find_Simple_Circuits(void)
{
	/********************************************************************************************
		Find_Simple_Circuits() finds each 'simple' circuit disjoint from the relators in the 
		Heegaard diagram and checks if it is primitive, a proper-power, or less than full rank.
	********************************************************************************************/
		
	register unsigned char 	*p,
							*q,
							*r;
							
	register unsigned int 	c,
							d,
							e,
							v,
							w;
							
	unsigned char 			*FacesVisited = NULL,
							*FacesVisitedList = NULL,
							InitialFace,
							NumPathsInCircuit,
							**PM = NULL,
							*PM_From = NULL,
							*PM_To = NULL,
							**PP = NULL,
							*PP_From = NULL,
							*PP_To = NULL,
							PossibleNewTerminalFace,
							*SRelator1 = NULL,
							TerminalFace,
							*T2 = NULL,
							tl,
							tr,
							V1,
							V2,
							V3,
							V4,
							x;
													
	int 					Big_Number = 50000,
							Error,
							ii,
							*PathsInCircuit = NULL,
							**P_From_Face = NULL,
							*pp,
							ss,
							SNumGenerators,
							SNumRelators;
							
	unsigned int 			edge,
							EL[6*VERTICES],
							ee,
							h,
							i,
							j,
							k,
							NumCircuitsFound,
							NumNotFullRank,
							NumPaths,
							vertex,
							vertexRS,
							vertexLE;
	
	long 					HSS = 0,
							length;
	
	unsigned long			max;
	
	Error = 0;
	NumPaths = 0;
	
	Fill_A(NumRelators);				
	if(ComputeValences_A()) return(0);
	Get_Matrix();
	for(i = 0; i < Vertices; i++) ZZ[i] = 0;
	if(Connected_(0,0) == FALSE) return(0);
	if(Sep_Pairs(0,0,1)) return(0);
	if(Planar(FALSE,TRUE) == TRUE) return(0);
	if(Diagram_Main()) return(0);

	PM = (unsigned char **) NewPtr(sizeof(long)*(NumEdges + 1));
	if(PM == NULL) Mem_Error();
	for(d = 1; d <= NumEdges; d++) PM[d] = NULL;	
	PP = (unsigned char **) NewPtr(sizeof(long)*(NumEdges + 1));
	if(PP == NULL) Mem_Error();
	for(d = 1; d <= NumEdges; d++) PP[d] = NULL;
	P_From_Face = (int **) NewPtr(sizeof(long)*(NumFaces + 1));
	if(P_From_Face == NULL) Mem_Error();
	for(d = 1; d <= NumFaces; d++) P_From_Face[d] = NULL;	
	PM_From = (unsigned char *) NewPtr(sizeof(long)*(NumFaces + 1));
	if(PM_From == NULL) Mem_Error();
	PP_From = (unsigned char *) NewPtr(sizeof(long)*(NumFaces + 1));
	if(PP_From == NULL) Mem_Error();
	PM_To = (unsigned char *) NewPtr(sizeof(long)*(NumFaces + 1));
	if(PM_To == NULL) Mem_Error();
	PP_To = (unsigned char *) NewPtr(sizeof(long)*(NumFaces + 1));
	if(PP_To == NULL) Mem_Error();
	PathsInCircuit = (int *) NewPtr(sizeof(int)*(NumFaces + 2));
	if(PathsInCircuit == NULL) Mem_Error();
	FacesVisitedList = (unsigned char *) NewPtr(sizeof(char)*(NumFaces+1));
	if(FacesVisitedList == NULL) Mem_Error();
	FacesVisited = (unsigned char *) NewPtr(sizeof(char)*(NumFaces+1));
	if(FacesVisited == NULL) Mem_Error();

	for(d = 1,max = 0L; d <= NumRelators; d++) if(LR[d] > max) max = LR[d];
	T2 = (unsigned char *) NewPtr(max + 2);
	if(T2 == NULL) Mem_Error();	
	for(d = 1; d <= 2*NumEdges; d++) EL[d] = d;
	NumPaths = 0;
					
	for(ss = 2*NumEdges; ss > 0; ss --)
		{		
		/**************************************************************************************
			Determine which edges of the original diagram are joined by the edge ss of the dual
			diagram.
		**************************************************************************************/	
		
		for(v = c = 0; c + VWG[v] < EL[ss]; v++) c += VWG[v];
		w = EL[ss] - c;
		for(c = ee = 0,vertexLE = d = FV[v]; c < w; c++)
			{
			ee += A[v][d];
			vertexLE = d;
			d = CO[v][d];
			}
		e = ee - 1;
		if(v & 1)
			{	
			e = OSA[v] - e;
			if(e >= V[v]) e -= V[v];
			if(e) e--;
			else e = V[v] - 1;
			}
		
		p = q = *DualRelators[(v >> 1) + 1] + e;
		q++;
		if(!*q) q = *DualRelators[(v >> 1) + 1];
		if(v & 1)
			{
			tr = *p;
			tl = *q;
			if(tr > 95) tr -= 32;
			else tr += 32;			
			if(tl > 95) tl -= 32;
			else tl += 32;			
			}
		else
			{
			tl = *p;
			tr = *q;
			}	
		
		length	= 0L;
		vertex  = v;
		V2 		= vertex;
		edge    = ee % V[v];
		e 		= edge;
		r		= T2;
		
		/**************************************************************************************
			Two relators have come from distinct vertices and are adjacent to each other at
			vertex v. Follow these two relators along until they diverge. The region where
			they are parallel is the "cancellation region".
		**************************************************************************************/	
		
		do
			{
			if(v & 1)
				{
				*r = (v >> 1) + 97;
				w = v - 1;
				}
			else
				{
				*r = (v >> 1) + 65;
				w = v + 1;
				}
			V4 = v;	
			r++;	
			e = OSA[v] - e;
			if(e >= V[v]) e -= V[v];
			length ++;
			if(length > max)
				{			
				Error = 2;
				goto END;
				}	
			v = FV[w];
			d = A[w][v];
			while(d <= e)
				{
				v = CO[w][v];
				d += A[w][v];
				}	
			if(e == (d - 1))	
				{
				/***********************************************************************
					If e = d - 1, then we have found the end of the cancellation region.
				***********************************************************************/
				
				*r++ = EOS;
				if(V4 & 1)
					V4 -= 1;
				else
					V4 += 1;
				vertexRS = v;
				
				/***********************************************************************
					Determine which edge of the dual diagram corresponds to the end of
					the cancellation region and delete it from the list of available
					edges.
				***********************************************************************/	
				
				for(d = 0,c = 1; d < w; d++) c += VWG[d];
				v = FV[w];
				d = A[w][v];
				while(d <= e)
					{
					c++;
					v = CO[w][v];
					d += A[w][v];
					}
				if(EL[c] == c)
					{
					ss --;
					v = EL[ss];
					EL[ss] = c;
					EL[c] = v;
					}	
				else
					{
					for(d = 1; d < ss && (EL[d] != c); d++) ;
					if(d < ss)
						{
						ss --;
						v = EL[ss];
						EL[ss] = c;
						EL[d] = v;
						}
					}
				if(vertexLE & 1)
					x = (vertexLE >> 1) + 97;
				else
					x = (vertexLE >> 1) + 65;
				V1 = vertexLE;
				V3 = vertexRS;
				
				HSS = r - T2;
				NumPaths ++;
				
				PP[NumPaths] = (unsigned char *) NewPtr(HSS);		
				if(PP[NumPaths] == NULL) Mem_Error();
				q = PP[NumPaths];	
				p = T2;
				while( (*q++ = *p++) ) ;

				PM[NumPaths] = (unsigned char *) NewPtr(HSS);		
				if(PM[NumPaths] == NULL) Mem_Error();
				Inverse(PP[NumPaths]);
				q = PM[NumPaths];	
				p = PP[NumPaths];
				while( (*q++ = *p++) ) ;
				Inverse(PP[NumPaths]);	

				for(i = 1; i <= NumFaces; i++)
					{
					FacesVisited[i] = FALSE;
					p = Face[i];
					j = 2;
					while(*p++ < VERTICES) j++;
					P_From_Face[i] = (int *)NewPtr(sizeof(int)*j);
					if(P_From_Face[i] == NULL) Mem_Error();
					P_From_Face[i][0] = 1;
					P_From_Face[i][j-1] = Big_Number;
					for(h = 1; h < j-1; h++) P_From_Face[i][h] = 0;	
					}
			
				
				/***********************************************************************************
							Identify the initial and terminal faces of this path.
				***********************************************************************************/
				for(h = 1,i = FALSE; h <= NumFaces; h++)
					{
					p = Face[h];
                	while((x = *p) < VERTICES) 
                		{
                		if(x == V1)
							{
							q = p;
							q++;
							if(*q == VERTICES) q = Face[h];
							if(*q == V2)
								{
								i = TRUE;
								break;
								}	
							}
						p++;
						}
                	if(i) break;	
					}
					
				for(j = 1,i = FALSE; j <= NumFaces; j++)
					{
					p = Face[j];
                	while((x = *p) < VERTICES) 
                		{
                		if(x == V3)
							{
							q = p;
							q++;
							if(*q == VERTICES) q = Face[j];
							if(*q == V4)
								{
								i = TRUE;
								break;
								}	
							}
						p++;
						}
                	if(i) break;	
					}

				PP_From[NumPaths] = h;
				PP_To[NumPaths]   = j;
				PM_From[NumPaths] = j;
				PM_To[NumPaths]   = h;
						
				break;
				}
			e = B[w][v] - e;
			}	
		while(v != vertex || e != edge);	
	}
	
	DisposePtr((char *) T2);
		
	for(ii = 1; ii <= NumPaths; ii++)
		{
		pp = P_From_Face[PP_From[ii]];
		while(*pp && *pp != Big_Number) pp++;
		*pp = ii;
		pp = P_From_Face[PP_To[ii]];
		while(*pp && *pp != Big_Number) pp++;
		*pp = -ii;	
		}	
		
	/* Save copies of Relators[1], NumGenerators, and NumRelators. */
	
	SNumGenerators = NumGenerators;
	SNumRelators   = NumRelators;
	
	if(Relators[1] != NULL)
		{
		SRelator1 = (unsigned char *) NewPtr(GetHandleSize((char **) Relators[1]));
		if(SRelator1 == NULL) Mem_Error();		
		p = *Relators[1];
		q = SRelator1;
		while((*q++ = *p++)) ;	
		}

	for(i = 1,NumCircuitsFound = NumNotFullRank = 0; i <= NumFaces; i++)
		{
		InitialFace = i;
		TerminalFace = i;
		NumPathsInCircuit = 0;
		for(j = 1; j <= InitialFace; j++) FacesVisited[j] = TRUE;
		FacesVisitedList[1] = InitialFace;
		while(1)
			{
			h = P_From_Face[TerminalFace][0];
			P_From_Face[TerminalFace][0] ++;
			ii = P_From_Face[TerminalFace][h];
			if(ii == Big_Number)
				{
				/*  We've reached a dead end. No new path from the current terminal face leads to a new face.
					Clear info for the current terminal face and reset the terminal face to the previous face. 	*/ 

				P_From_Face[TerminalFace][0] = 1;
				FacesVisited[TerminalFace] = FALSE;
				if(NumPathsInCircuit == 0) break; 
				NumPathsInCircuit --;
				TerminalFace = FacesVisitedList[NumPathsInCircuit + 1];		
				continue;
				}
			if(ii != Big_Number)
				{
				if(ii > 0) PossibleNewTerminalFace = PP_To[ii];
				if(ii < 0) PossibleNewTerminalFace = PM_To[-ii];
				if(PossibleNewTerminalFace == InitialFace)
					{
					if((NumPathsInCircuit == 0) && (ii < 0)) continue;
					if((NumPathsInCircuit > 0) && (ii == -PathsInCircuit[1])) continue;
					if((NumPathsInCircuit > 1) && (TerminalFace < FacesVisitedList[2])) continue;
					if((NumPathsInCircuit == 1) && (abs(ii) < abs(PathsInCircuit[1]))) continue;
					NumPathsInCircuit ++;
					PathsInCircuit[NumPathsInCircuit] = ii;
					NumCircuitsFound ++;										
					k = CHSP_Check_Simple_Circuits(NumCircuitsFound,PathsInCircuit,NumPathsInCircuit,PP,PM);	
					switch(k)
						{
						case 1:
						case 2:
						case 3:
							{
							NumNotFullRank ++;
				/*			printf(" Circuit Faces and Paths: ");
							for(j = 1; j <= NumPathsInCircuit; j++) 
								printf("F%d,P%d,",FacesVisitedList[j],PathsInCircuit[j]);
								printf("F%d", PossibleNewTerminalFace);						*/
							}
						}	
					NumPathsInCircuit --;
					continue;	
					}		
		
				if((PossibleNewTerminalFace > InitialFace) && (FacesVisited[PossibleNewTerminalFace] == FALSE))
					{
					TerminalFace = PossibleNewTerminalFace;
					NumPathsInCircuit ++;
					PathsInCircuit[NumPathsInCircuit] = ii;
					FacesVisited[TerminalFace] = TRUE;
					FacesVisitedList[NumPathsInCircuit + 1] = TerminalFace;
					}
				}			
			}	
		}
		
	/* Restore Relators[1], NumGenerators and NumRelators. */
	
	NumGenerators = SNumGenerators;
	NumRelators   = SNumRelators;
	
	if(SRelator1 != NULL)
		{
		if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
		Relators[1] = (unsigned char **) NewHandle(GetPtrSize((char *) SRelator1));
		if(Relators[1] == NULL) Mem_Error();		
		p = SRelator1;
		q = *Relators[1];
		while((*q++ = *p++)) ;	
		DisposePtr((unsigned char *) SRelator1);
		}
		
END:

	for(j = 1; j <= NumPaths; j++) 
		{
		if(PM[j] != NULL) DisposePtr((char *) PM[j]);
		if(PP[j] != NULL) DisposePtr((char *) PP[j]);		
		}

	for(j = 1; j <= NumFaces; j++) if(P_From_Face[j]) DisposePtr((int *) P_From_Face[j]);
	if(PM != NULL) 					DisposePtr((unsigned char **) PM);
	if(PP != NULL) 					DisposePtr((unsigned char **) PP);	
	if(P_From_Face != NULL) 		DisposePtr((int **) P_From_Face);
	if(PM_From != NULL) 			DisposePtr((unsigned char *)  PM_From);
	if(PP_From != NULL) 			DisposePtr((unsigned char *)  PP_From);
	if(PM_To != NULL) 				DisposePtr((unsigned char *)  PM_To);
	if(PP_To != NULL) 				DisposePtr((unsigned char *)  PP_To);
	if(PathsInCircuit != NULL) 		DisposePtr((int *) PathsInCircuit);
	if(FacesVisitedList != NULL) 	DisposePtr((unsigned char *) FacesVisitedList);
	if(FacesVisited != NULL) 		DisposePtr((unsigned char *) FacesVisited);

if(NumNotFullRank > 0) return(1);	
return(0);	
}

int CHSP_Check_Simple_Circuits(unsigned int NCF,int* My_PathsInCircuit,int NumPaths,unsigned char** MyPP,
	unsigned char** MyPM)
{
	unsigned char	*p,
					*q;
					
	int				i,
					j,
					k,
					SNumGenerators;
					
	unsigned int	C[125];					
	
	unsigned long	HS;	

	for(i = 1,HS = 1L; i <= NumPaths; i++)	
		{
		j   = abs(My_PathsInCircuit[i]);
		HS += GetPtrSize((char *) MyPP[j]);
		}
	HS -= NumPaths;	
	if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
	Relators[1] = (unsigned char **) NewHandle(HS);
	if(Relators[1] == NULL) Mem_Error();
	q = *Relators[1];	
	for(i = 1; i <= NumPaths; i++)
		{
		j = My_PathsInCircuit[i];
		if(j > 0)
			{
			p = MyPP[j];
			while( (*q++ = *p++) ) ;
			q--;
			}
		else
			{
			p = MyPM[-j];
			while( (*q++ = *p++) ) ;
			q--;
			}	
		}
	
	SNumGenerators = NumGenerators;
	NumRelators = 1;
	j = Find_Primitives(2);
	switch(j)
		{
		case 0:
		case 1: 
			/*******************************************************************************************
			 Check if Relators[1] is a defining relator, proper-power, full-rank or less than full-rank.
			********************************************************************************************/	
			for(i = 0; i < SNumGenerators; i++) C[i+'A'] = C[i+'a'] = 0;
			p = *Relators[1];
			while(*p)
				C[*p++]++;
			for(i = j = k = 0; i < SNumGenerators; i++,j+=2)
				{
				C[i+'A'] += C[i+'a'];
				C[j] = C[i+'A'];
				if(C[j] == 1) 
					{
					NumGenerators = SNumGenerators;
					return(1);
					}
				if(C[j]) k++;
				}
			if(k == 1) 
				{
				NumGenerators = SNumGenerators;
				return(2);				
				}
			if(k < SNumGenerators)
				{
				NumGenerators = SNumGenerators;
				return(3);				
				}
			if(k == SNumGenerators)
				{
				NumGenerators = SNumGenerators;
				return(0);				
				}			
			break;
		case TOO_LONG:
			break;
		}
		
	return(0);	
}

int Check_HS_Reps(int NumHSReps,int* HSRepL)
{
	unsigned char   *p,
                    *q;
    
    unsigned int    i,
                    k,
                    l,
                    m; 
    
    /******************************************************************************************
    				  Compare HS P1 with that saved in Copy_Of_Input[].
    ******************************************************************************************/
    
    CheckHSReps = FALSE;
    
 	for(i = 0; i < NumHSReps; i++)
		{		
		if(HSRepL[i] < 0) 
			{
			k = 1;
			continue;
			}
		if((i > 1) && Get_Next_Presentation_From_File()) return(1);
		if(HSRepL[i] >= INFINITE)
			{
			/* The canonical rep presentation is new. */
			l = HSRepL[i] - INFINITE;
			ReadPres = SUR_Num[l];
			if(NumRelators != NR[ReadPres]) return(1);						
			for(m = 1; m <= NumRelators; m++)
				{
				p = *SMGP[l][m];
				if(i == 1)
					q = *Copy_Of_Input[m];
				else
					q = *Relators[m];
				while(*p && (*q++ == *p++) ) ;
				if((*p == EOS) && (*q == EOS)) continue;
				else
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\n\nNote: Presentation HS 1 P %d does not match that in Input_Presentations!",k);
					return(1);
					}
				}
			k++;		
			}
		else 
			{
			ReadPres = HSRepL[i];
			if(NumRelators != NR[ReadPres]) return(1);			
			for(m = 1; m <= NumRelators; m++)
				{
				p = *SUR[ReadPres][m];
				if(i == 1)
					q = *Copy_Of_Input[m];
				else
					q = *Relators[m];				
				while(*p && (*q++ == *p++) ) ;
				if((*p == EOS) && (*q == EOS)) continue;
				else
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\n\nNote: Presentation HS 1 P %d does not match that in Input_Presentations!",k);
					return(1);
					}
				}
			k++;
			}
		}
	printf("\n\nThe new HS Reps match the original HS Reps in Input_Presentations!");
	return(0);
}

int Get_Next_Presentation_From_File()
{
    register unsigned char  *p,
                            *q,
                            t;                                     
    
    unsigned int            h,
                            i,
                            j;
                            
    long                    StrLength;

TOP:      

	/******************************************************************************************
	Look for the next nonempty line. This should be the identifier of the next Pres.
	******************************************************************************************/	
	
	do
		if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL) return(1);
	while(*Inst == '\n' || *Inst == '\r');

	p = Inst;
	while(1)
		{
		t = *p;
        if(t == '\n' || t == EOS || t == '\r')
            {
            *p = EOS;
            break;
            }
        p++;
        } 
        
    q = Inst;
    if((p - q) >= MAXLENGTH) goto TOP;
    
    /******************************************************************************************
    		Check that Inst is not just a string of '-'s or a string of '*'s.
    ******************************************************************************************/ 
    
    p = Inst;
    while(1)
    	{
    	if(*p == EOS) break;
    	if(*p != '-') break;
    	p++;
    	}
    if(*p == EOS) goto TOP;
    p = Inst;
    while(1)
    	{
    	if(*p == EOS) break;
    	if(*p != '*') break;
    	p++;
    	}
    if(*p == EOS) goto TOP;	
    q = Inst;  
    p = PresName;
    while((*p++ = *q++)) ;    
    
    /******************************************************************************************
    	Then look for the next following nonempty line. This should start the next Pres.
    ******************************************************************************************/   
    do
        if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL) return(2);
    while(*Inst == '\n' || *Inst == '\r');

    /******************************************************************************************
        Read in the relators, at one relator to each nonempty line, stripping off leading 
        spaces and tabs.
    ******************************************************************************************/
                            
    for(i = 1,NumRelators = 0; i <= MAXNUMRELATORS; i++)
        {
 		p = Inst;
        t = *p;
        h = 0;
        while(t == ' ' || t == '\t')
            {
            h++;
            p++;
            t = *p;
            }
        if(t == '\n' || t == EOS)
            {
            *p = EOS;
            break;
            }
        j = 0;       
        while( (t = *p) )
            {
            if(t == '\n' || t == ' ' || t == '\t')
                {
                *p = EOS;
                break;
                }
            j++;
            if(!isalpha(t))
            	{
				/******************************************************************************
					The current line is not a relator. Look for the next blank line.
				******************************************************************************/
				
				do
        			if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL) return(3);
    			while(*Inst != '\n' && *Inst != '\r');
    			
    			/******** Found a blank line. Go to the beginning. ***********/
    			
				goto TOP;
				}            
            p++;        
            }                                
  		if(*Inst== '\n') break;
        StrLength = strlen((char *) Inst);
        if(StrLength >= MAXLENGTH) return(4);
        if(StrLength == h) break;
        NumRelators ++;
        if(NumRelators == 1)
        	{
    		printf("\n\n------------------------------------");        	
        	printf("\n\n%s",PresName);
        	}
        if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
        Relators[i] = (unsigned char **) NewHandle(StrLength + 1 - h);
        if(Relators[i] == NULL) Mem_Error();
        LR[i] = StrLength - h;
        p = *Relators[i];    
  		q = Inst + h;
        while( (*p++ = *q++) ) ;
        if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL) break;         
        }
    if(NumRelators == MAXNUMRELATORS) return(5);
    
    Print_Relators(Relators, NumRelators);
    
    return(0);
}

int Is_IP_In_HS_Reps(int NumHSReps,int* HSRepL)
{
	unsigned char   Match,
					*p,
                    *q;
    
    unsigned int    i,
    				j,
                    k,
                    l,
                    m; 
    
    /******************************************************************************************
    				  See if the IP appears in the list of HS Reps.
    ******************************************************************************************/
    
    CheckHSReps = FALSE;
    
 	for(i = j = 0,Match = FALSE; i < NumHSReps; i++)
		{		
		if(HSRepL[i] < 0) 
			{
			k = 1;
			j++;
			continue;
			}
		if(HSRepL[i] >= INFINITE)
			{
			/* The canonical rep presentation is new. */
			l = HSRepL[i] - INFINITE;
			ReadPres = SUR_Num[l];
			if(CopyNumRelators != NR[ReadPres]) continue;						
			for(m = 1; m <= NumRelators; m++)
				{
				p = *SMGP[l][m];
				q = *Copy_Of_Input[m];
				while(*p && (*q++ == *p++) ) ;
				if((*p == EOS) && (*q == EOS)) continue;
				break;
				}
			if(m > NumRelators)	
				{
				Match = TRUE;
				break;
				}
			k++;		
			}
		else 
			{
			ReadPres = HSRepL[i];
			if(CopyNumRelators != NR[ReadPres]) continue;			
			for(m = 1; m <= NumRelators; m++)
				{
				p = *SUR[ReadPres][m];
				q = *Copy_Of_Input[m];
				while(*p && (*q++ == *p++) ) ;
				if((*p == EOS) && (*q == EOS)) continue;
				break;
				}
			if(m > NumRelators)	
				{
				Match = TRUE;
				break;
				}	
			k++;
			}
		}
		
	if(Match) printf("\n\nThe IP appears in the HS_Rep List as HS %d, P %d.",j,k);
	else printf("\n\n%s <-- Not a HS Rep!",PresName);
	if(Batch == 53 && H_Results != NULL) 
		{
		if(Match) fprintf(H_Results,"\n\n%s appears as HS %d, P %d",PresName,j,k);
		else fprintf(H_Results,"\n\n%s <-- Not a HS Rep!",PresName);
		}
	return(0);
}