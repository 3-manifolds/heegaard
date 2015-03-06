#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   31 Level_Transformations(int F1,int F2,int F3)
L  275 Get_Components(void)
L  322 Been_Seen(void)
L  353 Annulus(unsigned int V3,unsigned int V4,int TheComp,int F1)
L  511 Find_Path(int V1,int V2,int TheComp,int SepType)
L  692 Do_Aut_L(void)
L  802 Do_Aut_L_Long(int V1,int V2,int TheComp)
L 1089 Delete_Trivial_Generators(int Test)
L 1236 Level_Transformations_2(int F1,int F2,int F3,unsigned int MyNum)
L 1320 Do_Aut_L_Long_2(int V1,int V2,int TheComp)
L 1611 Slide_ValenceTwo_Comp(int TheComp,unsigned int VL,unsigned int VR,unsigned int MyNum)
L 1701 Level_Trans_Reset(unsigned int MyNum, unsigned int V3, unsigned int V4)
L 1741 Micro_Print_Level_Transformations_Reset(unsigned int MyNum)
L 1747 Count_Sep_Pairs(unsigned int Num_Saved_LPres)
L 1835 Is_Sep_Pair(unsigned int VL,unsigned int VR,unsigned int MyNum)
L 1873 Micro_Print_Level_Transformations(unsigned int TheComp,unsigned int V1,unsigned int V2,
	unsigned int Type)
L 1920 Slide_LComp(unsigned int VX, unsigned int VY, int TheComp, int SepType, unsigned int MyNum,
    int F1, int F2, int F3, int F5)
L 1994 Random_Sep_Pair(unsigned int WhichSLRPres)  	
********************************************************************************************/

#define MAX_NOT_NEW_PRES	50 

int		Found_L_Annulus;

int Level_Transformations(int F1,int F2,int F3)
{
	/******************************************************************************************
		This routine is called when the "reduced" Whitehead graph has a separating pair of
		vertices and hence does not have a unique embedding in the plane. Since the
		connectivity of the "reduced" Whitehead graph of a set of relators is not invariant
		under level transformations, it is sometimes possible to perform some length
		preserving i.e. level transformations on the relators and to obtain a modified set of
		relators whose associated Whitehead graph has no separating vertices. This routine
		calls itself recursively in a search for such a modified set of relators.
			Setting F1 FALSE suppresses printing of the message which the routine normally
		puts on the screen when it discovers an annulus. 
			If F2 is TRUE, when the routine finds the first new presentation without any pairs
		of separating vertices, it stops, saves this presentation, and returns 2.
			If F3 is TRUE, when the routine finds the first new presentation without any pairs
		of separating vertices, it stops and returns 2  -- without bothering to save the new
		presentation.
			If Level_Transformations() finds more than MAX_NOT_NEW_PRES consecutive
		presentations, all of which are on file, then the routine assumes it is retracing old
		ground and further searching will probably not yield anything new, so it stops and
		returns 6.
	******************************************************************************************/

	char			MyChar;
	
	unsigned char 	*p,
					*q;

	int 			k,
					MyNumSepComps,
					TheComp;
					
	unsigned int 	i,
					j,
					MyNum,
					Random,
					V3,
					V4,
					VL,
					VR,
					VLI,
					VRI;
	
	/******************************************************************************************
		Each invocation of Level_Transformations() "owns" a presentation. On entry, the
		routine first checks that the presentation being passed to it is distinct from those
		passed to all previous invocations. The routine Been_Seen() is called to check this.
		We also check at this point whether the previous invocation of Level_Transformations()
		has managed to get rid of all separating pairs of vertices in the graph. If so, we save
		this new modified presentation and return 1 as a signal of success. Otherwise we save
		a copy of our private presentation on the stack of presentations SLR[][] at the
		locations SLR[MyNum][1] ... SLR[MyNum][NumRelators] and start looking for level
		transformations.
	******************************************************************************************/	
	
	if(Num_Saved_LPres >= MAX_SAVED_LEVELS) return(7);
	k = Been_Seen();
	if(k < Num_Saved_LPres)
		{
		if(Micro_Print) printf("\n\nPresentation L%u is a duplicate of Presentation %d, and will be deleted.",
			Num_Saved_LPres + 1, k);
		return(0);
		}
	if((MyChar = mykbhit()) && !TestRealizability1)
		{
		if(MyChar == ' ')
			{
			Level_Interrupt = 1;
			return(6);
			}
		Level_Interrupt = 2;	
		}
	MyNum = Num_Saved_LPres;
	if(Num_Saved_LPres == 0) Found_L_Annulus = FALSE;
	Fill_A(NumRelators);
	Get_Matrix();	
	switch(Sep_Pairs(0,0,1))
		{
		case 0:
			if(F3) return(2);
			if(Num_Saved_LPres)
				{
				This_Pres = On_File();	
				if(This_Pres == NumFilled)
					{
					if((Batch == FALSE) && (BytesUsed > BytesAvailable))
						{
						if(UserSaidQuit) return(6);
						if( (UserSaidQuit = User_Says_Quit()) ) return(6);
						}
					NotNewPres = 0;
					if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);
					if(Dup_On_File < INFINITE)
						{
						This_Pres = Dup_On_File;
						if(Save_Pres(ReadPres,Dup_On_File,Length,1,70,1,0,0)) return(TOO_LONG);
						Mark_As_Duplicate(Dup_On_File);							
						}
					else
						{
						if(Save_Pres(ReadPres,0,Length,1,70,1,0,0)) return(TOO_LONG);							
						BDY[NumFilled - 1] = BDY[ReadPres];
						UDV[NumFilled - 1] = 0;
						if(F2) return(2);
						}
					}
				else
					if(++NotNewPres > MAX_NOT_NEW_PRES) return(6);
				}
			return(2);
		case 1:
			break;
		}

	/******************************************************************************************
			If a minimal length presentation P is realizable, then P and any presentation P'
		obtained from P by level-transformations must have a planar Whitehead graph, even if P
		and P' have separating pairs of vertices. 
			So we may be able to verify that P is not realizable by checking every P' obtained
		from P by level-transformations for planarity, even when we cannot find a presentation
		without any separating pairs of vertices.
	******************************************************************************************/

	if((TestRealizability1 || TestRealizability4) && Planar(TRUE,FALSE)) return(3);

	/******************************************************************************************
			At this point, Been_Seen() has checked that the incoming presentation in Relators[] 
		has a separating pair of vertices and is not on file in SLR[]. So save a copy of this
		incoming presentation in SLR[] for future reference.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		if(SLR[Num_Saved_LPres][i] != NULL) DisposeHandle((char **) SLR[Num_Saved_LPres][i]);
		SLR[Num_Saved_LPres][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(SLR[Num_Saved_LPres][i] == NULL) Mem_Error();
		q = *SLR[Num_Saved_LPres][i];	
		p = *Relators[i];
		while( (*q++ = *p++) ) ;
		}
	
	Num_Saved_LPres ++;
	
	/*******************************************************************************************
			Next, call Count_Sep_Pairs(). It counts the total number of separating pairs of 
		vertices in the reduced Whitehead graph RWG of the incoming presentation, and saves a
		list of these in the array SLR[Num__Saved_Pres - 1][0]. Then call Random_Sep-Pair() to
		select a random separating pair of vertices to process.
	*******************************************************************************************/
	
	j = Count_Sep_Pairs(Num_Saved_LPres);
	if(j == TOO_LONG) return(TOO_LONG);

	if(Random_Sep_Pair(MyNum)) return(0);
	BSV1[MyNum] = V1;
	BSV2[MyNum] = V2;
		
	/******************************************************************************************
		VL and VR are the current separating pair of vertices. We look for a component
		"TheComp", of the separation produced by deleting VL and VR, with the property that
		neither VL inverse or VR inverse are in TheComp. The objective is to "slide" the
		component TheComp around the Whitehead graph sliding first to the "left" and then
		(if necessary) to the "right" until either TheComp is attached to the main body of
		the graph at more than two distinct vertices, or we arrive at a vertex of attachment
		of TheComp with the property that its inverse vertex is a member of TheComp, or we
		manage to slide TheComp completely around the graph along some path and the
		presentation has returned to the original presentation. If the last case happens, then
		there is an annulus present. Otherwise we call Level_Transformations() recursively and
		pass it the new presentation.
	******************************************************************************************/
									
	while(1)
		{
		for(i = 0; i < Vertices; i++) XX[i] = 0;
		VL = V3 = BSV1[MyNum];
		VR = V4 = BSV2[MyNum];
		XX[VL] = VERTICES;
		XX[VR] = VERTICES;
		if(VL & 1)
			VLI = VL - 1;
		else
			VLI = VL + 1;		
		if(VR & 1)
			VRI = VR - 1;
		else
			VRI = VR + 1;		
		Get_Components();
		MyNumSepComps = NumSepComps;
		
		for(TheComp = 1; TheComp <= MyNumSepComps; TheComp++)
			{
			if(RandomizeSlides)	Random = abs(rand()) %2;
			else Random = 1;
			
			if(Random == 1) 
				{
				if(XX[VLI] != TheComp && XX[VRI] != TheComp)
					{
					j = Slide_LComp(VL,VR,TheComp,1,MyNum,F1,F2,F3);
					if(j > 1) return(j);
					if(Level_Trans_Reset(MyNum,V3,V4) == TOO_LONG) return(TOO_LONG);
					j = Slide_LComp(VR,VL,TheComp,1,MyNum,F1,F2,F3);
					if(j > 1) return(j);
					if(Level_Trans_Reset(MyNum,V3,V4) == TOO_LONG) return(TOO_LONG);					
					}
				}
			
			if(Random == 0)
				{
				if(XX[VLI] != TheComp && XX[VRI] != TheComp)
					{
					j = Slide_LComp(VR,VL,TheComp,1,MyNum,F1,F2,F3);
					if(j > 1) return(j);
					if(Level_Trans_Reset(MyNum,V3,V4) == TOO_LONG) return(TOO_LONG);
					j = Slide_LComp(VL,VR,TheComp,1,MyNum,F1,F2,F3);
					if(j > 1) return(j);
					if(Level_Trans_Reset(MyNum,V3,V4) == TOO_LONG) return(TOO_LONG);
					}
				}
			}	
		
		if(MyNumSepComps == 2) 
			{
			j = Level_Transformations_2(F1,F2,F3,MyNum);
			if(j > 1) return(j);
			}
		
		/*****************************************************************************************
			Ask Random_Sep_Pair() to select a random pair of unprocessed separating vertices 
			(V1,V2) from the set of separating pairs of vertices of the RWG of Presentation 
			MyNum. If all separating pairs of vertices have been processed, Random_Sep_Pair()
			returns 1, and we return(0). Otherwise, Random_Sep_Pair() returns 0, and we continue.
		******************************************************************************************/				
		
		if(Random_Sep_Pair(MyNum))
			{
			if(MyNum == 0)	Print_SLR(0,Found_L_Annulus);
			return(0);
			}
		BSV1[MyNum] = V1;
		BSV2[MyNum] = V2;	
		}
}

void Get_Components()
{
	/******************************************************************************************
		If the "reduced" Whitehead graph has a separating pair of vertices, this routine can
		be used to find all of the components resulting from the deletion of the separating
		vertices. The array XX[] is initialized by the calling routine, which sets all of the
		entries to zero, except for the two entries corresponding to the separating vertices.
		These entries are set to the value VERTICES. Upon return, the entries of XX[] which
		are equal to g are the vertices in the gth component of the graph.
	******************************************************************************************/	
	
	register unsigned int 	g,
							h,
							i,
							j,
							k;
							
	register unsigned int 	*p,
							*r;
	
	k = 2;
	g = 1;
	while(k < Vertices)
		{
		for(i = 0; i < Vertices && XX[i]; i++) ;
		XX[i] = g;
		if(++k >= Vertices) goto FOUND_COMPS;
		for(r = UpDate,*r = i,p = r + 1; r < p; r++)
			{
			i = *r;
			for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
				{
				if(i == j) continue;
				if(XX[j] == 0)
					{
					XX[j] = g;
					*p++ = j;
					if(++k >= Vertices) goto FOUND_COMPS;
					}
				}
			}
		if(k < Vertices) g++;		
		}
	FOUND_COMPS:	
	NumSepComps = g;					
}								

int Been_Seen()
{
	/******************************************************************************************
			This routine compares the presentation, which we are about to pass to
			Level_Transformations(), to all of the presentations passed in previous calls.
			If there is a duplication, then we do not need or want to make this call.
			If there is a match the routine returns the number of the matched presentation. 
			Otherwise it returns Num_Saved_LPres.
	******************************************************************************************/
		
	register int 	i,
					j;
	
	for(i = 0; i < Num_Saved_LPres; i++)
		{
		for(j = 1; j <= NumRelators; j++)
			{
			if(GetHandleSize((char **) Relators[j]) == GetHandleSize((char **) SLR[i][j]))
				{
				if(Compare_Str(*Relators[j],*SLR[i][j],
					GetHandleSize((char **) Relators[j]) - 1) == FALSE)
					break;
				}	
			else
				break;			
			}			
		if(j > NumRelators) break;							
		}
	return(i);		
}

int Annulus(unsigned int V3,unsigned int V4,int TheComp,int F1)
{
	/******************************************************************************************
		This routine is called by the subroutine Slide_LComp() of Level_Transformations() 
		when Slide_LComp() has determined that there is an annulus present. V3 and V4 are a 
		pair of vertices that separate the Heegaard diagram. TheComp gives us a set of 
		vertices which are "swallowed" by the annulus; while Temp8 contains a string which 
		represents the path on the Heegaard surface that serves as the center line of the 
		rest of the annulus.
	******************************************************************************************/
		
	unsigned char  *p,
					*ptr,
					*q,
					*r,
					x,
					y;
							
	unsigned int 	i,
					j;
							
	int				Separating;
	
	unsigned long	HS;
			
	if(V3 & 1)
		x = V3/2 + 97;
	else
		x = V3/2 + 65;
	if(V4 & 1)
		y = V4/2 + 97;
	else
		y = V4/2 + 65;	
	
	ptr = (unsigned char*) NewPtr((sizeof(char)*(GetHandleSize((char **) Temp8) + Vertices + 2)));
	if(ptr == NULL) Mem_Error();
	q = r = ptr;
	*q++ = x;
	*q++ = y;
	for(i = 0,Separating = TRUE; i < Vertices; i++)
		{
		if(XX[i] == TheComp)
			{
			if(i & 1)
				{
				if(XX[i-1] != TheComp) Separating = FALSE;
				y = i/2 + 97;
				}
			else
				{
				if(XX[i+1] != TheComp) Separating = FALSE;
				y = i/2 + 65;
				}
			*q++ = y;	
			}
		}
	*q++ = '@';
	p = *Temp8;
	while( (*q++ = *p++) ) ;
	
	if(Relators[0] != NULL) DisposeHandle((char **) Relators[0]);
	Relators[0] = (unsigned char **) NewHandle(q - r);	
	if(Relators[0] == NULL) Mem_Error();
	p = ptr;
	q = *Relators[0];	
	while( (*q++ = *p++) ) ;
	DisposePtr((unsigned char*) ptr);
	
	Canonical_Rewrite(Relators,TRUE,FALSE); 
	
	for(i = 0; i < NumFilled; i++) if(SURL[i] == Length && NG[i] == NumGenerators
		&& NR[i] == NumRelators)
		{
	 	for(j = 1; j <= NumRelators; j++)
	 		if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;
	 	if(j > NumRelators && Compare_Pres(i))
	 		{
	 		/******************************************************************************
	 			We have discovered that an annulus exists in the diagram of presentation i.
	 			Save the info about the annulus in the handle SUR[i][0].
	 		******************************************************************************/
	 		
	 		if(UDV[i] <= DONE || UDV[i] == SEP_PAIRS)
	 			{
	 			HS = GetHandleSize((char **) Relators[0]);
	 			if(SUR[i][0] != NULL) DisposeHandle((char **) SUR[i][0]);
				SUR[i][0] = (unsigned char **) NewHandle(HS);
				if(SUR[i][0] == NULL) Mem_Error();
				p = *Relators[0];
				q = *SUR[i][0];
				r = q;
				while( (*q++ = *p++) ) ;
				if((q-r) != HS) 
					{
					NumErrors ++;
					printf("\n\n6) Error in Presentation %u! |Relator[0]| = %lu, HS = %lu.",i + 1,q-r-1,HS);
					}
				BytesUsed += GetHandleSize((char **) Relators[0]);
				UDV[i] = ANNULUS_EXISTS;
				if(F1)
					{
					if(Separating)
						{
						if((Batch == 4 || Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE)
						printf("\n                    There is a separating annulus in diagram %u.",i + 1);
						NumSepAnnuli ++;
						}
					else
						{
						if((Batch == 4 || Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE)
						printf("\n                    There is a nonseparating annulus in diagram %u.",i + 1);
						NumNonSepAnnuli ++;
						}
					}
				}
	 		break;
	 		}
	 	}
	
	if(i == NumFilled)
		{
		if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);
		if(Save_Pres(ReadPres,0,Length,1,70,0,0,0)) return(TOO_LONG);
		HS = GetHandleSize((char **) Relators[0]);
		if(SUR[NumFilled - 1][0] != NULL) DisposeHandle((char **) SUR[NumFilled - 1][0]);
		SUR[NumFilled - 1][0] = (unsigned char **) NewHandle(HS);
		if(SUR[NumFilled - 1][0] == NULL) Mem_Error();
		p = *Relators[0];
		q = *SUR[NumFilled - 1][0];		
		r = q;
		while( (*q++ = *p++) ) ;
		if((q-r) != HS) 
			{
			NumErrors ++;
			printf("\n\n7) Error in Presentation %u! |Relator[0]| = %lu, HS = %lu.",NumFilled,q-r-1,HS);
			}
		BytesUsed += GetHandleSize((char **) Relators[0]);
		BDY[NumFilled - 1] = BDY[ReadPres];
		UDV[NumFilled - 1] = ANNULUS_EXISTS;
		if(F1)
			{
			if(Separating)
				{
				if((Batch == 4 || Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE)
				printf("\n                    There is a separating annulus in diagram %u.",NumFilled);
				NumSepAnnuli ++;
				}
			else
				{
				if((Batch == 4 || Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE)
				printf("\n                    There is a nonseparating annulus in diagram %u.",NumFilled);
				NumNonSepAnnuli ++;
				}
			}	
		}
	return(Separating);			
}

int Find_Path(int V1,int V2,int TheComp,int SepType)
{
	/******************************************************************************************
			Let S be the set of directed paths, in the Heegaard diagram, which originate at a
		vertex in TheComp and leave TheComp via vertex V1. Find_Path() finds the longest path
		P on the Heegaard surface such that:
			1) No vertex of P lies in TheComp.
			2) P's first vertex is V1.
			3) P is an initial segment of each member of S.
		NOTE: If each member of S returns to TheComp via vertex V2, then the Heegaard diagram
		contains an annulus, since each path in S is parallel to P, up to the point where P
		returns to TheComp. (The annulus in question swallows TheComp and follows P.)
	******************************************************************************************/
		
	register unsigned char	*p,
							*q,
							x,
							y,
							z;

	unsigned char 			T1[125],
							T2[125],
							**Temp;
								
	int						FP,
							i,
							j;

	long					max;
	
	unsigned long			HS;
								
	for(i = 65; i < 65 + NumGenerators; i++)
		{
		if(XX[(i << 1) - 130] == TheComp)
			T2[i] = 1;
		else
			T2[i] = 0;
		}
	for(i = 97; i < 97 + NumGenerators; i++)
		{
		if(XX[(i << 1) - 193] == TheComp)
			T2[i] = 1;
		else
			T2[i] = 0;
		}
	for(i = 65; i < 65 + NumGenerators; i++)
		{
		T1[i] = T2[i + 32];
		T1[i + 32] = T2[i];
		}

	if(SepType == 2)
		{
		for(i = 65; i < 65 + NumGenerators; i++)
			if(T1[i])
				T2[i] = 2;
			else
				T2[i] = 0;	
		for(i = 97; i < 97 + NumGenerators; i++)
			if(T1[i])
				T2[i] = 2;
			else
				T2[i] = 0;	
		if(V2 & 1)
			z = (V2 >> 1) + 97;
		else
			z = (V2 >> 1) + 65;
		T2[z] = 1;	
		}
	else
		{		
		for(i = 65; i < 65 + NumGenerators; i++)
			if((T2[i] == 0) && (T2[i + 32] == 1)) T2[i] = 2;
		for(i = 97; i < 97 + NumGenerators; i++)
			if((T2[i] == 0) && (T2[i - 32] == 1)) T2[i] = 2;
		}
		
	for(i = 1,max = 0L; i <= NumRelators; i++)	if(GetHandleSize((char **) Relators[i]) > max)
		max = GetHandleSize((char **) Relators[i]);
	
	if(Temp8 != NULL) DisposeHandle((char **) Temp8);
	Temp8 = (unsigned char **) NewHandle(max + 1);
	if(Temp8 == NULL) Mem_Error();
	
	if(V1 & 1)
		z = (V1 >> 1) + 97;
	else
		z = (V1 >> 1) + 65;
		
	/******************************************************************************************
								Find a maximal common path.
	******************************************************************************************/
	
	for(j = 0,FP = FALSE; j < 2; j++)
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		q = p + GetHandleSize((char **) Relators[i]) - 2;
		x = *q;
		while( (y = *p++) )
			{
			if((y == z) && T1[x])
				{
				if(FP)
					{
					for(q = *Temp8 + 1; *p && *q == *p; p++,q++) ;
					if(*p == EOS && *q)
						{
						for(p = *Relators[i]; *p && *q == *p; p++,q++) ;
						if(*q && T2[*p] != 1) *q = *p;	
						q++;	
						*q = EOS;				
						break;	/* Done with this relator. */
						}
					if(*q)
						{
						if(T2[*p] != 1) *q = *p;	
						q++;
						*q = EOS;	
						y = *p;
						p++;
						}
					else
						{
						q = p - 1;
						y = *q;
						}					
					}
				else
					{
					q = *Temp8;
					*q++ = z;
					while(*p && T2[*p] == 0) *q++ = *p++;
					if(*p == EOS)
						{
						p = *Relators[i];
						while(*p && T2[*p] == 0) *q++ = *p++;
						*q++ = *p;	
						*q = EOS;
						FP = TRUE;
						break;		/* Done with this relator. */				
						}
					*q++ = *p;	
					*q = EOS;
					FP = TRUE;
					y = *p;
					p++;
					}	
				}
			x = y;	
			}
		Inverse(*Relators[i]);
		}
	
	/******************************************************************************************
		The last char of Temp8 serves as a flag, and is not actually a part of the path P.
		The routine is set up so that T2[last char] == 1 iff there is an annulus present.
		We next strip this last char off of the end of Temp8, shorten Temp8 by 1, and test
		the last char.
	******************************************************************************************/
		
	for(q = *Temp8; *q; q++) ;
	q--;
	x = *q;
	*q = EOS;
	HS = q - *Temp8;
	
	if(Temp16 != NULL) DisposeHandle((char **) Temp16);
	Temp16 = (unsigned char **) NewHandle(HS + 1);
	if(Temp16 == NULL) Mem_Error();
	p = *Temp16;
	q = *Temp8;
	while( (*p++ = *q++) ) ;
	Temp = Temp8;
	Temp8 = Temp16;
	Temp16 = Temp;
	if(T2[x] == 1) return(1); 					/*	An annulus exists. */
	return(0);
}

int Do_Aut_L(void)
{		
	/******************************************************************************************
		This variant of Do_Aut() performs a sequence of level T-transformations on the
		presentation Relators[]. The specific sequence of T-transformations performed is
		determined by the string which Find_Path() returned in Temp8.
	******************************************************************************************/
	
	register unsigned char 	A,
							a,
							*p,
							*q,
							*r,
							x,
							y;
							
	register 				int i;
	
	unsigned char 			TX[125],
							TY[125];
							
	int						Vertex;
	
	unsigned int			Auts,
							j;
							
	long					HS;												
	
	Num_Level_Slides ++;
	Auts = GetHandleSize((char **) Temp8) - 1;
	
	for(j = 0; j < Auts; j++)
		{
		r = *Temp8 + j;
		if(*r < 97)
			{
			A = *r;
			Vertex = 2*A - 130;
			}
		else
			{
			A = *r - 32;
			Vertex = 2*A - 129;
			}	
		a = A + 32;
		
		for(i = 0; i < Vertices; i++) ZZ[i] = YY[i];
			ZZ[Vertex] = 1;		
		if(!(Vertex & 1)) for(i = 0; i < Vertices; i++) 
			{
			if(ZZ[i])
				ZZ[i] = 0;
			else
				ZZ[i] = 1;
			}	
		
		for(i = 0; i < NumGenerators; i++)
			if(ZZ[i << 1])
				TX[i+97] = TY[i+65] = FALSE;
			else
				TX[i+97] = TY[i+65] = TRUE;
		for(i = 0; i < NumGenerators; i++)
			if(ZZ[(i << 1) + 1])
				TX[i+65] = TY[i+97] = FALSE;
			else
				TX[i+65] = TY[i+97] = TRUE;
				
		for(i = 1; i <= NumRelators; i++) if(**Relators[i])
			{
			HS = GetHandleSize((char **) Relators[i]);
			if(HS > MAXLENGTH) return(TOO_LONG);
			if(Temp5 != NULL) DisposeHandle((char **) Temp5);
			Temp5 = (unsigned char **) NewHandle(2*HS);
			if(Temp5 == NULL) Mem_Error();
			p = *Temp5;	
			q = *Relators[i];
			x = *q++;
			while( (y = *q++) )
				{
				if(x != A && x != a) *p++ = x;
				if(TX[x] && !TY[y])
					*p++ = a;
				else
				if(!TX[x] && TY[y])
					*p++ = A;
				x = y;
				}
			if(x != A && x != a) *p++ = x;
			if(TX[x] && !TY[**Relators[i]])
				*p++ = a;
			else
			if(!TX[x] && TY[**Relators[i]])
				*p++ = A;
			*p = EOS;
			q = *Temp5;
			HS = p + 1 - q;
			if(HS > MAXLENGTH) return(TOO_LONG);			
			if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
            Relators[i] = (unsigned char **) NewHandle(HS);
            if(Relators[i] == NULL) Mem_Error();
            p = *Temp5;
            q = *Relators[i];
            while( (*q++ = *p++) ) ; 					
			}
		}
		
	TotalAuts += Auts;	
	return(NO_ERROR);				
}

int Do_Aut_L_Long(int V1,int V2,int TheComp)
{
	/******************************************************************************************
		Like Do_Aut_L(), this variant of Do_Aut() performs a sequence of level
		T-transformations on the presentation Relators[]. The specific sequence of
		T-transformations performed is determined by the string which Find_Path() returned in
		Temp8. If Level_Transformations() calls Do_Aut_L(), the number of consecutive
		level-transformations that Do_Aut_L() must perform is equal to the length of Temp8.
		If this length is greater than 5, (it is sometimes several thousand), then it is
		significantly more efficient to call Do_Aut_L_Long() instead of Do_Aut_L().
			Do_Aut_L_Long() works by first deleting certain appearances of Temp8 and its
		inverse from the relators, and then reinserting copies of Temp8 and its inverse
		at other locations in the relators.
	******************************************************************************************/

	register unsigned char	*p,
							*q,
							*r,
							w,
							x,
							y,
							z;
	
	unsigned char		T1[125],
						T2[125];
	
	int			i,
				j;
	
	unsigned long		Delete,
						Insert,
						HS,
						Li,
						L8;
	
	Num_Level_Slides ++;
	
	for(x = 65; x < 65 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 129] == TheComp)
			T1[x] = TRUE;
		else
			T1[x] = FALSE;	
		}
		
	for(x = 97; x < 97 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 194] == TheComp)
			T1[x] = TRUE;
		else
			T1[x] = FALSE;	
		}
	
	for(x = 65; x < 65 + NumGenerators; x++) T2[x] = T1[x + 32];
	for(x = 97; x < 97 + NumGenerators; x++) T2[x] = T1[x - 32];
	
	T2[0] =	T1[0] = T1[63] = T1[64] = T1[95] = T1[96] = EOS;
	
	for(x = 1; x <= NumRelators; x++) T1[x] = EOS;
		
	if(V1 & 1)
		z = (V1 >> 1) + 97;
	else
		z = (V1 >> 1) + 65;
	
	if(V2 & 1)
		w = (V2 >> 1) + 97;
	else
		w = (V2 >> 1) + 65;
		
	L8 = GetHandleSize((char **) Temp8) - 1;
		
	for(j = 0; j < 2; j++)
	for(i = 1; i <= NumRelators; i++)
		{
		Li = GetHandleSize((char **) Relators[i]);
		if(Li > L8 + 1)
			{
			Delete = 0L;
			q = *Relators[i] + Li - 1;
			p = q - L8;
			x = *p++;
			while( (y = *p) )
				{
				if((y == z) && T1[x])
					{
					*p = EOS;
					Delete = L8 + p - q;
					x = z;
					break;	
					}
				x = y;
				p++;	
				}
			if(Li > MAXLENGTH) return(TOO_LONG);
			if(Temp5 != NULL) DisposeHandle((char **) Temp5);
			Temp5 = (unsigned char **) NewHandle(Li);	
			if(Temp5 == NULL) Mem_Error();
			r = *Temp5;	
			p = *Relators[i] + Delete;
			if(Delete)
				{
				/*****************************************************************************
					The first Delete chars of Relators[i] are a proper terminal substring  of
					a copy of Temp8. If this copy of Temp8 ends at a vertex of The_Comp, it
					will need to be reinserted later. So, mark this location in Temp5 by
					inserting the char '64'. Otherwise, insert the char '63' in Temp5 to
					indicate that a deletion was made from this location but no reinsertion
					should take place.
				*****************************************************************************/
					
				if(T2[*p])
					*r++ = 64;
				else
					*r++ = 63;
				} 
			while( (y = *p) )
				{
				if((y == z) && T1[x])
					{
					p += L8;
					Delete++;
					if(*p == EOS)
						{
						if(T2[**Relators[i]])
							*r++ = 64;
						else
							*r++ = 63;	
						break;
						}
					if(T2[*p])
						*r++ = 64;
					else
						*r++ = 63;	
					x = EOS;	
					}
				else
					{
					*r++ = y;
					x = y;
					p++;
					}	
				}
			if(Delete)
				{
				T1[i] 		= TRUE;	
				*r++ 		= EOS;
				HS = r - *Temp5;
				if(HS > MAXLENGTH) return(TOO_LONG);				
				if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
				Relators[i] = (unsigned char **) NewHandle(HS);
				if(Relators[i] == NULL) Mem_Error();
				p = *Temp5;
				r = *Relators[i];
				while( (*r++ = *p++) ) ; 
				}	
			}
		Inverse(*Relators[i]);
		}
	
	/*****************************************************************************************
						Put a copy of the inverse of Temp8 in Temp9.
	*****************************************************************************************/
	
	HS = L8 + 1;
	if(HS > MAXLENGTH) return(TOO_LONG);
	if(Temp9 != NULL) DisposeHandle((char **) Temp9);
	Temp9 = (unsigned char **) NewHandle(HS);
	if(Temp9 == NULL) Mem_Error();
	q = *Temp9;
	p = *Temp8;
	while( (*q++ = *p++) ) ;
	Inverse(*Temp9);
	
	if(w < 97)
		z = w + 32;
	else
		z = w - 32;
		
	for(i = 1; i <= NumRelators; i++)
		{
		Li = GetHandleSize((char **) Relators[i]);
		
		/*************************************************************************************
			Determine how the length of Relator[i] will be changed by counting how many
			copies of the strings Temp8 and Temp9 need to be inserted into Relators[i].
		*************************************************************************************/
			
		Insert = 0L;
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while( (y = *p++) )
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64 || y == 96)
				{
				Insert ++;
				x = EOS;
				continue;
				}
			if((x == z) && T2[y])
				Insert ++;
			else	
			if((y == w) && T1[x])
				Insert ++;
			x = y;
			}
		
		if(Insert == 0L && T1[i] == EOS) continue;
			
		/*************************************************************************************
					Reserve some memory for the	new version of Relators[i].
		*************************************************************************************/
		
		HS = Li + Insert*L8;
		if(HS > MAXLENGTH) return(TOO_LONG);
		if(Temp5 != NULL) DisposeHandle((char **) Temp5);
		Temp5 = (unsigned char **) NewHandle(HS);	
		if(Temp5 == NULL) Mem_Error();
		r = *Temp5;
		
		/*************************************************************************************
									Rewrite Relators[i].
		*************************************************************************************/
			
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while( (y = *p++) )
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64)
				{
				q = *Temp8;
				while( (*r++ = *q++) ) ;
				r--;
				x = EOS;
				continue;	
				}
			if(y == 96)
				{
				q = *Temp9;
				while( (*r++ = *q++) ) ;
				r--;
				x = EOS;
				continue;	
				}
			if((x == z) && T2[y])
				{
				q = *Temp8;
				while( (*r++ = *q++) ) ;
				r--;	
				}
			else	
			if((y == w) && T1[x])
				{
				q = *Temp9;
				while( (*r++ = *q++) ) ;
				r--;
				}
			*r++ = y;	
			x = y;	
			}
		*r++ = EOS;
		HS = r - *Temp5;
		if(HS > MAXLENGTH) return(TOO_LONG); 		
		if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
		Relators[i] = (unsigned char **) NewHandle(HS);
		if(Relators[i] == NULL) Mem_Error();
		p = *Temp5;
		q = *Relators[i];
		while( (*q++ = *p++) ) ; 		
		}
	
	TotalAuts += L8;
	return(NO_ERROR);			
}

int Delete_Trivial_Generators(int Test)
{
	/***************************************************************************************
		This routine deletes a generator X from the relators, provided there is a relator
		of the form X, and the total number of appearances of X, together with its inverse,
		in the remaining relators is no more than one.
	***************************************************************************************/
	
	register unsigned char	*p,
							*q,
							s,
							t,
							x,
							y,
							z;
	
	unsigned char			**Temp;
							
	int			i,
				j,
				k,
				SNumRelators1,
				SNumRelators2,
				TheRelator;
							
	long			HS;						

	if(NumRelators == 1) return(FALSE);
			
	SNumRelators1 = NumRelators;

RETRY:
	SNumRelators2 = NumRelators;
												
	for(i = 1; i <= NumRelators; i++) if(GetHandleSize((char **) Relators[i]) == 2)
		{
		p = *Relators[i];
		x = *p;
		if(x > 95) x -= 32;
		y = x + 32;
		for(j = 1,k = 0; j <= NumRelators; j++)
			{
			if(i == j) continue;
			p = *Relators[j];
			while( (z = *p++) ) if(z == x || z == y)
				{
				k++;
				if(k > 1) goto TEST_NEXT_RELATOR;
				if(k == 1) TheRelator = j;
				}
			}
			
		if(k == 1)
			{
			HS = GetHandleSize((char **) Relators[TheRelator]);
			if(HS == 2)
				{
				/**************************************************************************
					If k = 1 and HS = 2, the presentation contains two relators of the
					form X, and deleting X would create an empty relator. We do not want
					to do this at this point.
				**************************************************************************/
						
				continue;
				}
			if(Test) return(TRUE);	
			
			if(Temp16 != NULL) DisposeHandle((char **) Temp16);
			Temp16 = (unsigned char **) NewHandle(HS - 1);
			if(Temp16 == NULL) Mem_Error();
			p = *Relators[TheRelator];
			q = *Temp16;
			while( (z = *p++) ) if(z != x && z != y) *q++ = z;
			*q = EOS;
			Temp = 					Relators[TheRelator];
			Relators[TheRelator] 	= Temp16;
			Temp16 					= Temp;
			Length -= 1;
			}
		
		if(Test) return(TRUE);	
		Length -= 1;	
		Temp 					= Relators[i];
		Relators[i] 			= Relators[NumRelators];
		Relators[NumRelators] 	= Temp;
		
		if(Micro_Print)
			{
			printf("\n\nRelator %d is trivial. ",i);
			printf("Deleted generator %c from the presentation.",x);
			if(i != NumRelators)
				printf("\nSwapped Relator %d with Relator %d and reduced the number of relators.",
					i,NumRelators);
			}
			
		/**********************************************************************************
					If x is not the largest generator, replace the largest generator
					throughout with x, and its inverse with y.								
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
						Update the variables NumGenerators, NumRelators and Vertices.
		**********************************************************************************/
		
		NumRelators --;
		NumGenerators --;
		Vertices -= 2;
		
		/**********************************************************************************
								Freely reduce the new relators.
		**********************************************************************************/
		
		if(Micro_Print) for(j = 1,Length = 0L; j <= NumRelators; j++)
			Length += GetHandleSize((char **) Relators[j]) - 1;
		if(Freely_Reduce() == TOO_LONG) return(TOO_LONG);
		Length = OrigLength;
		if(NumRelators == 1)
			{
			if(NumRelators == SNumRelators1) return(FALSE);
			return(TRUE);
			}
		TEST_NEXT_RELATOR:
		continue;	
		}
	if(NumRelators < SNumRelators2) goto RETRY;	
	if(NumRelators == SNumRelators1) return(FALSE);
	return(TRUE);	
}

int Level_Transformations_2(int F1,int F2,int F3,unsigned int MyNum)
{
	/******************************************************************************************
			Suppose {X,Y} is a separating pair of vertices of the "reduced" Whitehead graph G
		of a minimal length presentation, with x the inverse of X and y the inverse of Y.
		Such separating pairs of vertices are of two types. Say {X,Y} is a "Type 1" pair
		if Y = x, or if deleting X and Y from G yields a component C which contains both x
		and y. If no such component exists, say that {X,Y} is a "Type 2" pair.
			It is possible to show that if {X,Y} is a Type 2 pair, then the separation of G, 
		produced by deleting X and Y from G, has only two components, and there is no edge 
		joining X and Y in G.
			Level_Transformations() calls Level_Transformations_2() when the separation has
		two components.
			F1 and F2 are flags which are passed along by Level_Transformations_2() to any
		new invocations of Level_Transformations(). If F1 is FALSE, then printing of
		information about any annuli the routine finds is  suppressed. If F2 is TRUE, new
		invocations of Level_Transformations() will stop and return 2 when they find the
		first new presentation without any separating pairs of vertices.	
	******************************************************************************************/
	
	int 		k,
				TheComp;	
							
	unsigned int 	i,
					j,
					Random,
					VL,
					VR,
					VLI,
					VRI;
		
	/******************************************************************************************
			Let VL and VR be a pair of Type 2 separating vertices of G, and let TheComp be one
		of the two components obtained by deleting VL and VR from G. Let VLI be the inverse of
		VL and let VRI be the inverse of VR. If VRI lies in TheComp, let C be the subgraph
		consisting of TheComp together with VR. Then the objective is to "slide" C around the
		Whitehead graph by sliding to the "left" until either C is attached to the main body of
		the graph at more than two distinct vertices, or we arrive at a vertex of attachment V
		of C with the property that V inverse lies in C, or we manage to slide C completely
		around the graph along some path and the presentation has returned to the original
		presentation. If the last case happens, then there is an annulus present. Otherwise
		we call Level_Transformations() recursively and pass it the new presentation.
	******************************************************************************************/
									
	VL = BSV1[MyNum];
	VR = BSV2[MyNum];
	for(i = 0; i < Vertices; i++) XX[i] = 0;
	XX[VL] = VERTICES;
	XX[VR] = VERTICES;
	if(VL & 1)
		VLI = VL - 1;
	else
		VLI = VL + 1;		
	if(VR & 1)
		VRI = VR - 1;
	else
		VRI = VR + 1;		
	Get_Components();

	if(NumSepComps != 2) return(0);
	
	if(RandomizeSlides)	Random = abs(rand()) %2;
	else	Random = 1;
	
	for(k = 1; k <= 2; k++)
		{
		if(Random) TheComp = k;
		else TheComp = 3 - k;
		if(XX[VLI] != TheComp && XX[VRI] == TheComp)
			{
			if((j = Slide_LComp(VL,VR,TheComp,2,MyNum,F1,F2,F3)) > 1) return(j);
			if(Level_Trans_Reset(MyNum,VL,VR) == TOO_LONG) return(TOO_LONG);
			}
			
		if(XX[VLI] == TheComp && XX[VRI] != TheComp)			
			{
			if((j = Slide_LComp(VR,VL,TheComp,2,MyNum,F1,F2,F3)) > 1) return(j);
			if(Level_Trans_Reset(MyNum,VL,VR) == TOO_LONG) return(TOO_LONG);
			}
		}
			
	return(0);
}

int Do_Aut_L_Long_2(int V1,int V2,int TheComp)
{
	/******************************************************************************************
		Like Do_Aut_L(), this variant of Do_Aut() performs a sequence of level
		T-transformations on the presentation Relators[]. The specific sequence of
		T-transformations performed is determined by the string which Find_Path() returned
		in Temp8. If Level_Transformations_2() calls Do_Aut_L(), the number of consecutive
		level-transformations that Do_Aut_L() must perform is equal to the length of Temp8.
		If this length is greater than 5, (it is sometimes several thousand), then it is
		significantly more efficient to call Do_Aut_L_Long_2() instead of Do_Aut_L().
			Do_Aut_L_Long_2() works by first deleting certain appearances of Temp8 and its
		inverse from the relators, and then reinserting copies of Temp8 and its inverse
		at other locations in the relators.
	******************************************************************************************/

	register unsigned char	*p,
							*q,
							*r,
							w,
							x,
							y,
							z;
	
	unsigned char			T1[125],
							T2[125];
	
	int			i,
				j;
	
	unsigned long		Delete,
						Insert,
						Li,
						L8;
							
	Num_Level_Slides ++;						
	
	for(x = 65; x < 65 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 129] == TheComp)
			T1[x] = TRUE;
		else
			T1[x] = FALSE;	
		}
		
	for(x = 97; x < 97 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 194] == TheComp)
			T1[x] = TRUE;
		else
			T1[x] = FALSE;	
		}
	
	for(x = 65; x < 65 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 130] == TheComp)
			T2[x] = FALSE;
		else
			T2[x] = TRUE;	
		}
	for(x = 97; x < 97 + NumGenerators; x++)
		{
		if(XX[(x << 1) - 193] == TheComp)
			T2[x] = FALSE;
		else
			T2[x] = TRUE;	
		}
	
	T2[0] = T1[0] = T1[63] = T1[64] = T1[95] = T1[96] = EOS;
	
	for(x = 1; x <= NumRelators; x++) T1[x] = EOS;
		
	if(V1 & 1)
		z = (V1 >> 1) + 97;
	else
		z = (V1 >> 1) + 65;
	
	if(V2 & 1)
		w = (V2 >> 1) + 97;
	else
		w = (V2 >> 1) + 65;
		
	L8 = GetHandleSize((char **) Temp8) - 1;
		
	for(j = 0; j < 2; j++)
	for(i = 1; i <= NumRelators; i++)
		{
		Li = GetHandleSize((char **) Relators[i]);
		if(Li > L8 + 1)
			{
			Delete = 0L;
			q = *Relators[i] + Li - 1;
			p = q - L8;
			x = *p++;
			while( (y = *p) )
				{
				if((y == z) && T1[x])
					{
					*p = EOS;
					Delete = L8 + p - q;
					x = z;
					break;	
					}
				x = y;
				p++;	
				}
			if(Temp5 != NULL) DisposeHandle((char **) Temp5);
			Temp5 = (unsigned char **) NewHandle(Li);
			if(Temp5 == NULL) Mem_Error();
			r = *Temp5;
			p = *Relators[i] + Delete;
			if(Delete)
				{
				/*****************************************************************************
					The first Delete chars of Relators[i] are a proper terminal substring  of
					a copy of Temp8. If this copy of Temp8 ends at V2, it will need to be
					reinserted later. So, mark this location in Temp5 by inserting the char
					'64'. Otherwise, insert the char '63' in Temp5 to indicate that a
					deletion was made from this location but no reinsertion should take
					place.
				*****************************************************************************/
					
				if(*p == w)
					*r++ = 64;
				else
					*r++ = 63;
				} 
			while( (y = *p) )
				{
				if((y == z) && T1[x])
					{
					p += L8;
					Delete++;
					if(*p == EOS)
						{
						if(**Relators[i] == w)
							*r++ = 64;
						else
							*r++ = 63;	
						break;
						}
					if(*p == w)
						*r++ = 64;
					else
						*r++ = 63;	
					x = EOS;	
					}
				else
					{
					*r++ = y;
					x = y;
					p++;
					}	
				}
			if(Delete)
				{
				T1[i] 		= TRUE;	
				*r++ 		= EOS;				
				if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
				Relators[i] = (unsigned char **) NewHandle(r - *Temp5);
				if(Relators[i] == NULL) Mem_Error();
				r = *Temp5;
				p = *Relators[i];
				while( (*p++ = *r++) ) ; 
				}	
			}
		Inverse(*Relators[i]);
		}
	
	/*****************************************************************************************
						Put a copy of the inverse of Temp8 in Temp9.
	*****************************************************************************************/
	
	if(Temp9 != NULL) DisposeHandle((char **) Temp9);
	Temp9 = (unsigned char **) NewHandle(L8 + 1);	
	if(Temp9 == NULL) Mem_Error();
	q = *Temp9;
	p = *Temp8;
	while( (*q++ = *p++) ) ;
	Inverse(*Temp9);
	
	if(w < 97)
		z = w + 32;
	else
		z = w - 32;
	
	T1[0] = T1[63] = T1[64] = T1[95] = T1[96] = 1;
		
	for(i = 1; i <= NumRelators; i++)
		{
		Li = GetHandleSize((char **) Relators[i]);
		
		/*************************************************************************************
			Determine how the length of Relator[i] will be changed by counting how many
			copies of the strings Temp8 and Temp9 need to be inserted into Relators[i].
		*************************************************************************************/
			
		Insert = 0L;
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while( (y = *p++) )
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64 || y == 96)
				{
				Insert ++;
				x = EOS;
				continue;
				}
			if((x == z) && T2[y])
				Insert ++;
			else	
			if((y == w) && !T1[x])
				Insert ++;
			x = y;
			}
		
		if(Insert == 0L && T1[i] == EOS) continue;
			
		/*************************************************************************************
					Reserve some memory for the	new version of Relators[i].
		*************************************************************************************/
		
		if(Temp5 != NULL) DisposeHandle((char **) Temp5);
		Temp5 = (unsigned char **) NewHandle(Li + Insert*L8);	
		if(Temp5 == NULL) Mem_Error();
		r = *Temp5;
		
		/*************************************************************************************
									Rewrite Relators[i].
		*************************************************************************************/
			
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while( (y = *p++) )
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64)
				{
				q = *Temp8;
				while( (*r++ = *q++) ) ;
				r--;
				x = EOS;
				continue;	
				}
			if(y == 96)
				{
				q = *Temp9;
				while( (*r++ = *q++) ) ;
				r--;
				x = EOS;
				continue;	
				}
			if((x == z) && T2[y])
				{
				q = *Temp9;
				while( (*r++ = *q++) ) ;
				r--;	
				}
			else	
			if((y == w) && !T1[x])
				{
				q = *Temp8;
				while( (*r++ = *q++) ) ;
				r--;
				}
			*r++ = y;	
			x = y;	
			}
		*r++ = EOS;		
		if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
		Relators[i] = (unsigned char **) NewHandle(r - *Temp5);
		if(Relators[i] == NULL) Mem_Error();
		r = *Temp5;
		p = *Relators[i];
		while( (*p++ = *r++) ) ; 			
		}
	
	TotalAuts += L8;
	return(NO_ERROR);			
}

int Slide_ValenceTwo_Comp(int TheComp,unsigned int VL,unsigned int VR,unsigned int MyNum)
{
	char			xx,
					yy;
					
	unsigned char	x,
					j,
					*p,
					*q,
					**Temp;
					
	int				i;
	
	unsigned int	Width;
	
	/******************************************************************************************
							Check if each vertex in TheComp has valence two.
	******************************************************************************************/
							

	for(i = 0; i < Vertices; i++)
	if(XX[i] == TheComp && VWG[i] != 2) return(0);
		
	/******************************************************************************************
		Compute the Width of the annulus. This is the number of edges of the Whitehead graph
	which connect vertices of TheComp to VL (resp VR). (Note that because the presentation
	has been reduced to minimal length, the number of edges connecting vertices of TheComp
	to VL and the number of edges connecting vertices of TheComp to VR must be equal.)
	******************************************************************************************/	
	
	for(i = Width = 0; i < Vertices; i++) if(XX[i] == TheComp) Width += A[i][VL];
	
	/******************************************************************************************
		Check if there is a pair of vertices, say Vx and Vy, in the Whitehead graph, which are 
	connected only by edges of the annulus swallowing TheComp. This will be the case if and
	only if A[Vx][Vy] = Width. If this happens, and the pair (Vx,Vy) is not a pair of 
	separating vertices of the Whitehead graph, we will slide TheComp around the annulus until 
	it lies between Vx and Vy.												
	******************************************************************************************/
	
	p = *Temp8;
	if(*p == EOS) return(0);
	x = *p << 1;
	if(x < 194) x -= 129;
	else x -= 194;
	j = x;
	p++;
	while( (x = *p) )
		{
		x = x << 1;
		if(x < 194) x -= 130;
		else x -= 193;
		if((A[j][x] == Width) && (Is_Sep_Pair(j,x,MyNum) == FALSE))
			{
			*p = EOS;			
			if(Temp16 != NULL) DisposeHandle((char **) Temp16);
			Temp16 = (unsigned char **) NewHandle(p - *Temp8 + 1);
			if(Temp16 == NULL) Mem_Error();
			p = *Temp16;
			q = *Temp8;
			while( (*p++ = *q++) ) ;
			Temp = Temp8;
			Temp8 = Temp16;
			Temp16 = Temp;
			
			/* 	We have found a pair of vertices we are looking for.	*/
			
			if(Micro_Print)
				{
				if(j & 1)
					xx = j/2 + 97;
				else
					xx = j/2 + 65;
				if(x & 1)
					yy = x/2 + 97;
				else
					yy = x/2 + 65;     
				printf("\n\nTheComp lies in a 'Valence-Two Annulus'. Will slide TheComp around ");
				printf("this annulus so it lies between vertices %c and %c.",xx,yy);
				}
			return(1);
			}
		if(x & 1) j = x - 1;
		else j = x + 1;
		p++;
		}
		
	return(0);	
}				

int Level_Trans_Reset(unsigned int MyNum, unsigned int V3, unsigned int V4)
{
	/**********************************************************************************
		If a call to Level_Transformations() is unsuccessful, then we need to reset
		Relators[] and other globals that may have been changed before we try to slide
		TheComp in the other direction. This routine performs those restitutions.
	**********************************************************************************/	
	
	unsigned char 	*p,
					*q;
					
	unsigned int 	i,
					VL,
					VR;

	for(i = 1; i <= NumRelators; i++)
		{
		if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
		Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SLR[MyNum][i]));
		if(Relators[i] == NULL) Mem_Error();
		q = *Relators[i];
		p = *SLR[MyNum][i];
		while( (*q++ = *p++) ) ;
		}
		
	if(Micro_Print) Micro_Print_Level_Transformations_Reset(MyNum);	
	
	for(i = 0; i < Vertices; i++) XX[i] = 0;
	
	VL = V3;
	VR = V4;
	XX[VL] = VERTICES;
	XX[VR] = VERTICES;
	Fill_A(NumRelators);
	Get_Matrix();	
	Get_Components();
	
	return(0);
}

void Micro_Print_Level_Transformations_Reset(unsigned int MyNum)
{
    printf("\n\nReverting to looking for level-transformations of Presentation L%u:\n", MyNum + 1);
    Print_Relators(Relators,NumRelators);
}

unsigned int Count_Sep_Pairs(unsigned int Num_Saved_LPres)
{
	/**********************************************************************************
	  This routine calls Sep_Pairs() repeatedly in order to locate all pairs of 
	  separating vertices of the Reduced Whitehead Graph RWG of the presentation
	  Num_Saved_LPres -1. At termination, the array SLR[Num_Saved_LPres -1][0] 
	  contains a list of all separating pairs of vertices of the RWG. 
	  	The total number of separating pairs of vertices found is saved twice. Once in
	 Num_Sep_Vertex_Pairs[] at Num_Sep_Vertex_Pairs[Num_Saved_LPres -1]. And also in
	 UnUsed_Sep_Vertex_Pairs[] at UnUsed_Sep_Vertex_Pairs[Num_Saved_LPres -1]. (Later, 
	 the value in UnUsed_Sep_Vertex_Pairs[Num_Saved_LPres -1] is decremented each time 
	 a separating pair of vertices is selected for processing.)
	***********************************************************************************/
	
	char			x,
					y;
					
	unsigned char	*p,
					*q;
					
	unsigned int	i,
					Num_Sep_Pairs;
	
	Num_Sep_Pairs = 0;
	p = SepVertexList;
	
	if(Micro_Print && (Num_Saved_LPres == 1))
		printf("\nNote the initial presentation will be referred to as Presentation L1.\n");
		
	if(Sep_Pairs(0,0,1))
		{
		if(Micro_Print)
			{
			if(V1 & 1)
				x = V1/2 + 97;
			else
				x = V1/2 + 65;
			if(V2 & 1)
				y = V2/2 + 97;
			else
				y = V2/2 + 65;
			printf("\n(%c,%c)",x,y);
			}
		*p++ = V1;
		*p++ = V2;
		Num_Sep_Pairs ++;
		while(Sep_Pairs(V1,V2,0))
			{
			if(Micro_Print)
				{
				if(V1 & 1)
					x = V1/2 + 97;
				else
					x = V1/2 + 65;
				if(V2 & 1)
					y = V2/2 + 97;
				else
					y = V2/2 + 65;
				printf(", (%c,%c)",x,y);
				}
			*p++ = V1;
			*p++ = V2;	
			Num_Sep_Pairs ++;
			}
		}
	
	Num_Sep_Vertex_Pairs[Num_Saved_LPres - 1] = Num_Sep_Pairs;
	UnUsed_Sep_Vertex_Pairs[Num_Saved_LPres - 1] = Num_Sep_Pairs;
	
	if(SLR[Num_Saved_LPres - 1][0] != NULL) DisposeHandle((char **) SLR[Num_Saved_LPres - 1][0]);
	SLR[Num_Saved_LPres - 1][0] = (unsigned char **) NewHandle(sizeof(char)*(2*Num_Sep_Pairs + 1));
	if(SLR[Num_Saved_LPres - 1][0] == NULL) Mem_Error();
	q = *SLR[Num_Saved_LPres - 1][0];		
	for(i = 0,p = SepVertexList; i < 2*Num_Sep_Pairs; i++) *q++ = *p++;	
		*q = VERTICES;
		
	if(Micro_Print)
		{
		if(Num_Sep_Pairs == 1)
			printf(" is the only separating pair of vertices of the RWG of Presentation L%u.",Num_Saved_LPres);
		else 
			printf(" are the %u separating pairs of vertices of the RWG of Presentation L%u. ",
			Num_Sep_Pairs,Num_Saved_LPres);
		}
	
	return(Num_Sep_Pairs);	
}

int Is_Sep_Pair(unsigned int VL,unsigned int VR,unsigned int MyNum)
{
	/***************************************************************************************
		This routine checks whether the pair of vertices (VL,VR) is a separating pair of
	vertices of the Whitehead graph in SLR[MyNum][]. It scans the list of pairs of 
	separating vertices in SLR[MyNum][0] and returns TRUE if the pair (VL,VR) appears.
	Otherwise, it returns FALSE.
	***************************************************************************************/
	
	unsigned char	*p,
					*q;
	
	unsigned int	i;
	
	/***************************************************************************************
		Separating pairs are listed in SLR[MyNum][0] as ordered pairs (VL,VR) with VL < VR.
	Therefore, we swap VL and VR if necessary so VL < VR.	
	***************************************************************************************/
	
	if(VL > VR)
		{
		i 		= VL;
		VL 		= VR;
		VR 		= i;
		}
		
	p = *SLR[MyNum][0];	
	q = p + 1;
	while(*p < VERTICES)
		{
		if((*p == VL) && (*q == VR)) return(TRUE);
		p += 2;
		q += 2;
		}
	
	return(FALSE);
}

void Micro_Print_Level_Transformations(unsigned int TheComp,unsigned int V1,unsigned int V2,
	unsigned int Type)
{
    char            x,
    				y;
                                        
    unsigned char   *p;
    
    int             i;
    
    if(V1 & 1)
        x = V1/2 + 97;
    else
        x = V1/2 + 65;
    if(V2 & 1)
        y = V2/2 + 97;
    else
        y = V2/2 + 65;                
    printf("\n\nVertices %c and %c form a Type %u separating pair.",x,y,Type);
    
    if(Temp9 != NULL) DisposeHandle((char **) Temp9);
    Temp9 = (unsigned char **) NewHandle(2*VERTICES);
    if(Temp9 == NULL) Mem_Error();
    p = *Temp9;
    
    for(i = 0; i < Vertices; i++)
        {
        if(XX[i] == TheComp)
            {
            if(i & 1)
                *p++ = i/2 + 97;
            else
                *p++ = i/2 + 65;
            *p++ = ',';    
            }
        }
    p--;
    *p = EOS;
  
    printf("\nPerformed a level-transformation by sliding vertice(s):");
    printf("\n{%s}",*Temp9);
    printf("\nalong a path represented by: ");
    printf("%s",*Temp8);
    printf("\nto obtain Presentation L%u:\n",Num_Saved_LPres + 1);
    Print_Relators(Relators,NumRelators);
}

unsigned int Slide_LComp(unsigned int VX, unsigned int VY, int TheComp, int SepType, unsigned int MyNum,
    int F1, int F2, int F3)
{	
	unsigned int	i,
					j;

	switch(Find_Path(VX,VY,TheComp,SepType))
		{
		case 0:
			break;
		case TOO_LONG:
			return(TOO_LONG);
		case FATAL_ERROR:
			return(FATAL_ERROR);	
		case 1:
			if(Found_L_Annulus && !TestRealizability4) 
				{
				if(Micro_Print) printf("\n\nFound another annulus.");
				break;
				}		
			Found_L_Annulus = TRUE;	
			if(TestRealizability4) return(13);			
			if((SepType == 1) && (Slide_ValenceTwo_Comp(TheComp,VX,VY,MyNum))) break;
			if((SepType == 1) && F3) return(4);
			if(SepType == 2) XX[VY] = TheComp;
			switch(Annulus(VX,VY,TheComp,F1))
				{
				case 0:
					return(0);
				case 1:
					return(4);
				case TOO_LONG:
					return(TOO_LONG);	
				}
		}
		
	if(GetHandleSize((char **) Temp8) > 5)
		{
		if((SepType == 1) && (Do_Aut_L_Long(VX,VY,TheComp) == TOO_LONG)) return(TOO_LONG);
		if((SepType == 2) && (Do_Aut_L_Long_2(VX,VY,TheComp) == TOO_LONG)) return(TOO_LONG);
		}
	else
		{
		for(i = 0; i < Vertices; i++)
			{
			if(XX[i] == TheComp)
				YY[i] = 1;
			else
				YY[i] = 0;
			}

		if((SepType == 1) && (Do_Aut_L() == TOO_LONG)) return(TOO_LONG);	
		if(SepType == 2)
			{
			YY[VY] = 1;
			if(Do_Aut_L() == TOO_LONG) return(TOO_LONG);
			}
		}	
	if(Micro_Print) 
		{
		if(SepType == 1) Micro_Print_Level_Transformations(TheComp,VX,VY,1);	
		if(SepType == 2)
			{
			XX[VY] = TheComp;
			Micro_Print_Level_Transformations(TheComp,VX,VY,2);
			XX[VY] = VERTICES;
			}		
		}
	
	if((j = Level_Transformations(F1,F2,F3)) > 1) return(j);							
	
	return(0);	
}

int Random_Sep_Pair(unsigned int WhichSLRPres)
{
	char			x,
					y;
					
	unsigned char	*p;
	
	unsigned int	j,
					k,
					jj,
					kk;
					
	if((k = UnUsed_Sep_Vertex_Pairs[WhichSLRPres]) == 0) return(1);	
	UnUsed_Sep_Vertex_Pairs[WhichSLRPres] --;
	
	j = abs(rand()) % k;
	
	p = *SLR[WhichSLRPres][0];
	
	jj = 2*j;
	kk = 2*k - 2;
	
	V1 = p[jj];
	V2 = p[jj+1];
	
	p[jj] 	= p[kk];
	p[jj+1] = p[kk+1];
	
	p[kk] 	= V1;
	p[kk+1] = V2;
	
	if(Micro_Print) 
		{
		if(V1 & 1)
        	x = V1/2 + 97;
		else
			x = V1/2 + 65;
		if(V2 & 1)
			y = V2/2 + 97;
		else
			y = V2/2 + 65;     
		printf("\nSelected Pair %u = (%c,%c) of %u untried pairs from %u of the RWG of Presentation L%u.",
			j+1,x,y,k,Num_Sep_Vertex_Pairs[WhichSLRPres],WhichSLRPres + 1);
		}	
	
	return(0);
}
