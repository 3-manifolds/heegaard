#include "Heegaard.h"
#include "Heegaard_Dec.h"

#define MAX_NOT_NEW_PRES	10 
#define SAVE_LEVEL_PRES		TRUE

Level_Transformations(F1,F2,F3,F4,F5)
int 	F1,
		F2,
		F3,
		F4,
		F5;
{
	/******************************************************************************************
		This routine is called when the "reduced" Whitehead graph has a pair of separating
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
	
	unsigned char 	*p,
					*q;
							
	unsigned int 	i,
					VL,
					VR,
					VLI,
					VRI;
	
	char			MyChar;
							
	int 			MyNumSepComps,
					TheComp;
							
	unsigned int 	j,
					MyNum,
					NumReps,
					V3,
					V4;
	
	/******************************************************************************************
		Each invocation of Level_Transformations() "owns" a presentation. On entry, the
		routine first checks that the presentation being passed to it is distinct from those
		passed to all previous invocations. The routine Been_Seen() is called to check this.
		We also check at this point whether the previous invocation of Level_Transformations()
		has managed to get rid of all pairs of separating vertices in the graph. If so, we save
		this new modified presentation and return 1 as a signal of success. Otherwise we save
		a copy of our private presentation on the stack of presentations SLR[][] at the
		locations SLR[MyNum][1] ... SLR[MyNum][NumRelators] and start looking for level
		transformations.
	******************************************************************************************/	
	
	if(NumCalled >= MAX_SAVED_LEVELS) return(7);
/*	if(NumCalled >= 50) return(7);	*/
	if(Been_Seen() < NumCalled) return(0);
/*	j = Been_Seen();
	if(j < NumCalled)
		{
		printf("\nPres %u --> %u.",F4+2,j+2);
		return(0);
		}			*/
	if((MyChar = mykbhit()) && !TestRealizability1)
		{
		if(MyChar == ' ')
			{
			Level_Interrupt = 1;
			return(6);
			}
		Level_Interrupt = 2;	
		}
	MyNum = NumCalled;
	switch(Sep_Pairs(0,0))
		{
		case 0:
			if(F3) return(2);
			if(NumCalled)
				{
				This_Pres = On_File();	
				if(This_Pres == NumFilled)
					{
					if(BytesUsed > BytesAvailable)
						{
						if(UserSaidQuit) return(6);
						if(UserSaidQuit = User_Says_Quit()) return(6);
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
			return(1);	
		case 1:
			BSV1[MyNum] = V1;
			BSV2[MyNum] = V2;
			break;	
		case 2:
			return(4);	
		case TOO_LONG:
			return(TOO_LONG);
		}

	if(TestRealizability1 && Planar(TRUE,FALSE)) return(3);
	
	for(i = 1; i <= NumRelators; i++)
		{
		SLR[NumCalled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if((q = *SLR[NumCalled][i]) == NULL)
			{
			for(j = 1; j < i; j++) DisposeHandle((char **) SLR[NumCalled][j]);
			return(6);
			}
		p = *Relators[i];
		while(*q++ = *p++) ;
		}
	
/*	if(SAVE_LEVEL_PRES && Save_Pres(F4+1,0,Length,1,70,1,0,0)) return(TOO_LONG);	*/
				
	NumCalled ++;
		
	/******************************************************************************************
		VL and VR are the current pair of separating vertices. We look for a component
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
		for(TheComp = 1; TheComp <= MyNumSepComps; TheComp ++)
		if(XX[VLI] != TheComp && XX[VRI] != TheComp)
			{
			switch(Find_Path(VL,VR,TheComp,1,&NumReps))
				{
				case 0:
					break;
				case TOO_LONG:
					return(6);
				case FATAL_ERROR:
					return(FATAL_ERROR);	
				case 1:
					if(Slide_ValenceTwo_Comp(TheComp,VL,VR)) break;
					if(Delete_Trivial_Generators(TRUE))
						{
						if(F5)
							return(5);
						else
							goto SLIDE_RIGHT;
						}		
					if(F3) return(4);
					if(Annulus(V3,V4,TheComp,F1) == TOO_LONG) return(TOO_LONG);
					return(4);	
				}
			if(GetHandleSize((char **) Temp8) > 5L)
				{
				for(j = 0; j < NumReps; j++)
					if(Do_Aut_L_Long(VL,VR,TheComp) == TOO_LONG) return(TOO_LONG);
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
				for(j = 0; j < NumReps; j++) if(Do_Aut_L() == TOO_LONG)
					return(TOO_LONG);
				}		
			if(Micro_Print) Micro_Print_Level_Transformations(TheComp,VL,VR,1,NumReps);	
			Fill_A(NumRelators);
			Get_Matrix();								
			j = Level_Transformations(F1,F2,F3,MyNum,F5);
			if(j > 1) return(j);
					
			/**********************************************************************************
				If a call to Level_Transformations() was unsuccessful, then we need to reset
				Relators[] and other globals that may have been changed before we try to slide
				TheComp to the "right". The following code is really a duplication of the
				previous code except that the roles of VL and VR reversed so that we are now
				sliding TheComp in the other direction.
			**********************************************************************************/	

SLIDE_RIGHT:			
			for(i = 1; i <= NumRelators; i++)
				{
				ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SLR[MyNum][i]));
				if((q = *Relators[i]) == NULL) return(TOO_LONG);
				p = *SLR[MyNum][i];
				while(*q++ = *p++) ;
				}
			if(Micro_Print) Micro_Print_Level_Transformations_Reset();			
			for(i = 0; i < Vertices; i++) XX[i] = 0;
			VL = V3;
			VR = V4;
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
			Fill_A(NumRelators);
			Get_Matrix();	
			Get_Components();
			switch(Find_Path(VR,VL,TheComp,1,&NumReps))
				{
				case 0:
					break;
				case TOO_LONG:
					return(TOO_LONG);
				case FATAL_ERROR:
					return(FATAL_ERROR);					
				case 1:
					if(Slide_ValenceTwo_Comp(TheComp,VL,VR)) break;
					if(Delete_Trivial_Generators(TRUE))
						{
						if(F5)
							return(5);
						else
							goto NEXT_COMP;
						}					
					if(F3) return(4);
					if(Annulus(V3,V4,TheComp,F1) == TOO_LONG) return(TOO_LONG);
					return(4);	
				}
			if(GetHandleSize((char **) Temp8) > 5L)
				{
				for(j = 0; j < NumReps; j++)
					if(Do_Aut_L_Long(VR,VL,TheComp) == TOO_LONG) return(TOO_LONG);
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
				for(j = 0; j < NumReps; j++) if(Do_Aut_L() == TOO_LONG)
					return(TOO_LONG);			
				}
			if(Micro_Print) Micro_Print_Level_Transformations(TheComp,VR,VL,1,NumReps);	
			Fill_A(NumRelators);
			Get_Matrix();								
			j = Level_Transformations(F1,F2,F3,MyNum,F5);
			if(j > 1) return(j);
				
			/**********************************************************************************
				If a call to Level_Transformations() was unsuccessful, then we need to reset
				Relators[] and other globals that may have been changed before we try to slide
				other components of the initial separation around.
			**********************************************************************************/	

NEXT_COMP:			
			for(i = 1; i <= NumRelators; i++)
				{
				ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SLR[MyNum][i]));
				if((q = *Relators[i]) == NULL) return(TOO_LONG);
				p = *SLR[MyNum][i];
				while(*q++ = *p++) ;
				}
			if(Micro_Print) Micro_Print_Level_Transformations_Reset();			
			for(i = 0; i < Vertices; i++) XX[i] = 0;
			VL = V3;
			VR = V4;
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
			Fill_A(NumRelators);
			Get_Matrix();	
			Get_Components();
			}
		
		if(MyNumSepComps == 2)
			{
			j = Level_Transformations_2(F1,F2,F3,MyNum,F5);
			if(j > 1) return(j);
			}
		
		/**************************************************************************************
			When Sep_Pairs(i,j) is called, it starts searching for separating pairs of
			vertices in the graph beginning with the ordered pair of vertices following the
			pair (i,j) --(We always have j > i.)-- If Sep_Pairs(i,j) finds a pair of
			separating vertices, it returns that pair in the globals V1 and V2. So when the
			routine Level_Transformations() has dealt with the pair of separating vertices
			saved in BSV1[MyNum] and BSV2[MyNum], we call Sep_Pairs((BSV1[MyNum],BSV2[MyNum])
			to see if the graph has any other pairs of separating vertices. If there are none,
			we return(0), otherwise we update BSV1[MyNum] and BSV2[MyNum] and go to the top
			of this loop.
		**************************************************************************************/				
		
		if(Sep_Pairs(BSV1[MyNum],BSV2[MyNum]) == FALSE || MajorVert == 2)
			return(0);
		BSV1[MyNum] = V1;
		BSV2[MyNum] = V2;	
		}
}

Get_Components()
{
	register unsigned int g,h,i,j,k;
	register unsigned int *p,*r;
	
	/******************************************************************************************
		If the "reduced" Whitehead graph has a pair of separating vertices, this routine can
		be used to find all of the components resulting from the deletion of the separating
		vertices. The array XX[] is initialized by the calling routine, which sets all of the
		entries to zero except for the two entries corresponding to the separating vertices.
		These entries are set to the value VERTICES. Upon return, the entries of XX[] which
		are equal to g are the vertices in the gth component of the graph.
	******************************************************************************************/	
	
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

Been_Seen()
{
	/******************************************************************************************
			This routine compares the presentation, which we are about to pass to
			Level_Transformations(), to all of the presentations passed in previous calls.
			If there is a duplication, then we do not need or want to make this call.
			Otherwise it could be a long wait. If there is a match the routine returns the
			number of the matched presentation. Otherwise it returns NumCalled.
	******************************************************************************************/
		
	register int i,j;
	
	for(j = 1; j <= NumRelators; j++) HLock((char **) Relators[j]);
	for(i = 0; i < NumCalled; i++)
		{
		for(j = 1; j <= NumRelators; j++)
			{
			if(GetHandleSize((char **) Relators[j]) == GetHandleSize((char **) SLR[i][j]))
				{
				HLock((char **) SLR[i][j]);
				if(Compare_Str(*Relators[j],*SLR[i][j],
					GetHandleSize((char **) Relators[j]) - 1) == FALSE)
					{
					HUnlock((char **) SLR[i][j]);
					break;
					}
				HUnlock((char **) SLR[i][j]);
				}	
			else
				break;			
			}			
		if(j > NumRelators) break;							
		}
	for(j = 1; j <= NumRelators; j++) HUnlock((char **) Relators[j]);			
	return(i);		
}

Annulus(V3,V4,TheComp,F1)
unsigned int 	V3,
				V4;
int 			TheComp,
				F1;
{
	/******************************************************************************************
		This routine is called by Level_Transformations() when Level_Transformations() has
		determined that there is an annulus present. V3 and V4 are a pair of vertices that
		separate the Heegaard diagram. TheComp gives us a set of vertices which are 
		"swallowed" by the annulus; while Temp8 contains a string which represents the path
		on the Heegaard surface which serves as the center line of the rest of the annulus.
	******************************************************************************************/
		
	register unsigned char  *p,
							*q,
							x,
							y;
							
	register unsigned int 	i,
							j;
			
	if(V3 & 1)
		x = V3/2 + 97;
	else
		x = V3/2 + 65;
	if(V4 & 1)
		y = V4/2 + 97;
	else
		y = V4/2 + 65;	
	ReallocateHandle((char **) Relators[0],GetHandleSize((char **) Temp8) + Vertices + 2);
	if((p = *Relators[0]) == NULL) return(TOO_LONG);
	*p++ = x;
	*p++ = y;
	for(i = 0; i < Vertices; i++)
		{
		if(XX[i] == TheComp)
			{
			if(i & 1)
				y = i/2 + 97;
			else
				y = i/2 + 65;
			*p++ = y;	
			}
		}
	*p++ = '@';
	q = *Temp8;
	while(*p++ = *q++) ;
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
				SUR[i][0] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[0]));
				if((p = *SUR[i][0]) == NULL) return(TOO_LONG);
				q = *Relators[0];
				while(*p++ = *q++) ;
				BytesUsed += GetHandleSize((char **) Relators[0]);
				UDV[i] = ANNULUS_EXISTS;
				if(F1) printf("\n					There is an annulus in diagram %u.",i + 1);
				}
	 		break;
	 		}
	 	}						
	if(i == NumFilled)
		{
		if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);
		if(Save_Pres(ReadPres,0,Length,1,70,0,0,0)) return(TOO_LONG);		
		SUR[NumFilled - 1][0] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[0]));
		if((p = *SUR[NumFilled - 1][0]) == NULL) return(TOO_LONG);
		q = *Relators[0];
		while(*p++ = *q++) ;
		BytesUsed += GetHandleSize((char **) Relators[0]);
		BDY[NumFilled - 1] = BDY[ReadPres];
		UDV[NumFilled - 1] = ANNULUS_EXISTS;
		if(F1) printf("\n					There is an annulus in diagram %u.",NumFilled);
		}
	return(NO_ERROR);			
}

int Find_Path(int V1,int V2,int TheComp,int Type,unsigned int *NumRepsPtr)
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
							T2[125];
								
	int						FP,
							i,
							j;

	long					max;
	
	unsigned long			HS,
							L,
							L1;
								
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

	if(Type == 2)
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
	
	ReallocateHandle((char **) Temp8,max + 1);	
	if(*Temp8 == NULL) return(TOO_LONG);
	
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
		while(y = *p++)
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
		HLock((char **) Relators[i]);
		Inverse(*Relators[i]);
		HUnlock((char **) Relators[i]);			
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
	if(T2[x] == 1)
		{
		SetHandleSize((char **) Temp8,HS + 1);
		*NumRepsPtr = 1;		
		return(1); 					/*	An annulus exists. */
		}
	if(V2 & 1)
		y = (V2 >> 1) + 65;
	else
		y = (V2 >> 1) + 97;
	q--;
	x = *q;			
	if(x != y || Type == 2)
		{
		SetHandleSize((char **) Temp8,HS + 1);
		*NumRepsPtr = 1;		
		return(0);
		}
	
	/******************************************************************************************
		If execution gets to this point, the path in Temp8 corresponds to sliding TheComp
		around a circuit that brings TheComp back to again lie between V2 and V1. We will
		next look for the smallest proper power of *Temp8 which starts at a vertex of TheComp
		and ends at a vertex which is not part of TheComp. This corresponds to the situation
		where TheComp can be pushed around a spiral path which represents a proper power of
		Temp8.
	******************************************************************************************/	
	
	for(j = 0,L1 = BIG_NUMBER; j < 2; j++)
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		q = p + GetHandleSize((char **) Relators[i]) - 2;
		x = *q;
		if(L1 > 1L) while(y = *p++)
			{
			if((y == z) && T1[x])
				{
				for(L = 0L,q = *Temp8 + 1; *p ; p++, q++)
					{
					if(*q == EOS)
						{
						q = *Temp8;
						L++;	
						}
					if(*q != *p)
						{
						if(T2[*p] == 1) break;
						if(L < L1)
							L1 = L;
						break;		/* Done with this relator. */
						}
					}
				if(*p == EOS)
					{
					for(p = *Relators[i]; *p; p++, q++)
						{
						if(*q == EOS)
							{
							q = *Temp8;
							L++;	
							}
						if(*q != *p)
							{
							if(T2[*p] == 1) break; /* Disregard paths which end in TheComp. */
							if(L < L1) L1 = L;
							break;	/* Done with this relator. */
							}
						}
					break;	
					}						
				q = p - 1;
				y = *q;
				}
			x = y;	
			}
		HLock((char **) Relators[i]);
		Inverse(*Relators[i]);
		HUnlock((char **) Relators[i]);			
		}
		
	if(L1 == 1L || L1 == BIG_NUMBER)
		{
		SetHandleSize((char **) Temp8,HS + 1);
		*NumRepsPtr = 1;
		return(0);		
		}
	
	SetHandleSize((char **) Temp8,HS + 1);
	*NumRepsPtr = L1;
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
							TY[125],
							**Temp;
							
	int						Vertex;
	
	unsigned int			j,
							Auts;
							
	long					HS;												
	
	Auts = GetHandleSize((char **) Temp8) - 1;
	TotalAuts += Auts;
	
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
			ReallocateHandle((char **) Temp5,2*HS);
			if((p = *Temp5) == NULL) return(TOO_LONG);		
			q = *Relators[i];
			x = *q++;
			while(y = *q++)
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
			SetHandleSize((char **) Temp5,HS);
			Temp = Relators[i];
			Relators[i] = Temp5;
			Temp5 = Temp;					
			}
		}
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
	
	unsigned char 			T1[125],
							T2[125],
							**Temp;
	
	int						i,
							j;
	
	unsigned long			Delete,
							Insert,
							HS,
							Li,
							L8;
	
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
			while(y = *p)
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
			ReallocateHandle((char **) Temp5,Li);
			if((r = *Temp5) == NULL) return(TOO_LONG);	
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
			while(y = *p)
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
				SetHandleSize((char **) Temp5,HS);
				Temp 		= Relators[i];
				Relators[i] = Temp5;
				Temp5 		= Temp;
				}	
			}
		HLock((char **) Relators[i]);
		Inverse(*Relators[i]);
		HUnlock((char **) Relators[i]);
		}
	
	/*****************************************************************************************
						Put a copy of the inverse of Temp8 in Temp9.
	*****************************************************************************************/
	
	HS = L8 + 1;
	if(HS > MAXLENGTH) return(TOO_LONG);
	ReallocateHandle((char **) Temp9,HS);
	if((q = *Temp9) == NULL) return(TOO_LONG);
	p = *Temp8;
	while(*q++ = *p++) ;
	HLock((char **) Temp9);
	Inverse(*Temp9);
	HUnlock((char **) Temp9);
	
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
		while(y = *p++)
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
		ReallocateHandle((char **) Temp5,HS);
		if((r = *Temp5) == NULL) return(TOO_LONG);
		
		/*************************************************************************************
									Rewrite Relators[i].
		*************************************************************************************/
			
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while(y = *p++)
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64)
				{
				q = *Temp8;
				while(*r++ = *q++) ;
				r--;
				x = EOS;
				continue;	
				}
			if(y == 96)
				{
				q = *Temp9;
				while(*r++ = *q++) ;
				r--;
				x = EOS;
				continue;	
				}
			if((x == z) && T2[y])
				{
				q = *Temp8;
				while(*r++ = *q++) ;
				r--;	
				}
			else	
			if((y == w) && T1[x])
				{
				q = *Temp9;
				while(*r++ = *q++) ;
				r--;
				}
			*r++ = y;	
			x = y;	
			}
		*r++ = EOS;
		HS = r - *Temp5;
		if(HS > MAXLENGTH) return(TOO_LONG); 
		SetHandleSize((char **) Temp5,HS);
		Temp 		= Relators[i];
		Relators[i] = Temp5;
		Temp5 		= Temp;			
		}
	
	TotalAuts += L8;
	return(NO_ERROR);			
}

Delete_Trivial_Generators(Test)
int		Test;
{
	register unsigned char	*p,
							*q,
							s,
							t,
							x,
							y,
							z;
	
	unsigned char			**Temp;
							
	int						i,
							j,
							k,
							SNumRelators1,
							SNumRelators2,
							TheRelator;
							
	long					HS;						
	
	/***************************************************************************************
		This routine deletes a generator X from the relators, provided there is a relator
		of the form X, and the total number of appearances of X, together with its inverse,
		in the remaining relators is no more than one.
	***************************************************************************************/

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
			while(z = *p++) if(z == x || z == y)
				{
				k++;
				if(k > 1) goto TEST_NEXT_RELATOR;
				if(k == 1) TheRelator = j;
				}
			}
			
		if(k == 1)
			{
			HS = GetHandleSize((char **) Relators[TheRelator]);
			if(HS == 2L)
				{
				/**************************************************************************
					If k = 1 and HS = 2, the presentation contains two relators of the
					form X, and deleting X would create an empty relator. We do not want
					to do this at this point.
				**************************************************************************/
						
				continue;
				}
			if(Test) return(TRUE);		
			p = *Relators[TheRelator];
			q = p;
			while(z = *p++) if(z != x && z != y) *q++ = z;
			*q = EOS;
			SetHandleSize((char **) Relators[TheRelator],HS - 1);
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
				printf("\nSwapped Relator %d with Relator %d and reduced the number of relators.",i,NumRelators);
			if(Micro_Print_F)
				{
				fprintf(myout,"\n\nRelator %d is trivial. ",i);
				fprintf(myout,"Deleted generator %c from the presentation.",x);
				if(i != NumRelators)
					fprintf(myout,"\nSwapped Relator %d with Relator %d and reduced the number of relators.",i,NumRelators);
				}
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
				while(z = *p)
					{
					if(z == s) *p = x;
					if(z == t) *p = y;
					p++;
					}
				}
			if(Micro_Print)
				{
				printf("\n\nRewrote the relators by replacing %c with %c.",s,x);
				if(Micro_Print_F)
					fprintf(myout,"\n\nRewrote the relators by replacing %c with %c.",s,x);
				}	
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

int Level_Transformations_2(int F1,int F2,int F3,unsigned int MyNum,int F5)
{
	/******************************************************************************************
			Suppose {X,Y} is a pair of separating vertices of the "reduced" Whitehead graph G
		of a minimal length presentation, with x the inverse of X and y the inverse of Y.
		Such separating pairs of vertices are of two types. Say {X,Y} is a "Type 1" pair
		if Y = x, or if deleting X and Y from G yields a component C which contains both x
		and y. If no such component exists, say that {X,Y} is a "Type 2" pair.
			We mention that if {X,Y} is a Type 2 pair, then the separation of G, produced by
		deleting X and Y from G, has only two components, and there is no edge joining X and
		Y in G.
			Level_Transformations() calls Level_Transformations_2() when the separation has
		two components.
			F1 and F2 are flags which are passed along by Level_Transformations_2() to any
		new invocations of Level_Transformations(). If F1 is FALSE, then printing of
		information about any annuli the routine finds is  suppressed. If F2 is TRUE, new
		invocations of Level_Transformations() will stop and return 2 when they find the
		first new presentation without any pairs of separating vertices.	
	******************************************************************************************/
	
	unsigned char 	*p,
					*q;
							
	unsigned int 	i,
					VL,
					VR,
					VLI,
					VRI;
							
	int 			TheComp;
							
	unsigned int 	j,
					NumReps;
		
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
	if(NumSepComps != 2 && !TestRealizability1) return(0);
	if(NumSepComps == 2)
	for(TheComp = 1; TheComp <= 2; TheComp ++)
		{
		if(XX[VLI] != TheComp && XX[VRI] == TheComp)
			{
			switch(Find_Path(VL,VR,TheComp,2,&NumReps))
				{
				case 0:
					break;
				case TOO_LONG:
					return(6);
				case FATAL_ERROR:
					return(FATAL_ERROR);						
				case 1:
					if(Delete_Trivial_Generators(TRUE))
						{
						if(F5)
							return(5);
						else
							goto SLIDE_RIGHT;
						}		
					XX[VR] = TheComp;
					if(Annulus(VL,VR,TheComp,F1) == TOO_LONG) return(TOO_LONG);
					return(4);	
				}
			if(GetHandleSize((char **) Temp8) > 5L)
				j = Do_Aut_L_Long_2(VL,VR,TheComp);
			else
				{
				for(i = 0; i < Vertices; i++)
					{
					if(XX[i] == TheComp)
						YY[i] = 1;
					else
						YY[i] = 0;
					}
				YY[VR] = 1;	
				j = Do_Aut_L();
				}					
			if(j == TOO_LONG) return(TOO_LONG);
			if(Micro_Print)
				{
				XX[VR] = TheComp;
				Micro_Print_Level_Transformations(TheComp,VL,VR,2,NumReps);
				XX[VR] = VERTICES;
				}
			Fill_A(NumRelators);
			Get_Matrix();								
			j = Level_Transformations(F1,F2,F3,MyNum,F5);
			if(j > 1) return(j);
					
			/******************************************************************************
				If a call to Level_Transformations() was unsuccessful, then we need to
				reset Relators[] and other globals that may have been changed before we
				try to slide the other component together with VL to the "right".
			******************************************************************************/	

SLIDE_RIGHT:			
			for(i = 1; i <= NumRelators; i++)
				{
				ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SLR[MyNum][i]));
				if((q = *Relators[i]) == NULL) return(TOO_LONG);
				p = *SLR[MyNum][i];
				while(*q++ = *p++) ;
				}
			if(Micro_Print) Micro_Print_Level_Transformations_Reset();
			for(i = 0; i < Vertices; i++) XX[i] = 0;
			XX[VL] = VERTICES;
			XX[VR] = VERTICES;
			Fill_A(NumRelators);
			Get_Matrix();	
			Get_Components();
			}
			
		if(XX[VLI] == TheComp && XX[VRI] != TheComp)			
			{
			switch(Find_Path(VR,VL,TheComp,2,&NumReps))
				{
				case 0:
					break;
				case TOO_LONG:
					return(TOO_LONG);
				case FATAL_ERROR:
					return(FATAL_ERROR);						
				case 1:
					if(Delete_Trivial_Generators(TRUE))
						{
						if(F5)
							return(5);
						else
							goto NEXT_COMP;
						}		
					XX[VL] = TheComp;
					if(Annulus(VR,VL,TheComp,F1) == TOO_LONG) return(TOO_LONG);
					return(4);	
				}
			if(GetHandleSize((char **) Temp8) > 5L)
				j = Do_Aut_L_Long_2(VR,VL,TheComp);
			else
				{
				for(i = 0; i < Vertices; i++)
					{
					if(XX[i] == TheComp)
						YY[i] = 1;
					else
						YY[i] = 0;
					}
				YY[VL] = 1;					
				j = Do_Aut_L();
				}
			if(j == TOO_LONG) return(TOO_LONG);
			if(Micro_Print)
				{
				XX[VL] = TheComp;
				Micro_Print_Level_Transformations(TheComp,VR,VL,2,NumReps);
				XX[VL] = VERTICES;
				}
			Fill_A(NumRelators);
			Get_Matrix();								
			j = Level_Transformations(F1,F2,F3,MyNum,F5);
			if(j > 1) return(j);
				
			/******************************************************************************
				If a call to Level_Transformations() was unsuccessful, then we need to
				reset Relators[] and other globals that may have been changed before we
				try to slide the other component of the separation around.
			******************************************************************************/	

NEXT_COMP:			
			for(i = 1; i <= NumRelators; i++)
				{
				ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SLR[MyNum][i]));
				if((q = *Relators[i]) == NULL) return(TOO_LONG);
				p = *SLR[MyNum][i];
				while(*q++ = *p++) ;
				}
			if(Micro_Print) Micro_Print_Level_Transformations_Reset();	
			for(i = 0; i < Vertices; i++) XX[i] = 0;
			XX[VL] = VERTICES;
			XX[VR] = VERTICES;
			Fill_A(NumRelators);
			Get_Matrix();	
			Get_Components();
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
	
	unsigned char 			T1[125],
							T2[125],
							**Temp;
	
	int						i,
							j;
	
	unsigned long			Delete,
							Insert,
							Li,
							L8;
	
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
			while(y = *p)
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
			ReallocateHandle((char **) Temp5,Li);
			if((r = *Temp5) == NULL) return(TOO_LONG);	
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
			while(y = *p)
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
				SetHandleSize((char **) Temp5,r - *Temp5);
				Temp 		= Relators[i];
				Relators[i] = Temp5;
				Temp5 		= Temp;
				}	
			}
		HLock((char **) Relators[i]);
		Inverse(*Relators[i]);
		HUnlock((char **) Relators[i]);
		}
	
	/*****************************************************************************************
						Put a copy of the inverse of Temp8 in Temp9.
	*****************************************************************************************/
		
	ReallocateHandle((char **) Temp9,L8 + 1);
	if((q = *Temp9) == NULL) return(TOO_LONG);
	p = *Temp8;
	while(*q++ = *p++) ;
	HLock((char **) Temp9);
	Inverse(*Temp9);
	HUnlock((char **) Temp9);
	
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
		while(y = *p++)
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
				
		ReallocateHandle((char **) Temp5,Li + Insert*L8);
		if((r = *Temp5) == NULL) return(TOO_LONG);
		
		/*************************************************************************************
									Rewrite Relators[i].
		*************************************************************************************/
			
		p = *Relators[i];
		q = p + Li - 2;
		x = *q;
		while(y = *p++)
			{
			if(y == 63 || y == 95)
				{
				x = EOS;
				continue;
				}
			if(y == 64)
				{
				q = *Temp8;
				while(*r++ = *q++) ;
				r--;
				x = EOS;
				continue;	
				}
			if(y == 96)
				{
				q = *Temp9;
				while(*r++ = *q++) ;
				r--;
				x = EOS;
				continue;	
				}
			if((x == z) && T2[y])
				{
				q = *Temp9;
				while(*r++ = *q++) ;
				r--;	
				}
			else	
			if((y == w) && !T1[x])
				{
				q = *Temp8;
				while(*r++ = *q++) ;
				r--;
				}
			*r++ = y;	
			x = y;	
			}
		*r++ = EOS;
		SetHandleSize((char **) Temp5,r - *Temp5);
		Temp 		= Relators[i];
		Relators[i] = Temp5;
		Temp5 		= Temp;			
		}
	
	TotalAuts += L8;
	return(NO_ERROR);			
}

int Slide_ValenceTwo_Comp(int TheComp,unsigned int VL,unsigned int VR)
{
int				i;

unsigned char	x,
				j,
				*p;

	/******************************************************************************************
							Check if each vertex in TheComp has valence two.
	******************************************************************************************/
							

	for(i = 0; i < Vertices; i++)
	if(XX[i] == TheComp && VWG[i] != 2) return(0);

	/******************************************************************************************
			Check that TheComp is joined to ech of VL and VR with only one edge.
	******************************************************************************************/
			

	for(i = 0; i < Vertices; i++)
	if(XX[i] == TheComp && (A[i][VL] > 1 || A[i][VR] > 1)) return(0);
	
	/******************************************************************************************
	Check if there is a pair of vertices in the Whitehead graph which are joined only
	by an edge of the annulus which swallows TheComp. If there is, we will slide TheComp
	until it joins these two vertices.												
	******************************************************************************************/
	
	
	p = *Temp8;
	if(*p == EOS) return(0);
	x = *p << 1;
	if(x < 194) x -= 129;
	else x -= 194;
	j = x;
	p++;
	while(x = *p)
		{
		x = x << 1;
		if(x < 194) x -= 130;
		else x -= 193;
		if(A[j][x] == 1)	/* 	We have found the edge we are looking for.	*/
			{
			*p = EOS;
			SetHandleSize((char **) Temp8,p - *Temp8 + 1);
			return(1);
			}
		if(x & 1) j = x - 1;
		else j = x + 1;
		p++;
		}
	return(0);	
}				
