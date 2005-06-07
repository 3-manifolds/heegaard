#include "Heegaard.h"
#include "Heegaard_Dec.h"																																											 

unsigned long 	length;

unsigned int Valence_Two(F1)
int		F1;
{
	/******************************************************************************************
		Here is a brief description of this routine which attempts to resolve ambiguities
	about how edges should be identified, in the Heegaard diagram, for a pair of vertices
	both of which have valence two in the reduced Whitehead graph of a presentation.
	
		Suppose P is a minimal length presentation with Whitehead graph WG, and reduced
	Whitehead graph RWG. Suppose that RWG has no pairs of separating vertices, so that RWG
	and WG have unique embeddings in the plane. Suppose that both vertex X of RWG and the
	the inverse vertex x of X in RWG, have valence two in RWG. And, suppose that the program
	has determined that there are two possible ways to identify the edges of WG, which meet
	X with the edges of WG, which meet x. Let B be one of the bands of parallel edges in WG,
	which meet vertex X, and suppose the edges in B do not join vertex X to vertex x. We note
	that, since there is an ambiguity about how these identifications should be made, B
	contains at least two edges.
		The program has computed the two possible ways that the edges meeting vertex X can
	be identified with the edges meeting vertex x, and has stored these possibilities in two
	arrays, which we will call id1[] and id2[]. Next, the program assumes that the
	identification in id1[] is correct. Using the information in id1[], the program can 
	trace the left-hand edge of B, and the right-hand edge of B, backward over vertex X, and
	it can follow the paths taken by the left-hand and right-hand edges of B until they leave
	vertex x, and split apart, with each path going to a distinct vertex of WG. (If this
	"splitting" did not occur, then there would not be any ambiguity about how edges should
	be identified at X and x.)
		Suppose the splitting takes place so that one path leaves vertex x and goes to vertex
	a of WG, while the other path leaves vertex x and goes to vertex b of WG. The program now
	traces these paths in the opposite direction, starting from vertex a in one case, and
	starting from vertex b, in the other. The program continues to trace these two paths
	through the diagram, using known information about identifications at each new vertex
	reached, until one of the following occurs:
		
		1) 	The two paths traced out have each become longer than the longest relator in P. In
			this case, P is not realizable, and the program reports that fact.
		
		2)	The two paths diverge, and proceed to distinct vertices of WG.
		
		3)	The two paths have reached a vertex say Y, which is one of a pair of vertices of
			valence two in RWG, and the identification of the edges of WG meeting Y with the
			edges of WG meeting y is still ambiguous.
			
		Suppose 2) occurs, and suppose that the path which originated at vertex a goes to
	vertex C of WG, while the path that originated at vertex b goes to vertex D of WG. Thus,
	if P is realizable, and the identification in id1[] is correct for the pair {X,x}, then
	subwords of the form AWC, and BWD must occur in the relators of P, where W is the "word"
	traced out by the two paths up to the point where they diverged. While, if the
	identification in id2[] is correct for the pair {X,x}, then	subwords of the form BWC, and
	AWD must occur in the relators of P, where W is the "word" traced out by the two paths up
	to the point where they diverged.
		Now it is not hard to show that, under these circumstances, if P is realizable, then
	at least two and at most 3 of the the 4 subwords  {AWC, BWD, AWD, BWC} can occur in the
	relators of P. Thus, by scanning the relators in P, the program can determine which of
	these possibilities occurs, and then it can determine which identification, of the edges
	meeting X with the edges meeting x, is correct. This reduces the number of pairs of
	vertices with ambiguous identifications, and the program can now use this information to
	try to resolve ambiguous identifications of other pairs of vertices.


								VALENCE-TWO-ANNULI IN DIAGRAMS.
								
		If case 3) occurs, the program marks the band B as having been traced, and looks for
	another band of edges which it has not traced. If the program cannot find any band of
	edges, which it has not traced, and which leads to case 1), or case 2), then the diagram
	contains a cycle of valence two vertices, which form a "Valence-Two-Annulus". In this
	situation, there is an annulus A in the underlying handlebody H, such that: ¶A swallows
	an ambiguous vertex of valence two, then ¶A follows the left and right hand edges of the
	band of parallel edges, which leave this vertex, to another ambiguous vertex of valence
	two, swallows that vertex, and continues, in this manner, until it returns to the original
	vertex of valence two.
		If a Valence-Two-Annulus exists, then the program can arbitrarily choose either of
	the two possible ways to identify the edges meeting one pair of valence-two vertices,
	provided at least one vertex of the pair lies in the cycle which constitutes the annulus,
	and then this choice determines the identifications of all other pairs of ambiguous
	valence-two vertices, for which at least one vertex of the pair lies in the cycle. By
	repeating this step, the program can eventually determine the identifications of all
	ambiguous pairs of valence-two vertices.
		There are a couple of degenerate situations which can occur. One possibility is that
	all of the vertices in RWG are of valence-two, and the annulus in question has swallowed
	everything. In this case, we can ignore the annulus and proceed. Another possibility is
	that the annulus in question has swallowed all but two vertices of RWG. In this case, the
	omitted vertices form a pair of inverse vertices, and again the annulus does not cause
	a problem, so we can ignore it and proceed. (Observe that in this case, ¶A represents a
	power of a free generator of P.)
		
	******************************************************************************************/		
		
	register unsigned int 	i,
							j;
	
	int						AnnulusExists;
							
	unsigned long			max;
	
	AnnulusExists = FALSE;
	
	for(i = 0; i < Vertices; i++) if(GV2[i >> 1])
		{
		GV2L[i] = FV[i];
		GV2R[i] = CO[i][FV[i]];
		}
	for(i = 1,max = 0L; i <= NumRelators; i++) if(LR[i] > max) max = LR[i];
	max += 4;
	ReallocateHandle((char **) Temp9,max);
		if(*Temp9 == NULL) return(TOO_LONG);
	ReallocateHandle((char **) Temp10,max);
		if(*Temp10 == NULL) return(TOO_LONG);	
	while(NGV2)
		{
		for(i = 0, j = NGV2; i < Vertices; i++)
			{
			if(GV2[i >> 1] && !GV2[GV2L[i] >> 1])
				switch(Resolve(i,FV[i],max))
					{
					case 0:
						GV2L[i] = Stopper;
						break;
					case 1:
						/**********************************************************************
										The presentation is not realizable.
						**********************************************************************/	
						
						return(FATAL_ERROR);
					case 2:
						NGV2 --;
						GV2[i >> 1] = FALSE;
						Test_Sub_Str(i >> 1);
						break;	
					}
			if(GV2[i >> 1] && !GV2[GV2R[i] >> 1])
				switch(Resolve(i,CO[i][FV[i]],max))
					{
					case 0:
						GV2R[i] = Stopper;
						break;
					case 1:
						/**********************************************************************
										The presentation is not realizable.
						**********************************************************************/
							
						return(FATAL_ERROR);
					case 2:
						NGV2 --;
						GV2[i >> 1] = FALSE;
						Test_Sub_Str(i >> 1);
						break;		
					}
			}	
		if(j == NGV2)
			{
			/**********************************************************************************
				If j = NGV2, then the program did not manage to resolve any of the ambiguous
				vertices of valence two. Call Valence_Two_Annulus() to see whether there is a
				nondegenerate annulus present.
			**********************************************************************************/
			
			if(F1) return(NO_ERROR);
			if(NGV2 < NumGenerators) switch(Valence_Two_Annulus())
				{
				case NO_ERROR:
					break;
				case TOO_LONG:
					return(TOO_LONG);
				case TRUE:
					AnnulusExists = TRUE;
					break;	
				}	
			for(i = 0; i < NumGenerators; i++) if(GV2[i])
				{
				NGV2 --;
				GV2[i] = FALSE;
				break;
				}
			}				
		}
	if(AnnulusExists) return(V2_ANNULUS_EXISTS);		
	return(NO_ERROR);			
}

Resolve(i,j,max)
unsigned int	i,
				j;				
unsigned long	max;
{
	/******************************************************************************************
		This routine is called by Valence_Two(). It traces out two paths and creates two
		strings which correspond to these paths. In the case where these paths run parallel 
		to each other and then diverge, we get strings of the form: A........B and C........D
		which are of equal length and differ only in their initial and terminal letters.
		Temp9 and Temp10 are handles to these two strings.
	******************************************************************************************/
		
	register unsigned char 	*p,
							*q;
							
	register unsigned int 	e1,
							e2,
							v,
							w,
							d;
							
	unsigned char 			x;

	/******************************************************************************************
								Backup until we get a "split".
	******************************************************************************************/
	
	p = *Temp9;
	q = *Temp10;
	if(j == FV[i])
		{
		e1 = 0;
		e2 = A[i][j] - 1;
		}
	else
		{
		e1 = A[i][FV[i]];
		e2 = V[i] - 1;
		}
	if(e1 == e2) return(1);
	
	/******************************************************************************************
		The case where e1 = e2 is an error!!!. However, if this occurs this routine will "hang"
		unless we do something, so we "fake" success, return and let somebody else deliver
		the bad news.
	******************************************************************************************/
	
	v = i;
	length = 0L;
	max -= 2;
	while(1)
		{
		e1 = OSA[v] - e1;
		if(e1 >= V[v])
			e1 -= V[v];
		e2 = OSA[v] - e2;
		if(e2 >= V[v])
			e2 -= V[v];
		if(v & 1)
			w = v - 1;
		else
			w = v + 1;	
		v = FV[w];
		d = A[w][v];
		if(++length > max) return(1);
		if(e1 <= e2)
			{
			while(d <= e1)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			if(d <= e2)
				{
				if(v & 1)
					x = (v >> 1) + 65;
				else
					x = (v >> 1) + 97;
				*p++ = x;
				while(d <= e2)
					{
					v = CO[w][v];
					d += A[w][v];
					}
				if(v & 1)
					x = (v >> 1) + 65;
				else
					x = (v >> 1) + 97;
				*q++ = x;			
				break;
				}
			}
		else
			{
			while(d <= e2)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			if(d <= e1)
				{
				if(v & 1)
					x = (v >> 1) + 65;
				else
					x = (v >> 1) + 97;
				*q++ = x;
				while(d <= e1)
					{
					v = CO[w][v];
					d += A[w][v];
					}
				if(v & 1)
					x = (v >> 1) + 65;
				else
					x = (v >> 1) + 97;
				*p++ = x;			
				break;
				}
			}	
		e1 = B[w][v] - e1;
		e2 = B[w][v] - e2;							
		}
		
	/******************************************************************************************
	 						Go forward until we get to vertex j.
	 *****************************************************************************************/
	 
	v = w;
	length = 1L;
	while(1)
		{
		if(v == j) break;
		if(v & 1)
			{
			x = (v >> 1) + 97;
			*p++ = x;
			*q++ = x;
			w = v - 1;
			}
		else
			{
			x = (v >> 1) + 65;
			*p++ = x;
			*q++ = x;
			w = v + 1;
			}
		if(++length > max) return(1);	
		e1 = OSA[v] - e1;
		if(e1 >= V[v])
			e1 -= V[v];
		e2 = OSA[v] - e2;
		if(e2 >= V[v])
			e2 -= V[v];
		v = FV[w];
		d = A[w][v];
		if(e1 <= e2)
			{
			while(d <= e1)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			}
		else
			{
			while(d <= e2)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			}	
		e1 = B[w][v] - e1;
		e2 = B[w][v] - e2;							
		}
		
	/******************************************************************************************
	 	Continue forward until we get a split or we reach an ambiguous vertex of valence two,
	 					or we run out of space in Temp9 and Temp10.
	******************************************************************************************/
		
	while(1)	
		{
		if(GV2[v >> 1])
			{
			/**********************************************************************************
				We have reached an ambiguous vertex of valence two. Label this vertex as the
				Stopper and return.
			**********************************************************************************/
				
			Stopper = v;
			return(0);
			}
		if(v & 1)
			{
			x = (v >> 1) + 97;
			*p++ = x;
			*q++ = x;
			w = v - 1;
			}
		else
			{
			x = (v >> 1) + 65;
			*p++ = x;
			*q++ = x;
			w = v + 1;
			}
		e1 = OSA[v] - e1;
		if(e1 >= V[v])
			e1 -= V[v];
		e2 = OSA[v] - e2;
		if(e2 >= V[v])
			e2 -= V[v];	
		if(++length > max)
			{
			/**********************************************************************************
				We have followed this pair of parallel paths until their length has become
				more than two characters longer than the length of the longest relator in the
				presentation. This is impossible if the presentation is realizable. Return 1
				whereupon Valence_Two() will return FATAL_ERROR.
			**********************************************************************************/	
			
			*p = EOS;
			*q = EOS;
			return(1);
			}		
		v = FV[w];
		d = A[w][v];
		if(e1 <= e2)
			{
			while(d <= e1)
				{
				v = CO[w][v];
				d += A[w][v];
				}	
			if(d <= e2)
				{
				/******************************************************************************
					The two paths we have been following have finally diverged and gone to
					different vertices. Determine what these two vertices are, record them in
					Temp9 and Temp10 and return.
				******************************************************************************/
					
				if(v & 1)
					x = (v >> 1) + 97;
				else
					x = (v >> 1) + 65;
				*p++ = x;
				while(d <= e2)
					{
					v = CO[w][v];
					d += A[w][v];
					}
				if(v & 1)
					x = (v >> 1) + 97;
				else
					x = (v >> 1) + 65;
				*q++ = x;
				*p = EOS;
				*q = EOS;
				length ++;		
				return(2);
				}	
			}
		else
			{
			while(d <= e2)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			if(d <= e1)
				{
				/******************************************************************************
					The two paths we have been following have finally diverged and gone to
					different vertices. Determine what these two vertices are, record them in
					Temp9 and Temp10 and return.
				******************************************************************************/
					
				if(v & 1)
					x = (v >> 1) + 97;
				else
					x = (v >> 1) + 65;
				*q++ = x;
				while(d <= e1)
					{
					v = CO[w][v];
					d += A[w][v];
					}
				if(v & 1)
					x = (v >> 1) + 97;
				else
					x = (v >> 1) + 65;
				*p++ = x;
				*p = EOS;
				*q = EOS;
				length ++;							
				return(2);
				}	
			}	
		e1 = B[w][v] - e1;
		e2 = B[w][v] - e2;							
		}	
}

Test_Sub_Str(j)
unsigned int j;
{
	/*****************************************************************************************
		Test_Sub_Str(j) is called to determine which of the two possible alternative values
		of "offset" is correct that in OSA[2j] or that in OSB[2j]. It is not hard to see that
		if the presentation in Relators[] is realizable, (and there is no annulus present)
		then only one of these alternatives can be correct. If both of the strings Temp9
		and Temp10 appear as substrings of proper powers of the relators or of proper powers
		of cyclic conjugates of the relators, then OSA[2j] must be correct. Otherwise OSB[2j]
		must be correct. So this routine calls Sub_Str() to determine whether Temp9 and
		Temp10 appear as substrings of the relators and then makes the appropriate
		determination.
	*****************************************************************************************/
		
	register int 	i;
	
	unsigned int 	t;
	
	int 			Found_Temp9,
					Found_Temp10;
	
	Found_Temp9 = FALSE;
	for(i = 1; i <= NumRelators; i++)
		if((length <= LR[i] + 2) && Sub_Str(*Temp9,*Relators[i]))
			{
			Found_Temp9 = TRUE;
			break;
			}		
	if(Found_Temp9 == FALSE)
		{
		HLock((char **) Temp9);
		Inverse(*Temp9);
		HUnlock((char **) Temp9);
		for(i = 1; i <= NumRelators; i++)
			if((length <= LR[i] + 2) && Sub_Str(*Temp9,*Relators[i]))
				{
				Found_Temp9 = TRUE;
				break;
				}
		}	
	if(Found_Temp9 == FALSE)
		{
		i = 2*j;
		t = OSA[i];
		OSA[i] = OSB[i];
		OSB[i] = t;
		t = OSA[i+1];
		OSA[i+1] = OSB[i+1];
		OSB[i+1] = t;
		return;
		}
		
	Found_Temp10 = FALSE;
	for(i = 1; i <= NumRelators; i++)
		if((length <= LR[i] + 2) && Sub_Str(*Temp10,*Relators[i]))
			{
			Found_Temp10 = TRUE;
			break;
			}		
	if(Found_Temp10 == FALSE)
		{
		HLock((char **) Temp10);
		Inverse(*Temp10);
		HUnlock((char **) Temp10);
		for(i = 1; i <= NumRelators; i++)
			if((length <= LR[i] + 2) && Sub_Str(*Temp10,*Relators[i]))
				{
				Found_Temp10 = TRUE;
				break;
				}
		}
	if(Found_Temp10 == FALSE)
		{
		i = 2*j;
		t = OSA[i];
		OSA[i] = OSB[i];
		OSB[i] = t;
		t = OSA[i+1];
		OSA[i+1] = OSB[i+1];
		OSB[i+1] = t;
		return;
		}
	return;								
}

Sub_Str(S1,S2)
unsigned char 	*S1,
				*S2;
{	
	/******************************************************************************************
		Given two strings S1 and S2, this routine returns TRUE if string S1 is a substring
		of a proper power of a cyclic conjugate of string S2 and otherwise returns FALSE.
	******************************************************************************************/	
	
	register unsigned char 	*p,
							*q,
							*r,
							x;
							
	x = *S1;
	for(r = S2; *r; r++) if(*r == x)
		{
		p = S1;
		q = r;
		while(*p)
			{
			if(*q == EOS) q = S2;
			if(*p != *q) break;
			p++;
			q++;	
			}
		if(*p == EOS) return(TRUE);	
		}
	return(FALSE);
}


Valence_Two_Annulus()
{
	/******************************************************************************************
		Valence_Two_Annulus() is called by Valence_Two() when Valence_Two() has determined that
		an annulus exists. This routine is supposed to determine which vertices of valence two
		are swallowed by the annulus and to also determine what curve on the Heegaard surface
		the centerline of the annulus follows.
	******************************************************************************************/
		
	register unsigned char 	*ptr1,
							*ptr2,
							*q,
							*r;
							
	register unsigned int 	d,
							e,
							t,
							v,
							w;
	
	int						Error;
							
	unsigned int 			i,
							j;						
							
	long 					HS,
							length1,
							length2;
	
	/******************************************************************************************
		First find a vertex of valence two to start things off. Then proceed from that vertex
		and trace out the path of the annulus. Any vertex of valence two that we encounter is
		swallowed by the annulus. The other vertices that we meet are part of the path of the
		centerline of the annulus. The routine traces out the path of the annulus twice. The
		first run is done only to determine how long the path of the centerline is. We are
		doing this because otherwise it seems we would have to allocate a handle or pointer
		large enough to hold a string of length equal to the sum of the lengths of all the
		relators in order to guarantee that it would hold the string that represents the 
		centerline of the annulus. By making a dry run first we can avoid having to ask for 
		such a large chunk of memory.
	******************************************************************************************/
	
	Error = FALSE;
	for(i = 0; i < Vertices && !GV2[i >> 1]; i++) ;
	w = i;
	v = FV[w];
	length1 = length2 = 0L;
	while(1)
		{
		if(GV2[v >> 1])
			{
			length1 ++;
			if(v == i) break;
			t = v;
			v = CO[v][w];
			w = t;
			}
		else
			{
			if(GV2[w >> 1])
				{
				e = 0;
				t = FV[v];
				while(t != w)
					{
					e += A[v][t];
					t = CO[v][t];
					}
				}	
			if(v & 1)
				w = v - 1;
			else
				w = v + 1;
			e = OSA[v] - e;
			length2 ++;
			if(e >= V[v]) e -= V[v];
			v = FV[w];
			d = A[w][v];
			while(d <= e)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			e = B[w][v] - e;
			}
		}
		
	/*****************************************************************************************
		If length1 >= Vertices - 2, then all but a single inverse pair of vertices are
		swallowed by this annulus. Since this is a trivial situation, we quit and return.
	*****************************************************************************************/
		
	if(length1 >= (long)(Vertices - 2)) return(NO_ERROR);
	
	if(TestRealizability3) return(TRUE);
				
	ptr1 = (unsigned char *) NewPtr(length1 + 1);
	if((q = ptr1) == NULL) return(TOO_LONG);	
	ptr2 = (unsigned char *) NewPtr(length2 + 1);
	if((r = ptr2) == NULL)
		{
		DisposePtr((char *) ptr1);
		return(TOO_LONG);
		}
	w = i;
	v = FV[w];
	while(1)
		{
		if(GV2[v >> 1])
			{
			if(v & 1)
				*q++ = (v >> 1) + 97;
			else
				*q++ = (v >> 1) + 65;
			if(v == i) break;
			t = v;
			v = CO[v][w];
			w = t;
			}
		else
			{
			if(GV2[w >> 1])
				{
				e = 0;
				t = FV[v];
				while(t != w)	
					{
					e += A[v][t];
					t = CO[v][t];
					}
				}			
			if(v & 1)
				{
				*r++ = (v >> 1) + 97;
				w = v - 1;
				}
			else
				{
				*r++ = (v >> 1) + 65;
				w = v + 1;
				}
			e = OSA[v] - e;
			if(e >= V[v]) e -= V[v];
			v = FV[w];
			d = A[w][v];
			while(d <= e)
				{
				v = CO[w][v];
				d += A[w][v];
				}
			e = B[w][v] - e;
			}
		}
	*r = EOS;
	*q = EOS;
	
	/******************************************************************************************
		Since we are going to call Canonical_Rewrite() on the relators in Relators[], and that
		routine modifies the relators, we will get an error when the program returns from
		Valence_Two(), and the program finds that the relators realized by the Heegaard diagram
		are not the relators in Relators[]. To avoid this, we need to save a copy of the
		relators and restore things before we return.
	******************************************************************************************/
		
	for(i = 1; i <= NumRelators; i++)
		{
		ReallocateHandle((char **) DelRelators[i],GetHandleSize((char **) Relators[i]));
		if((r = *DelRelators[i]) == NULL)
			{
			Error = TRUE;
			goto END;
			}
		q = *Relators[i];
		while(*r++ = *q++) ;
		}
	ReallocateHandle((char **) Relators[0],length1 + length2 + 2);
	if((r = *Relators[0]) == NULL)
		{
		Error = TRUE;
		goto END;
		}
	q = ptr1;
	while(*r++ = *q++) ;
	r--;
	*r++ = '@';
	q = ptr2;
	while(*r++ = *q++) ;
	Canonical_Rewrite(Relators,TRUE,FALSE);
	
	/******************************************************************************************
			Run through our list of presentations to see if we already have this one.
	******************************************************************************************/
		
	for(i = 0; i < NumFilled; i++)
 	if(SURL[i] == Length  
		&& NG[i] == NumGenerators
		&& NR[i] == NumRelators)
		{
	 	for(j = 1; j <= NumRelators; j++)
	 		if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;
	 	if(j > NumRelators && Compare_Pres(i))
	 		{
	 		if(UDV[i] <= DONE)
	 			{
				SUR[i][0] = (unsigned char **) NewHandle(length1 + length2 + 2);
				if((r = *SUR[i][0]) == NULL)
					{
					Error = TRUE;
					goto END;
					}	
				q = *Relators[0];
				while(*r++ = *q++) ;
				BytesUsed += GetHandleSize((char **) SUR[i][0]);
				UDV[i] = V2_ANNULUS_EXISTS;
				if(!DrawingDiagrams)
				printf("\n					There is a 'valence-two' annulus in diagram %u.",
				i + 1);
				}
	 		break;
	 		}
	 	}						
	if(i == NumFilled && !DrawingDiagrams)
		{
		if((NumFilled >= MAX_SAVED_PRES - 3) || Save_Pres(ReadPres,0,Length,1,80,0,0,0))
			{
			Error = TRUE;
			goto END;
			}		
		SUR[NumFilled - 1][0] = (unsigned char **) NewHandle(length1 + length2 + 2);
		if((r = *SUR[NumFilled - 1][0]) == NULL)
			{
			Error = TRUE;
			goto END;
			}
		q = *Relators[0];
		while(*r++ = *q++) ;
		BytesUsed += GetHandleSize((char **) SUR[NumFilled - 1][0]);
		BDY[NumFilled - 1] = BDY[ReadPres];
		UDV[NumFilled - 1] = V2_ANNULUS_EXISTS;
		printf("\n					There is a 'valence-two' annulus in diagram %u.",
		NumFilled);
		}
		
	/******************************************************************************************
							Restore the relators and the array LR[].
	******************************************************************************************/
			
	for(i = 1; i <= NumRelators; i++)
		{
		HS = GetHandleSize((char **) DelRelators[i]);
		ReallocateHandle((char **) Relators[i],HS);
		if((r = *Relators[i]) == NULL)
			{
			Error = TRUE;
			goto END;
			}
		LR[i] = HS - 1;		
		q = *DelRelators[i];
		while(*r++ = *q++) ;
		}
END:			
	ReallocateHandle((char **) Relators[0],4L);	
	DisposePtr((char *) ptr1);
	DisposePtr((char *) ptr2);
	if(Error) return(TOO_LONG);
	return(TRUE);			
}

Missing_Gen()
{
	/*****************************************************************************************
		Missing_Gen() is called when Find_Flow_A() has returned TRUE. This happens when there
		is a generator which no longer appears in the relators. In this case the manifold 
		has some handles: (either S1 X D2's or S1 X S2's). This routine is supposed to "split
		off" the presentations of these handles. It creates two new presentations: one 
		consists of empty relators corresponding to the handles, and the other consists
		of the remaining rewritten nonempty relators. These two new presentations are then
		treated as two new "summands" of the manifold.
	*****************************************************************************************/
		
	unsigned char 	*p,
					*q;
						
	int 			i,
					j,
					MissingGens,
					MissingRels,
					SaveCS,
					SNumGenerators,
					SNumRelators;
	
	unsigned int	SaveUDV;								
			
	SNumGenerators = NumGenerators;
	SNumRelators = NumRelators;
	Rewrite_Input();
	MissingGens = SNumGenerators - NumGenerators;
	MissingRels = SNumRelators - NumRelators;
	
	if(MissingGens)
		{
		for(i = 1, Length = 0L; i <= NumRelators; i++)
			Length += GetHandleSize((char **) Relators[i]);
		Length -= NumRelators;			
				
 		/**************************************************************************************
 			If the program already has as many summands as it can handle, flag any other
 			presentations corresponding to this summand so that we will quit processing them.
 		**************************************************************************************/
		
		if(Length && TotalComp > MAXNUMCOMPONENTS - 3)
			{
			Mark_As_Found_Elsewhere(CurrentComp);
			SysBeep(5);
			printf("\n\nStopping because the program cannot deal with any more summands. Sorry!");	
			return(TOO_MANY_COMPONENTS);
			}
		
		if(Length == 0L)
			{
			NumGenerators = SNumGenerators;
			NumRelators = SNumRelators;
			if(Micro_Print)
				{
				printf("\n\nThe presentation has reduced to an empty presentation.");
				if(Micro_Print_F)
					fprintf(myout,"\n\nThe presentation has reduced to an empty presentation.");
				}
			if((NumFilled >= MAX_SAVED_PRES - 3)
				|| Save_Pres(ReadPres,0,Length,1,30,0,0,0)) return(TOO_LONG);
			BDY[NumFilled - 1] = BDY[ReadPres];
			if(BDY[ReadPres] == FALSE)
				{
				UDV[NumFilled - 1] = S1_X_S2;
				Mark_As_Found_Elsewhere(CurrentComp);
				CS[CurrentComp] = 2;
				if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
					{
					CBC[CurrentComp][0] = 1;
					CBC[CurrentComp][1] = BDRY_UNKNOWN;
					}								
				}
			else	
				UDV[NumFilled - 1] = DONE;
			return(NO_ERROR);	
			}
			
		Canonical_Rewrite(Relators,FALSE,FALSE);
	
		SaveCS = CS[CurrentComp];
		SaveUDV = UDV[ReadPres];
		switch(SaveUDV)
			{
			case SPLIT:
			case ANNULUS_EXISTS:
			case V2_ANNULUS_EXISTS:
				return(NO_ERROR);
			}
		if(!CS[CurrentComp]) CS[CurrentComp] = TRUE;

		if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);
										
		TotalComp 							++;
		UDV[ReadPres] 						= SPLIT;
		NCS[ReadPres] 						= 2;				
		Daughters[ReadPres] 				= NumFilled;
		BDY[NumFilled]						= BDY[ReadPres];
		ComponentNum[NumFilled] 			= TotalComp;	
		ER[NumFilled]  						= -2;
		FR[NumFilled]						= ReadPres;
		MLC[TotalComp][NumGenerators] 		= Length;
		NG[NumFilled] 						= NumGenerators;
		NR[NumFilled] 						= NumRelators;
		PRIM[NumFilled] 					= 30;
		SURL[NumFilled] 					= Length;
		UDV[NumFilled] 						= 0;
		TP[NumFilled] 						= NumRelators;
		OnStack 						    += 2*NumGenerators;
		
		for(i = 1; i <= NumRelators; i++)
			{
			SUR[NumFilled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));			
			if((q = *SUR[NumFilled][i]) == NULL)
				{
				for(j = 1; j < i; j++) DisposeHandle((char **) SUR[NumFilled][j]);
				MLC[TotalComp][NumGenerators] 	= BIG_NUMBER;
				UDV[ReadPres] 					= SaveUDV;
				CS[CurrentComp] 				= SaveCS;
				OnStack 						-= 2*NumGenerators;
				TotalComp 						--;
				return(TOO_LONG);
				}
			p = *Relators[i];
			while(*q++ = *p++) ;									
			}
			
		BytesUsed += Length;
									
		NumFilled ++;
		SaveMinima = TRUE;
		
		if(Micro_Print)
			{
			if(MissingGens == 1)
				{
				printf("\n\nA generator no longer appears in the relators.");
				printf("\nThe manifold therefore has a 'handle' of some kind.");
				if(Micro_Print_F)
					{
					fprintf(myout,"\n\nA generator no longer appears in the relators.");
					fprintf(myout,"\nThe manifold therefore has a 'handle' of some kind.");
					}
				}
			else
				{
				printf("\n\n%d generators no longer appear in the relators.",MissingGens);
				printf("\nThe manifold therefore has %d 'handles' of some kind.",MissingGens);
				if(Micro_Print_F)
					{
					fprintf(myout,"\n\n%d generators no longer appear in the relators.",MissingGens);
					fprintf(myout,"\nThe manifold therefore has a 'handle' of some kind.");
					}
				}
			printf("\nThe presentation is currently:\n");
			Print_Relators(Relators,NumRelators,stdout);
			printf("\n\nFound Presentation %u and Presentation %u.\n",
				NumFilled,NumFilled + 1);	
			if(Micro_Print_F)
				{
				fprintf(myout,"\nThe presentation is currently:\n");			
				Print_Relators(Relators,NumRelators,myout);	
				fprintf(myout,"\n\nFound Presentation %u and Presentation %u.\n",
					NumFilled,NumFilled + 1);
				}
			}	
		fprintf(stdout,"\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
		fprintf(stdout,"Gen%3d  Rel%3d  Length%6lu  From%6d  MG",
			NumGenerators,NumRelators,Length,ReadPres + 1);
		
		TotalComp 					++;
		BDY[NumFilled] 				= 3;
		ComponentNum[NumFilled] 	= TotalComp;			
		ER[NumFilled]  				= 0;
		FR[NumFilled]				= ReadPres;
		MLC[TotalComp][MissingGens] = 0L;	
		NG[NumFilled] 				= MissingGens;
		NR[NumFilled] 				= MissingRels;
		PRIM[NumFilled]             = 30;
		SURL[NumFilled] 			= 0L;
		TP[NumFilled]				= 0;
		
		for(i = 1; i <= MissingRels; i++)
			{
			SUR[NumFilled][i] = (unsigned char **) NewHandle(1L);			
			q = *SUR[NumFilled][i];
			*q = EOS;									
			}
		
		if(BDY[ReadPres] == FALSE)
			{
			BDY[NumFilled] = FALSE;
			UDV[NumFilled] = S1_X_S2;
			CS[TotalComp] = 2;
			CBC[TotalComp - 1][0] = 1;
			CBC[TotalComp - 1][1] = BDRY_UNKNOWN;
			CBC[TotalComp][0] = 1;
			CBC[TotalComp][1] = BDRY_UNKNOWN;
			}
		else
			{		
			UDV[NumFilled] = S1_X_X2;
			CS[TotalComp] = 3;
			}	
		NumFilled ++;
		SaveMinima = TRUE;
						
		fprintf(stdout,"\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
		fprintf(stdout,"Gen%3d  Rel%3d  Length%6lu  From%6d  MG",
			MissingGens,MissingRels,0L,ReadPres + 1);
		
		NumFilled -= 2;
		i = On_File();
		if(Dup_On_File < INFINITE)
			{
			UDV[NumFilled] = DUPLICATE;
			Daughters[NumFilled] = Dup_On_File;
			}
		else if(i < NumFilled)
			{
			UDV[NumFilled] = DUPLICATE;
			Daughters[NumFilled] = i;
			}
		NumFilled += 2;
		}
		
	return(NO_ERROR);																			
}	

int Empty_Relator_D(void)
{
	/******************************************************************************************
		Empty_Relator_D() is called when the dual presentation of the presentation of a closed
		manifold contains an empty relator. The routine "splits" off a presentation
		corresponding to the set of empty relators, decrements the number of relators and
		calls Save_Pres() to save the new presentation.
	******************************************************************************************/
			
	register int 	i,
					j;
	
	unsigned char 	**Temp,
					*q;
	
	int 			NumEmptyRelators,
					SaveCS;
	
	unsigned int	SaveUDV;

	unsigned long	SaveMyLength;
	
	/******************************************************************************************
		First, take care of the case where the presentation has reduced to the empty
									presentation.
	******************************************************************************************/
	
	if(Length == 0L)
		{
		Missing_Gen();
		PRIM[NumFilled - 1] = 40;
		return(NO_ERROR);
		}
										
	/******************************************************************************************
		Next, check whether the presentation which "split" is on file. If it is, we don't
								want to save another copy.
	******************************************************************************************/	

	for(i = 1; i <= NumRelators; i++)
		{ 
		Temp 			= Relators[i];
		Relators[i] 	= DualRelators[i];
		DualRelators[i] = Temp;	
		}
	
	SaveMyLength = Length;
	Length 		 = SLength;
	
	switch(Splitting_Pres_On_File(40,2))
		{
		case TOO_LONG:
			return(TOO_LONG);
		case TOO_MANY_COMPONENTS:
			return(TOO_MANY_COMPONENTS);
		case TRUE:
			return(TRUE);	
		case NO_ERROR:
			break;
		}
	
	for(i = 1; i <= NumRelators; i++)
		{ 
		Temp 			= Relators[i];
		Relators[i] 	= DualRelators[i];
		DualRelators[i] = Temp;	
		}
	
	Length					= SaveMyLength;
	SaveUDV 				= UDV[ReadPres];
	SaveCS					= CS[CurrentComp];	
	switch(SaveUDV)
		{
		case SPLIT:
		case ANNULUS_EXISTS:
		case V2_ANNULUS_EXISTS:
			return(NO_ERROR);
		}
		
	NumEmptyRelators = 0;
	for(i = NumRelators; i > 0; i--) if(GetHandleSize((char **) Relators[i]) == 1L)
		{
		NumEmptyRelators ++;
		if(i < NumRelators)
			{
			Temp = Relators[NumRelators];
			Relators[NumRelators] = Relators[i];
			Relators[i] = Temp;
			}
		NumRelators --;
		}
	
	if(Micro_Print)
		{
		printf("\n\nThe current dual presentation contains %d empty relator(s):\n",
			NumEmptyRelators);
		if(Micro_Print_F)	
			fprintf(myout,"\n\nThe current dual presentation contains %d empty relator(s):\n",
				NumEmptyRelators);
		NumRelators += NumEmptyRelators;	
		Print_Relators(Relators,NumRelators,stdout);
		if(Micro_Print_F)
			Print_Relators(Relators,NumRelators,myout);
		NumRelators -= NumEmptyRelators;
		}
	
	switch(Find_Flow_A(NORMAL,FALSE))
		{
		case TOO_LONG:
			UDV[ReadPres] 	= SaveUDV;
			CS[CurrentComp] = SaveCS;
			ReadPres = SReadPres;
			return(TOO_LONG);
		case 1:
			NumRelators += NumEmptyRelators;
			switch(Missing_Gen())
				{
				case TOO_LONG:
					UDV[ReadPres] 	= SaveUDV;
					CS[CurrentComp] = SaveCS;
					ReadPres = SReadPres;
					return(TOO_LONG);
				case TOO_MANY_COMPONENTS:
					UDV[ReadPres] 	= SaveUDV;
					CS[CurrentComp] = SaveCS;
					ReadPres = SReadPres;
					return(TOO_MANY_COMPONENTS);
				case NO_ERROR:
					PRIM[NumFilled - 2] = 40;
					PRIM[NumFilled - 1] = 40;
					ReadPres = SReadPres;
					return(FALSE);
				}	
		case NO_ERROR:
			ReadPres = SReadPres;
			return(TRUE);
		}
}

int Empty_Relator_BS(void)
{
	/******************************************************************************************
		Empty_Relator_BS() is called if a bandsum has created an empty relator. The routine
		"simplifies" the presentation corresponding to the set of nonempty relators, decrements
		the number of relators and calls Save_Pres() to save the new presentation.
	******************************************************************************************/


	Delete_Dups();
	NumRelators -= NumEmptyRels;
			
	switch(Find_Flow_A(NORMAL,FALSE))
		{
		case TOO_LONG:
			return(TOO_LONG);
		case 1:
			/**********************************************************************************
				This case should not occur. Assuming that the incoming presentation has
				minimal length, full rank and a connected Whitehead graph, then the only
				possible type of empty relator that can be generated by a bandsum, is one in
				which the empty relator contracts on the Heegaard surface. For suppose that
				R3 is a bandsum of R1 and R2, and R3 is empty. Then R3 must be trivial. Hence
				R1 and R2 are parallel. But then any cutting disk for the handlebody H, which
				is disjoint from the curves in P - R2 + R3, is also disjoint from P.
			**********************************************************************************/
			
			return(TRUE);	
		default:
			break;
		}
			
	if(On_File() == NumFilled)
		{
		if(Micro_Print)
			{
			printf("\n\nDeleted an empty relator created by a bandsum. The presentation is currently:\n");
			Print_Relators(Relators,NumRelators,stdout);
			if(Micro_Print_F)
				{
				fprintf(myout,"\n\nDeleted an empty relator created by a bandsum. The presentation is currently:\n");	
				Print_Relators(Relators,NumRelators,myout);
				}
			}
		if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);	
		if(Dup_On_File < INFINITE)
			{
 			if(Save_Pres(ReadPres,Dup_On_File,Length,1,13,1,0,0)) return(1);								
 			Mark_As_Duplicate(Dup_On_File);
 			}
 		else
 			{
			if(Save_Pres(ReadPres,0,Length,1,13,1,0,0)) return(1);		 						
			BDY[NumFilled - 1] = BDY[ReadPres];
			UDV[NumFilled - 1] = 0;
			}	
		}
		
	return(NO_ERROR);
}
