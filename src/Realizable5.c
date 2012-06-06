#include "Heegaard.h"
#include "Heegaard_Dec.h"

int New_Relator(int F1)
{
	/******************************************************************************************
		New_Relator() looks for a pair of relators that it can bandsum together to form a 
		new relator. It makes use of data about the Heegaard diagram found by the routines in
		Diagram_Main().
			The routine only looks for and performs certain restricted kinds of bandsums,
		namely those which correspond to bandsumming two relators together along a band which
		corresponds to an edge of the dual Heegaard diagram. These edges correspond to "corners"
		of the major faces of the original Heegaard diagram. Note that forming bandsums in 
		this restricted fashion is the only way to avoid forming "waves" in the diagram.
	******************************************************************************************/
		
	register unsigned char 	*p,
							*q;
							
	register unsigned int 	c,
							d,
							e,
							v,
							w;
							
	unsigned char 			s,
							t,
							x,
							y,
							sx,
							sy;
													
	int 					Depth,
							ss,
							uu;
							
	unsigned int 			edge,
							edgeLE,
							edgeLS,
							edgeRE,
							edgeRS,
							EL[6*VERTICES],
							ee,
							eLE,
							eLS,
							eRE,
							eRS,
							vertex,
							vertexLS,
							vertexRS,
							vLS,
							vRS,
							vrt;
	
	long 					Diff,
							LOR[MAXNUMRELATORS + 1],
							length,
							slength;
	
	unsigned long			max;
	
	for(d = 1,max = 0L; d <= NumRelators; d++) if(LR[d] > max) max = LR[d];
	if(SRError == 3) for(d = 1; d <= NumRelators; d++)
		LOR[d] = GetHandleSize((char **) OutRelators[d]);
	Minimum = BIG_NUMBER;
	
	/******************************************************************************************
		Set the parameter Depth. The number of edges along which the routine will form
		bandsums is determined by this parameter.
	******************************************************************************************/	
	
	if(!GoingUp)
		Depth = 2*NumEdges;
	else
		Depth = 3*(2*NumGenerators - UDV[ReadPres]);
			
	for(d = 1; d <= 2*NumEdges; d++) EL[d] = d;
	uu = 0;
	for(ss = 2*NumEdges; ss > 0; ss --)
		{
		/**************************************************************************************
			Choose a random edge of the dual diagram along which to perform a bandsum.
			Remove the chosen edge from the list of available edges so that we won't try it
			again.
		**************************************************************************************/	
		
		c = abs(rand()) % ss;
		c++;
		v = EL[ss];
		EL[ss] = EL[c];
		EL[c] = v;
		
		/**************************************************************************************
			Determine which edges of the original diagram the chosen edge, of the dual
			diagram, joins.
		**************************************************************************************/	
		
		for(v = c = 0; c + VWG[v] < EL[ss]; v++) c += VWG[v];
		w = EL[ss] - c;
		for(c = ee = 0,d = FV[v]; c < w; c++)
			{
			ee += A[v][d];
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
			
		/**************************************************************************************
				Check whether this edge we are trying to bandsum along has its ends on
				distinct relators. If it does, we can go ahead and form a bandsum.
		**************************************************************************************/		
		
		p = q = *DualRelators[(v >> 1) + 1] + e;
		q++;
		if(! *q) q = *DualRelators[(v >> 1) + 1];
		uu ++;
		if((abs(*p - *q) != 32) && (*p != *q)) 	
			{
			x 		= *p;
			y 		= *q;
			length	= 0L;
			vertex  = v;
			edgeRE 	= ee;
			edgeLE 	= edgeRE - 1;
			edgeRE  = edgeRE % V[v];	
			edge 	= edgeRE;
			e 		= edge;
			
			/*******************************************************************************
				Two distinct relators have come from distinct vertices and are adjacent to 
				each other at vertex v. Follow these two relators along until they diverge.
				The region where they are parallel is the "cancellation region".
			*******************************************************************************/	
			
			do
				{
				if(v & 1)
					w = v - 1;
				else
					w = v + 1;
				e = OSA[v] - e;
				if(e >= V[v]) e -= V[v];
				length ++;
				if(length > max)
					{
					/***********************************************************************
						This is an error that should not occur. But if the impossible
						happens, we fail gracefully by returning and pretending we could not
						find the bandsum.
					***********************************************************************/
					
					printf("\nlength > max in a bandsum!");
					Minimum = BIG_NUMBER;
					return;
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
					
					edgeRS = B[w][v] - e;
					vertexRS = v;
					v = CO[w][v];
					d = d % V[w];
					edgeLS = B[w][v] - d;
					vertexLS = v;
					
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
						uu ++;
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
							uu ++;
							}
						}		
					break;
					}
				e = B[w][v] - e;
				}	
			while(v != vertex || e != edge);
			
			/*******************************************************************************
				Increment the number of bandsums examined and use the information about the
				length of the cancellation region to determine how the total length of the
				relators would be changed by performing this bandsum and replacing one of
				the relators with the result.
			*******************************************************************************/							
			
			Band_Sums ++;
			for(d = 0; d < 2; d++)
				{
				if(d == 0)
					{
					s = x;
					t = y;
					}
				else
					{
					s = y;
					t = x;
					}
					
				/***************************************************************************
						Compute Diff, which is the amount by which the length of the
						presentation would be changed by performing this bandsum and
						replacement. If Diff < Minimum and !GoingUp or Diff > 0, replace
						the current saved set of parameters with the set of parameters we
						have just found.
				***************************************************************************/			
				
				if(s < 'a')
					Diff = LOR[s - 64] - (length << 1) - 1;
				else	
					Diff = LOR[s - 96] - (length << 1) - 1;
				if(Diff < Minimum)
					{
					if(!GoingUp || Diff > 0L)
						{
						Minimum = Diff;
						sx 		= x;
						sy 		= y;
						eRS 	= edgeRS;
						eLS 	= edgeLS;
						eRE 	= edgeRE;
						eLE 	= edgeLE;
						vRS 	= vertexRS;
						vLS 	= vertexLS;
						vrt 	= vertex;
						slength = length;
						if(t < 'a')
							Word1 = t - 64;
						else	
							Word1 = t - 96;
						if(s < 'a')
							Word2 = s - 64;
						else	
							Word2 = s - 96;	
						}
					}
				else
				if(Diff == Minimum)
					{
					if((!GoingUp || Diff > 0L) && abs(rand()) % 2)
						{
						sx 		= x;
						sy 		= y;
						eRS 	= edgeRS;
						eLS 	= edgeLS;
						eRE 	= edgeRE;
						eLE 	= edgeLE;
						vRS 	= vertexRS;
						vLS 	= vertexLS;
						vrt 	= vertex;
						slength = length;
						if(t < 'a')
							Word1 = t - 64;
						else	
							Word1 = t - 96;
						if(s < 'a')
							Word2 = s - 64;
						else	
							Word2 = s - 96;	
						}
					}
				}		
			}		
		if(uu >= Depth || (!GoingUp && Minimum < 0L && (uu << 2) >= Depth))
			{
			/**********************************************************************************
				If uu >= Depth, we have looked at bandsums along at least Depth number of 
				randomly chosen edges. If GoingUp and Minimum < BIG_NUMBER, call Long_Mult() to
				perform our chosen bandsum. Otherwise return.
					If !GoingUp, there is a bandsum which reduces the length of the
				presentation, and we have examined at least one-fourth of the possible
				bandsums, call Long_Mult(). 
			**********************************************************************************/		
			
			if(Minimum == BIG_NUMBER || !F1) return(NO_ERROR);	
			Long_Mult(sx,sy,eRS,eLS,eRE,eLE,vRS,vLS,vrt,slength,max);		
			return(NO_ERROR);
			}
		}
		
	if(Minimum < BIG_NUMBER && F1)
		Long_Mult(sx,sy,eRS,eLS,eRE,eLE,vRS,vLS,vrt,slength,max);
		return(NO_ERROR);			
}

int Long_Mult(unsigned char x,unsigned char y,unsigned int edgeRS,unsigned int edgeLS,unsigned int edgeRE,
				unsigned int edgeLE,unsigned int vertexRS,unsigned int vertexLS,unsigned int vertex,long slength,
				unsigned long max)
{
	/******************************************************************************************
				This routine is called by New_Relator() to perform a bandsum.
	******************************************************************************************/
		
	register unsigned char 	*p,
							*r;
							
	register unsigned int 	d,
							e,
							v,
							w;
							
	unsigned long 			HS;
	
	HS = max + 2;
	if(HS > MAXLENGTH)
		{
		Minimum = BIG_NUMBER;
		return(NO_ERROR);
		}
	ReallocateHandle((char **) Temp1,HS);
	if((r = *Temp1) == NULL)
		{
		Minimum = BIG_NUMBER;
		return(NO_ERROR);	
		}
	HS = 2*max + 2;
	if(HS > MAXLENGTH)
		{
		Minimum = BIG_NUMBER;
		ReallocateHandle((char **) Temp1,4L);
		return(NO_ERROR);
		}		
	ReallocateHandle((char **) Temp2,HS);
	if(*Temp2 == NULL)
		{
		Minimum = BIG_NUMBER;
		ReallocateHandle((char **) Temp1,4L);
		return(NO_ERROR);	
		}		
	e = edgeLS;
	v = vertexLS;
	
	/******************************************************************************************
		Start at edge e of vertex v and record the word formed by this path, in the string
		*Temp1, until we arrive at edge edgeLE of vertex vertex.
	******************************************************************************************/	
	
	while(v != vertex || e != edgeLE)
		{
		if(v & 1)
			{
			*r = (v >> 1) + 'a';
			w = v - 1;
			}
		else
			{
			*r = (v >> 1) + 'A';
			w = v + 1;
			}
		e = OSA[v] - e;
		if(e >= V[v]) e -= V[v];
		r++;
		v = FV[w];
		d = A[w][v];
		while(d <= e)
			{
			v = CO[w][v];
			d += A[w][v];
			}
		e = B[w][v] - e;
		}
	*r = EOS;
	if((r = *Temp2) == NULL)
		{
		Minimum = BIG_NUMBER;
		ReallocateHandle((char **) Temp1,4L);
		return(NO_ERROR);	
		}
	v = vertexRS;
	e = edgeRS;
	
	/******************************************************************************************
		Next, start at edge e of vertex v and record the word formed by this path, in the
		string *Temp2, until we arrive at edge edgeRE of vertex vertex.
	******************************************************************************************/	
	
	while(v != vertex || e != edgeRE)
		{
		if(v & 1)
			{
			*r = (v >> 1) + 'a';
			w = v - 1;
			}
		else
			{
			*r = (v >> 1) + 'A';
			w = v + 1;
			}
		e = OSA[v] - e;
		if(e >= V[v]) e -= V[v];
		r++;
		v = FV[w];
		d = A[w][v];
		while(d <= e)
			{
			v = CO[w][v];
			d += A[w][v];
			}
		e = B[w][v] - e;
		}
	*r = EOS;
	
	/******************************************************************************************
					Invert *Temp2 and concatenate *Temp1 and *Temp2.
	******************************************************************************************/
			
	Inverse_Nr(*Temp2);  		
	p = *Temp1;
	while(*r++ = *p++) ;
	
	/******************************************************************************************
		If vertex = vertexLS and edgeLS = edgeLE or vertex = vertexRS and edgeRS = edgeRE,
		then one of the two relators just banded together is a subword of the other. 
		In this case, after deleting the shorter relator from the longer relator, the new
		relator may not be freely reduced. So we have to call Freely_Reduce_Nr() to freely
		reduce this relator.
	******************************************************************************************/	
	
	if((vertex == vertexLS) && (edgeLS == edgeLE))
		{
		Freely_Reduce_Nr();
		if(Minimum == BIG_NUMBER) return(NO_ERROR);
		HS = GetHandleSize((char **) Temp2);	
		}	
	else
		{	 
		if((vertex == vertexRS) && (edgeRS == edgeRE))
			{
			Freely_Reduce_Nr();
			if(Minimum == BIG_NUMBER) return(NO_ERROR);
			HS = GetHandleSize((char **) Temp2);	
			}	
		else
			{	
			if(x > 96) x -= 32;
			if(y > 96) y -= 32;	   
			HS = GetHandleSize((char **) OutRelators[x - 64])
				   + GetHandleSize((char **) OutRelators[y - 64]) - (slength << 1) - 1;
			if(HS > MAXLENGTH)
				{
				Minimum = BIG_NUMBER;
				return(NO_ERROR);
				}		   	   	   	   
			SetHandleSize((char **) Temp2,HS);	
			}
		}
	if(HS == 1L) EmtyRel = TRUE;
	if(HS > MAXLENGTH) Minimum = BIG_NUMBER;			
	return(NO_ERROR);				 																														
}

Inverse_Nr(p)
register unsigned char *p;
{
	/******************************************************************************************
		This routine transforms the string pointed to by p into its inverse by first replacing
		each char in the string with its inverse and then writing the string backwards.
	******************************************************************************************/
	
	register unsigned char *q,
							x;
							
	q = p;
	while(*q) 
		{
		if(*q < 95) *q += 32;
		else *q -= 32;
		q++;
		}
	q--;
	while(p < q)
		{
		x = *q;
		*q = *p;
		*p = x;
		p++;
		q--;
		}					
}

Freely_Reduce_Nr()
{		
	/******************************************************************************************
					This routine freely reduces the new relator found in Temp2.
	******************************************************************************************/
	
	register unsigned char 	*p,
							*q;				
							
	register char 			x;
	
	ReallocateHandle((char **) Temp12,GetHandleSize((char **) Temp2) + 2);
	if((p = *Temp12) == NULL)
		{
		Minimum = BIG_NUMBER;
		return;
		}	
	q = *Temp2;
	*p = '@';
	while(*q)
		{
		x = *p - *q;
		if(x == 32 || x == -32)
			p--;
		else 
			{
			p++;
			*p = *q;
			}	
		q++;	
		}
	q = *Temp12;
	q++;
	if(p > q) while(*p - *q == 32 || *p - *q == -32)
		{
		p--;
		q++;
		}	
	p++;
	*p = EOS;		
	SetHandleSize((char **) Temp2,p + 1 - q);
	if((p = *Temp2) == NULL)
		{
		Minimum = BIG_NUMBER;
		return;
		}
	while(*p++ = *q++) ;						
}		

Canonical_Rewrite(MyRelators,Annulus_Exists,ReportSymmetries)
unsigned char 	***MyRelators;

int 			Annulus_Exists,
				ReportSymmetries;
{
	/******************************************************************************************
		The program needs some means of determining whether two presentations correspond to 
		isomorphic Heegaard diagrams.
			Canonical_Rewrite() rewrites presentations in a canonical form so that two
		presentations have isomorphic Heegaard diagrams only if the two rewritten
		presentations are identical. (The converse holds provided there are no pairs of
		separating vertices or annuli.) The idea is to put a presentation into a canonical
		form under the action of the group of substitutions generated by permutations of the
		generators and replacements of generators by their inverses together with permutations
		of the relators and replacements of relators by cyclic conjugates or cyclic conjugates
		of their inverses.
		 	To get a canonical form, we impose the following ordering on presentations.
		Given presentations P = { R1, , , ,Rn} and P' = { R1', , ,Rn'} with the same
		number of relators, say that P precedes P' if there exists an integer i such that the
		length of Ri is greater than the length of Ri'. While if Ri and Ri' have the
		same length, say that Ri precedes Ri' if Ri precedes Ri' in the lexicographic
		ordering induced by the ordering A,B,C,..Z,a,b,c,...z of the generators and their
		inverses.
			Canonical_Rewrite() then works essentially as follows. First, the routine sorts
		the relators according to length. Then the routine maintains a stack of current optimal 
		substitutions which minimize the presentation up to relators of length say j. If better
		substitutions are found, the stack is truncated, while if new substitutions are found
		which also minimize the presentation in the ordering, these substitutions are added to
		the stack. At termination, the stack holds a set of optimal substitutions. We may then
		rewrite the relators using any of the set of substitutions on the stack. 
			Finally, we finish with a call to Final_Rewrite() which sorts the relators and
		replaces rewritten relators with cyclic conjugates of themselves or of their inverses
		in order to put the relators in their final canonical form.
	******************************************************************************************/
					 
	register unsigned char  *p12,
							*p22,
							*s,
							x,
							y;
							
	register int 			i,
							j;
							
	register unsigned int	k;							
							
	unsigned char 			Changed[VERTICES/2],
							Gen1,
							Gen2,
							*p,
							*p11,
							*p21,
							*Str,
							*Str1,
							*Str2,
							*T1,
							*T2,
							**Temp;
	
	int 					CurrentR1,
							CurrentR2,
							Depth,
							F2,
							Invert,
							Match,
							NumChanges,
							NumDups,
							Mem_Flag,
							SaveR2,
							Stack_Ptr,
							UpdateStr1;
							
	unsigned int			Exp,
							Max_Exp;						
													
	long 					length;
	
	/******************************************************************************************
		Sort the relators of the presentation according to length. Set length equal to the
		length of the longest relator.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++) LR[i] = GetHandleSize((char **) MyRelators[i]) - 1;
	LR[0] = 10000000L;
	for(i = 2; i <= NumRelators; i++)
		{
		j = i;
		length = LR[i];
		if(LR[j-1] < length)
			{
			Temp = MyRelators[i];
			while(LR[j-1] < length)
				{
				LR[j] = LR[j-1];
				MyRelators[j] = MyRelators[j-1];
				j--;
				}
			LR[j] = length;
			MyRelators[j] = Temp;
			}
		}		

	length = LR[1];
	
	/******************************************************************************************
			Get some memory for the strings Str, Str1, Str2 and for the array T2[].
	******************************************************************************************/
		
	ReallocateHandle((char **) Temp12,length + 2);
	ReallocateHandle((char **) Temp13,length + 1);
	ReallocateHandle((char **) Temp14,length + 1);
	ReallocateHandle((char **) Temp15,125L);
	
	if((Str  = *Temp12) == NULL) return(TOO_LONG);
	if((Str1 = *Temp13) == NULL) return(TOO_LONG);
	if((Str2 = *Temp14) == NULL) return(TOO_LONG);
	if((T2   = *Temp15) == NULL) return(TOO_LONG);
	
	HLock((char **) Temp12);
	HLock((char **) Temp13);
	HLock((char **) Temp14);
	HLock((char **) Temp15);
	
	/******************************************************************************************
		Get a handle RWR[0] of size 125 bytes. Make T1 point to the first of these bytes
		and initialize T1[]. T1[0] gives the current number of generators for which
		substitutions have been determined. T1[1] -> T1[NumRelators] are flags
		used to keep track of which relators have been used. One relator is used each time
		depth is incremented. T1[A] -> T1[A + Numgenerators] and T1[a] ->
		T1[a + NumGenerators] are used to hold info about substitutions of the generators.
	******************************************************************************************/
		
	RWR[0] = (unsigned char **) NewHandle(125L);
	if((T1 = *RWR[0]) == NULL)
		{
		HUnlock((char **) Temp12);
		HUnlock((char **) Temp13);
		HUnlock((char **) Temp14);
		HUnlock((char **) Temp15);
		return(TOO_LONG);
		}
	for(i = 'A' + NumGenerators - 1; i >= 'A'; i--) T1[i] = EOS;
	for(i = 'a' + NumGenerators - 1; i >= 'a'; i--) T1[i] = EOS;
	for(i = NumRelators; i >= 0; i--) T1[i] = EOS;
	T2[0] = 1;
	
	/******************************************************************************************
		Set Mem_Flag FALSE, Mem_Flag will become TRUE only if the number of potential
		substitutions on the stack of potential substitutions exceeds the preset limit
		MAXNUMDUPS.
	******************************************************************************************/	
	
	Mem_Flag = FALSE;
	
	/******************************************************************************************
			Numdups keeps track of the number of handles that are on the stack.
	******************************************************************************************/
		
	for(Depth = 1,NumDups = 0,Max_Exp = 0; Depth <= NumRelators; Depth ++)
		{
		length = LR[Depth];
		if(length <= 0L) break;
		
		/**************************************************************************************
			Make Stack_Ptr point to the 'substitution' on the top of the stack. Lock the
			corresponding handle, RWR[Stack_Ptr], and make T1 point to its data. Copy the
			data about substitutions for the generators into T2[].
		**************************************************************************************/
			
		Stack_Ptr = NumDups;		
		HLock((char **) RWR[Stack_Ptr]);
		T1 = *RWR[Stack_Ptr];
		for(i = 'A'; i < 'A' + NumGenerators; i++) T2[i] = T1[i];
		for(i = 'a'; i < 'a' + NumGenerators; i++) T2[i] = T1[i];
		
		/**************************************************************************************
			Scan through the relators and find the first relator which is unused and whose
			length equals the current value of length. Set CurrentR2 equal to i if the first
			relator satisfying these requirements is relator i. Make a copy of this relator
			in the string Str2.
		**************************************************************************************/

		for(i = 1; i <= NumRelators; i++)
			if(!T1[i] && LR[i] == length) break;
		CurrentR2 = i;	
		p22 = *MyRelators[CurrentR2];	
			
		p12 = Str2;
		while(*p12++ = *p22++) ;
		
		/**************************************************************************************
			Initialize Gen2, NumChanges, UpdateStr1, Invert, Match, and the pointers p21, p22
			and s. Set *s equal to 125, which is greater than the ASCII decimal value of any
			alphabetic character.
		**************************************************************************************/
			
		Gen2 = T1[0];	/* 	Gen2 gives the current number of generators for which we have
							determined substitutions.					    				 */
		NumChanges = 0;	/*	Numchanges records the additional number of generators for which
							we have determined tentative substitutions.						 */					
		UpdateStr1 = TRUE;	/*	UpdateStr1 is set to True whenever the relator CurrentR2 is
								changed or inverted.										 */
		Invert = TRUE;	/*	Invert is a flag which tells the program whether it has examined
							the inverse of the relator CurrentR2.							 */
		p21 = p22 = Str2;
		s = Str;
		*s = 125;
		k = length;
		while(1)
			{
			/**********************************************************************************
				Set y eaual to the character pointed to by p22. If we know what we want to
				substitute for this character replace y by its image under the substitution.
				Otherwise.....
			**********************************************************************************/
				
			y = *p22;
			if(T2[y])
				y = T2[y];
			else
				{
				/******************************************************************************
					Find the first generator that has not been used as a substitute and
					use it as a substitute for y.
					Then update the entry of T2 corresponding to the inverse of 'y' to reflect
					this substitution. Replace y with its image under this substitution and
					increment the number of generators for which we have determined
					substitutions. Keep track of the changes made in the array Changed[],
					because we may want to undo them later.
				******************************************************************************/
					
				T2[y] = Gen2 + 'A';
				if(y > 'Z')
					{
					T2[y - 32] = Gen2 + 'a';
					Changed[NumChanges++] = y - 32;
					}
				else
					{
					T2[y + 32] = Gen2 + 'a';
					Changed[NumChanges++] = y;
					}
				y = T2[y];
				Gen2 ++;
				}
					
			x = *s;
			if(x == EOS) 
				{
				/******************************************************************************
					If x = EOS which is binary 0, we need to "translate" the next character in
					the string Str1. If no substitute for this character has been determined
					we choose one now.
				******************************************************************************/
					
				p12++;
				if(*p12 == EOS) p12 = Str1;
				x = *p12;
				if(T1[x])
					x = T1[x];
				else
					{
					T1[x] = Gen1 + 'A';
					if(x > 'Z')
						T1[x - 32] = Gen1 + 'a';
					else
						T1[x + 32] = Gen1 + 'a';
					x = T1[x];
					Gen1 ++;
					}
				*s++ = x;
				*s = EOS;		
				}
			else
				s++;
								
			if(x == y)
				{
				p22 ++;
				if(*p22 == EOS) p22 = Str2;
				if(--k == 0)
					{
					j = TRUE;
					if(length > 1L)
						{
						/**********************************************************************
							If length <= 1, this "match" comes from a permutation of trivial
							relators of length 1. The program chooses to treat such matches as
							trivial and ignores them.
						**********************************************************************/
							
						NumDups ++;
						if(NumDups < MAXNUMDUPS)
							{
							/******************************************************************
								A "match" has occured. The program has found more than one way
								to get the same optimal image string at this current depth.
								Save the relevant data in a new handle. Increment the number of
								handles on the stack i.e. NumDups. Note however, that Stack_Ptr
								is not changed.
							******************************************************************/
								
							RWR[NumDups] = (unsigned char **)NewHandle(125L);
							p = p12;
							p12 = *RWR[NumDups];
							if(p12 == NULL)
								{
								NumDups --;
								Mem_Flag = TRUE;
								p12 = p;
								printf("\nOut of memory set aside for potential symmetries. Sorry!");
								goto _END;	
								}
							for(i = 'A'; i < 'A' + NumGenerators; i++) p12[i] = T2[i];
							for(i = 'a'; i < 'a' + NumGenerators; i++) p12[i] = T2[i];
							for(i = 1; i <= NumRelators; i++) p12[i] = T1[i];
							p12[CurrentR2] = Depth;
							p12[0] = Gen2;
							p12 = p;
							}
						else
							{
							NumDups --;
							Mem_Flag = TRUE;
							printf("\nOut of memory set aside for potential symmetries. Sorry!");
							goto _END;		
							}
						}					
					}
				else
					j = FALSE;
				}
			else
				j = TRUE;
								
			if(j)
				{
				if(x > y)			
_UPDATE:			{		
					/**************************************************************************
						If x > y, the current image string in Str is not optimal. If there
						were some new handles created on the basis of "matches" with this
						now known to be non-optimal image string, delete those new handles.
						If UpdateStr1 is TRUE, replace Str1 with a copy of Str2, set
						UpdateStr1 FALSE, and set CurrentR1 equal to CurrentR2. Then update
						the pointers p11, p12 and s. Change the last char of the image string
						from x to y. Finally, copy the data about the new substitution from
						T2 into T1.
					**************************************************************************/	
					
					if(NumDups > Stack_Ptr)
						{
						for(i = NumDups; i > Stack_Ptr; i--) DisposeHandle((char **) RWR[i]);
						NumDups = Stack_Ptr;
						}
					if(UpdateStr1)
						{
						UpdateStr1 = FALSE;
						CurrentR1 = CurrentR2;
						p = p22;
						p12 = Str1;
						p22 = Str2;
						while(*p12++ = *p22++) ;
						p22 = p;
						}
					p11 = Str1 + (p21 - Str2);	
	 				p12 = Str1 + (p22 - Str2);
					s --;
					*s++ = y;
					*s = EOS;
					for(i = 'A'; i < 'A' + NumGenerators; i++) T1[i] = T2[i];
					for(i = 'a'; i < 'a' + NumGenerators; i++) T1[i] = T2[i];
					T1[0] = Gen2;
					Gen1 = Gen2;
					}					
				
				/******************************************************************************
					Delete any changes made to T2, and Gen2, and reset NumChanges to 0.
					Then look for a new position from which to start reading Str2.
					Observe that if Str2 contains a subword of the form XZZZZZZ... with
					p21 pointing to the initital Z, then making p21 point to any other Z
					in this substring cannot lead to a new optimal image for Str2. Hence, we
					need only check the initial letter of any such substring as a starting
					point.	
				******************************************************************************/
				
				if(NumChanges)
					{	
					for(i = 0; i < NumChanges; i++)
						T2[Changed[i]] = T2[Changed[i] + 32] = EOS;
					Gen2 -= NumChanges;	
					NumChanges = 0;
					}
				p22 = p21;		
				if(Gen2)
					{
					y = *p22;
					while(y == *p22) p22++;
					y = *Str;
					if(Gen2 + 'A' > y)
						while((x = T2[*p22]) > y || x == EOS) p22++;
					else	
						while(T2[*p22] > y) p22++;
					}
				else
					{
					/*************************************************************************
							Locate the appearances of maximal exponents.
					*************************************************************************/	
					y = *p22;
					for(Exp = 0; y == *p22; Exp++,p22++) ;
					if(Exp > Max_Exp) Max_Exp = Exp;
					if(y = *p22) while(1)
						{
						for(Exp = 0; y == *p22; Exp++,p22++) ;
						if(Exp >= Max_Exp)
							{
							Max_Exp = Exp;
							p22 -= Exp;
							break;
							}
						if((y = *p22) == EOS)
							{
							p22--;
							y = *p22++;
							for(p21 = p22 - Exp,p = p22,p22 = Str2; y == *p22; Exp++,p22++) ;
							if(Exp >= Max_Exp)
								{
								Max_Exp = Exp;
								p22 = p21;
								}
							else
								p22 = p;		
							break;
							}			
						}
					}	
				p21 = p22;
				k = length;
				
				if(*p21 == EOS)
					{
					/**************************************************************************
						We have reached the end of Str2. If Str2 has not been inverted
						replace Str2 with its inverse, and set the Invert flag FALSE.
						Otherwise look for an unused relator with its length equal to the
						current value of length.
					**************************************************************************/	
					
					if(Invert)
						{
						Invert = FALSE;
						Inverse(Str2);
						}
					else
						{
						Invert = TRUE;
						CurrentR2 ++;
						F2 = FALSE;
						if(CurrentR2 <= NumRelators) for(i = CurrentR2; i <= NumRelators; i++)
							{
							if(LR[i] < length) break;
							if(!T1[i])
								{
								F2 = TRUE;
								CurrentR2 = i;
								p22 = *MyRelators[i];
								break;
								}
							}				
						if(F2)
							{
							/******************************************************************
								Since F2 has been changed from FALSE to TRUE, we found an
								unused relator with its length equal to the current
								value of length. Make a copy of this relator in Str2.
							******************************************************************/	
							
							p = p12;
							p12 = Str2;
							while(*p12++ = *p22++) ;
							p12 = p;
							}
						else
							{
							/******************************************************************
								The program may have been able to determine that the image of
								Str1 starting from the char *p11 is optimal without scanning 
								all of Str1 and without determining substitutions for all of
								the generators that appear in Str1. So finish reading Str1
								adding new substitutions to T1 if necessary. Then update the
								number of generators for which we know substitutions and flag
								the relator CurrentR1 as having been used.
							******************************************************************/
									
							if(Gen1 < NumGenerators)
								{
								p12 = p11;
								while(x = *p12++) if(!T1[x])
									{
									T1[x] = Gen1 + 'A';
									if(x > 'Z')
										T1[x - 32] = Gen1 + 'a';
									else
										T1[x + 32] = Gen1 + 'a';
									if(++Gen1 == NumGenerators) break;
									}
								if(Gen1 < NumGenerators)
									{			
									p12 = Str1;
									y = *p11;
									*p11 = EOS;
									while(x = *p12++) if(!T1[x])
										{
										T1[x] = Gen1 + 'A';
										if(x > 'Z')
											T1[x - 32] = Gen1 + 'a';
										else
											T1[x + 32] = Gen1 + 'a';
										if(++Gen1 == NumGenerators) break;
										}
									*p11 = y;	
									}	
								}	
									
							T1[0] = Gen1;
							T1[CurrentR1] = Depth;
							
							/******************************************************************
								Use the table of substitutions to translate Str1 into Str.
								The image in Str then becomes the initial optimal image of a
								relator at the current value of depth. Then break out of this
								while loop and process any other handles on the stack.
							******************************************************************/
								
							p22 = Str;
							p12 = p11;
							while(x = *p12++) *p22++ = T1[x];
							y = *p11;
							*p11 = EOS;
							p12 = Str1;
							while(x = *p12++) *p22++ = T1[x];
							*p22 = EOS;
							*p11 = y;
							
							break;
							}	
						}	
					UpdateStr1 = TRUE;
					p21 = Str2;
					}
				p22 = p21;
				s = Str;		
				}									
			}	
		HUnlock((char **) RWR[Stack_Ptr]);
		
		/**************************************************************************************
			Next, process any other handles which were on the stack when depth was
			incremented. Note that this code for processing the remaining handles on the stack
			differs from the code for processing the handle on the top of the stack in a few
			ways. For example, we no longer need Str1 and its associated pointers p11 and
			p12. What happens when a "match" occurs is also different.
		**************************************************************************************/
		
		Stack_Ptr --;
		for( ; Stack_Ptr >= 0; Stack_Ptr--)
			{
			HLock((char **) RWR[Stack_Ptr]);
			T1 = *RWR[Stack_Ptr];
			for(i = 'A'; i < 'A' + NumGenerators; i++) T2[i] = T1[i];
			for(i = 'a'; i < 'a' + NumGenerators; i++) T2[i] = T1[i];
			
			/**********************************************************************************
				Find the first relator which has not been used and which has its length
				equal to the current value of length. Make a copy of this relator in Str2.
			**********************************************************************************/
				
			for(i = 1; i <= NumRelators; i++)
				if(!T1[i] && LR[i] == length) break;
			CurrentR2 = i;	
			p22 = *MyRelators[CurrentR2];	
				
			p12 = Str2;
			while(*p12++ = *p22++) ;
			
			Gen1 = Gen2 = T1[0];
			NumChanges = 0;
			UpdateStr1 = TRUE;
			Invert = TRUE;
			Match = 0;
			p21 = p22 = Str2;
			s = Str;
			k = length;
			while(1)
				{
				y = *p22;
				if(T2[y])
					y = T2[y];
				else
					{
					T2[y] = Gen2 + 'A';
					if(y > 'Z')
						{
						T2[y - 32] = Gen2 + 'a';
						Changed[NumChanges++] = y - 32;
						}
					else
						{
						T2[y + 32] = Gen2 + 'a';
						Changed[NumChanges++] = y;
						}
					y = T2[y];
					Gen2 ++;
					}
						
				x = *s++;
				if(x == y)
					{
					p22 ++;
					if(*p22 == EOS) p22 = Str2;
					if(--k == 0)
						{
						j = TRUE;
						Match ++;
						if(Match == 1)
							{
							/******************************************************************
								The first "match" that occurs allows the handle RWR[Stack_Ptr]
								to survive this round at the current value of depth.
								Update the data in RWR[Stack_Ptr]. Record the fact that
								CurrentR2 was the relator used, at this depth, to get the
								optimal image, by setting T1[CurrentR2] to Depth. Note that
								if more matches occur we will need to know that T1[CurrentR2]
								was just marked as used. So set SaveR2 equal to CurrentR2.
							******************************************************************/
							
							if(Gen2 > Gen1)
								{
								T1[0] = Gen2;	
								for(i = 'A'; i < 'A' + NumGenerators; i++) T1[i] = T2[i];
								for(i = 'a'; i < 'a' + NumGenerators; i++) T1[i] = T2[i];
								}
							T1[CurrentR2] = Depth;
							SaveR2 = CurrentR2;	
							}
						else
						if(length > 1L)
							{
							NumDups ++;
							if(NumDups < MAXNUMDUPS)
								{
								/**************************************************************
									An additional match has been found, and there is room on
									the stack to save it. Here is where we need SaveR2. We must
									set p[SaveR2] FALSE before setting p[CurrentR2] to Depth.
									Otherwise, if CurrentR2 is no longer equal to SaveR2, then
									the data in this new handle would indicate that two
									relators had been used at this depth. Bad! 
								**************************************************************/
									
								RWR[NumDups] = (unsigned char **)NewHandle(125L);
								p12 = *RWR[NumDups];
								if(p12 == NULL)
									{
									NumDups --;
									Mem_Flag = TRUE;
									printf("\nOut of memory set aside for potential symmetries. Sorry!");
									goto _END;	
									}
								for(i = 'A'; i < 'A' + NumGenerators; i++) p12[i] = T2[i];
								for(i = 'a'; i < 'a' + NumGenerators; i++) p12[i] = T2[i];
								for(i = 1; i <= NumRelators; i++) p12[i] = T1[i];
								p12[SaveR2] = EOS;
								p12[CurrentR2] = Depth;
								p12[0] = Gen2;
								}
							else
								{
								NumDups --;
								Mem_Flag = TRUE;
								printf("\nOut of memory set aside for potential symmetries. Sorry!");
								goto _END;		
								}
							}					
						}
					else
						j = FALSE;
					}
				else
					j = TRUE;
									
				if(j)
					{
					if(x > y) 
						{
						/**********************************************************************
							If x > y, then the image string in Str, which we thought was the
							optimal image string, at this depth, is not optimal. It has just
							been bettered. Set UpdateStr1 to TRUE and jump to _UPDATE. We will
							then truncate the stack, treat RWR[Stack_Ptr] as the new top of
							the stack and look for a new optimal image string for this depth.
						**********************************************************************/	
							
						UpdateStr1 = TRUE;
						if(Match) T1[SaveR2] = EOS;
						goto _UPDATE;
						}
											
					/**************************************************************************
						Since j is TRUE and x > y is false, we either have x < y, or we have
						read all the way around Str2 and p22 = p21. In either case we want to
						delete any changes made to T2, and Gen2, and reset NumChanges to 0.
						Then we look for a new position from which to start reading Str2.
						Observe that if Str2 contains a subword of the form XZZZZZZ... with
						p21 pointing to the initital Z, then making p21 point to any other Z
						in this substring cannot lead to an optimal image for Str2. Hence, we
						need only check the initial letter of any such substring as a starting
						point.	
					**************************************************************************/
					
					if(NumChanges)
						{
						for(i = 0; i < NumChanges; i++)
							T2[Changed[i]] = T2[Changed[i] + 32] = EOS;
						Gen2 -= NumChanges;	
						NumChanges = 0;
						}
					p12 = p21;		
					y = *p12;
					while(y == *p12) p12++;
					y = *Str;
					if(Gen2 + 'A' > y)
						while((x = T2[*p12]) > y || x == EOS) p12++;
					else	
						while(T2[*p12] > y) p12++;
					p21 = p12;
					k = length;
					
					if(*p21 == EOS)
						{
						/**********************************************************************
							We have reached the end of Str2. If Str2 has not been inverted
							replace Str2 with its inverse, and set the Invert flag FALSE.
							Otherwise look for an unused relator with its length equal to the
							current value of length.
						**********************************************************************/	
						
						if(Invert)
							{
							Invert = FALSE;
							Inverse(Str2);
							}
						else
							{
							Invert = TRUE;
							CurrentR2 ++;
							F2 = FALSE;
							if(CurrentR2 <= NumRelators) for(i = CurrentR2; i <= NumRelators; i++)
								{
								if(LR[i] < length) break;
								if(!T1[i])
									{
									F2 = TRUE;
									CurrentR2 = i;
									p22 = *MyRelators[i];
									break;
									}
								}				
							if(F2)
								{
								/**************************************************************
									Since F2 has been changed from FALSE to TRUE, we found an
									unused relator with its length equal to the current value
									of length. Make a copy of this relator in Str2.
								**************************************************************/	
								
								p12 = Str2;
								while(*p12++ = *p22++) ;
								}	
							else
								{
								/**************************************************************
									No unused relators with length equal to the current
									value of length exist. That means we are done processing 
									the handle at RWR[Stack_Ptr]. If no match has been found,
									then the substitution given by the data in RWR[Stack_Ptr]
									is no longer competitive. So swap the handle RWR[Stack_Ptr]
									with the new current top of stack, namely RWR[NumDups],
									dispose of the handle RWR[NumDups] and decrement the number
									of things on the stack. Then break out of the while loop
									and process the next handle on the stack.
										Note that if we haven't deleted RWR[Stack_Ptr], then
									the data in RWR[Stack_Ptr] is current because it was
									updated when we found a match.
								**************************************************************/	
									 
								if(Match == 0)
									{
									Temp = RWR[Stack_Ptr];
									RWR[Stack_Ptr] = RWR[NumDups];
									RWR[NumDups] = Temp;
									HUnlock((char **) RWR[NumDups]);
									DisposeHandle((char **) RWR[NumDups]);
									NumDups --;
									}
								break;
								}	
							}	
						p21 = Str2;
						}
					p22 = p21;
					s = Str;		
					}									
				}
			HUnlock((char **) RWR[Stack_Ptr]);
			}
		if(NumDups == 0 && **RWR[0] == NumGenerators) break;		
		}									
	
	/******************************************************************************************
				If we aren't reporting symmetries, rewrite the Relators.
	******************************************************************************************/
		
	if(ReportSymmetries == FALSE)
		{
		T1 = *RWR[0];
		if(Compute_Stabilizers)
			{
			printf("\n        ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)	
				printf("%c",T1[i]);	
			fprintf(myout,"\n        ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)	
				fprintf(myout,"%c",T1[i]);
			}
		if(Micro_Print)
			{
			printf("\n\nRewrote the presentation using the substitution: ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)	
				printf("%c",T1[i]);
			}		
		for(i = 1; i <= NumRelators; i++)
			{
			s = *MyRelators[i];
			while(x = *s) *s++ = T1[x];
			}
		if(Annulus_Exists)
			{
			T1['@'] = '@';
			s = *Relators[0];
			while(x = *s) *s++ = T1[x];
			}			
		}
	else
		{
		/**************************************************************************************
				Report symmetries. But first, filter out duplicate substitutions.
		**************************************************************************************/
		for(i = NumDups; i > 0; i--)
		for(j = i-1; j >= 0; j--)
			{
			p12 = *RWR[i] + 65;
			p22 = *RWR[j] + 65;
			for(k = NumGenerators; k > 0 && *p12++ == *p22++; k--) ;
			if(k == 0)
				{
				Temp = RWR[i];
				RWR[i] = RWR[NumDups];
				RWR[NumDups] = Temp;
				DisposeHandle((char **) RWR[NumDups]);
				NumDups --;
				break;
				}
			}
							
		if(NumDups > 0)
			{
			if(Mem_Flag == FALSE)
				{
				if(WhichInput == -1)
					{
					if(NumDups == 1)
						{
						printf("\nThe following substitution induces an automorphism of the initial presentation.\n");
						if(NoReport == FALSE)
						fprintf(myout,"\nThe following substitution induces an automorphism of the initial presentation.\n");
						}
					else
						{
						printf("\nThe following substitutions induce automorphisms of the initial presentation.\n");
						if(NoReport == FALSE)
						fprintf(myout,"\nThe following substitutions induce automorphisms of the initial presentation.\n");
						}
					}
				else
					{	
					if(NumDups == 1)
						{
						printf("\nThe following substitution induces an automorphism of presentation %d:",
						WhichInput + 1);
						if(NoReport == FALSE)
						fprintf(myout,"\nThe following substitution induces an automorphism of presentation %d:",
						WhichInput + 1);
						}
					else
						{
						printf("\nThe following substitutions induce automorphisms of presentation %d:",
						WhichInput + 1);
						if(NoReport == FALSE)
						fprintf(myout,"\nThe following substitutions induce automorphisms of presentation %d:",
						WhichInput + 1);
						}
					}	
				HLock((char **) RWR[0]);
				T1 = *RWR[0];
				for(i = 1; i <= NumDups; i++)
					{
					HLock((char **) RWR[i]);
					Report_Symmetries(*RWR[i],i,NumDups);
					HUnlock((char **) RWR[i]);
					}
				HUnlock((char **) RWR[0]);
				}
			else
				if(NoReport == FALSE)
					{
					fprintf(myout,"\nWe ran out of memory set aside for potential substitutions");
					fprintf(myout," for presentation %d",WhichInput + 1);
					}		
			}
		}	
_END:					
	HUnlock((char **) Temp12);
	HUnlock((char **) Temp13);
	HUnlock((char **) Temp14);
	HUnlock((char **) Temp15);		
	for(i = 0; i <= NumDups; i++) DisposeHandle((char **) RWR[i]);
	
	if(Mem_Flag) return(TOO_LONG);
	
	/******************************************************************************************
				If we aren't reporting symmetries, call Final_Rewrite().
	******************************************************************************************/	
	
	if(ReportSymmetries == FALSE && Final_Rewrite(MyRelators) == TOO_LONG)
		return(TOO_LONG);
	
	return(NO_ERROR);			
}		

Final_Rewrite(MyRelators)
unsigned char	***MyRelators;
{
	/******************************************************************************************
		Final_Rewrite() is called by Canonical_Rewrite(). It proceeds through the relators
		and replaces each relator with a cyclic conjugate of itself or a cyclic conjugate of
		its inverse which comes first in our ordering. Then if there are relators of equal
		length, it sorts these relators to put them in order.
			Together, the routines Canonical_Rewrite() and Final_Rewrite() guarantee that two
		presentations have isomorphic Heegaard diagrams iff the presentations are identical.
		Thus determining whether two presentations give rise to isomorphic Heegaard diagrams
		becomes trivial, both for the program and the user.
	******************************************************************************************/
		
	register unsigned char 	*p11,
							*p12,
							*p22,
							x,
							y;
							
	register unsigned int 	j,
							k,
							length;
													
	unsigned char 			*p21,
							**Relator,
							**Temp;
	
	for(j = 1; j <= NumRelators; j++)
		{
		Relator = MyRelators[j];	
		
		length = LR[j] + 1;	
		if(length <= 2)
			{
			p11 = *Relator;
			x = *p11;
			if(x && x > 'Z') *p11 -= 32;
			continue;
			}
	
		/*************************************************************************************
		   Find a cyclic conjugate of the current Relator which is lexicographically minimal.
		*************************************************************************************/	
		
		HLock((char **) Relator);
		p11 = p12 = p21 = p22 = *Relator;
		p21++;
		p22++;
		k = length;
		while(1)
			{
			if(--k == 0) break;
			x = *p12++;
			y = *p22++;
			if(y == x)
				{
				if(*p12 == EOS) p12 = *Relator;
				if(*p22 == EOS) p22 = *Relator;	
				}
			else
				{
				if(y < x) p11 = p21;
				p12 = p21;
				y = *p12;
				while(*p12 == y) p12++;
				y = *p11;
				while(*p12 > y) p12++;
				if(*p12 == EOS) break;
				p21 = p12;
				p22 = p12;
				p12 = p11;
				k = length;	
				}		
			}
			
		/*************************************************************************************
			Save a copy of this lexicographically minimal cyclic conjugate in Temp6.
		*************************************************************************************/		
		
		ReallocateHandle((char **) Temp6,length);
		if((p22 = *Temp6) == NULL)
			{
			HUnlock((char **) Relator);
			return(TOO_LONG);
			}
		p12 = p11;
		while(*p22++ = *p12++) ;
		p22--;
		p12 = *Relator;
		x = *p11;
		*p11 = EOS;
		while(*p22++ = *p12++) ;
		*p11 = x;
		
		/*************************************************************************************
			Invert the current Relator and find a cyclic conjugate of the inverse of the
			current Relator which is lexicographically minimal.
		*************************************************************************************/	
		
		Inverse(*Relator);
		p11 = p12 = p21 = p22 = *Relator;
		p21++;
		p22++;
		k = length;
		while(1)
			{
			if(--k == 0) break;
			x = *p12++;
			y = *p22++;
			if(y == x)
				{
				if(*p12 == EOS) p12 = *Relator;
				if(*p22 == EOS) p22 = *Relator;	
				}
			else
				{
				if(y < x) p11 = p21;
				p12 = p21;
				y = *p12;
				while(*p12 == y) p12++;
				y = *p11;
				while(*p12 > y) p12++;
				if(*p12 == EOS) break;
				p21 = p12;
				p22 = p12;
				p12 = p11;
				k = length;	
				}	
			}
				
		/*************************************************************************************
			Save a copy of this lexicographically minimal cyclic conjugate in Temp7.
		*************************************************************************************/			
		
		ReallocateHandle((char **) Temp7,length);
		if((p22 = *Temp7) == NULL)
			{
			HUnlock((char **) Relator);
			return(TOO_LONG);
			}
		p12 = p11;
		while(*p22++ = *p12++) ;
		p22--;
		p12 = *Relator;
		x = *p11;
		*p11 = EOS;
		while(*p22++ = *p12++) ;
		*p11 = x;
		
		/*************************************************************************************
			Compare the strings we saved in Temp6 and Temp7. Replace the current Relator
			with whichever string is lexicographically minimal.
		*************************************************************************************/	
		
		HUnlock((char **) Relator);
		p22 = *Temp6 + length - 1;
		*p22 = 125;
		for(p11 = *Temp6,p12 = *Temp7; *p11 == *p12; p11++,p12++) ;
		*p22 = EOS;
		if(*p11 <= *p12)
			{
			Temp = Relator;
			Relator = Temp6;
			Temp6 = Temp;
			}
		else
			{
			Temp = Relator;
			Relator = Temp7;
			Temp7 = Temp;
			}
						
		MyRelators[j] = Relator;
				
		}
		
	/******************************************************************************************
		If there are relators with the same length, sort these relators to put them in
		lexicographic order.
	******************************************************************************************/	
	
	for(j = 1; j < NumRelators; j++)
		{
		length = LR[j];
		if(length > 0L)
		for(k = j+1; k <= NumRelators && LR[k] == length; k++)
			{
			p22 = *MyRelators[j] + length;
			*p22 = 125;
			for(p11 = *MyRelators[j],p12 = *MyRelators[k]; *p11 == *p12; p11++,p12++) ;
			*p22 = EOS;
			if(*p11 > *p12)
				{
				Temp = MyRelators[j];
				MyRelators[j] = MyRelators[k];
				MyRelators[k] = Temp;
				}	
			}	
		}
	
	return(NO_ERROR);		
}			
	
Find_Symmetries(Flag)
int		Flag;
{
	/******************************************************************************************
			This routine calls Canonical_Rewrite() for each nontrivial presentation that the
		program found and saved in SUR[][]. If Canonical_Rewrite() terminates with more than
		one entry on its stack of optimal substitutions, then there exists a substitution
		which induces a symmetry of the presentation.
			If there are symmetries present in a presentation, Canonical_Rewrite() calls 
		Report_Symmetries() which prints them out.
			The routine will also report symmetries which are present in the initial set of 
		input relators.	
	******************************************************************************************/
		
	register unsigned char 	*p,
							*q;
							
	int 					i,
							SNoReport;
	
	SNoReport = NoReport;
	NoReport = FALSE;
	printf("\n\n    Looking for symmetries. . .\n");
	NumSymmetries = 0;
	
	if(Modified_Init_Pres || Flag)
		{
		WhichInput = -1;
		NumRelators = CopyNumRelators;
		NumGenerators = CopyNumGenerators;
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));				
			if((q = *Relators[i]) == NULL) return(TOO_LONG);
			p = *Copy_Of_Input[i];
			while(*q++ = *p++) ;
			}
		Canonical_Rewrite(Relators,FALSE,TRUE);
		}
	
	if(Flag == FALSE)
		{	
		for(WhichInput = 0; WhichInput < NumFilled; WhichInput++)
		if(SURL[WhichInput] != 0L && UDV[WhichInput] != DUPLICATE)
			{			
			NumRelators = NR[WhichInput];
			NumGenerators = NG[WhichInput];
			for(i = 1; i <= NumRelators; i++)
				{
				ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[WhichInput][i]));				
				if((q = *Relators[i]) == NULL) return(TOO_LONG);
				p = *SUR[WhichInput][i];
				while(*q++ = *p++) ;
				}
			Canonical_Rewrite(Relators,FALSE,TRUE);	
			}
		}	
		
	if(NumSymmetries == 0)
		{
		printf("\n\n    No non-trivial symmetries exist.");
		fprintf(myout,"\n\n    No non-trivial symmetries exist.");
		}
	NoReport = SNoReport;	
	if(NoReport == TRUE && NumSymmetries)
		{
		printf("\n\n    Info about these symmetries has been printed in 'Heegaard_Results'.");
		printf("\n\nALSO SAVE THE CURRENT PRESENTATIONS TO 'Heegaard_Results' ?  HIT 'y' OR 'n'.");
		GET_RESPONSE1:
		switch(WaitkbHit())
			{
			case 'y':
				Report(Band_Sums,NumDiagrams,OnStack,0,1,1,0,0,1,0);
				break;
			case 'n':
				break;
			default:
				SysBeep(5);
				goto GET_RESPONSE1;
			}	
		}
	return(NO_ERROR);				
}	

void Report_Symmetries(unsigned char *T2,int j,int NumDups)		
{
	/******************************************************************************************
		This routine is called by Canonical_Rewrite() when we are reporting symmetries of 
		presentations to the user.
	******************************************************************************************/
		
	unsigned char 	*T;
	
	int 			i;

	T = (unsigned char *) NewPtr(100L);
	NumSymmetries ++;
	
	for(i = 'A'; i < 'A' + NumGenerators; i++)
		{
		if(T2[i])
			T[i] = T2[i];
		else
			T[i] = i;	
		}
		
	if(NumDups == 1)
		{
		printf("\n       ");
		for(i = 'A'; i < 'A' + NumGenerators; i++)
			printf("%2c",i);
		printf("\n       ");
		for(i = 'A'; i < 'A' + NumGenerators; i++)	
			printf("%2c",T[i]);
		if(NoReport == FALSE)
			{		
			fprintf(myout,"\n       ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)
				fprintf(myout,"%2c",i);
			fprintf(myout,"\n       ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)	
				fprintf(myout,"%2c",T[i]);
			}
		}
	else
		{
		printf("\n%4d)  ",j);
		for(i = 'A'; i < 'A' + NumGenerators; i++)
			printf("%2c",i);
		printf("\n       ");
		for(i = 'A'; i < 'A' + NumGenerators; i++)	
			printf("%2c",T[i]);
		if(NoReport == FALSE)
			{		
			fprintf(myout,"\n%4d)  ",j);
			for(i = 'A'; i < 'A' + NumGenerators; i++)
				fprintf(myout,"%2c",i);
			fprintf(myout,"\n       ");
			for(i = 'A'; i < 'A' + NumGenerators; i++)	
				fprintf(myout,"%2c",T[i]);
			}
		}	
	DisposePtr((char *) T);	
}


void Find_Cancellation_Paths(void)
{
	/******************************************************************************************
		Find_Cancellation_Paths() finds those paths which join the major faces of the
		Heegaard diagram and are disjoint from the relators.
	******************************************************************************************/
		
	register unsigned char 	*p,
							*q,
							*r;
							
	register unsigned int 	c,
							d,
							e,
							v,
							w;
							
	unsigned char 			x,
							y,
							**SP,
							sx,
							sy,
							tl,
							tr;
													
	int 					Error,
							ss;
							
	unsigned int 			edge,
							EL[6*VERTICES],
							ee,
							h,
							i,
							j,
							k,
							NumPaths,
							vertex,
							vertexLS,
							vertexRS,
							vertexLE,
							vertexRE;
	
	long 					HS,
							LOR[MAXNUMRELATORS + 1],
							length;
	
	unsigned long			max;
	
	Error = 0;
	NumPaths = 0;
	
	r = (unsigned char *) NewPtr(100L);		
	printf("\n\nENTER A PRESENTATION FROM 0 TO %u FOR WHICH YOU WANT PATHS AND HIT 'return'.     ",NumFilled);
	for(i = j = 0; j < NumFilled; j++) if(SURL[j] == 0L) i ++;
	if(i)
		{
		if(i == 1)
			{
			for(i = 0; i < NumFilled && SURL[i]; i++) ;
			i++;
			printf("\n\nExcept for presentation %d which is a presentation of S1 X S2 or S1 X D2 and is empty.     ",i);
			}
		else
			{
			j = 0;
			j += printf("\n\nExcept for presentations: ");
			for(h = 0,k = 1; h < NumFilled; h++) if(SURL[h] == 0L)
				{
				h++;
				j += printf("{%d,",h);	
				break;
				}
			for( ; h < NumFilled; h++) if(SURL[h] == 0L)
				{
				if(++k < i)
					j += printf("%d,",h+1);
				else
					j += printf("%d}.",h+1);	
				if(j > 80)
					{
					j = 0;
					printf("\n");
					}
				}	
			printf("\nThese are presentations of S1 X S2 (s) or S1 X D2 (s) and are empty.     ");
			}
		}
	GET_RESPONSE1:
	WhichInput = -1;		
	ReadString((char *)r, GetPtrSize(r));
	sscanf((char *) r,"%d",&WhichInput);	
	if(WhichInput == 0)
		{
		printf("\n\nFinding paths of the original presentation.\n");
		NumRelators = CopyNumRelators;
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));				
			p = *Copy_Of_Input[i];
			if((q = *Relators[i]) == NULL)
				{
				DisposePtr((char *) r);
				printf("\n\n	Memory Error. Sorry!");
				Error = 1;
				goto END;
				}
			while(*q++ = *p++) ;
			}
		}
	else
		{
		if(WhichInput < 1 || WhichInput > NumFilled || SURL[WhichInput-1] == 0L)
			{
			SysBeep(5);
			goto GET_RESPONSE1;
			}	
		printf("\n\nFinding paths of presentation %d.\n",WhichInput);

		WhichInput --;
		NumRelators = NR[WhichInput];
		NumGenerators = NG[WhichInput];
		Vertices = 2*NumGenerators;
		Length = SURL[WhichInput];
		Saved_Vertices = 0;
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[WhichInput][i]));				
			p = *SUR[WhichInput][i];
			if((q = *Relators[i]) == NULL)
				{
				DisposePtr((char *) r);
				printf("\n\n	Memory Error. Sorry!");
				Error = 1;
				goto END;
				}			
			while(*q++ = *p++) ;
			}
		}
	DisposePtr((char *) r);
		
	if(Find_Flow_A(NORMAL,FALSE))
		{
		printf("\n\nUnable to find paths. Sorry!");
		Error = 1;
		goto END;
		}
	if(Whitehead_Graph())
		{
		printf("\n\nUnable to find paths. Sorry!");
		Error = 1;
		goto END;
		}
	
	for(d = 1,max = 0L; d <= NumRelators; d++) if(LR[d] > max) max = LR[d];
	if(SRError == 3) for(d = 1; d <= NumRelators; d++)
		LOR[d] = GetHandleSize((char **) OutRelators[d]);
	
	ReallocateHandle((char **) Temp1,max + 2);
	if(*Temp1 == NULL)
		{
		printf("\n\nNot enough memory available. Sorry!");
		Error = 1;
		goto END;
		}
				
	for(d = 1; d <= 2*NumEdges; d++) EL[d] = d;
	NumPaths = 0;
	SP = (unsigned char **) NewPtr(sizeof(long)*(NumEdges + 1));
	if(SP == NULL)
		{
		printf("\n\nNot enough memory available. Sorry!");
		Error = 1;
		goto END;
		}
	
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
		vertexRE = CO[v][vertexLE];	
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
		edge    = ee % V[v];
		e 		= edge;
		r		= *Temp1;
		
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
			r++;	
			e = OSA[v] - e;
			if(e >= V[v]) e -= V[v];
			length ++;
			if(length > max)
				{			
				printf("\n\nError in finding paths. Sorry!");
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
				vertexRS = v;
				vertexLS = CO[w][v];
				
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
				if(vertexRE & 1)
					y = (vertexRE >> 1) + 97;
				else
					y = (vertexRE >> 1) + 65;
				if(vertexLS & 1)
					sx = (vertexLS >> 1) + 97;
				else
					sx = (vertexLS >> 1) + 65;
				if(vertexRS & 1)
					sy = (vertexRS >> 1) + 97;
				else
					sy = (vertexRS >> 1) + 65;
				
				HS = r - *Temp1;
				NumPaths ++;				
				SP[NumPaths] = (unsigned char *) NewPtr(HS);		
				if((q = SP[NumPaths]) == NULL)
					{
					NumPaths --;
					printf("\n\nNot enough memory available. Sorry!");
					Error = 2;
					goto END;
					}
				p = *Temp1;
				while(*q++ = *p++) ;	
									
				printf("\n %3u) %c %c %c %c %c %c: %s",NumPaths,x,y,sx,sy,tl,tr,*Temp1);
				fprintf(myout,"\n %3u) %c %c %c %c %c %c: %s",NumPaths,x,y,sx,sy,tl,tr,*Temp1);			
				break;
				}
			e = B[w][v] - e;
			}	
		while(v != vertex || e != edge);	
	}

END:	
	printf("\n\n HIT ANY KEY TO CONTINUE.");
	WaitkbHit();

	if(Error == 0) CP_Concatanate_Paths(SP,NumPaths);	

	for(j = 1; j <= NumPaths; j++) DisposePtr((char *) SP[j]);
	if(Error != 1) DisposePtr((char *) SP);
}

void CP_Concatanate_Paths(unsigned char **MySP,int NumPaths)
{
	char			*r;
	
	unsigned char	*p,
					*q;
					
	int				i,
					j,
					NewEntry,
					NumEntries,
					PathList[200];
	
	unsigned long	HS;

	ANOTHER_PATH:	
	printf("\n\nPlease enter the numbers of the paths you wish to concatenate.");
	printf("\nTerminate path entry by entering '0'.\n");
	NumEntries = 0;
		
	GET_RESPONSE1:
	r = (char *) NewPtr(100);		
	ReadString((char *)r, GetPtrSize(r));	
	sscanf((char *) r,"%d",&NewEntry);
	DisposePtr(r);
	if(NewEntry)
		{
		if(NewEntry < - NumPaths || NewEntry > NumPaths)
			{
			SysBeep(5);
			printf("\nThe last entry was out of bounds. Please reenter it.");
			goto GET_RESPONSE1;
			}
		PathList[NumEntries++] = NewEntry;
		goto GET_RESPONSE1;
		}
	
	if(NumEntries == 0)
		{
		GET_RESPONSE2:
		printf("\n\nTry another path? Hit 'y' or 'n'.");
		switch(WaitkbHit())
			{
			case 'y':
				goto ANOTHER_PATH;
			case 'n':
				return;
			}			
		}

	for(i = 0,HS = 1L; i < NumEntries; i++)	
		{
		j= abs(PathList[i]);
		HS += GetPtrSize((char *) MySP[j]) - 1;
		}
	ReallocateHandle((char **) Relators[1],HS);
	HLock((char **) Relators[1]);
	if((q = *Relators[1]) == NULL)
		{
		printf("\nMemory Error. Sorry!");
		HUnlock((char **) Relators[1]);
		return;
		}
	for(i = 0; i < NumEntries; i++)
		{
		j = PathList[i];
		if(j > 0)
			{
			p = MySP[j];
			while(*q++ = *p++) ;
			q--;
			}
		else
			{
			j = -j;
			Inverse(MySP[j]);
			p = MySP[j];
			while(*q++ = *p++) ;
			q--;
			Inverse(MySP[j]);
			}	
		}
	HUnlock((char **) Relators[1]);
				
	j = CP_Find_Primitives();
	switch(j)
		{
		case 0:
			printf("\n This path has a good graph.");
			break;
		case 1:
			printf("\n The graph of this path has a cut vertex.");
			break;
		case 2:
			printf("\n The graph of this path is not connected.");
			break;
		case TOO_LONG:
			printf("\n Memory error. Sorry!");
			break;
		}
	switch(j)
		{
		case 1:
		case 2:
			printf("\n Augment path by adding more edges? Hit 'y' or 'n'.");
			GET_RESPONSE3:
			switch(WaitkbHit())
				{
				case 'y':
					goto GET_RESPONSE1;
					break;
				case 'n':
					goto ANOTHER_PATH;
					break;
				}	
		}		
	
	goto ANOTHER_PATH;				
}

CP_Fill_AA()
{
	/******************************************************************************************
		This routine takes Relators[1] and examines each pair of consecutive letters that
		appear in the relator. For each such pair of letters, it adds the appropriate edge
		to the array AA[][].
	******************************************************************************************/
	
	register unsigned char 	i,
							j,
							*p,
							x;
							
	register unsigned int	*q;							
	
	for(i = 0; i < Vertices; i++)
	for(j = 0,q = AA[i]; j < Vertices; j++,q++) *q = 0;
	
	p = *Relators[1];
	if(*p == EOS) return;
	x = *p << 1;
	if(x < 194) x -= 129;
	else x -= 194;
	i = x;
	p++;
	while(x = *p)
		{
		x = x << 1;
		if(x < 194) x -= 130;
		else x -= 193;
		AA[i][x]++;
		if(x & 1) i = x - 1;
		else i = x + 1;
		p++;
		}
		
	for(i = 0; i < Vertices - 1; i++)
	for(j = i + 1; j < Vertices; j++)
		{
		AA[i][j] += AA[j][i];
		AA[j][i] = AA[i][j];
		}						
}

int CP_Find_Primitives(void)
{
	/******************************************************************************************
		This procedure is used to determine whether the relator Relators[1] is a primitive or a
		proper power of a primitive. It exploits the fact that if Relators[1] is a primitive or
		a proper power of a primitive, then the Whitehead graph of Relators[1] is either not
		connected or has a separating vertex. (This is assuming, as we may, that more than one
		generator appears in Relators[1].) If the parameter flag = 2, then the routine is
		being called by CheckPrimitivity() to check whether Relators[1] is a primitive or
		proper power of a free generator. If the parameter flag = 1, then this routine is
		being called by Reduce_Genus() under circumstances where we know that a primitive
		exists. In this case, we want to apply enough automorphisms to the relators to reduce
		Relators[1] to a defining relator and then call Defining_Relator() to make the
		substitution which deletes appearances of Relators[1] in the other relators.
		If the parameter flag = 3, this routine is being called by Reduce_Genus when we
		know that Relators[1] is a proper power. As in the case when the parameter flag
		= 1, we want to perform enough automorphisms to reduce Relators[1] to a proper
		power and then make the appropriate substitutions in the remaining relators.
	******************************************************************************************/
	
	register unsigned char	*p;
		
	register unsigned int 	h,
							i,
							j,
							k,
							m,
							*q;
							
	unsigned int 			C[125],
							VG[VERTICES],
							root;
							
	unsigned long			SLength1,
							SLength2;

	/******************************************************************************************
						Check whether Relators[1] is a defining relator.
	******************************************************************************************/		
	
	for(i = 0; i < NumGenerators; i++) C[i+'A'] = C[i+'a'] = 0;
	p = *Relators[1];
	while(*p)
		C[*p++]++;
	for(i = j = 0; i < NumGenerators; i++,j+=2)
		{
		C[i+'A'] += C[i+'a'];
		C[j] = C[i+'A'];
		if(C[j] == 1) return(TRUE);
		}

	SLength1 = GetHandleSize((char **) Relators[1]);
		
	while(1)
		{
		CP_Fill_AA();
		for(i = 0; i < Vertices; i++)
			{
			for(j = k = 0; j < Vertices; j++) if(AA[i][j])
				{
				AJ3[i][k] = j;
				k ++;
				}
			AJ3[i][k] = VERTICES;
			}
			
		/**************************************************************************************
		 						Check whether the graph is connected.
		 *************************************************************************************/
		 
		for(i = 0; i < Vertices; i++) ZZ[i] = 0;
		
		if(Connected_AJ3(0,0) == FALSE) return(2);
				
		/**********************************************************************************
			The graph is connected. Use vertex '0' as the root of a depth-first-search of
			the graph.
		**********************************************************************************/	
		
		for(j = 0; j < Vertices; j++) Number[j] = 0;
		NumVert 	= 1;
		Number[0] 	= 1;
		Lowpt[0] 	= 1;
		Father[0] 	= 0;
		VG[0] 		= 0;
		root		= 0;
		for(q = UpDate,*q = root; q >= UpDate; q--)
			{
			NEW_VERT:
			i = *q;
			for(k = VG[i]; (j = AJ3[i][k]) < VERTICES; k++)
				{
				if(Number[j] == 0)
					{
					NumVert 	++;
					Number[j] 	= NumVert;
					Lowpt[j] 	= NumVert;
					Father[j] 	= i;
					VG[j] 		= 0;
					VG[i] 		= k + 1;
					q 			++;
					*q 			= j;
					goto 		NEW_VERT;		
					}
				if(j != Father[i] && Number[j] < Lowpt[i]) Lowpt[i] = Number[j];		
				}
			if(Lowpt[i] < Lowpt[Father[i]]) Lowpt[Father[i]] = Lowpt[i];
			}
		
		/**********************************************************************************
			Check whether the root of the depth-first-search tree is a cut vertex. This 
			will be the case iff the root has more than one son.
		**********************************************************************************/
		
		for(i = k = 0; (j = AJ3[root][i]) < VERTICES; i++)
			if(Father[j] == root && ++k > 1) break;
		
		/**********************************************************************************
			If k > 1, the root of the depth-first-search tree is a cut vertex.
			Check whether the root is a non-trivial cut vertex.		
		**********************************************************************************/
		
		if(k > 1)
			{
			for(i = 0; i < Vertices; i++) ZZ[i] = 0;
			ZZ[0] = 3;
			Connected_AJ3(1,1);
			for(i = k = 0; (h = AJ3[root][i]) < VERTICES; i++)
				if(Father[h] == root && ++k > 2) goto DO_AUT;
			for(i = k = 0; i < Vertices; i++) if(ZZ[i] == 1 && ++k > 1) break;
			}
		
		/*********************************************************************************
			If k <= 1, the root of the depth-first search tree is a trivial cut vertex.
			Look for other cut vertices.
		**********************************************************************************/
					
		if(k <= 1) for(j = 0; j < Vertices; j++)
			{
			k = Father[j];
			if(Lowpt[j] >= Number[k] && k != root)
				{
				/*************************************************************************
					The vertex Father[j] is a cut vertex, check whether it is non-trivial.
				*************************************************************************/
				
				for(i = 0; i < Vertices; i++) ZZ[i] = 0;
				ZZ[k] = 3;
				if(k & 1)
					Connected_AJ3(k-1,1);
				else
					Connected_AJ3(k+1,1);
				for(i = m = 0; (h = AJ3[k][i]) < VERTICES; i++)
					if(Father[h] == k && ++m > 1) goto DO_AUT;		
				for(i = k = 0; i < Vertices; i++) if(ZZ[i] == 1 && ++k > 1) goto DO_AUT;
				}
			}	
		if(j == Vertices) return(FALSE);
		
		/**********************************************************************************
			If j = Vertices, then there is no "cut vertex" and Relators[1] is not a
			primitive or a proper power of a generator. So we return FALSE. Otherwise
			the vertex Father[j] is a "cut vertex" and there is an automorphism which
			reduces the length of Relators[1]. Perform the appropriate automorphism
			and, if necessary, return to the top of this loop.
		**********************************************************************************/
			
		DO_AUT:
		return(TRUE);	/* Remove this line in order to perform automorphisms. */
		i = Father[j];
		if(i & 1)
			{
			ZZ[i] = 0;
			for(j = 0; j < Vertices; j++)
				{
				if(ZZ[j] == 1)
					ZZ[j] = 0;
				else
				if(ZZ[j] == 0)
					ZZ[j] = 1;
				}	

			if(CP_Do_Aut(i-1) == TOO_LONG) return(TOO_LONG);
				
			SLength2 = GetHandleSize((char **) Relators[1]);
			C[i-1] -= SLength1 - SLength2;			
			if(C[i-1] == 1) return(TRUE);
			SLength1 = SLength2;		
			}
		else
			{
			ZZ[i] = 0;
				
			if(CP_Do_Aut(i) == TOO_LONG) return(TOO_LONG);
				
			SLength2 = GetHandleSize((char **) Relators[1]);
			C[i] -= SLength1 - SLength2;			
			if(C[i] == 1) return(TRUE);
			SLength1 = SLength2;			
			}	
		}	
}

int CP_Do_Aut(unsigned int Source)
{		
	/******************************************************************************************
		This routine performs the automorphism(s) determined by the routines Find_Flow_A()
		Find_Primitives() and Level_Transformations(). The exact automorphism to be performed
		is specified by the entries of the global array ZZ[]. The routine performs a T-
		transformation corresponding to a splitting of the vertices of the Whitehead graph
		into two subsets: those vertices accessible from the "sink" and the remaining vertices.
		If 'A' is the source vertex and 'a' is the sink vertex. and X is a generator, (X != A),
		then the T-transformation acts on X according to the following table.
			1)	If 'X' not accessible and 'x' not accessible, then X: --> AXa.
			2)  If 'X' not accessible and 'x'     accessible, then X: --> AX.
			3)  If 'X'     accessible and 'x'     accessible, then X: --> X.
			4)  If 'X'     accessible and 'x' not accessible, then X: --> Xa.
		Perhaps a better way to describe Do_Aut, is to think of it as acting on oriented edges
		of the Whitehead graph. If E is an edge and both ends of E are accessible or both
		ends of E are nonaccessible we do nothing. Otherwise either an "A" or an "a" is 
		inserted in the current relator depending upon the orientation of E.	
		Going in to the routine, the total number of appearances of "A"(s) and "a"(s) in the
		relators was equal to the valence of vertex 'A'. After the T-transformation, the total
		number of appearances of "A"(s) and "a"(s) in the relators has been changed to a 
		value equal to the flow from vertex 'A' to vertex 'a' in the Whitehead graph. The
		appearances of the remaining generators in the relators are unchanged. We also note 
		the pleasant property of T-transformations that if the original relators are
		freely reduced, then the modified relators are also freely reduced.	
	******************************************************************************************/
	
	register unsigned char 	A,
							a,
							*p,
							*q,
							x,
							y;
							
	register 				int i;
	
	unsigned char 			TX[125],
							TY[125],
							**Temp;
	
	unsigned long			HS;						
			
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
	A = ((Source >> 1) + 65);
	a = A + 32;
			
	HS = GetHandleSize((char **) Relators[1]);
	if(HS > MAXLENGTH) return(TOO_LONG);
	ReallocateHandle((char **) Temp5,HS);
	if((p = *Temp5) == NULL) return(TOO_LONG);
	q = *Relators[1];
	if(*q == EOS)
		return(0);
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
	*p = EOS;
	q = *Temp5;
	HS = p + 1 - q;
	if(HS > MAXLENGTH) return(TOO_LONG);
	SetHandleSize((char **) Temp5,HS);
	Temp = Relators[1];
	Relators[1] = Temp5;
	Temp5 = Temp;
		
	TotalAuts ++;	
	return(0);				
}	
