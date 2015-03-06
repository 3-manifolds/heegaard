#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   21 Find_Level_Transformations(int F1,int F2)
L  132 Find_Level_Flow(unsigned int Source)
L  215 Init_SCdfsR(void)
L  258 SCdfsR(unsigned char w)
L  304 Init_TCR(void)
L  355 TCdfsR(unsigned char v)
L  384 Get_LTs(unsigned int Source,int F1,int F2)
L  451 Do_LT(unsigned int Source,int F1,int F2)
********************************************************************************************/

unsigned char	NumVisited,
				SCompNum,
				Stptr;
				
unsigned int	SEL_Num;				
					
int Find_Level_Transformations(int F1,int F2)
{
	/******************************************************************************************
		This routine can be used to find the complete orbit of a minimal length presentation P
		under level-transformations, provided P has a connected Whitehead graph, and the orbit 
		of P has no more than MAX_SAVED_PRES members. (Currently MAX_SAVED_PRES = 50,000.)
			The routine Init_Find_Level_Transformations() should be called before this routine
		is called. Init_Find_Level_Transformations() checks if the presentation P has minimal 
		length, and also checks if the Whitehead graph of P is connected.
	******************************************************************************************/
	
	unsigned char	i,
					j,
					k,
					*p,
					*q;
					
	int				Return_Value;
	
	unsigned int	Source;
					
    /******************************************************************************************
        Copy the presentation which we want level transformations for into Relators[].
    ******************************************************************************************/
    
    if(F2 == 2)
    	{
		for(i = 1; i <= NumRelators; i++)
			{
			if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
			Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SLP[ReadPres][i]));
			if(Relators[i] == NULL) Mem_Error(); 
			q = *Relators[i];   
			p = *SLP[ReadPres][i];
			while( (*q++ = *p++) ) ;
			} 	
    	}
    else
    	{	        
	    for(i = 1; i <= NumRelators; i++)
			{
			if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
			Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][i]));
			if(Relators[i] == NULL) Mem_Error();
			q = *Relators[i];    
			p = *SUR[ReadPres][i];
			while( (*q++ = *p++) ) ;
			}	
		}
        
    /******************************************************************************************
        Call Fill_A(NumRelators) and ComputeValences_A() to initialize the arrays A[][],
        and VA[]. Then set up AJ3[][] for use by Find_Level_Flow(), Init_SCdfsR(), and 
        SCdfsR() in this file.
    ******************************************************************************************/
    
    Fill_A(NumRelators);
    ComputeValences_A();        
    for(i = 0; i < Vertices; i++)
        {
        for(j = k = 0; j < Vertices; j++) if(A[i][j])
            {
            if(i == j) continue;
            AJ3[i][k] = j;
            k++;
            }
        AJ3[i][k] = VERTICES;
        }  
        
    /******************************************************************************************
            If F2 = 2, we check whether the Whitehead graph of this presentation is connected.
            It is important to do this because Find_Level_Transformations() will not work
            correctly if the Whitehead graph of the presentation is not connected!!
    ******************************************************************************************/
    
    if((NumFilled == 1) && (F2 == 2))
    	{                         
		for(i = 0; i < Vertices; i++) ZZ[i] = 0;
		if(Connected_AJ3(0,0) == FALSE) return(NOT_CONNECTED);	 
        }            
        
    /******************************************************************************************
            Find the level T-transformations corresponding to each generator in turn.
    ******************************************************************************************/
  
  	Source = 0;
  	if(NumFilled == 1) Left[0] = Right[0] = INFINITE;
  	
  	while(1)
  		{
  		Find_Level_Flow(Source);
  		Init_SCdfsR();
  		if(SCompNum > 3)
  			{			
  			Init_TCR();
  			if( (Return_Value = Get_LTs(Source,F1,F2)) ) return(Return_Value);
  			}
  		Source += 2;
  		if(Source >= Vertices) break;
  		
		/**************************************************************************************
		  If Source < Vertices, we need to restore A[][] before processing the next generator.
		**************************************************************************************/
		
  		for(i = 0; i < Vertices; i++)
  		for(k = 0; (j = AJ3[i][k]) < VERTICES; k++) A[i][j] = A[j][i] = (A[i][j] + A[j][i]) >> 1;
  		}
  		
	return(0);
}

int Find_Level_Flow(unsigned int Source)
{
    /******************************************************************************************
        This routine is a variant of Find_Flow_A(). However, unlike Find_Flow_A(), this
        routine only finds a maximal flow in the Whitehead graph between a single pair of
        vertices, namely the source and sink passed to it. (Note that we reverse the roles of
        Source and Sink here because it makes things slightly more convenient for following
        routines.)
    ******************************************************************************************/
            
    register unsigned int   i,
    					    j,
    			    		max,
			    			min,
			    			*p,
			    			*q,
			    			*r;
                            
    unsigned int            Flow,
    			    		MaxFlow,
    			    		S[VERTICES],
    			    		Sink;
            
	Sink = Source;
	Source = Sink + 1;
    MaxFlow = VA[Source/2];
    if(MaxFlow == 0) return(2);
    Flow = A[Source][Sink];
    A[Source][Sink] = 0;
    A[Sink][Source] *= 2;                
    while(1)
        {        
        for(i = 0,p = ZZ,q = InQueue; i < Vertices; i++,p++,q++) *p = *q = 0;
        ZZ[Sink] = INFINITE;
        InQueue[Source] = TRUE;
        for(r = UpDate,*r = Sink,p = r + 1; r < p; r++) 
            {
            i = *r;
            InQueue[i] = FALSE;
            max = ZZ[i];
            if(max > ZZ[Source]) for(q = AJ3[i]; (j = *q) < VERTICES; q++)
                {
                if(max > ZZ[j])
                    {
                    if(A[j][i] < max)
                        min = A[j][i];
                    else
                        min = max;    
                    if(min > ZZ[j] && min > ZZ[Source])    
                        {
                        ZZ[j] = min;
                        S[j] = i;
                        if(!InQueue[j])
                            {
                            InQueue[j] = TRUE;
                            *p++ = j;
                            }
                        }
                    }
                }
            }                                                    
        max = ZZ[Source];
        Flow += max;
        if(max)
            {
            i = Source;
            while(i != Sink)
                {
                j = S[i];
                A[i][j] -= max;                                                                
                A[j][i] += max;                                    
                i = j;
                }        
            }
        else
            {
            if(Flow < MaxFlow) return(1);    /* The presentation does not have minimal length. */    
            return(0);    
            }
        if(Flow >= MaxFlow)    return(0);                                                
        }                                                                                                            
}
       			
void Init_SCdfsR(void)
{
	/********************************************************************************
	  This routine initializes Tarjan's Depth-First-Search Strong-Component routine.
	*********************************************************************************/
	
	unsigned char	i,
					j,
					k;
					
	SEL_Num = 0;
	Stptr = 0;
	SCompNum = 0;
	NumVisited = 1;	

	for(i = 0; i < Vertices; i++) Num[i] = 0;
	
	/*********************************************************************************
		 Later, we process the graph G', whose vertices are the strong componets of 
	the Whitehead Graph WG, as modified by the maximal flow from source to sink. 
		The only edges of G' arise from edges of the Whitehead graph WG that are 
	saturated by the maximal flow from source to sink in WG. We save the endpoints of 
	each of these saturated edges in the lists SatEdgeList1[] and SatEdgeList2[]. 
	**********************************************************************************/
			
	for(i = 0; i < Vertices; i++)
	for(k = 0; (j = AJ3[i][k]) < VERTICES ; k++) if(A[i][j] && !A[j][i])
		{
		SatEdgeList1[SEL_Num] = i;
		SatEdgeList2[SEL_Num] = j;
		SEL_Num++;
		}	
	
	/******************************************************************************** 
	  Make repeated calls to SCdfsR() until Num[v] is positive for each vertex v
	of the Whitehead Graph WG. Thus SCdfsR() is called a total of Vertices times,
	once for each vertex of WG. The Stong Component to which v belongs is left in 
	SComp[v], while the total number of strong components found is left in SCompNum.
	*********************************************************************************/
	
	for(i = 0; i < Vertices; i++) if(!Num[i]) SCdfsR(i);	
}				
				
void SCdfsR(unsigned char w)	
{
	/****************************************************************************** 
	  This routine uses Tarjan's Depth-First-Search Strong Component algorithm.
	 (The code is an adaptation of code in Sedgewick's Algorithms in C, Part 5.) 
	******************************************************************************/
			
	unsigned char 	k,
					min,
					v;
	
	Num[w] 			= NumVisited++;
	Low[w] 			= Num[w];
	min 			= Low[w];
	LStack[Stptr++] = w;
	
	for(k = 0; (v = AJ3[w][k]) < VERTICES ; k++) if(A[w][v])
		{
		if(!Num[v]) SCdfsR(v);
		if(Low[v] < min) min = Low[v];
		}
	if(min < Low[w])
		{
		Low[w] = min;
		return;
		}
		
	/******************************************************************************** 
		If min == Low[w], we have located a strong component. The vertices on the
	stack down to w comprise the component. If v is a vertex of this component, we
	remove v from the stack, set SComp[v] = SCompNum, and add the value Vertices 
	to Low[v]. (Adding Vertices to Low[v] forms a sentinel which prevents the 
	value of Low[v] from interfering with the location of later strong components.) 
	*********************************************************************************/	
	
	do
		{
		SComp[(v = LStack[--Stptr])] = SCompNum;
		Low[v] += Vertices;
		}
	while (LStack[Stptr] != w);
	SCompNum++;
	
	return;
}			

void Init_TCR(void) 	
{
	/********************************************************************************  
	 This routine sets up the depth-first-search Transitive Closure routine TCdfsR(). 
	*********************************************************************************/
					
	unsigned char 	i,
					j,
					k;
		
	/*********************************************
	  Clear the arrays Num[], TC[][], and MM[][]. 
	**********************************************/
	
	for(i = 0; i < SCompNum; i++)
		{
		Num[i] = 0;
		for(j = 0; j < SCompNum; j++)
			{
			TC[i][j] = 0;
			MM[i][j] = 0;
			}
		}
		
	/***********************************************************************************
	  Add the edges whose endpoints are SComp[SatEdgeList1[]] and SComp[SatEdgeList2[]]
	to MM[][]. Then make MM[][] pretty by eliminating loops.
	***********************************************************************************/
	
	for(k = 0; k < SEL_Num; k++)	
		{
		i = SComp[SatEdgeList1[k]];
		j = SComp[SatEdgeList2[k]];
		MM[i][j] = 1;
		}
	
	for(i = 0; i < SCompNum; i++) MM[i][i] = 0;
		
	NumVisited = 1;
	
	/**************************************************************************** 
	  Make repeated calls to TCdfsR() until Num[v] is positive for each vertex v
	of ???. Then set each diagonal element of TC[][] equal to 1 so that each 
	vertex of ??? lies in its transititive closure.
	*****************************************************************************/
	
	for(i = 0; i < SCompNum; i++) if(!Num[i]) TCdfsR(i);
	
	for(i = 0; i < SCompNum; i++) TC[i][i] = 1;	
}

void TCdfsR(unsigned char v)	
{				
	/*****************************************************************************
		Let G' be the directed acyclic graph or "DAG", whose vertices are the 
	strong componets of the WHG G as modified by the maximal fow from Sink to 
	Source. The edges of G' are the edges of the WHG G that are saturated with 
	flow. The endpoints of these edges have been saved in the arrays 
	SatEdgeList1[] and SatEdgeList2[].
	 	This routine uses depth-first-search to locate the Transitive Closure 
	of each vertex in G'. Note TCdfsR() is called a total of SCompNum times, 
	once for each vertex of ???.
	  (The code is an adaptation of code in Sedgewick's Algorithms in C, Part 5.)
	******************************************************************************/	
	
	unsigned char 	i,
					u;
					
	Num[v] = NumVisited++;
	for(u = 0; u < SCompNum; u++) if(MM[v][u])
		{
		TC[v][u] = 1;
		if(Num[u] > Num[v]) continue;
		if(!Num[u]) TCdfsR(u);
		for(i = 0; i < SCompNum; i++) if(TC[u][i] == 1)
			TC[v][i] = 1;
		}
	return;					
}

int Get_LTs(unsigned int Source,int F1,int F2)
{
	unsigned char	i,
					j,
					k,
					Stptr,
					TopOfStack;
					
	int				Return_Value;				
	
	for(i = 0; i < SCompNum; i++) 
		{
		SinkSet[i] = 0;
		OnLStack[i] = 0;
		}
		
	Stptr = 0;
	TopOfStack = 0;
	
	while(1)
		{
		/*****************************************************************************
		  Look for the first strong component whose number i is greater than or equal
		  to the current value of TopOfStack which also satisfies:
			1) Component i is not in the SinkSet[].
			2) The transitive closure of i does not contain any vertex in OnLStack[]. 
		******************************************************************************/
		
		for(i = TopOfStack; i < SCompNum && !SinkSet[i]; i++)
			{
			for(j = 0; j < SCompNum; j++) if(TC[i][j] && OnLStack[j]) break;
			if(j == SCompNum) break;
			}
		if(i < SCompNum)
			{
			OnLStack[i] = 1;
			LStack[Stptr++] = i;
			TopOfStack = i + 1;
			for(j = 0; j < SCompNum; j++) SinkSet[j] += TC[i][j];
			
			/*************************************************************************
				Count the current number k of strong components in the SinkSet, and
			check if 1 < k < (SCompnum - 1). (There is a nontrivial Whitehead 
			automorphism, with its acting generator corresponding to the current 
			source vertex, which yields a level T-transformation iff k satisfies  
								1 < k < (SCompnum - 1).)
			*************************************************************************/
			
			for(j = k = 0; j < SCompNum; j++) if(SinkSet[j]) k++;
			if(k > 1 && (k + 1) < SCompNum)
				{ 
				if( (Return_Value = Do_LT(Source, F1, F2)) ) 
				    return(Return_Value);
				}
			}
		else
			{
			if(Stptr == 0) break;
			TopOfStack = LStack[--Stptr];
			OnLStack[TopOfStack] = 0;
			for(j = 0; j < SCompNum; j++) SinkSet[j] -= TC[TopOfStack][j];
			TopOfStack++;
			}	
		}
	return(0);	
}

int Do_LT(unsigned int Source,int F1,int F2)
{
	unsigned char	i,
					j,
					*p,
					*q;
					
	unsigned int	k;				
	
	for(i = 0; i < Vertices; i++)
	for(j = 0; j < SCompNum; j++)
		{
		if(OnLStack[j] && TC[j][SComp[i]])
			{
			ZZ[i] = 1;
			break;
			}
		ZZ[i] = 0;	
		}
		
	/**************************************************************************
		Call Do_Aut() to perform the level transformation specified in ZZ[].
		Increment Num_Level_Transformations.
	**************************************************************************/
	
	if(Micro_Print)
		printf("\n\nPerformed the following level-transformation:");

	if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);

	Num_Level_Transformations ++;
	
	if(Micro_Print)
		{
		printf("\n\nThe presentation is currently:\n");
		Print_Relators(Relators,NumRelators);
		}		
	
	/**************************************************************************
		Call In_File() or On_File() to see whether this new presentation is
		already on file. If not, save the new presentation.
	**************************************************************************/
	
	switch(F2)
		{
		case 0:
			{
			Length = SURL[ReadPres];
			k = On_File();
			if(k == NumFilled)
				{
				if(Dup_On_File < INFINITE)
					 {
					 if(Save_Pres(ReadPres,Dup_On_File,Length,1,75,0,0,0)) return(5);
					 Mark_As_Duplicate(Dup_On_File);                                 
					 }
				 else
					 {
					 if(Save_Pres(ReadPres,0,Length,1,75,0,0,0)) return(5);
					 UDV[NumFilled - 1] = 0;
					 BDY[NumFilled - 1] = BDY[ReadPres];
					 }
				 if(NumFilled >= MAX_SAVED_PRES - 3) return(FULL_HOUSE);
			 
				 /*****************************************************************
					If F1 is TRUE, check if the presentation contains a relator
					of length <= 2.
				 *****************************************************************/
			 
				 if(F1)
					 {
					 for(i = 1; i <= NumRelators; i++)
						 {
						 if(GetHandleSize((char **) Relators[i]) <= 3) return(3);
						 }
					 }
				}
			break;
			}
		case 1:
			{
			Length = SURL[ReadPres];
			k = In_File();
			if(k == TOO_LONG) return(5);
			if((k == NumFilled - 1) && (NumFilled >= MAX_SAVED_PRES - 3)) return(FULL_HOUSE);
			break;
			}
		case 2:
			{
			Length = SLength;
			k = In_File2(FALSE, Relators);
			if(k == TOO_LONG) return(5);
			if((k == NumFilled - 1) && (NumFilled >= MAX_SAVED_PRES - 3)) return(FULL_HOUSE);
			break;
			}
		}

	/**************************************************************************
					Restore the original presentation.
		This is necessary because Do_Aut changes the presentation in
		Relators[] into its image under the level transformation and we do
		not want to compose level transformations at this point. In particular,
		the data in the arrays used above will be invalid unless Relators[] 
		has been restored.			
	**************************************************************************/
	
	if(F2 == 2)
		{
		for(k = 1; k <= NumRelators; k++)
			{
			if(Relators[k] != NULL) DisposeHandle((char **) Relators[k]);
			Relators[k] = (unsigned char **) NewHandle(GetHandleSize((char **) SLP[ReadPres][k]));
			if(Relators[k] == NULL) Mem_Error();
			q = *Relators[k];	
			p = *SLP[ReadPres][k];
			while( (*q++ = *p++) ) ;
			}
		
		}
	else
		{	
		for(k = 1; k <= NumRelators; k++)
			{
			if(Relators[k] != NULL) DisposeHandle((char **) Relators[k]);
			Relators[k] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][k]));
			if(Relators[k] == NULL) Mem_Error();
			q = *Relators[k];	
			p = *SUR[ReadPres][k];
			while( (*q++ = *p++) ) ;
			}
		}						
		
	return(0);
}
