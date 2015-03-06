#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   10 Init_Find_Level_Transformations(int Print)
L   82 Test_LT_For_Pseudo_Min(void)
L  437 In_File(void)
********************************************************************************************/
   
int Init_Find_Level_Transformations(int Print)
{
    /******************************************************************************************
        This routine should be called before Find_Level_Transformations() is called. It
        checks whether the presentation of interest has minimal length and has a connected
        Whitehead graph.
    ******************************************************************************************/
                 
    register unsigned char     *p,
    						   *q;   
    						   
    int                        i;						   
    
    /******************************************************************************************
        Copy the presentation which we want level transformations for into Relators[].
    ******************************************************************************************/
        
    for(i = 1; i <= NumRelators; i++)
        {
        if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
        Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][i]));
        if(Relators[i] == NULL) Mem_Error();
        q = *Relators[i];           
        p = *SUR[ReadPres][i];
        while( (*q++ = *p++) ) ;
        }
    Length = SURL[ReadPres];
        
    /******************************************************************************************
                    Check whether this presentation has minimal length.
        If this presentation does not have minimal length, we could still look 
        for level transformations, but it does not seem very interesting to do so.            
    ******************************************************************************************/
    
    if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG)
        {
        printf("\n\nPresentation %d is too long!",ReadPres + 1);
        return(TOO_LONG);
        }
    if(Automorphisms)
        {
        if(Print)
            {
            printf("\n\nPresentation %d does not have minimal length.",ReadPres + 1);
            printf("\n\nHeegaard will only find level transformations for presentations ");
            printf("that have minimal length.");
            }
        return(1);
        }
                        
    /******************************************************************************************
                    Check whether the Whitehead graph of this presentation is connected.
            It is important to do this because Find_Level_Transformations() will not work
            correctly if the Whitehead graph of the presentation is not connected!!
    ******************************************************************************************/
                                
    for(i = 0; i < Vertices; i++) ZZ[i] = 0;
    if(Connected_AJ3(0,0) == FALSE)
        {
        if(Print)
            {
            printf("\n\nPresentation %d does not have a connected Whitehead graph.",
            ReadPres + 1);
            printf("\n\nHeegaard will only find level transformations for presentations ");
            printf("with connected graphs!");
            }
        return(2);
        }  
    return(0);
}        


void Test_LT_For_Pseudo_Min()
{
    /******************************************************************************************
        This is a rather crude routine which takes the set of presentations in the orbit of a
        presentation under level transformations, and checks them for separating pairs of
        vertices, reducibility of their dual diagrams and pseudo-minimality. It does
        essentially no error checking!!!!
    ******************************************************************************************/
        
    unsigned char   *AA,
    		    	*p,
    		    	*q,
    		    	**Temp;
    
    int             Realizable;
                    
    unsigned int    Flag,
    		    	h,
    		    	i,
		    		j,
		    		k,
		    		Num_Reducible,
		    		Num_Sep_Vert,
		    		Num_Pseudo_Min,
		    		Num_Standard;
                    
    unsigned long   SLength;                
                    
    unsigned int Whitehead_Graph();
    
    DrawingDiagrams = TRUE;
    TestRealizability1 = TRUE;
    Realizable = -1;
    Num_Pseudo_Min = Num_Reducible = Num_Sep_Vert = Num_Standard = 0;
    if(Batch != 7) printf("\n");
    
    AA = (unsigned char*) NewPtr((sizeof(char)*(NumFilled+1)));
	if(AA == NULL) Mem_Error();
        
    for(ReadPres = 0; ReadPres < NumFilled; ReadPres++)
        {     
        /**************************************************************************************
                    Copy the presentation that we want to test into Relators[].
        **************************************************************************************/
        
        NumGenerators = NG[ReadPres];
        NumRelators = NR[ReadPres];
        Vertices = 2*NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
            Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[ReadPres][i]));
            if(Relators[i] == NULL) Mem_Error();
            q = *Relators[i];           
            p = *SUR[ReadPres][i];
            while( (*q++ = *p++) ) ;
            }
        WhichInput = ReadPres;        
        Fill_A(NumRelators);
        Saved_Vertices = 0;
        Flag = Whitehead_Graph();
        if(Flag == 0)
            {
            switch(Realizable)
                {
                case -1:
                    Realizable = TRUE;
                    break;
                case 0:
                    Realizable = 3;
                    break;
                case 1:
                    break;
                case 2:
                    Realizable = 3;
                    break;
                case 3:
                    break;
                }
            AA[ReadPres] = 1;   
            if(NumGenerators == NumRelators)
                {
                SLength = SURL[ReadPres];
                for(i = 1; i <= NumRelators; i++)
                    { 
                    Temp                = Relators[i];
                    Relators[i]         = DualRelators[i];
                    DualRelators[i]     = Temp;    
                    }
                if(Freely_Reduce() == TOO_LONG) continue;
                Length = OrigLength;
                Find_Flow_A(NORMAL,FALSE);
                if(Length == SLength)
                    {
                    Num_Pseudo_Min ++;
                    AA[ReadPres] = 2;
                    }
                else
                    {
		    		Num_Reducible ++;
                    AA[ReadPres] = 3;
                    if(Length == NumRelators)
                        {
                        Num_Standard ++;
                        AA[ReadPres] = 4;
                        }
                    }
                }
            }    
        else
            {
            AA[ReadPres] = 1;
            switch(Flag)
                {
                case NON_PLANAR:
                case FATAL_ERROR:
                    switch(Realizable)
                        {
                        case -1:
                            Realizable = FALSE;
                            break;
                        case 0:
                            break;
                        case 1:
                            Realizable = 3;
                            break;
                        case 2:
                            Realizable = 3;
                            break;
                        case 3:
                            break;
                        }
                    break;
                case SEP_PAIRS:
                    AA[ReadPres] = 0;
                    UDV[ReadPres] = SEP_PAIRS;
                    if(V1 & 1)
                        LSP[ReadPres] = V1/2 + 97;
                    else
                        LSP[ReadPres] = V1/2 + 65;
                    if(V2 & 1)
                        LSQ[ReadPres] = V2/2 + 97;
                    else
                        LSQ[ReadPres] = V2/2 + 65;
		    		Num_Sep_Vert ++;
                    break;
                case NON_UNIQUE_1:
                case NON_UNIQUE_2:
                case NON_UNIQUE_3:
                case NON_UNIQUE_4:
                case V2_ANNULUS_EXISTS:
                    if(Realizable < 0)
                        {
                        TestRealizability2 = TRUE;
                        if(Whitehead_Graph() == NO_ERROR) Realizable = 2;
                        TestRealizability2 = FALSE;
                        }
                    break;
                default:
                    switch(UDV[ReadPres])
                        {
                        case NON_UNIQUE_1:
                        case NON_UNIQUE_2:
                        case NON_UNIQUE_3:
                        case NON_UNIQUE_4:
                            if(Realizable < 0)
                                {
                                TestRealizability2 = TRUE;
                                if(Whitehead_Graph() == NO_ERROR) Realizable = 2;
                                TestRealizability2 = FALSE;
                                }
                            break;
                        }
                    break;                
                }
            }
        continue;            
        }

    DrawingDiagrams = FALSE;
    TestRealizability1 = FALSE;
     
    if(Batch == 7)
    	{
    	printf("|WOSepVert| %u.",NumFilled - Num_Sep_Vert);
    	if(H_Results != NULL) fprintf(H_Results,"|WOSepVert| %u.",NumFilled - Num_Sep_Vert);
    	DisposePtr((unsigned char*) AA);
    	return;
    	}
          
    if(Num_Sep_Vert < NumFilled)
        {
        j = k = 0;
        j += printf("\n\nPresentations with no separating pairs of vertices:    ");
        for(h = 0; h < NumFilled; h++) if(AA[h])
            {
            if(++k < NumFilled - Num_Sep_Vert)
                j += printf("{%3d,",h+1);
            else
                j += printf("{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i])
            {
            if(++k < NumFilled - Num_Sep_Vert)
                j += printf("%3d,",i+1);
            else
                j += printf("%3d}.",i+1);    
            if(j > 83)
                {
                j = 0;
                printf("\n");
                }
            }                   
        }    
    
    if(Num_Pseudo_Min)
        {
        j = k = 0;
        j += printf("\n\nPresentations which are pseudo-minimal:    ");
        for(h = 0; h < NumFilled; h++) if(AA[h] == 2)
            {
            if(++k < Num_Pseudo_Min)
                j += printf("{%3d,",h+1);
            else
                j += printf("{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i] == 2)
            {
            if(++k < Num_Pseudo_Min)
                j += printf("%3d,",i+1);
            else
                j += printf("%3d}.",i+1);    
            if(j > 83)
                {
                j = 0;
                printf("\n");
                }
            }          
        }
	
    if(Num_Reducible)
        {
        j = k = 0;
        j += printf("\n\nPresentations which are reducible:     ");
        for(h = 0; h < NumFilled; h++) if(AA[h] >= 3)
            {
            if(++k < Num_Reducible)
                j += printf("{%3d,",h+1);
            else
                j += printf("{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i] >= 3)
            {
            if(++k < Num_Reducible)
                j += printf("%3d,",i+1);
            else
                j += printf("%3d}.",i+1);
            if(j > 83)
                {
                j = 0;
                printf("\n");
                }
            }
        if(Num_Standard)
            {
            j = k = 0;
            j += printf("\n\nPresentations whose duals are Standard:     ");
            for(h = 0; h < NumFilled; h++) if(AA[h] == 4)
                {
                if(++k < Num_Standard)
                    j += printf("{%3d,",h+1);
                else
                    j += printf("{%3d}.",h+1);    
                break;
                }
            for(i = h + 1; i < NumFilled; i++) if(AA[i] == 4)
                {
                if(++k < Num_Standard)
                    j += printf("%3d,",i+1);
                else
                    j += printf("%3d}.",i+1);
                if(j > 83)
                    {
                    j = 0;
                    printf("\n");
                    }
                }
            }
        }
    
    printf("\n\nTotal presentations with no separating pairs of vertices:  %u.",NumFilled - Num_Sep_Vert);
    printf("\n\nTotal pseudo-minimal presentations:  %u.",Num_Pseudo_Min);
    printf("  (Ignore this number if the manifold is not closed.)");
    printf("\n\nTotal reducible presentations:  %u.",Num_Reducible);
    printf("  (Ignore this number if the manifold is not closed.)");
    
    switch(Realizable)
        {
        case -1:
            printf("\n\nUnable to determine if any presentations in the orbit are realizable.");
            break;
        case 0:
            printf("\n\nNone of the presentations in the orbit are realizable.");
            break;
        case 1:
            printf("\n\nAll presentations in the orbit are realizable.");
            break;
        case 2:
            printf("\n\nUnable to determine if any presentations in the orbit are realizable.");
            break;
        case 3:
            printf("\n\nSome kind of error occured.");
            printf("\n\nHeegaard found both realizable and unrealizable presentations in the orbit!!!??");
            break;            
        }
         
		printf("\n\nScroll back if necessary to see these lists.");    
	
		printf("\n\nPLEASE SELECT ONE OF THE FOLLOWING ALTERNATIVES."); 
		printf("\n    0) Show no presentations in the orbit.");
		printf("\n    1) Show only presentations in the orbit which have no separating pairs of vertices.");
		printf("\n    2) Show only presentations in the orbit which are pseudo-minimal."); 
		printf("\n    3) Show only presentations in the orbit which are reducible.");
		printf("\n    4) Show all presentations in the orbit.");
		printf("\nHIT 0,1,2,3,OR 4.        ");
		GET_RESPONSE4:    
		switch(WaitkbHit())
			{
			case '0':
				break;
			case '1':
				if(Num_Sep_Vert < NumFilled)
					Report(0,0,0,0,0,1,0,1,1,AA);
				break;
			case '2':
				if(Num_Pseudo_Min)
					Report(0,0,0,0,0,1,0,2,1,AA);
				break;
			case '3':
				if(Num_Reducible)
					Report(0,0,0,0,0,1,0,3,1,AA);
				break;
			case '4':
				Report(0,0,0,0,0,0,0,0,1,NULL);
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE4;                    
			}        
    DisposePtr((unsigned char*) AA);              
}

unsigned int In_File(void)
{    
    register unsigned char  *p,
    						*q;
    
    unsigned char           *r;
                            
    int                     i,
    						Result;
    
    unsigned int            Node;
    
    Size                    HSP,
    						HSQ;     
    
    Canonical_Rewrite(Relators,FALSE,FALSE);
    Node = 0;
    while(1)
        {
        for(i = 1,Result = 0; i <= NumRelators; i++)
            {
            HSP = LR[i] + 1;
            HSQ = GetHandleSize((char **) SUR[Node][i]);
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
            r = *Relators[i] + LR[i];
            *r = 125;
            for(p = *Relators[i],q = *SUR[Node][i];    *p == *q; p++,q++) ;
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
                    if(Compute_Stabilizers)
                         printf("  %d -> %d",ReadPres + 1,NumFilled + 1);
                    if(Save_Pres(ReadPres,0,Length,1,75,0,3,0)) return(TOO_LONG);
                    UDV[NumFilled - 1] = 0;
                    BDY[NumFilled - 1] = BDY[ReadPres];
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
                SURNumX[Node] ++;     
                return(Node);
            case -1:
                if(Right[Node] == INFINITE)
                    {
                    if(Compute_Stabilizers)
                         printf("  %d -> %d",ReadPres + 1,NumFilled + 1);
                    if(Save_Pres(ReadPres,0,Length,1,75,0,3,0)) return(TOO_LONG);
                    UDV[NumFilled - 1] = 0;
                    BDY[NumFilled - 1] = BDY[ReadPres];
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
