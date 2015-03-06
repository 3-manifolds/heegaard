
#include "Heegaard.h"
#include <ctype.h>
#include <string.h>

#ifndef MAC
#include <termios.h>
#include <stdio.h>
struct termios normal_termio;
#endif

/****************************** function prototypes *****************************************
L  334 main(int argv, char **argc)
L  923 Rewrite_Input(void)
L 1058 Delete_Dups(void)
L 1141 Compare_Input_Pres(void)
L 1169 Compare_Pres(int k)
L 1197 Compare_Dual_Pres(int k)
L 1224 Get_Initial_Diagram(PrintFlag)
L 1993 Save_Pres(unsigned int From,unsigned int Daut,unsigned long Len,int F1,int F2,int F3,unsigned char F4,char F5)		
L 2221 Non_Unique_Initial_Diagram(void)
L 2347 On_File(void)
L 2450 Init_G_Variables(void)
L 2529 Delete_Old_Presentations(void)
L 2549 ReRun_A_Presentation(void)
L 2798 Get_Presentation_From_File(void)
L 2929 Get_Presentation_From_KeyBoard(void)
L 2996Check_Realizability_Of_The_Initial_Presentation(void)
L 3022 Display_Diagram_Of_The_Initial_Presentation(void)
L 3117 Get_Simplification_Parameters_From_User(int Flag1,int Flag2)
L 3278 Turn_Micro_Print_On(void)
L 3320 Find_Level_Transformations_Of_The_Initial_Presentation(void)
L 3496 Reduce_The_Initial_Presentation_To_Minimal_Length(void)
L 3563 Initial_Realizability_Check(void)
L 3593 Print_Realizability(int Del_Only_Triv_Rel, unsigned int WhichPres)
L 3621 Realization_Warning(void)
L 3628 char mykbhit(void)
L 3661 char WaitkbHit(void)
********************************************************************************************/

FILE 
    *H_Results,
    *Gvizdata,
    *input_relators;

char
	Check_If_HS_Rep,
	CheckHSReps,
    *ER,
    FoundSF,
    FoundFiniteSF;
    
unsigned char
	Batch = 0,
	BPrintAnnulusData,
	BPrintNotConnectedData,
	B5PrintBdryComps,
	B5PrintDualRelators,
	B5PrintPaths,
	B5TestSimpleCircuits,
	B10B11ConSumPres,
	B10B11Finite,
	B10B11HSReps,
	B10B11Recognized,	
	B14B15PrintPres,
	*CBC[MAXNUMCOMPONENTS],
	*BCF,
	*BCWG,
	BreadthFirstSearch,
	*CO[VERTICES],
	**Copy_Of_Input[MAXNUMRELATORS + 1],
	**Copy_Of_Rel_1[MAXNUMRELATORS + 1],
	**Copy_Of_Rel_2[MAXNUMRELATORS + 1],			
	*CS,
	**CD_Surgery_Rel[MAXNUMRELATORS + 1],
	*DD,
	*DeletedEdgePtr,
	*DeletedEdges,
	**DelRelators[MAXNUMRELATORS + 1],
	DepthFirstSearch,
	**DualRelators[MAXNUMRELATORS + 1],
	**Exp_Surgery_Rel[MAXNUMRELATORS + 1],	
	*Face[2*VERTICES],
	*FV,
	*GBC,
	GoingDown,
	GoingUp,
	GoingUpDown,
	*Inst,
	**KorLRelators[MAXNUMRELATORS + 1],
	*MM[VERTICES],
	*N1H,
	*NCS,
	*NS1XD2,
	*NS1XS2,	
	**OutRelators[MAXNUMRELATORS + 1],
	*PresName,
	*QPM,			
	**Relators[MAXNUMRELATORS + 1],
	**TopOfChain[MAXNUMRELATORS + 1],
	*RWR[MAXNUMDUPS],
	***SLR[MAX_SAVED_LEVELS],
	***SMGP[MAX_MIN_GEN_PRES],
	****SUR,
	*T[(VERTICES)/2],
	**Temp1,
	**Temp2,
	**Temp3,
	**Temp4,
	**Temp5,
	**Temp6,
	**Temp7,
	**Temp8,
	**Temp9,
	**Temp10,
	**Temp11,
	**Temp12,
	**Temp13,
	**Temp14,
	**Temp15,
	**Temp16,
	*TP,
	*Low,
	*MM[VERTICES],
	*Num, 
	*OnLStack,
	*SatEdgeList1,
	*SatEdgeList2,
	*SComp,
	*SepVertexList,
	*SinkSet,
	*LStack,
	*TC[VERTICES],
	TestRealizability1,
    TestRealizability2,
    TestRealizability3,
    TestRealizability4,
    **WirtingerL[MAXNUMRELATORS + 1],
    **WirtingerM[MAXNUMRELATORS + 1];
    
int 
    *BDY,
    BdryData,
    Betti_Number,
    Boundary,
    Compute_Stabilizers,
    Connected,
    CopyNumGenerators,
    CopyNumRelators,
    Count,
    CurrentComp,
    Delete_Only_Short_Primitives,
    Did_Cutting_Disk_Surgery,    
    Did_Exponent_Surgery,
    Do_Not_Reduce_Genus,
    DrawingDiagrams,
    EmtyRel,
    Find_All_Min_Pres,
    *Flags,
    FormBandsums,
    FoundPower,
    FoundPrimitive,
    *GB[VERTICES],
    HS_Rep_NumGens,
    *InPS,
    Input,
    Knot_Or_Link,
    Level_Interrupt,
    MajorVert,
    Micro_Print,
    Micro_Print_F,
    Modified_Init_Pres,
    *NG,
    NG_TOC,
    NGV2,
    NonPlanar,    
    NoReport,
    NotConnectedError,
    *NR,
    NR_TOC,
    NumBdryComps,
    NumComps,
    NumCuttingDiskSurgeryRel,
    NumDelRelators,
    NumEdges,
    NumEmptyRels,
    NumExp_Sur_Rel,    
    NumFaces,
    NumGenerators,
    NumKnot_Or_Link_Rel,
    *NumLoops,   
    NumRelators,
    NumSepComps,
    NumTimes,
    Num_Torsion_Values,
    OnlyReducingBandsums,
    *PRIM,
    RandomizeSlides,
    ReadPres,
    *SaveBdry,
    Save_Init_Pres,
    SaveMinima,
    Saved_Vertices,
    SepPairs,
    **SFSols,
    SFSolV[20],    
    SRError,
    SReadPres,
    SSSReadPres,
    *SURNumX,    
    TotalComp,
    UserSaidQuit,
    Vertices,
    WhichInput,
    *X,
    *Y;    
                             
unsigned int
    *A[VERTICES],
    *AA[VERTICES],
    *AJ1[VERTICES],
    *AJ2[VERTICES],
    *AJ3[VERTICES],
    *AT,
    *B[VERTICES],
    *Bdry,
    *BeenChecked,
    *BSV1,
    *BSV2,
    *ComponentNum,
    CouldNotRemove,
    *Daughters,
    *DF,
    *DRA[2*MAXNUMRELATORS],
    Dup_On_File,
    *ED[(VERTICES)/2],
    *EXL[(VERTICES)/2],
    *EXP[(VERTICES)/2],
    *EXR[(VERTICES)/2],
    *Father,
    *FR,
    From_BANDSUM,
    From_DUALIZE,
    *GV2,
    *GV2L,
    *GV2R,
    *HegSplNum,
    *HegSplNxt,
    *InDisk,
    *InQueue,
    *IV,
    *Left,
    LensSpace,    
    *Lowpt,
    MaxLength,
    Mergers,
    MyMaxSavedPres,
    *NEBC,    
    *NEX[(VERTICES)/2],
    *NFBC,    
    NotNewPres,
    *NRBC,    
    *Number,
    NumFilled,
    NumRealizable,
    NumRelTooLong,   
    Num_Saved_LPres,
    *Num_Sep_Vertex_Pairs,
    NumNonSepAnnuli,
    NumNotConnected,
    NumPresExamined,
    NumSepAnnuli,
    NumSymmetries,
    NumVert,
    OnStack,
    *OSA,
    *OSB,
    *PG,
    *Right,
    Start_Level_Search,
    Starting_Pres,
    Stopper,
    *SUR_Num,
    *SV,
    This_Pres,
    *TV,
    *UDV,
    *UnUsed_Sep_Vertex_Pairs,
    *UpDate,
    *V,
    *VA,
    *VWG,
    V1,
    V2,
    Word1,
    Word2,
    *XX,
    *YY,     
    *ZZ,
    *zz;

long
    Band_Sums,
    MaxTotalLength,
    Minimum,
    NumDiagrams,
    Recip_P,
    Recip_Q;

unsigned long 
    Automorphisms,
    BytesAvailable,
    BytesUsed,
    HS_Rep_Length,
    Length,
    *LR,
    *LSP,
    *LSQ,
    MinLength,
    *MLC[MAXNUMCOMPONENTS],
    NumDualized,
    NumErrors,
    Num_Level_Slides,
    Num_Level_Transformations,
    OrigLength,
    P,
    Q,
    SLength,
    SNum_Level_Slides,
    *SURL,
    TOCLength,
    TotalAuts;

int main(int argv, char **argc)
{
    char            FoundAlexPoly,
        			FoundMeridianReps;
                
    unsigned char 	*p,
                   	*q;
                            
    unsigned int    i,
    				j;
    
    long            Scratch;

#ifndef MAC
    struct termios newtermio;

    tcgetattr(0, &normal_termio);
    newtermio = normal_termio;
    cfmakeraw(&newtermio);
    newtermio.c_iflag = ~(IGNBRK|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|IXON);
    newtermio.c_oflag = OPOST|ONLCR;
    newtermio.c_lflag = ~(FLUSHO|ECHONL|ICANON|ISIG|IEXTEN);
    tcsetattr(0, 0, &newtermio);
#endif
    
    if( (i = Do_Initialization()) )
        {
        printf("\n\nHeegaard was unable to allocate memory for all of its data structures.");
        printf("\nTry increasing the amount of memory which Heegaard can use.");
        printf("\nFailure at allocation number %u.",i);
        
#ifndef MAC
	tcsetattr(0, 0, &normal_termio);
#endif
        return(0);
        }
    
    Input = INITIAL_PRES;
 
    printf("\n\n                                 HEEGAARD");
    printf("\n                               BY JOHN BERGE");
    printf("\n                             jberge@charter.net");
    printf("\n                                   2/5/15\n");
    printf("\n A PROGRAM FOR STUDYING 3-MANIFOLDS VIA PRESENTATIONS AND HEEGAARD DIAGRAMS.\n");
    printf("\n        Copyright 1995-2015 by John Berge, released under GNU GPLv2+.");
 
    if((Gvizdata = fopen("Heegaard_Diagrams.dot","w+")) == NULL)
 		printf("\n\nUnable to create file Gvizdata used by Graphviz() to display Heegaard diagrams");						
 	else
 		{
 		Heegaard_Splash_Screen();
		printf("\n\nNote! Open the file 'Heegaard_Diagrams.dot' in Graphviz() to see Heegaard's Heegaard diagrams.");
		}
		
     printf("\n\n	Note!	Hit the 'space-bar' to interrupt long computations.");
     printf("\n            Hitting 's' during long computations should provide a status report.");  
     
_BEGIN:
    if(Input == BEGIN)
        {
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;
        printf("\n");
        #ifdef DEBUGGING
            printf("\nHIT 'b' TO PRINT DEBUGGING INFO INTO THE FILE SIMPLIFY_H_ERRORS.");
        #endif
        printf("\nHIT 'd' TO SEE DATA FOR THE HEEGAARD DIAGRAMS.");
        printf("\nHIT 'n' TO ENTER A NEW PRESENTATION.");
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
        if(Knot_Or_Link || Did_Exponent_Surgery || Did_Cutting_Disk_Surgery)
            printf("\nHIT 'r' TO TRY ANOTHER SURGERY ON A PREVIOUS PRESENTATION, OR TO RERUN A PRESENTATION.");        
        else
            printf("\nHIT 'r' TO RERUN A PRESENTATION.");
        printf("\nHIT 's' TO LOOK FOR SYMMETRIES.");
        printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS NOW IN MEMORY.");
        printf("\nHIT 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY BY SUMMAND NUMBER,");
        printf("\n	NUMGENERATORS, NUMRELATORS, LENGTH AND 'LEXICOGRAPHIC' ORDER.");
        printf("\n	(NOTE THAT AFTER SORTING ONE CAN FIND CANONICAL ORBIT REPS, ETC.)");

		printf("\n");
        GET_RESPONSE1:    
        switch(WaitkbHit())
            {
            case 'b':
                #ifdef DEBUGGING
                    DrawingDiagrams = FALSE;
                    Debug();
                    Input = BEGIN;
                    goto _BEGIN;
                #else
                    if(Batch == FALSE) SysBeep(5);
                    DrawingDiagrams = FALSE;
                    Input = BEGIN;
                    goto GET_RESPONSE1;                    
                #endif
                
            case 'd':
                Display_Diagrams();
                Input = BEGIN;
                goto _BEGIN;
                   
            case 'n': 
                Input = INITIAL_PRES;
                DrawingDiagrams = FALSE;
                Did_Exponent_Surgery = FALSE;
                Did_Cutting_Disk_Surgery = FALSE;
                Delete_Old_Presentations();
                break;
                            
            case 'q':
                goto _QUIT;
                            
            case 'r':
                if(ReRun_A_Presentation())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }        
                break;
                
            case 's':
                Find_Symmetries(FALSE);
                Input = BEGIN;
                goto _BEGIN;
                
            case 'v':
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,NULL);
                Input = BEGIN;
                goto _BEGIN;
                
            case 'w':
                printf("\n\n     Sorting presentations. . . .");
                Sort_Presentations_In_Memory(1);
                Input = BEGIN;
                goto _BEGIN;
                                            
            default:
                if(Batch == FALSE) SysBeep(5);
                DrawingDiagrams = FALSE;
                Input = BEGIN;
                goto GET_RESPONSE1;                                
            }            
        } 
           
    if(Input == INITIAL_PRES)
        {
        PRINT_INPUT_OPTIONS: 
        printf("\n\nHIT 'B' TO DO SOME BATCH PROCESSING.");           
        printf("\nHIT 'f' IF THE PRESENTATION WILL COME FROM THE FILE 'Input_Presentations'.");
        printf("\nHIT 'k' IF THE PRESENTATION WILL BE ENTERED FROM THE KEYBOARD.");
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
		printf("\n");
        GET_RESPONSE2:        
        switch(WaitkbHit())
            {
            case 'B': BatchProcessing();
            	DrawingDiagrams = FALSE;
                Did_Exponent_Surgery = FALSE;
                Did_Cutting_Disk_Surgery = FALSE;
            	Batch = 0;
            	Input = INITIAL_PRES;
            	goto PRINT_INPUT_OPTIONS;
            case 'f':
                if(Get_Presentation_From_File()) goto PRINT_INPUT_OPTIONS;
                break;     
            case 'k':
                if(Get_Presentation_From_KeyBoard()) goto PRINT_INPUT_OPTIONS;
                break;
            case 'q':
            	printf("\n\n");
                goto _QUIT;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE2;
            }
            
        /**************************************************************************************
            Call Wirtinger() and check whether this presentation seems to meet enough criteria
            to be the Wirtinger presentation of a knot or link and, if so, give the user the
            option of performing Dehn-fillings.
        **************************************************************************************/
        
        Knot_Or_Link = FALSE;        
        if( (j = Wirtinger()) )
        	{
        	if(j == 6) 
        		{
        		printf("\n\n	This presentation resembles a Wirtinger presentation of a knot or link.");
        		printf("\n\n	However, it does not seem to have the right number of relators!!");
        		printf("\n\n(Heegaard expects a Wirtinger presentation obtained from an n crossing knot ");
        		printf("or link projection to have n relators---with one relator corresponding to ");
        		printf("each crossing of the projection.)");
        		}
        	}
        if(Knot_Or_Link && (j != 17))
            {
            printf("\n\nThe surgered presentation is:\n");
            Print_Relators(Relators,NumRelators);            
            }                    
        }
    if(Input == RERUN) Input = INITIAL_PRES;        
    if(Input == INITIAL_PRES)
		_STABILIZE_RETURN:
		_EXPONENT_SURGERY_RETURN:
		_CUTTING_DISK_SURGERY_RETURN:    
        {
        /**************************************************************************************
            Echo the initial relators to the output so we will have a copy of them. Then call
            Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
            which serves as the initial presentation for Heegaard. 
        **************************************************************************************/    
        
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            Scratch += GetHandleSize((char **) Relators[i]);
        Scratch -= NumRelators;
        printf("\n\nThis presentation has length %ld ",Scratch);        									
        if(Freely_Reduce() == TOO_LONG)
            {
            printf("\n\nThis presentation is too long!!");
            if(Batch == FALSE) SysBeep(5);
            Input = BEGIN;
            goto _BEGIN;
            }
        if(Scratch > OrigLength)
            {
            printf("\nand freely reduces to the following presentation of length %lu.\n",
                OrigLength);
            Print_Relators(Relators,NumRelators);
            Scratch = OrigLength;
            }
        else
            printf("and is freely reduced.");
        Micro_Print = TRUE;
        Micro_Print_F = TRUE;    
        if(Rewrite_Input())
            {
            printf("\n\nThere must be at least one non-empty relator!!");
            if(Batch == FALSE) SysBeep(5);
            Micro_Print = FALSE;
            Micro_Print_F = FALSE;
            Input = BEGIN;
            goto _BEGIN;
            }
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;    

        /**************************************************************************************
            Save a copy of the initial set of relators in case we want to refer to them later.
        **************************************************************************************/
        
        CopyNumRelators     = NumRelators;
        CopyNumGenerators   = NumGenerators;
        HS_Rep_Length 		= OrigLength;
        HS_Rep_NumGens		= NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
            Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
            if(Copy_Of_Input[i] == NULL) Mem_Error();
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while( (*p++ = *q++) ) ;
            }
                                        
        /**************************************************************************************
                Call Init_G_Variables() to initialize some global variables. Then call
                Canonical_Rewrite() to rewrite the presentation in canonical form.
        **************************************************************************************/
        
        Init_G_Variables();
        Length = Scratch;
        printf("\n\nGen %d, Rel %d\n",NumGenerators,NumRelators);   
        Micro_Print = TRUE;
        Micro_Print_F = TRUE;
          					
        if(Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG)
            {
            printf("\n\nThis presentation has too many symmetries!!");
            if(Batch == FALSE) SysBeep(5);
            Micro_Print = FALSE;
            Micro_Print_F = FALSE;
            }
        Micro_Print = FALSE;
        Micro_Print_F = FALSE; 
        					                         
        if(Compare_Input_Pres() == FALSE)
            {
            if(Batch == 0) printf("\n\n");
            printf(" The rewritten initial presentation is:");
            if(Batch == 0) printf("\n");
            Print_Relators(Relators,NumRelators);
            printf("\n");
            }
        					
        /**************************************************************************************
                        Update the copy of the initial set of relators.
        **************************************************************************************/
        
        CopyNumRelators       = NumRelators;
        CopyNumGenerators     = NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
            Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
	        if(Copy_Of_Input[i] == NULL) Mem_Error();    
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while( (*p++ = *q++) ) ;
            }
    
        /**************************************************************************************
                                Present the user with some options.
        **************************************************************************************/
        
        FoundAlexPoly = FALSE;
        FoundMeridianReps = FALSE;                               
    _OPTIONS:
        printf("\n");
        if(NumGenerators == 2 && NumRelators == 1)
        	{
        	if(FoundAlexPoly == FALSE)
        		printf("\nHit 'a' TO COMPUTE THE ALEXANDER POLYNOMIAL OF R1.");
        	if(FoundMeridianReps == FALSE)
        		{	
        		printf("\nHIT 'b' TO FIND 'meridian' REPS M1 & M2 OF H[R]. THIS ALLOWS ONE TO CHECK IF R IS A KNOT RELATOR.");
        		printf("\nN.B.'Meridian' reps depend only on the relator R and don't depend on embedding in S^3.");
        		}
        	}
        printf("\nHIT 'c' TO CHECK REALIZABILITY OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'C' TO CHECK IF THE INITIAL PRESENTATION IS A \042HS REP\042. (Heegaard will stop and alert the user if");  
        printf("\n    a sequence of handle-slides of the initial presentation P yields a presentation P' with |P'| < |P|.)");
        printf("\nHIT 'd' TO SEE DATA FOR THE HEEGAARD DIAGRAM OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'e' TO TRY EXPONENT SURGERY ON THE INITIAL PRESENTATION.");
        printf("\nHIT 'h' TO FIND THE INTEGRAL FIRST HOMOLOGY OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'k' TO TRY CUTTING DISK SURGERY ON THE INITIAL PRESENTATION.");
        printf("\nHIT 'l' TO FIND LEVEL TRANSFORMATIONS OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'm' TO LOOK FOR ALL MINIMAL PRESENTATIONS USING DEPTH-FIRST SEARCH.");
        printf("\nHIT 'M' TO LOOK FOR ALL MINIMAL PRESENTATIONS USING BREADTH-FIRST SEARCH.");       
        printf("\nHIT 'n' TO ENTER A NEW PRESENTATION.");
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
        printf("\nHIT 'r' TO REDUCE AND SIMPLIFY THE INITIAL PRESENTATION USING DEPTH-FIRST SEARCH.");
        printf("\nHIT 'R' TO REDUCE AND SIMPLIFY THE INITIAL PRESENTATION USING BREADTH-FIRST SEARCH.");        
        printf("\nHIT 's' TO FIND SYMMETRIES OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'u' TO STABILIZE THE PRESENTATION WHILE PRESERVING REALIZABILITY.");
        printf("\nHIT 'v' TO REVIEW THE INITIAL PRESENTATION.");
        printf("\nHIT 'x' TO SIMPLIFY A PRESENTATION BY SUCCESSIVELY DELETING PRIMITIVES WITHOUT CHECKING REALIZABILITY.");
        printf("\nHIT 'X' TO FIND THE PRESENTATIONS OBTAINED BY DELETING PRIMITIVES FROM JUST THE INITIAL PRESENTATION.");        
        printf("\nHIT 'z' TO REDUCE THE INITIAL PRESENTATION TO MINIMAL LENGTH.\n");
        GET_RESPONSE3:        
        switch(WaitkbHit())
            {
            case 'a':
            	if(FoundAlexPoly == TRUE) 
            		{
            		printf("\n\nBeen Done! Choose a different option from the previous list.");
            		goto GET_RESPONSE3;
            		}
            	if(NumGenerators != 2 || NumRelators != 1) goto GET_RESPONSE3;
            	FoundAlexPoly = TRUE;
            	AlexanderPolynomial(*Relators[1]);
				Length = GetHandleSize((char **) Relators[1]) - 1;
				Fill_A(1);
				Get_Matrix();            	
            	printf("\n\nCHOOSE ANOTHER OPTION FROM THE PREVIOUS LIST.");
            	goto GET_RESPONSE3;
            
            case 'b':
            	if(FoundMeridianReps == TRUE)
            		{
            		printf("\n\nBeen Done! Choose a different option from the previous list.");
            		goto GET_RESPONSE3;
            		}
            	if(NumGenerators != 2 || NumRelators != 1) goto GET_RESPONSE3;
            	FoundMeridianReps = TRUE;
            	if(Is_Knot_Relator())
            		{
            		printf("\n\nCHOOSE ANOTHER OPTION FROM THE PREVIOUS LIST.");
            		Length = GetHandleSize((char **) Relators[1]) - 1;
            		Fill_A(1);
            		Get_Matrix();
            		goto GET_RESPONSE3;
            		}
            	printf("\n\n	Note: If simplifying < A,B | M1,M2 > yields S^3, R 1 is a knot relator.");
            	printf("\n  If Heegaard produces a pseudo-minimal PM diagram != S^3, R 1 is not a knot relator."); 
       			Input = RERUN;
       			goto _BEGIN;
            	
            case 'c':
                if(Check_Realizability_Of_The_Initial_Presentation())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                break;
 
			case 'C':
				Check_If_HS_Rep 	= TRUE;
				Get_Simplification_Parameters_From_User(FALSE,TRUE);
                DepthFirstSearch	= FALSE;
                BreadthFirstSearch	= TRUE;
                Find_All_Min_Pres 	= FALSE; 
                CheckHSReps 		= TRUE;               
                goto _GET_INITIAL_DIAGRAM;
                               
            case 'd':
                Display_Diagram_Of_The_Initial_Presentation();
                break;
                
            case 'e':
                if(Try_Exponent_Surgery()) 
                	{
                	printf("\n\n	Did not do Exponent Surgery. Choose another option from the list above.");
                	goto GET_RESPONSE3;
                	}
                Did_Exponent_Surgery = TRUE;
                Input = INITIAL_PRES;
                goto _EXPONENT_SURGERY_RETURN;
                    
            case 'h':
                printf("\n\nComputing the integral first homology of the initial presentation . . .");
                Compute_Homology();
                printf("\n\nChoose another option from the list above.");
                goto GET_RESPONSE3;
            
            case 'k':
                if(Try_Cutting_Disk_Surgery())
                	{
                	printf("\n\n	Did not do Cutting_Disk_Surgery. Choose another option from the list above.");
                	goto GET_RESPONSE3;
                	}
                Input = INITIAL_PRES;
                goto _CUTTING_DISK_SURGERY_RETURN;
                
            case 'I':
            case 'l':
                if(Find_Level_Transformations_Of_The_Initial_Presentation())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                break;
                            
            case 'm':
                Get_Simplification_Parameters_From_User(FALSE,TRUE);
                DepthFirstSearch	= TRUE;
                BreadthFirstSearch	= FALSE;
                Find_All_Min_Pres 	= TRUE;
                goto _GET_INITIAL_DIAGRAM; 
                  
            case 'M':
                Get_Simplification_Parameters_From_User(FALSE,TRUE);
                DepthFirstSearch	= FALSE;
                BreadthFirstSearch	= TRUE;
                Find_All_Min_Pres 	= TRUE;
                goto _GET_INITIAL_DIAGRAM; 
                             
            case 'n':
                Input = INITIAL_PRES;
                goto _BEGIN;
                            
            case 'q':
                goto _QUIT;    
            
            case 'r':
                Get_Simplification_Parameters_From_User(FALSE,TRUE);
                DepthFirstSearch	= TRUE;
                BreadthFirstSearch	= FALSE;
                Find_All_Min_Pres 	= FALSE;                
                goto _GET_INITIAL_DIAGRAM;
                
             case 'R':
                Get_Simplification_Parameters_From_User(FALSE,TRUE);
                DepthFirstSearch	= FALSE;
                BreadthFirstSearch	= TRUE;
                Find_All_Min_Pres 	= FALSE;
                goto _GET_INITIAL_DIAGRAM;          
            
            case 's':
                printf("\n\n    This is the initial presentation:\n");
                Print_Relators(Relators,NumRelators);
                i = NoReport;
                NoReport = FALSE;
                Find_Symmetries(TRUE);
                NoReport = i;
                break;
                
            case 'u':
            	if(Stabilize())
            		{
            		printf("\n\n	Did not stabilize. Choose another option from the list above.");
            		goto GET_RESPONSE3;
            		}
            	Input = INITIAL_PRES;
            	goto _STABILIZE_RETURN;    
            
            case 'v':
                printf("\n\n    This is the initial presentation:\n");
                Print_Relators(Relators,NumRelators);    
                goto _OPTIONS;
                
            case 'x':
            	Just_Delete_Primitives(FALSE); 
                Input = BEGIN;
                goto _BEGIN;
                
            case 'X':
            	Just_Delete_Primitives(TRUE); 
                Input = BEGIN;
                goto _BEGIN; 
                                   
            case 'z':
                if(Reduce_The_Initial_Presentation_To_Minimal_Length())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                break;
                
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE3;
            }
        
        NumRelators		= CopyNumRelators;
        NumGenerators 	= CopyNumGenerators;
        Vertices 		= 2*NumGenerators;
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            {
            if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
            Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i]));
            if(Relators[i] == NULL) Mem_Error();     
            Scratch += GetHandleSize((char **) Copy_Of_Input[i]);                
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while( (*q++ = *p++) ) ;
            }
        Scratch -= NumRelators;                    
        Delete_Old_Presentations();
        Init_G_Variables();
        Length = Scratch;    
        goto _OPTIONS;
            
    
        /**************************************************************************************
            Call Get_Initial_Diagram(TRUE) to see whether this presentation is realizable.
        **************************************************************************************/

 _GET_INITIAL_DIAGRAM:                        
        switch(Get_Initial_Diagram(TRUE))
            {
            case 0:
                SaveMinima = TRUE;
                Input = RESET;
                break;
            case 1:
                TestRealizability1 = FALSE;
                TestRealizability2 = FALSE;
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,NULL);
                Input = BEGIN;
                goto _BEGIN;
            case 2:
                TestRealizability1 = FALSE;
                TestRealizability2 = FALSE;
                NoReport = FALSE;
                Input = BEGIN;
                goto _BEGIN;                    
            case REDUCE_GENUS:
                Input = REDUCE_GENUS;
                break;
            }        
        }
    
    TestRealizability1 = FALSE;
    TestRealizability2 = FALSE;
            
    Get_Diagrams();
    Input = BEGIN;
    goto _BEGIN;
    
_QUIT:    
#ifndef MAC
    tcsetattr(0, 0, &normal_termio);
#endif
    return(0);            
}        


int Rewrite_Input()
{
    /******************************************************************************************
        This routine is called after input of the initial set of relators. It reads through
        the relators to determine which generators appear and, if necessary, it rewrites the
        relators so that generators appear consecutively starting with A. It returns 1 if the
        relators are empty and otherwise returns 0.
    ******************************************************************************************/
    
    register unsigned char     *p,
                            	s,
                            	t,
                            	x,
                            	y,
                            	z;
                            
    register unsigned int     	i,
                            	j,
                            	k;
                            
    int                        	h,
                            	m;                        
                                                    
    unsigned char             **Temp;
    
    unsigned int             TT[125];
    
    /*****************************************************************************************
        Locate and count the number of nonempty relators. The number of nonempty relators is
        equal to k. TT[i] is set to true if Relators[i] is nonempty.
    *****************************************************************************************/
    
    for(i = 1,k = 0; i <= NumRelators; i++)
        {
        p = *Relators[i];
        if(*p)
            {
            k++;
            TT[i] = 1;
            }
        else
            TT[i] = 0;
        }
    TT[++i] = 0;
    
    /******************************************************************************************
        If k = 0, the presentation Heegaard was given has reduced to the empty presentation.
        Return(1);
    ******************************************************************************************/
        
    if(k == 0)
        {
        NumGenerators = 0;
        return(1);        
        }    
        
    /******************************************************************************************
        Rearrange the order of the relators so that all the nonempty relators come first and
        all empty relators come last as i runs from 1 to NumRelators. Set NumRelators equal 
        to the number of nonempty relators.
    ******************************************************************************************/    
            
    if(k < NumRelators)
        {        
        for(i = 1, j = NumRelators;  ;  )
            {                     
            while(TT[i]) i++;
            while(!TT[j]) j--;
            if(i >= j) break;
            Temp = Relators[i];
            Relators[i] = Relators[j];
            Relators[j] = Temp;
            TT[i] ++;
            TT[j] = 0;
            }
        NumRelators = k;    
        }
        
    /*****************************************************************************************
        Scan the relators and count the number of generators that appear. This number
        determines the initial Heegaard genus and determines the size of the arrays needed.
    *****************************************************************************************/
                    
    for(i = 65; i < 125; i++) TT[i] = 0;
    for(i = 1; i <= NumRelators; i++)
        {
        p = *Relators[i];
        while( (z = *p++) ) TT[z]++;
        }
    NumGenerators = 0;
    for(i = 65, j = 97; i < 91; i++,j++) if(TT[i] || TT[j])
        {
        NumGenerators ++;
        TT[i]++;
        }
    Vertices = 2*NumGenerators;
    
    /******************************************************************************************
        Rewrite the relators so that generators appear in consecutive alphabetical order 
        starting with A and continuing without omissions to the last generator. 
        This keeps arrays as small as possible and simplifies some other things as well.
    ******************************************************************************************/    
    
    for(h = m = 0,i = 65,j = 90; ; )
        {
        while(TT[i]) i++;
        while(!TT[j]) j--;
        if(i >= j) break;
        x = i;
        y = x + 32;
        s = j;
        t = j + 32;
        h ++;
        if(Micro_Print)
            {
            if(h == 1)
                printf("\n\nMade the following replacements of generators in the presentation:\n");
            printf("%d) %c -> %c ",h,s,x);  
            }
        for(k = 1; k <= NumRelators; k++)
            {
            p = *Relators[k];
            while( (z = *p) )
                {
                if(z == s) *p = x;
                if(z == t) *p = y;
                p++;
                }
            }
        TT[i]++;
        TT[j] = 0;
        }    
    return(0);
}

int Delete_Dups(void)
{
    /******************************************************************************************
            This routine checks for and counts relators in the input presentation which
            are cyclic conjugates of another relator or cyclic conjugates of the inverse 
            of another relator. It returns the number of distinct nonempty relators.
    ******************************************************************************************/
    
    register unsigned char  **Temp;
        
    register int            i,
                            j;
    
    int                     SNumRelators;
    
    unsigned long           HS;
    
    SNumRelators = NumRelators;
    NumEmptyRels = 0;
    
    /*****************************************************************************************
                        First, find and remove any empty relators.
    *****************************************************************************************/
    
    for(i = 1; i < SNumRelators; i++)
        {
        HS = GetHandleSize((char **) Relators[i]);
        if(HS == 1)                        
            {
            SNumRelators--;
            NumEmptyRels ++;
                
            /*********************************************************************************
            Relators[i] is empty, find the last nonempty relator and swap it with Relators[i].
            *********************************************************************************/
                
            for(j = SNumRelators + 1; j > i; j--)
                {
                if(GetHandleSize((char **) Relators[j]) > 1)
                    {
                    Temp = Relators[i];
                    Relators[i] = Relators[j];
                    Relators[j] = Temp;
                    break;
                    }
                SNumRelators--;
                NumEmptyRels++;    
                }
            }
        }
    
    /*****************************************************************************************
        After any empty relators have been removed, find and delete any relators which
        are cyclic conjugates of another relator or of the inverse of another relator.            
    *****************************************************************************************/
    
    for(i = 1; i < SNumRelators; i++)    
        {
        HS = GetHandleSize((char **) Relators[i]);
        for(j = SNumRelators; j > i; j--)
            if(HS == GetHandleSize((char **) Relators[j])
                 && Compare_Str(*Relators[i],*Relators[j],HS - 1))
                    {
                    Temp = Relators[j];
                    Relators[j] = Relators[SNumRelators];
                    Relators[SNumRelators] = Temp;
                    SNumRelators--;
                    }                
        Inverse(*Relators[i]);
        for(j = SNumRelators; j > i; j--)        
            if(HS == GetHandleSize((char **) Relators[j])
                && Compare_Str(*Relators[i],*Relators[j],HS - 1))
                    {
                    Temp = Relators[j];
                    Relators[j] = Relators[SNumRelators];
                    Relators[SNumRelators] = Temp;
                    SNumRelators--;
                    }        
        Inverse(*Relators[i]);                    
        }    
    return(SNumRelators);                
}

int Compare_Input_Pres()
{
    /******************************************************************************************
                This routine compares the presentation in Copy_Of_Input[] and the
            presentation given by Relators[]. It returns TRUE if they are identical and
            FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char 	*p,
                            *q;
                            
    register int             i;
    
    if(NumRelators != CopyNumRelators) return(FALSE);
    for(i = 1; i <= NumRelators; i++)
        {
        p = *Relators[i];
        q = *Copy_Of_Input[i];
        while(*p && (*p == *q))
            {
            p++;
            q++;
            }
        if(*p || *q) return(FALSE);
        }        
    return(TRUE);            
}

int Compare_Pres(k)
register int k;
{
    /******************************************************************************************
                This routine compares the presentation in SUR[k][] and the
            presentation given by Relators[]. It returns TRUE if they are identical and
            FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char 	*p,
                            *q;
                            
    register int             i;
    
    for(i = 1; i <= NumRelators; i++)
        {
        p = *Relators[i];
        q = *SUR[k][i];
        while(*p && (*p == *q))
            {
            p++;
            q++;
            }
        if(*p || *q) return(FALSE);
        }        
    return(TRUE);            
}

int Compare_Dual_Pres(k)
register int k;
{
    /******************************************************************************************
            This routine compares the presentation in SUR[k][] and the presentation given by
            DualRelators[]. It returns TRUE if they are identical and FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q;
                            
    register int             i;
    
    for(i = 1; i <= NumRelators; i++)
        {
        p = *DualRelators[i];
        q = *SUR[k][i];
        while(*p && (*p == *q))
            {
            p++;
            q++;
            }
        if(*p || *q) return(FALSE);
        }        
    return(TRUE);                
}

int Get_Initial_Diagram(PrintFlag)
int        PrintFlag;
{
    /******************************************************************************************
                        Get_Initial_Diagram() is supposed to determine 
            if the initial presentation is realizable by a Heegaard diagaram.
    ******************************************************************************************/
        
    register unsigned char  *p,
                            *q;
    
    register int            i,
                            j;
                            
    unsigned int            Flag1,
                            Flag2,
                            Flag3,
                            Flag4,
                            SNumFilled;                                                
    
    int                     Del_Only_Triv_Rel,
                            DistinctNonEmpty,
                            FirstPass,
                            NumReTries1,
                            NumReTries2,
                            SMicro_Print,
                            SMicro_Print_F,
                            SNumRelators;
                            
    unsigned long           HS;                        
    
    unsigned int 			Whitehead_Graph();
    unsigned int 			Reduce_Genus();        
    
    Input 				= NORMAL;
    TestRealizability1 	= TRUE;
    TestRealizability2 	= FALSE;
    Del_Only_Triv_Rel  	= TRUE;
    NumReTries1 		= 0;
    NumReTries2 		= 0;
    FirstPass   		= TRUE;
    BdryData 			= TRUE;

    if(NumGenerators > 2 && Initial_Realizability_Check()) return(2);
        
_RESTART:
    if(mykbhit() == ' ' || Level_Interrupt)
        {
        Level_Interrupt = FALSE;
        SMicro_Print = Micro_Print;
        SMicro_Print_F = Micro_Print_F;
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;
        LIST_OPTIONS2:
        printf("\n\nHIT 't' TO TERMINATE THIS RUN.");
        printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS.");
        printf("\nHIT 'x' TO TOGGLE MICRO_PRINTING ON AND OFF.");
        printf("\nHIT 'r' TO RESUME RUNNING THIS EXAMPLE.");
        GET_RESPONSE2:        
        switch(WaitkbHit())
            {
            case 'r':
                break;
            case 't':
                return(1);
            case 'v':
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,NULL);
                goto LIST_OPTIONS2;
            case 'x':
            	if(SMicro_Print)
					{
					printf("\nMicro_Printing IS ON. TURN Micro_Printing OFF ? HIT 'y' OR 'n'.\n");
					GET_RESPONSE1:
					switch(WaitkbHit())
						{
						case 'y':	
							SMicro_Print 	= FALSE;
							SMicro_Print_F 	= FALSE;
							break;
						case 'n':
							SMicro_Print 	= TRUE;
							SMicro_Print_F 	= TRUE;
							break;
						default:
							if(Batch == FALSE) SysBeep(5);
							goto GET_RESPONSE1;
						}
					}
				else
					{
					printf("\nTURN Micro_Printing ON ? HIT 'y' OR 'n'.\n");
					GET_RESPONSE3:					
					switch(WaitkbHit())
						{
						case 'y':
							SMicro_Print 	= TRUE;
							SMicro_Print_F 	= TRUE;			
							break;
						case 'n':
							SMicro_Print 	= FALSE;
							SMicro_Print_F 	= FALSE;
							break;
						default:
							if(Batch == FALSE) SysBeep(5);
							goto GET_RESPONSE3;
						}
					}	
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE2;
            }
        NoReport = TRUE;    
        Micro_Print = SMicro_Print;
        Micro_Print_F = SMicro_Print_F;    
        printf("\n\n     Resumed running. . . .\n");                    
        }    
        
    SNumRelators = NumRelators;
    
    if(!FirstPass && !Do_Not_Reduce_Genus) switch(Delete_Trivial_Generators(FALSE))
        {
        case 0:
            break;
        case 1:
            if(!Micro_Print)
                {
                if(SNumRelators == NumRelators + 1)
                    printf("\n\n                    Deleted a trivial generator.\n");
                else    
                    printf("\n\n                    Deleted %d trivial generators.\n",
                        SNumRelators - NumRelators);
                }
            break;
        case TOO_LONG:
			printf("\n\n                    This presentation may be too long!");
            return(1);
        }
    
    FirstPass = FALSE;
    
    if(NumGenerators == 1)
        {
        DistinctNonEmpty = Delete_Dups();
        if(On_File() == NumFilled)
            {
            if(DistinctNonEmpty > 1)
                {
                Fatal_Error();                
                return(1);            
                }
            Print_Realizability(Del_Only_Triv_Rel,NumFilled);                    
            if(DistinctNonEmpty == 0)
                {
                if(Save_Pres(ReadPres,0,Length,1,2,0,0,0)) return(1);
                UDV[NumFilled - 1] = S1_X_X2;
                Mark_As_Found_Elsewhere(CurrentComp);    
                return(1);    
                }
            CBC[CurrentComp][0] = NumRelators - NumEmptyRels;
            CBC[CurrentComp][1] = BDRY_UNKNOWN;
            HS = GetHandleSize((char **) Relators[1]) - 1;    
            if(HS == 1)
                {
                if(Save_Pres(ReadPres,0,Length,1,2,0,0,0)) return(1);                    
                UDV[NumFilled - 1] = THREE_SPHERE;
                BDY[NumFilled - 1] = FALSE;
                Mark_As_Found_Elsewhere(CurrentComp);
                return(1);
                }
            Canonical_Rewrite(Relators,FALSE,FALSE);
            if(Save_Pres(ReadPres,0,Length,1,2,0,0,0)) return(1);
            BDY[NumFilled - 1] = FALSE;
            UDV[NumFilled - 1] = GENERIC_LENS_SPACE;
            LSP[NumFilled - 1] = HS;
            LSQ[NumFilled - 1] = 1;        
            return(1);
            }
        return(1);    
        }
        
    if(Length == 0)
        {
        printf("\n\nThe presentation has reduced to a presentation of length zero.");
        Print_Realizability(Del_Only_Triv_Rel,NumFilled);
        if(On_File() == NumFilled)
        Save_Pres(ReadPres,0,Length,1,2,0,0,0);
        return(1);
        }    
    
    Flag4 = Find_Flow_A(NORMAL,FALSE);
    
    if(Micro_Print)
        {
        if(Automorphisms)
            {
            printf("\n\n%lu automorphism(s) reduced the length to %lu.",
                Automorphisms,Length);
            printf("\n\nThe presentation is currently:\n");
            Print_Relators(Relators,NumRelators);
            }
        else
            printf("\n\nThe current set of relators has minimal length of %lu.",Length);       
        }
            
    if(Flag4 == TOO_LONG)
        {
		printf("\n\n                    This presentation may be too long!");        
        return(1);
        }                
    if(Flag4 == 1)
        {
        Modified_Init_Pres = TRUE;
        if(NumRelators == 1)
            {
            printf("\n\n                    This relator is not of full rank.");
            printf("\n\n                    Reducing the number of generators. . .\n");
            }
        else
            {    
            printf("\n\n                    The initial relators are not of full rank.");
            printf("\n\n                    Reducing the number of generators. . .\n");
            }
                
        if(NumFilled == 0)
            {
            j = NumRelators;
            NumRelators = CopyNumRelators;
            for(i = 1, Length = 0L; i <= NumRelators; i++)
                Length += GetHandleSize((char **) Copy_Of_Input[i]);
            Length -= NumRelators;
            Canonical_Rewrite(Copy_Of_Input,FALSE,FALSE);
            if(Save_Pres(ReadPres,0,Length,2,2,1,0,0)) return(1);
            UDV[NumFilled - 1] = 0;
            BDY[NumFilled - 1] = BDY[ReadPres];
            NumRelators = j;
            }
        else
            {
            This_Pres = On_File();
            if(This_Pres == NumFilled)
                {
                if(Dup_On_File < INFINITE)
                    {
                    This_Pres = Dup_On_File;
                    if(Save_Pres(ReadPres,Dup_On_File,Length,1,2,1,0,0)) return(1);
                    Mark_As_Duplicate(Dup_On_File);
                    BDY[NumFilled - 1] = BDY[ReadPres];
                    ReadPres = This_Pres;
                    CurrentComp = ComponentNum[Dup_On_File];        
                    }
                else
                    {
                    if(Save_Pres(ReadPres,0,Length,1,2,1,0,0)) return(1);
                    UDV[NumFilled - 1] = 0;
                    BDY[NumFilled - 1] = BDY[ReadPres];
                    ReadPres = This_Pres;
                    }
                }    
            }    
        
        SNumFilled = NumFilled;
        if(Missing_Gen()) return(1);
        if(NumFilled == SNumFilled)
            {
            FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
            NumReTries1 ++;
            if(NumReTries1 > 3)
                {
                printf("\n\n                    Unable to determine whether the presentation is realizable.");
                printf("\n\nTry computing the orbit of this presentation under level-transformations. Then check members of");
        		printf("\nthe orbit for the absence of separating pairs of vertices. If such presentations exist, Heegaard");
        		printf("\ncan determine whether this presentation is realizable.");
                return(1);                
                }
            printf("\n\n                    Retrying the initial presentation.");
            NumRelators  		= CopyNumRelators;
            NumGenerators     	= CopyNumGenerators;
            Vertices         	= 2*NumGenerators;
            Del_Only_Triv_Rel  	= TRUE;
            FirstPass         	= TRUE;
            for(i = 1,Length = 0L; i <= NumRelators; i++)
                {
                if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
                Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i]));
				if(Relators[i] == NULL) Mem_Error();
                Length += GetHandleSize((char **) Copy_Of_Input[i]);                
                p = *Copy_Of_Input[i];
                if(Relators[i] == NULL)
                    {
                    printf("\n\n    Memory Error. Sorry!");
                    return(1);
                    }
                q = *Relators[i]; 
                while( (*q++ = *p++) ) ;
                }
            Length -= NumRelators;
            ReadPres = 0;
            CurrentComp = 1;
            goto _RESTART;
            }
                
        if(NumGenerators == 1)
            {
            DistinctNonEmpty = Delete_Dups();
            if(DistinctNonEmpty > 1)
                {
                Fatal_Error();                
                return(1);            
                }
            Print_Realizability(Del_Only_Triv_Rel,NumFilled - 1);                    
            if(DistinctNonEmpty == 0)
                {
                UDV[NumFilled - 2] = S1_X_X2;
                Mark_As_Found_Elsewhere(ComponentNum[NumFilled - 2]);
                return(1);    
                }
            CBC[ComponentNum[NumFilled - 2]][0] = NumRelators - NumEmptyRels;
            CBC[ComponentNum[NumFilled - 2]][1] = BDRY_UNKNOWN;
            
            HS = GetHandleSize((char **) Relators[1]) - 1;        
            if(HS == 1)
                {
                BDY[NumFilled - 2] = FALSE;
                BDY[NumFilled - 1] = TRUE;
                UDV[NumFilled - 2] = THREE_SPHERE;
                Mark_As_Found_Elsewhere(ComponentNum[NumFilled - 2]);            
                return(1);
                }
                
            BDY[NumFilled - 2] = FALSE;
            BDY[NumFilled - 1] = TRUE;
            UDV[NumFilled - 2] = GENERIC_LENS_SPACE;
            LSP[NumFilled - 2] = HS;
            LSQ[NumFilled - 2] = 1;                            
            return(1);
            }

        CurrentComp = ComponentNum[NumFilled - 2];
        ReadPres = This_Pres = NumFilled - 2;
        FirstPass = TRUE;
        printf("\n\n                    Looking for diagrams of summand %d. . .\n",
            CurrentComp);
        goto _RESTART;    
        }
        
    if(Automorphisms) Modified_Init_Pres = TRUE;    
    This_Pres = On_File();
    if(Micro_Print)
        {
        printf("\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
        Print_Relators(Relators,NumRelators);
        }                                            
    if(This_Pres == NumFilled)
        {
        if(Dup_On_File < INFINITE)
            {
            This_Pres = Dup_On_File;
            if(Save_Pres(ReadPres,Dup_On_File,Length,1,2,1,0,0)) return(1);
            Mark_As_Duplicate(Dup_On_File);
            BDY[NumFilled - 1] = BDY[ReadPres];
            ReadPres = This_Pres;
            CurrentComp = ComponentNum[Dup_On_File];        
            }
        else
            {
            if(Save_Pres(ReadPres,0,Length,1,2,1,0,0)) return(1);
            UDV[NumFilled - 1] = 0;
            BDY[NumFilled - 1] = BDY[ReadPres];
            ReadPres = This_Pres;
            }
        }
    Save_Init_Pres = This_Pres;    
    Fill_A(NumRelators);
    Input = NORMAL;
    Saved_Vertices = 0;
    NumDiagrams ++;
    
    if( (Flag1 = Whitehead_Graph()) ) switch(Flag1)
        {
        case NON_PLANAR:
            printf("\n\n                    The Whitehead graph is nonplanar.");
        case FATAL_ERROR:
            Fatal_Error();        
            return(1);
        case TOO_LONG:
			printf("\n\n                    This presentation may be too long!");
            return(1);
        case TOO_MANY_COMPONENTS:
            return(1);        
        case NON_UNIQUE_1:
        case NON_UNIQUE_2:
        case NON_UNIQUE_3:
        case NON_UNIQUE_4:
            UDV[This_Pres] = Flag1;
        case V2_ANNULUS_EXISTS:    
            if(Non_Unique_Initial_Diagram())
            	{
            	if(Batch == 4 || Batch == 10 || Batch == 11 || Batch == 53)
                return(1);
                }
            return(REDUCE_GENUS);
        case REDUCE_GENUS:
            goto _RESTART;    
        case NOT_CONNECTED:
            /**********************************************************************************
                The Whitehead graph corresponding to the initial presentation is not
                connected. Hence the corresponding Heegaard diagram is reducible.
            **********************************************************************************/
            
            if(!Del_Only_Triv_Rel)
                Realization_Warning();
            printf("\n\n                    The Whitehead graph is not connected.");
            printf("\n\n                    Please check each summand separately.");
            BdryData = FALSE;
            return(1);
                        
        /**************************************************************************************
                If Heegaard gets to this point, then Heegaard cannot find the diagram
            corresponding to the initial presentation. Before we give up, we let Heegaard
            try to find some level transformations which might change the presentation into
            one for which Heegaard can find the corresponding diagram. If Heegaard can't
            find any level transformations that work,then we look for ways to reduce the genus.
        **************************************************************************************/    
                
        case SEP_PAIRS:
            printf("\n\n                    The Whitehead graph has a separating pair of vertices.\n");
            UDV[This_Pres] = SEP_PAIRS;
            if(V1 & 1)
                LSP[This_Pres] = V1/2 + 97;
            else
                LSP[This_Pres] = V1/2 + 65;
            if(V2 & 1)
                LSQ[This_Pres] = V2/2 + 97;
            else
                LSQ[This_Pres] = V2/2 + 65;                
            Num_Saved_LPres = 0;
            NotNewPres = 0;
            SNum_Level_Slides = Num_Level_Slides;
            ReadPres = This_Pres; 
            switch(Flag3 = Level_Transformations(TRUE,TRUE,FALSE))
                {
                case 0:
                    break;
                case 1:
                    break;
                case 2:
                    Fill_A(NumRelators);
                    Saved_Vertices = 0;
                    NumDiagrams ++;
                    ReadPres = This_Pres = NumFilled - 1;
                    Flag2 = Whitehead_Graph();
                    switch(Flag2)
                        {
                        case NON_PLANAR:
                            printf("\n\n                    The Whitehead graph is nonplanar.");
                        case FATAL_ERROR:
                            Fatal_Error();        
                            return(1);
                        case TOO_LONG:
							printf("\n\n                    This presentation may be too long!");
                            return(1);
                        case NON_UNIQUE_1:
                        case NON_UNIQUE_2:
                        case NON_UNIQUE_3:
                        case NON_UNIQUE_4:
                            UDV[This_Pres] = Flag2;
                        case V2_ANNULUS_EXISTS: 
                        	if(Batch != 4 && Batch != 10 && Batch != 11) return(1); 
                            if(Non_Unique_Initial_Diagram())
                                return(1);
                            return(REDUCE_GENUS);    
                        case NO_ERROR:
                            printf("\n\n	After examining %u diagram(s) and performing %lu Sep-Vert-Slide(s):"
                            	,Num_Saved_LPres, Num_Level_Slides - SNum_Level_Slides);
                            Print_Realizability(Del_Only_Triv_Rel,ReadPres + 1);
                            if(UDV[This_Pres] == V2_ANNULUS_EXISTS)
                                printf("\n\n                    However, the realization is not unique because an annulus exists.");    
                            if(Micro_Print)
                                {
                                printf("\n\nStarted with Presentation %d:\n",ReadPres + 1);
                                Print_Relators(Relators,NumRelators); 
                                }                                
                            goto DELETE_REDUNDANT;
                        break;    
                        }        
                    break;
                case 3:
                    printf("\n\n                    The Whitehead graph is nonplanar.");
                    Fatal_Error();        
                    return(1);
                case 4:
                    break;
                case 5:
                    printf("\n\n          After some level transformations, ");
                    printf("the presentation contains a trivial generator.");
                    if(!Do_Not_Reduce_Genus)
                        {
                        goto _RESTART;
                        }
                    break;    
                case 6:
                    break;
                case 7:
                    printf("\n\n            Out of memory for level transformations.");
                    break;    
                case TOO_LONG:
                    break;
                case FATAL_ERROR:
                    Fatal_Error();
                    return(1);  
                 default:
                 	break;
                }
            if(Num_Saved_LPres) for(i = 1; i <= NumRelators; i++)
                {
                if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
                Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SLR[0][i]));
                if(Relators[i] == NULL) Mem_Error();
                p = *Relators[i];    
                q = *SLR[0][i];
                while( (*p++ = *q++) ) ;
                } 
            if(Micro_Print)
                {
                printf("\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
                Print_Relators(Relators,NumRelators);
                }    
            if(Flag3 != 5)
                {                
                printf("\n\n                    Unable to remove the separating pair of vertices by Sep_Vert_Slides.");
                if((Batch != 4 && Batch != 10 && Batch != 11) && Do_Not_Reduce_Genus) return(1);
                if(NumReTries1 <= 3 && !Do_Not_Reduce_Genus)
                	{
                	if(Batch == FALSE)
                		{
						printf("\n\nHeegaard would like to try to reduce the genus of the presentation. However, while deleting primitives\n");
						printf("preserves fundamental groups, there is a small probability genus reduction will convert an unrealizable\n");
						printf("into a realizable presentation.\n");
						printf("The alternative is to abort, rerun the initial presentation, reduce the initial presentation to minimal-\n");
						printf("length, have Heegaard compute the full orbit of the initial presentation under level-transformations, and\n");
						printf("then to test the members of the orbit for the absence of separating pairs of vertices.");
						printf("\n\n                    Hit 'p' to allow Heegaard to proceed. Hit 'q' to abort.");
                    	}
                    GET_RESPONSE5:
                    if(Batch == FALSE)        
                    switch(WaitkbHit())
						{
						case 'p':		
							break;
						case 'q':
							return(1);
						default:
							if(Batch == FALSE) SysBeep(5);
							goto GET_RESPONSE5;
                    	}
                    }
                }
            
        default:
        	
        	if(NumRelators > 1 && NumReTries1 <= 3 && !Do_Not_Reduce_Genus)
				{
				if(Batch) printf("\n\nTrying to delete primitives from the presentation.");
				
				/**********************************************************************************
					Call Reduce_Genus() to see if there are primitives which can be removed.
				**********************************************************************************/
				
				SNumRelators = NumRelators;
				SReadPres = This_Pres;                
				switch(Reduce_Genus(NORMAL,TRUE,TRUE))
					{
					case NO_ERROR:
						break;
					case FATAL_ERROR:
						Fatal_Error();
						return(1);
					case TOO_LONG:
					case CAN_NOT_DELETE:
						if(EmtyRel) break;
						printf("\n\n                    Unable to delete a primitive.");                 
						FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
						NumReTries1 ++;
						printf("\n\n                    Retrying the initial presentation.");
						NumRelators         = CopyNumRelators;
						NumGenerators       = CopyNumGenerators;
						Vertices            = 2*NumGenerators;
						Del_Only_Triv_Rel   = TRUE;
						FirstPass           = TRUE;
						for(i = 1,Length = 0L; i <= NumRelators; i++)
							{
							if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
							Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i]));
							if(Relators[i] == NULL) Mem_Error();
							Length += GetHandleSize((char **) Copy_Of_Input[i]);                
							p = *Copy_Of_Input[i];
							if(Relators[i] == NULL)
								{
								printf("\n\n    Memory Error. Sorry!");
								return(1);
								}
							q = *Relators[i];    
							while( (*q++ = *p++) ) ;
							}
						Length -= NumRelators;
						ReadPres = 0;
						CurrentComp = 1;
						goto _RESTART;
					}
					
				if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
					{
					Del_Only_Triv_Rel = FALSE;
					if(LensSpace)
						{
						if(SNumRelators > NumRelators)
							Realization_Warning();
						printf("\n\n                    This manifold is a lens space.");
						return(REDUCE_GENUS);
						}
					if(FoundPrimitive)
						{
						FoundPrimitive = FALSE;
						printf("\n\n                    Heegaard found relator(s) which are primitive");
						printf("\n                    and deleted their consequences from the presentation.");
						}
					if(FoundPower)
						{
						FoundPower = FALSE;
						printf("\n\n                    Heegaard found relator(s) which are proper powers");
						printf("\n                    and deleted their consequences from the presentation.");
						}
					if(EmtyRel)
						{
						if(SNumRelators > NumRelators)
							Realization_Warning();
						EmtyRel = FALSE;
						printf("\n\n                    This manifold is a connected sum.");
						printf("\n                    Please check each summand separately.\n");
						return(1);
						}
					printf("\n");
					goto _RESTART;
					}
				}
            
        	if(Delete_Only_Short_Primitives)
				{
				ReadPres = This_Pres = Save_Init_Pres;
				if(Get_Relators_From_SUR(ReadPres))
					{
					printf("\n\n    Memory Error. Sorry!");
					return(1);
					}
				Num_Saved_LPres = 0;
				NotNewPres = 0;
				SNum_Level_Slides = Num_Level_Slides;
				if(Micro_Print)
					{
					printf("\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
					Print_Relators(Relators,NumRelators);
					}
				Fill_A(NumRelators);
				Get_Matrix();
				switch(Flag3 = Level_Transformations(TRUE,FALSE,FALSE))
					{
					case 3:
						printf("\n\n                    The Whitehead graph is nonplanar.");
						Fatal_Error();        
						return(1);                                        
					case 5:
						if(Micro_Print)
							{
							printf("\n\nThe presentation is currently:\n");
							Print_Relators(Relators,NumRelators);   
							}
						printf("\n\n	After examining %u diagram(s) and performing %lu Sep-Vert-Slide(s):"
							,Num_Saved_LPres, Num_Level_Slides - SNum_Level_Slides);
						printf("the presentation contains a trivial generator.");                  
						goto _RESTART;
					case 7:
						printf("\n\n            Out of memory for level transformations.");
						break;                        
					case FATAL_ERROR:
						Fatal_Error();
						return(1);                
					default:
						break;                    
					}
				}
        
        	if(!Delete_Only_Short_Primitives && NumReTries2 < 3 && Flag1 == SEP_PAIRS && Flag3 == 4)
				{
				NumReTries2 ++;    
				DistinctNonEmpty = Delete_Dups();
				j = NumRelators - DistinctNonEmpty;
				if(j)
					{
					Del_Only_Triv_Rel = FALSE;
					FirstPass = TRUE;
					NumReTries1 = 0;
					printf("\n\n                    Unable to remove the annulus which is present.");
					if(j == 1)
						{
						if(Batch == FALSE) SysBeep(5);
						printf("\n\n    Modifying the initial presentation by deleting a duplicated or empty relator.\n");
						printf("\n    NOTE: This modification may change the homeomorphism type of M.\n");
						}                
					else
						{
						if(Batch == FALSE) SysBeep(5);
						printf("\n\n    Modifying the initial presentation by deleting %d duplicated or empty relators.\n", j);
						printf("\n    NOTE: This modification may change the homeomorphism type of M.\n");
						}
					NumRelators = DistinctNonEmpty;
					for(i = 1,Length = 0L; i <= NumRelators; i++)
						Length += GetHandleSize((char **) Relators[i]);
					Length -= NumRelators;        
					printf("\n\nThe modified presentation is:\n");
					Print_Relators(Relators,NumRelators);
					goto _RESTART;
					}
				}
                    
        printf("\n\n                    Unable to determine whether the presentation is realizable.");
        printf("\n\nTry computing the orbit of this presentation under level-transformations. Then check members of");
        printf("\nthe orbit for the absence of separating pairs of vertices. If such presentations exist, Heegaard");
        printf("\ncan determine whether this presentation is realizable.");
        return(1);    
        }
        
    /******************************************************************************************
            Otherwise, if Heegaard gets to this point, then Heegaard found the diagram
        corresponding to the initial presentation. So we are in business.
    ******************************************************************************************/    
    
    Print_Realizability(Del_Only_Triv_Rel,Save_Init_Pres + 1);
    if(UDV[This_Pres] == V2_ANNULUS_EXISTS)
        printf("\n\n                    However, the realization is not unique because an annulus exists.");    
        
    DELETE_REDUNDANT:
    	
    if(PrintFlag)            
    printf("\n\n                    Looking for other diagrams. . . .\n");	
        
    if(NumGenerators > 1)
        {
        Get_Bdry_Comps(TRUE,FALSE,This_Pres);
        if(BCWG[1] == BDRY_UNKNOWN)
            {
            BDY[This_Pres] = FALSE;
            Boundary = FALSE;
            }
        else
            {
            BDY[This_Pres] = TRUE;
            Boundary = TRUE;
            }
        for(i = 0; (CBC[CurrentComp][i] = BCWG[i]) < BDRY_UNKNOWN; i++) ;
        if(CS[ComponentNum[This_Pres] + 1] == 3) MG_Bdry_Comp_Data(This_Pres);    
        if(BCWG[0] > 1 || (BCWG[0] && NumBdryComps > BCWG[0]))
            Delete_Redundant_Relators();
        for(i = 0; i < NumFilled - 1; i++) ER[i] = -1;
        }

    return(0);            
}

int Save_Pres(unsigned int From,unsigned int Daut,unsigned long Len,int F1,int F2,int F3,unsigned char F4,char F5)
{
    /******************************************************************************************
        Save_Pres() is called when Heegaard has determined that a presentation should be
        put on file for further processing. It saves a copy of the presentation in the array
        SUR[][] and saves some other data about the presentation in the arrays listed below.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q,
                            *r;
                            
    unsigned char           ***MyRelators;
                    
    register int            i,
                            j;
                    
    int                     Error = FALSE;
        
    unsigned long			HS;
                
    switch(F1)
        {
        case 0:
            MyRelators = DualRelators;
            break;
        case 1:
            MyRelators = Relators;
            break;
        case 2:
            MyRelators = Copy_Of_Input;
            break;
        }

	if(F4 == 1)
		{
		if(SUR[NumFilled][0] != NULL) DisposeHandle((char **) SUR[NumFilled][0]);
		SUR[NumFilled][0] = (unsigned char **) NewHandle(sizeof(char)*(NumRelators + 1));
		if(SUR[NumFilled][0] == NULL) Mem_Error();
        p = *SUR[NumFilled][0];    
        *p++ = NumRelators;
        for(i = 1; i <= NumRelators; i++)  *p++ = 1;
		}
		
    for(i = 1; i <= NumRelators; i++)
        {
        HS = GetHandleSize((char **) MyRelators[i]);
        if(SUR[NumFilled][i] != NULL) DisposeHandle((char **) SUR[NumFilled][i]);
        SUR[NumFilled][i] = (unsigned char **) NewHandle(HS);            
        if(SUR[NumFilled][i] == NULL) Mem_Error();
        p = *MyRelators[i];
        q = *SUR[NumFilled][i];    
		r = q;
		while( (*q++ = *p++) ) ;  
		if((q-r) != HS) printf("\n\n1) Error in Presentation %u! |Relator[%d]| = %lu, HS = %lu.",NumFilled + 1,i,q-r-1,HS);
        }    

    ComponentNum[NumFilled] 	= CurrentComp;
    Daughters[NumFilled]    	= Daut;
    ER[NumFilled]            	= F5;
    FR[NumFilled]            	= From;
    NG[NumFilled]             	= NumGenerators;
    NR[NumFilled]             	= NumRelators;
    PRIM[NumFilled]         	= F2;
    SURL[NumFilled]         	= Len;
    SURNumX[NumFilled]			= 1;
    if(NumFilled && NG[ReadPres] == NumGenerators && ComponentNum[ReadPres] == CurrentComp)
    	{
    	HegSplNum[NumFilled] = HegSplNum[ReadPres];
    	j = HegSplNxt[ReadPres];
    	HegSplNxt[NumFilled] = j;
    	HegSplNxt[ReadPres] = NumFilled;
    	Mergers ++;
    	}
    else
    	{
    	HegSplNum[NumFilled] = NumFilled;
    	HegSplNxt[NumFilled] = NumFilled;
    	}
    NumLoops[NumFilled] 	 = 0;		
    BytesUsed                += Len;
    
    if(ER[From] < 0 && ComponentNum[From] == CurrentComp)
        ER[NumFilled] = ER[From];
    
    if(F3)
        TP[NumFilled] = NumRelators;
    else
        TP[NumFilled] = FALSE;            
    
    if(Len < MLC[CurrentComp][NumGenerators]) 
        MLC[CurrentComp][NumGenerators] = Len;

    NumFilled ++;
    SaveMinima = TRUE;
    
    switch(F4)
    	{
    	case 0:
    		OnStack += 2*NumGenerators;
    		break;
    	case 1:
    	case 2:
    		OnStack += NumRelators;
    		break;
    	case 3:
    		OnStack = NumFilled - ReadPres;
    	}
        
    if(Micro_Print)
        printf("\n\nSaved the current presentation as: Presentation %u.\n",NumFilled);
    
    if(Batch == FALSE)
    	{    
    	printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,CurrentComp);
    	printf("Gen%3d  Rel%3d  Length%6lu  From%6u  ",NumGenerators,NumRelators,Len,From + 1);
    	}
    
    if(Batch == FALSE)
    switch(F2)
		{
		case 1:
			printf("DR");
			break;
		case 2:
			printf("IP");
			break;
		case 3:
			printf("LS");
			break;
		case 4:
			printf("1G");
			break;
		case 5:
			printf("S3");
			break;
		case 6:
			printf("FP");
			break;
		case 7:
			printf("BC");
			break;
		case 8:
		case 108:
			printf("PM");
			break;
		case 9:
			printf("PM");
			break;
		case 10:
			printf("BC");
			break;
		case 11:
			printf("CF");
			break;
		case 12:
			printf("ER");
			break;
		case 13:
			printf("Er");
			break;            
		case 20:
			printf("NC");
			break;
		case 30:
			printf("MG");
			break;
		case 60:
			printf("PP");
			break;    
		case 70:
			printf("Lt");
			break;
		case 75:
			printf("LT");
			break;    
		case 80:
			printf("A2");
			break;    
		default:
			break;
		}
		
	if(Check_If_HS_Rep)
		{
		if(Len < HS_Rep_Length && NumGenerators <= HS_Rep_NumGens)
			{
			if(NumGenerators == HS_Rep_NumGens)
				{
				printf("\n\n		The Initial Presentation %s is not a Heegaard Splitting Rep!",PresName); 
				printf("\n	Pres %u lies on the same Heegaard surface and has length %lu < %lu !",
					NumFilled,Len,HS_Rep_Length);
				if(Batch == 4 && H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep!",PresName);	
				}
			if(NumGenerators < HS_Rep_NumGens)
				{
				printf("\n\n		The Initial Presentation %s is Stabilized and hence not a Heegaard Splitting Rep!",PresName); 
				printf("\n	Pres %u lies on a Heegaard surface of genus %d < the initial genus %d !",
					NumFilled,NumGenerators,HS_Rep_NumGens);
				if(Batch == 4 && H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep!",PresName);	
				}			
			if(Batch == FALSE)
				{
				printf("\n\nShould Heegaard Continue? Hit 'y' or 'n'.");
				Check_If_HS_Rep = CheckHSReps = FALSE;
				GET_RESPONSE1:			
				switch(WaitkbHit())
					{
					case 'y':
						break;
					case 'n':				
						Check_If_HS_Rep = 7;
						break;
					default: goto GET_RESPONSE1;
					}
				}
			else
				Check_If_HS_Rep = 7;		
			}
		}
    if(Error)
        {
        printf("\n    Memory Error. Sorry!");
        return(TOO_LONG);
        }
    return(NO_ERROR);        
}

int Non_Unique_Initial_Diagram(void)
{
    /******************************************************************************************
        This routine is called when Heegaard has determined that if the initial presentation
        is realizable by a Heegaard diagram, then the realization is not unique. The routine
        determines whether we have a presentation that is realizable, but not uniquely so, or
        we have a presentation that is not realizable at all. It returns 0 if the presentation
        is realizable and returns 1 if it is not realizable.
    ******************************************************************************************/
    
    int            	i;
                    
    unsigned int    AnnulusExists,
                    Flag1;                
    
    unsigned int 	Whitehead_Graph();
    unsigned int 	Reduce_Genus();
    
    /******************************************************************************************
                            Check whether a valence-two annulus exists.
    ******************************************************************************************/
    
    AnnulusExists         	= FALSE;    
    TestRealizability2     	= TRUE;
    TestRealizability3     	= TRUE;
    DrawingDiagrams     	= FALSE;
    
    Flag1 = Whitehead_Graph();
    
    TestRealizability3     	= FALSE;
    
    if(Flag1 == V2_ANNULUS_EXISTS)
        {
        Flag1 = Whitehead_Graph();
        AnnulusExists = TRUE;
        }
        
    switch(Flag1)
        {
        case NO_ERROR:
            printf("\n\n                    Presentation %d is realizable.",This_Pres + 1);
            printf("\n\n                    But the realization is not unique.");
            if(Batch) return(NO_ERROR);
            printf("\n\n                    Hit 'p' to have Heegaard investigate one possible realization.");
            printf("\n\n                    Hit any other key to abort further investigation.");          
            switch(WaitkbHit())
				{
				case 'p': break;
				default: return(NO_ERROR);
				}
					
            if(AnnulusExists)
                {
                printf("\n\n    NOTE: Heegaard is investigating one diagram that realizes the initial presentation.");
                printf("\n          Subsequent results may not apply to all realizations.");
                }
                
            if(NumGenerators > 1)
                {
                Get_Bdry_Comps(TRUE,FALSE,This_Pres);
                if(BCWG[1] == BDRY_UNKNOWN)
                    {
                    BDY[This_Pres] = FALSE;
                    Boundary = FALSE;
                    }
                else
                    {
                    BDY[This_Pres] = TRUE;
                    Boundary = TRUE;
                    }
                for(i = 0; (CBC[CurrentComp][i] = BCWG[i]) < BDRY_UNKNOWN; i++) ;
                if(CS[ComponentNum[This_Pres] + 1] == 3) MG_Bdry_Comp_Data(This_Pres);    
                if(BCWG[0] > 1 || (BCWG[0] && NumBdryComps > BCWG[0]))
                    Delete_Redundant_Relators();
                for(i = 0; i < NumFilled - 1; i++) ER[i] = -1;
                }                
            
            SReadPres = This_Pres;
            switch(Reduce_Genus(NORMAL,TRUE,TRUE))
                {
                case NO_ERROR:
                    break;
                case FATAL_ERROR:
                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                    Fatal_Error();
                    return(1);
                case TOO_LONG:
                case CAN_NOT_DELETE:
                	if(EmtyRel) break;
                    printf("\n\n                    Unable to delete a primitive.");                 
                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                    return(1);        
                }
            if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
                {
                if(!AnnulusExists)
                    {
                    printf("\n\nNOTE: Heegaard is investigating one diagram that realizes the initial presentation.");
                    printf("\n      Subsequent results may not apply to all realizations.");
                    }
                if(FoundPrimitive)
                    {
                    printf("\n\n                    Heegaard found a relator which is a primitive");
                    printf("\n                    and deleted its consequences from the presentation.");
                    }
                if(FoundPower)
                    {
                    printf("\n\n                    Heegaard found a relator which is a proper power");
                    printf("\n                    and deleted its consequences from the presentation.");
                    }
                if(LensSpace)
                    printf("\n\n                    This manifold is a lens space.");
                if(EmtyRel)
                    printf("\n\n                    This manifold is a connected sum.");
                printf("\n");
                return(0);
                }    
            return(1);
        case FATAL_ERROR:
            printf("\n\n                    The initial presentation is not realizable.\n");
            Print_Relators(Relators,NumRelators);
            return(1);
        }
    return(NO_ERROR);            
}

unsigned int On_File(void)
{
    /******************************************************************************************
        On_File() is called to check whether the presentation in Relators[] is already on file
        in SUR[][]. If it is on file, it may either be a duplicate of another presentation
        from the CurrentComp or a duplicate of a presentation from some other summand. The
        global Dup_On_File gives the presentation duplicated in the case that the current
        presentation is a duplicate of a presentation from another summand.
            Otherwise, if the presentation is not currently on file, the routine returns
        the value NumFilled.
    ******************************************************************************************/    
        
    register int  	i,
    				ii,
                   	j,
                   	jj,
                   	k;
    
    if(Length > 0L)    Canonical_Rewrite(Relators,FALSE,FALSE);        
/*     for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
     if(SURL[i] == Length  
        && NG[i] == NumGenerators
        && NR[i] == NumRelators)
        {
         for(j = 1; j <= NumRelators; j++)
             if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;
         if(j > NumRelators && Compare_Pres(i))
             {
             if(ComponentNum[i] != CurrentComp)
                 {
                 if(Dup_On_File == INFINITE) Dup_On_File = i;
                 }
             else
                 break;
             }    
         }							*/
 
 /*	Modified 2/18/15 to just check for duplications within a given component number. */  
       
    for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
     if(ComponentNum[i] == CurrentComp && SURL[i] == Length && NG[i] == NumGenerators && NR[i] == NumRelators)
        {
         for(j = 1; j <= NumRelators; j++)
             if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;
         if(j > NumRelators && Compare_Pres(i)) break;    
         }    
     
     if(Micro_Print)
         {
         if(i < NumFilled)
             printf("\n\nA copy of the current presentation is already on file as presentation %d.",
                 i+1);
         else if(Dup_On_File < INFINITE)
             printf("\n\nThe current presentation is a duplicate of presentation %d of summand %d.",
                 i+1,ComponentNum[Dup_On_File]);
         }  
         
     if(i < NumFilled) 
     	{
     	if(SSSReadPres == i) NumLoops[i] ++;
     	else SURNumX[i] ++; 
		if(NG[ReadPres] == NG[i] && ComponentNum[i] == CurrentComp)
			{
			if(HegSplNum[i] < HegSplNum[ReadPres])
				{
				ii = HegSplNum[i];
				jj = HegSplNum[ReadPres];
				k = ReadPres;
				while(1)
					{
					if(HegSplNum[k] == ii) break;
					HegSplNum[k] = ii;
					k = HegSplNxt[k];
					}
				ii = HegSplNxt[i];
				jj = HegSplNxt[ReadPres];
				HegSplNxt[ReadPres] = ii;
				HegSplNxt[i] = jj;	
				Mergers ++;
				}
			if(HegSplNum[i] > HegSplNum[ReadPres])
				{
				ii = HegSplNum[i];
				jj = HegSplNum[ReadPres];
				k = i;
				while(1)
					{
					if(HegSplNum[k] == jj) break;
					HegSplNum[k] = jj;
					k = HegSplNxt[k];
					}
				ii = HegSplNxt[i];
				jj = HegSplNxt[ReadPres];
				HegSplNxt[ReadPres] = ii;
				HegSplNxt[i] = jj;
				Mergers ++;
				}			
			}        
        }
          
     return(i);                    
}

void Init_G_Variables()
{
    /******************************************************************************************
        This routine initializes some of the global variables used by Heegaard.
    ******************************************************************************************/
        
    register int    i,
                    j;
                    
    Band_Sums                       = 0L;
    BDY[0]                          = 2;
    BdryData                        = FALSE;
    BytesAvailable                  = BYTES_AVAILABLE;
    BytesUsed                       = 0L;
    Check_If_HS_Rep					= FALSE;
    CheckHSReps						= FALSE;
    Compute_Stabilizers             = FALSE;
    Count                           = 0;
    CouldNotRemove					= 0;
    CurrentComp                     = 1;
    Delete_Only_Short_Primitives    = FALSE;
    Do_Not_Reduce_Genus             = TRUE;
    DrawingDiagrams                 = FALSE;
    EmtyRel                         = FALSE;
    Find_All_Min_Pres               = FALSE;
    FormBandsums                    = FALSE;
    FoundPower                      = FALSE;
    FoundPrimitive                  = FALSE;
    From_BANDSUM                    = 0;
    From_DUALIZE                    = 0;
    Length                          = 0L;
    LensSpace                       = FALSE;
    Level_Interrupt                 = FALSE;
    Mergers							= 0;
    Micro_Print                     = FALSE;
    Micro_Print_F                   = FALSE;
    Minimum                         = 0L;
    Modified_Init_Pres              = FALSE;
    NoReport                        = TRUE;
    NumDiagrams                     = 0L;
    NumDualized                     = 0L;
    NumErrors						= 0L;
    NumFilled                       = 0;
    Num_Level_Transformations       = 0L;
    Num_Level_Slides				= 0l;
    NumRelTooLong					= 0;
    NumNonSepAnnuli					= 0;
    NumNotConnected					= 0;
    Num_Saved_LPres					= 0;
    NumSepAnnuli					= 0;
    NumTimes                        = 0;
    OnStack                         = 0;
    RandomizeSlides					= TRUE;
    ReadPres                        = 0;
    SaveMinima                      = FALSE;
    Saved_Vertices                  = 0;
    Start_Level_Search              = 0;
    Starting_Pres                   = 0;
    TestRealizability1				= FALSE;
    TestRealizability2				= FALSE;
    TestRealizability3				= FALSE;
    TestRealizability4				= FALSE;
    TotalAuts                       = 0L;
    TotalComp                       = 1;
    UserSaidQuit                    = FALSE;
    Vertices                        = 2*NumGenerators;
    
    for(i = 0; i < MAX_SAVED_PRES; i++)
        {
        SURL[i] = 0L;
        QPM[i] = 0;
        }
    for(i = 0; i <  MAXNUMCOMPONENTS; i++)
    for(j = 1; j <= MAXNUMGENERATORS; j++) MLC[i][j] = BIG_NUMBER;
    for(i = 0; i <  MAXNUMCOMPONENTS; i++) CBC[i][0] = BDRY_UNKNOWN;
    for(i = 0; i <= MAXNUMCOMPONENTS; i++) CS[i] = EOS;
    
}

void Delete_Old_Presentations(void)
{   
	unsigned int    i,
                    j;
                              
    for(i = 0; i < NumFilled; i++)
        {
        if(UDV[i] == ANNULUS_EXISTS || UDV[i] == V2_ANNULUS_EXISTS) UDV[i] = PRIM[i] = 0;
        for(j = 0; j <= NR[i]; j++) if(SUR[i][j] != NULL)
        	{
        	DisposeHandle((char **) SUR[i][j]);
        	SUR[i][j] = NULL;
        	}
        }
    NumFilled         = 0;
    BytesAvailable    = BYTES_AVAILABLE;
    BytesUsed         = 0L;
    UserSaidQuit      = FALSE;             
}

int ReRun_A_Presentation()
{
    unsigned char   *p,
                    *q,
                    *ptr = NULL,
                    revnum[10];
    
    unsigned int    PresNum;
    
    unsigned int    h,
                    i,
                    j,
                    k;

    DrawingDiagrams = FALSE;

    if(NoReport)
        {
        printf("\n\n     HIT 'v' TO REVIEW THESE PRESENTATIONS.");
        printf("\n     OR HIT ANY OTHER KEY TO CONTINUE.");        
        switch(WaitkbHit())
            {
            case 'v':
                REVIEW:
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,NULL);
                printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
                GET_RESPONSE2:                
                switch(WaitkbHit())
                    {
                    case 'y':
                        goto REVIEW;
                    case 'n':
                        break;
                    default:
                        if(Batch == FALSE) SysBeep(5);
                        goto GET_RESPONSE2;
                    }
                break;
            default:
                break;        
            }
        }
    
    if(Did_Cutting_Disk_Surgery)
        {
        printf("\n\nRETRIEVE THE LAST PRESENTATION USED FOR CUTTING DISK SURGERY ?  HIT 'y' OR 'n'.    ");
        GET_RESPONSE5:        
        switch(WaitkbHit())
            {
            case 'y':
                NumRelators = NumCuttingDiskSurgeryRel;
                for(i = 1; i <= NumRelators; i++)
                    {
                    if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
                    Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) CD_Surgery_Rel[i]));
                    if(Relators[i] == NULL) Mem_Error();
                    p = *CD_Surgery_Rel[i];
                    q = *Relators[i];    
                    while( (*q++ = *p++) ) ;
                    }
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators);
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                Input = RERUN;
                Did_Cutting_Disk_Surgery = FALSE;
                return(0);
            case 'n':
                Did_Cutting_Disk_Surgery = FALSE;
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE5;
            }    
        }    
    
    if(Did_Exponent_Surgery)
        {
        printf("\n\nRETRIEVE THE LAST PRESENTATION USED FOR EXPONENT SURGERY ?  HIT 'y' OR 'n'.    ");
        GET_RESPONSE4:      
        switch(WaitkbHit())
            {
            case 'y':
                NumRelators = NumExp_Sur_Rel;
                for(i = 1; i <= NumRelators; i++)
                    {
                    if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
                    Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Exp_Surgery_Rel[i]));
                    if(Relators[i] == NULL) Mem_Error();
                    p = *Exp_Surgery_Rel[i];
                    q = *Relators[i];    
                    while( (*q++ = *p++) ) ;
                    }
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators);
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                Input = RERUN;
                Did_Exponent_Surgery = FALSE;
                return(0);
            case 'n':
                Did_Exponent_Surgery = FALSE;
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE4;
            }    
        }
        
    if(Knot_Or_Link)
        {
        printf("\n\nRETRIEVE THE LAST WIRTINGER PRESENTATION OF A KNOT OR LINK ?  HIT 'y' OR 'n'.    ");
        GET_RESPONSE3:       
        switch(WaitkbHit())
            {
            case 'y':
                for(i = 1; i <= NumKnot_Or_Link_Rel; i++)
                    {
                    if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
                    Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) KorLRelators[i]));
                    p = *KorLRelators[i];
                    if(Relators[i] == NULL) Mem_Error();
                    q = *Relators[i];    
                    while( (*q++ = *p++) ) ;
                    }
                NumRelators = NumKnot_Or_Link_Rel;   
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators); 
                j = Wirtinger();                                 
                if(j && (j != 17))
                    {
                    printf("\nSurgery is not possible. Sorry!");
                    Knot_Or_Link = FALSE;
                    break;
                    }        
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                printf("\nThe surgered presentation is:\n");
                Print_Relators(Relators,NumRelators);
                Input = RERUN;
                return(0);        
            case 'n':
                Knot_Or_Link = FALSE;
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE3;
            }
        }

    ptr = (unsigned char *) NewPtr(100); 
    if(ptr == NULL) Mem_Error(); 
    printf("\n\nENTER A PRESENTATION FROM 0 TO %u THAT YOU WANT RERUN AND HIT 'return'.     ",NumFilled);
    for(i = j = 0; j < NumFilled; j++) if(SURL[j] == 0) i ++;
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
            for(h = 0,k = 1; h < NumFilled; h++) if(SURL[h] == 0)
                {
                h++;
                j += printf("{%d,",h);    
                break;
                }
            for( ; h < NumFilled; h++) if(SURL[h] == 0)
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
    ReadString((char *)ptr, GetPtrSize(ptr));
    sscanf((char *) ptr,"%d",&WhichInput);    
    if(WhichInput == 0)
        {
        printf("\n\nRerunning the original presentation.\n");
        NumRelators = CopyNumRelators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
            Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i])); 
            if(Relators[i] == NULL) Mem_Error();
            p = *Copy_Of_Input[i];
            q = *Relators[i];    
            while( (*q++ = *p++) ) ;
            }
        }
    else
        {
        if(WhichInput < 1 || WhichInput > NumFilled || SURL[WhichInput-1] == 0)
            {
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE1;
            }    
        printf("\n\nRerunning presentation %d.\n",WhichInput);
        PresNum = WhichInput;
        sprintf((char*)PresName,"A rerun of Presentation  ");
        p = PresName + 24;
        i = 0;
        while (PresNum)
            {
            revnum[i++] = PresNum - (PresNum/10)*10;
            PresNum /= 10;
            }    
        for ( ; i ;i--) *p++ = revnum[i-1] + '0';
        *p++ = '.';
        *p = EOS;
        WhichInput --;
        NumRelators = NR[WhichInput];
        for(i = 1; i <= NumRelators; i++)
            {
            if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
            Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) SUR[WhichInput][i]));
            if(Relators[i] == NULL) Mem_Error();
            p = *SUR[WhichInput][i];
            q = *Relators[i];              
            while( (*q++ = *p++) ) ;
            }
        }
    DisposePtr((char *) ptr); 
    Delete_Old_Presentations();
    NoReport = TRUE;
    
    Print_Relators(Relators,NumRelators);
    Input = RERUN;
    return(0);
}

int Get_Presentation_From_File()
{
    register unsigned char  *p,
                            *q,
                            t;
    
    unsigned int            h,
                            i,
                            j;
                            
    long                    StrLength;                                                
    
    if((input_relators = fopen("Input_Presentations","r+")) == NULL)
        {
        if(Batch == FALSE) SysBeep(5);
        printf("\nUnable to open the file 'Input_Presentations'.");
        printf("\nPlease locate the file 'Input_Presentations', make sure it is closed,");
        printf("\nand place it in the parent folder of the folder containing Heegaard.");
        return(1);
        }
     
    printf("\n\nPLEASE INDICATE WHICH PRESENTATION YOU WANT USED BY ENTERING ITS IDENTIFIER:\n\n    ");        
    ReadString((char *)PresName, GetPtrSize(PresName));

    /******************************************************************************************
                    Look for a presentation with the given identifier.
    ******************************************************************************************/    
    
    fseek(input_relators,0L,0);
    do
        {
        if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL)
            {
            printf("\n\n	Unable to find this presentation in the file 'Input_Presentations'.");
            printf("\nCheck that the identifier you entered appears in 'Input_Presentations':");
            printf("\n1) Preceded by the same number of spaces and tabs as appear in your identifier."); 
            printf("\n2) Followed in 'Input_Presentations' by a blank space.");
            if(Batch == FALSE) SysBeep(5);
            return(1);
            }
   		p = Inst;
        q = PresName;
        while(*p && *q && *p == *q)
            {
            p++;
            q++;
            }
        if(*p == '\n' || *p == ' ' || *p == '\t') *p = EOS;        
        }
    while(*p != *q);
    
    /******************************************************************************************
                             Look for the next non-empty line.
    ******************************************************************************************/
    
    do
        if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL)
            {
            printf("\n\nUnable to find any non-trivial relators in this presentation!");
            if(Batch == FALSE) SysBeep(5);
            return(1);
            }
    while(*Inst == '\n' || *Inst == '\r');

    /******************************************************************************************
        Read in the relators, at one relator to each non-empty line, stripping off leading
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
        printf("\n R%3u)   ",i); 
        j = 0;       
        while( (t = *p) )
            {
            if(t == '\n' || t == ' ' || t == '\t')
                {
                *p = EOS;
                break;
                }     
            putc(t,stdout);
            j++;
            if(!isalpha(t))
                {
                if(Batch == FALSE) SysBeep(5);
                printf("\nPlease check relator %2u. Char number %u is %c = %d.\n",i,j,t,t);
                printf("\nA relator can only contain upper and lower case letters!\n");
                printf("\nPlease make sure each relator is on its own line!");
                printf("\nAnd the presentation is terminated by an empty line!");
                return(1);
                }
            p++;        
            }                                
  		if(*Inst== '\n') break;
        StrLength = strlen((char *) Inst);
        if(StrLength >= MAXLENGTH)
            {
            if(Batch == FALSE) SysBeep(5);
            printf("\nRelator %d is too long!!",NumRelators + 1);
            return(1);
            }
        if(StrLength == h) break;
        NumRelators ++;
        if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
        Relators[i] = (unsigned char **) NewHandle(StrLength + 1 - h);
        if(Relators[i] == NULL) Mem_Error();
        LR[i] = StrLength - h;
        p = *Relators[i];    
  		q = Inst + h;
        while( (*p++ = *q++) ) ;
        if(fgets((char *) Inst,MAXLENGTH,input_relators) == NULL) break;         
        }
    if(NumRelators == MAXNUMRELATORS)
        printf("\n\nHeegaard will not accept any more than %d relators.",NumRelators);
    return(0);
}

int Get_Presentation_From_KeyBoard()
{
    unsigned char   *p,
                    *q,
                    *r,
                    t;
    
    int             Flag1;
    
    unsigned int    i;
                            
    unsigned long   StrLength;                                                
    
    printf("\n\nPLEASE ENTER THE RELATORS.\n");
    printf("    Note 1: Use 'Input_Presentations' to enter presentations having individual relators longer than 1024.\n");
    printf("    Note 2: Hit 'return' twice to terminate entry of the presentation. \n\n");
    for(i = 1,NumRelators = 0; i <= MAXNUMRELATORS; i++)
        {
        do
            {
            r = Inst;
            p = r;
            printf("R%3u)   ",i);
            ReadString((char *)Inst, GetPtrSize(Inst));
            t = *Inst;
            if(t == '\r' || t == '\n')
                {
                *p = *r = EOS;
                break;
                }

            Flag1 = FALSE;
            while( (t = *p++) )
                {
                if(!isalpha(t))
                    {
                    Flag1 = TRUE;
                    if(Batch == FALSE) SysBeep(5);
                    printf("\nA relator can contain only upper and lower case letters!\n");
                    if(t == ' ')
                        printf("\nNo blanks or trailing spaces!\n");
                    printf("\nPLEASE ENTER RELATOR %2u.\n",i);
                    break;
                    }
                }                                
            }
        while(Flag1);
        StrLength = strlen((char *) r);
        if(StrLength == 0) break;
        NumRelators ++;
        if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
        Relators[i] = (unsigned char **) NewHandle(StrLength + 1);
        if(Relators[i] == NULL) Mem_Error();
        LR[i] = StrLength;
        p = *Relators[i];    
        q = r;
        while( (*p++ = *q++) ) ;            
        }
    if(NumRelators == MAXNUMRELATORS)
        printf("\n\nHeegaard will not accept more than %d relators.",NumRelators);        
    printf("\n\nPlease enter a name by which Heegaard can refer to this presentation,");
    printf("\nand then hit 'return'.");
    printf("\n\nSAVE THIS PRESENTATION AS: ");
    ReadString((char *)PresName, GetPtrSize(PresName));
    return(0);
}

int Check_Realizability_Of_The_Initial_Presentation()
{
    Get_Simplification_Parameters_From_User(FALSE,FALSE);
    if(Get_Initial_Diagram(FALSE) == 2)
        NoReport = FALSE;
    else    
        {
        Report(0,NumDiagrams,0,0,0,0,1,0,1,NULL);
        }

	if(Batch != 3)   
    printf("\n\nDISCARD ALL BUT THE INITIAL PRESENTATION ?  HIT 'y' OR 'n'.");
    if(Batch == 3) return(0);
    GET_RESPONSE1:    
    switch(WaitkbHit())
        {
        case 'y':
            return(0);
        case 'n':
            return(1);
        default:
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE1;
        }        
}

void Display_Diagram_Of_The_Initial_Presentation(void)
{
    int             Flag1,
                    i,
                    SNumGenerators;
    
    unsigned long  	Scratch;
    
    DrawingDiagrams = TRUE;
    if(NumFilled)
        WhichInput = NumFilled - 1;
    else
        {
        WhichInput = 0;
        UDV[WhichInput] = 0;
        }
    Flag1 = FALSE;
    switch(Find_Flow_A(NORMAL,FALSE))
        {
        case 1:
            Flag1 = TRUE;
            SNumGenerators = NumGenerators;
            Rewrite_Input();
            break;
        case TOO_LONG:
            printf("\n\n     This presentation may be too long! Unable to display its diagram. Sorry!");
            DrawingDiagrams = FALSE;
            return;    
        }
    if(Length == 0)
        {
        printf("\n\nThe presentation has reduced to the empty presentation.");
        printf("\nThere is no diagram to display.");
        DrawingDiagrams = FALSE;
        return;
        }    
    if(Automorphisms)
        {
        if(Batch == FALSE) SysBeep(5);
        printf("\n\n                    NOTE!");
        printf("\n\nThe initial presentation does not have minimal length. The diagram you see will");
        printf("\ncorrespond to the following rewritten minimal length version of the initial presentation.");
        if(Flag1)
            {
            if(SNumGenerators - NumGenerators == 1)
                printf("\n    Note: a generator disappeared from the original presentation.");
            else
                printf("\n    Note: %d generators disappeared from the original presentation.",
                SNumGenerators - NumGenerators);
            printf("\n          The missing generator(s) will not appear in the Heegaard diagram.");    
            }  
        Vertices = 2*NumGenerators;      
        Canonical_Rewrite(Relators,FALSE,FALSE);
        Fill_A(NumRelators);
        Saved_Vertices = 0;
        printf("\n\nThe minimal length presentation is:\n");
        Print_Relators(Relators,NumRelators);   
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            Scratch += GetHandleSize((char **) Relators[i]);
        Scratch -= NumRelators;
        printf("\nThis presentation has length %ld ",Scratch);
        if(Batch == 0)
        	{
        	printf("\n\n    HIT ANY KEY TO SEE DATA FOR THE HEEGAARD DIAGRAM.");
        	WaitkbHit();
        	} 
        }
    Get_Matrix();    
    ComputeValences_A();  
    Check_Connected();
    SepPairs = Sep_Pairs(0,0,1);
    if(SepPairs == TOO_LONG)
        {
        printf("\n\n     This presentation may be too long! Unable to display its diagram. Sorry!");
        DrawingDiagrams = FALSE;
        return;
        }    
    if(SepPairs)
        {
        if(V1 & 1)
            LSP[WhichInput] = V1/2 + 97;
        else
            LSP[WhichInput] = V1/2 + 65;
        if(V2 & 1)
            LSQ[WhichInput] = V2/2 + 97;
        else
            LSQ[WhichInput] = V2/2 + 65;
        if(UDV[WhichInput] <= DONE)    
            UDV[WhichInput] = SEP_PAIRS;            
        }
    NonPlanar = Planar(FALSE,TRUE);
    Print_Graph(FALSE,0,0,0);
    DrawingDiagrams = FALSE;
}

void Get_Simplification_Parameters_From_User(int Flag1,int Flag2)
{
	if(Check_If_HS_Rep)
    	{
		FormBandsums 					= TRUE;
		OnlyReducingBandsums 			= FALSE; 
		Delete_Only_Short_Primitives 	= FALSE;
        Do_Not_Reduce_Genus 			= TRUE;
        return; 	
    	}
    if(Flag2 == FALSE) OnlyReducingBandsums = FormBandsums = FALSE;
    if(Flag1)
        {
        if(BreadthFirstSearch)
        	printf("\nHeegaard is using Breadth-First search.\n");
        if(DepthFirstSearch)
          	printf("\nHeegaard is using Depth-First search.\n");      
        if(FormBandsums)
            {
            if(OnlyReducingBandsums)
                printf("\nHeegaard is currently forming only length reducing bandsums.\n");            
            else
                printf("\nHeegaard is currently forming all possible bandsums.\n");
            }
        else
            printf("\nHeegaard is not currently forming any bandsums.\n");
        }
    if(Flag2)
        {
        printf("\nCREATE NEW DIAGRAMS BY FORMING BANDSUMS ? HIT 'y' OR 'n'.\n");
        GET_RESPONSE1:        
        switch(WaitkbHit())
            {
            case 'y':
                printf("\n    HIT 'a' TO FORM ALL POSSIBLE BANDSUMS.");
                printf("\n    HIT 'r' TO FORM ONLY LENGTH REDUCING BANDSUMS.");
                printf("\n    HIT 'n' TO FORM NO BANDSUMS.    ");
                GET_RESPONSE2:                
                switch(WaitkbHit())
                    {
                    case 'a':
                        FormBandsums 			= TRUE;
                        OnlyReducingBandsums 	= FALSE;
                        break;
                    case 'r':
                        FormBandsums 			= TRUE;
                        OnlyReducingBandsums 	= TRUE;
                        break;
                    case 'n':
                        OnlyReducingBandsums = FormBandsums = FALSE;
                        break;    
                    default:
                        if(Batch == FALSE) SysBeep(5);
                        goto GET_RESPONSE2;
                    }
                printf("\n");    
                break;        
            case 'n':
                OnlyReducingBandsums = FormBandsums = FALSE;                
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE1;    
            }    
        }
    if(Flag1)
        {
        if(!Do_Not_Reduce_Genus && !Delete_Only_Short_Primitives)
            {
            printf("\nHeegaard is currently trying to delete all primitive relators.\n");
            printf("\nCONTINUE TO DELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
            }
        else
            {
            if(Delete_Only_Short_Primitives && !Do_Not_Reduce_Genus)
            printf("\nHeegaard is currently deleting only primitive relators of length 1 and 2.\n");
            else
            printf("\nHeegaard is not currently attempting to delete any primitive relators.\n");
            printf("\nCHANGE TO DELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
            }
        }
    else
    	if(Batch != 3)
    		{
			printf("\nDELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
  
			printf("\n***************************************************************************");
			printf("\n        Be somewhat careful when choosing to delete primitives!");
			printf("\n    Although deleting primitives always preserves fundamental groups,");
			printf("\nsimplifying a presentation P, whose realizability has not been verified, by");
			printf("\ndeleting primitives, may falsely suggest P is realizable. For example,");
			printf("\ndeleting C from P = < A,B,C | AAACBBC, C > yields P' = < A,B | AAABB >,");
			printf("\nand P' is realizable, but P is not.");
			printf("\n***************************************************************************\n");
			}
    if(Batch == 3)
    	{
 		Delete_Only_Short_Primitives 	= FALSE;
        Do_Not_Reduce_Genus 			= TRUE;
        return;   	
    	}
    GET_RESPONSE3:    
    switch(WaitkbHit())
        {
        case 'y':
            Do_Not_Reduce_Genus 			= FALSE;
            Delete_Only_Short_Primitives 	= FALSE;
            break;
        case 'n':    
            printf("\nDELETE PRIMITIVE RELATORS OF LENGTH 1 AND LENGTH 2 ? HIT 'y' OR 'n'.\n");
            GET_RESPONSE4:            
            switch(WaitkbHit())
                {
                case 'y':
                    Delete_Only_Short_Primitives 	= TRUE;
                    Do_Not_Reduce_Genus 			= FALSE;
                    break;
                case 'n':
                    Delete_Only_Short_Primitives 	= FALSE;
                    Do_Not_Reduce_Genus 			= TRUE;
                    break;
                default:
                    if(Batch == FALSE) SysBeep(5);
                    goto GET_RESPONSE4;
                }
            break;
        default:
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE3;
        }
        
    if(!Flag1) Turn_Micro_Print_On();
    
    if(Flag1)
        {
        if(Find_All_Min_Pres)
            {
            printf("\nHeegaard is currently trying to make broad use of level-transformations.\n");
            printf("\nCONTINUE TO MAKE BROAD USE OF LEVEL-TRANSFORMATIONS ? HIT 'y' OR 'n'.\n");
            }
        else
            {
            printf("\nHeegaard is currently trying to make only limited use of level-transformations.\n");
            printf("\nCHANGE STRATEGY TO MAKE BROAD USE OF LEVEL-TRANSFORMATIONS ? HIT 'y' OR 'n'.\n");
            }
        GET_RESPONSE5:       
        switch(WaitkbHit())
            {
            case 'y':
                Find_All_Min_Pres = TRUE;
                break;
            case 'n':
                Find_All_Min_Pres = FALSE;
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE5;
            }
    	}		   
}

void Turn_Micro_Print_On(void)
{
	if(Micro_Print)
		{
		printf("\nMicro_Printing IS ON. TURN Micro_Printing OFF ? HIT 'y' OR 'n'.\n");
		GET_RESPONSE1:		
		switch(WaitkbHit())
			{
			case 'y':	
				Micro_Print 	= FALSE;
				Micro_Print_F 	= FALSE;
				break;
			case 'n':
				Micro_Print 	= TRUE;
				Micro_Print_F 	= TRUE;
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE1;
			}
		}
	else
		{
		printf("\nTURN Micro_Printing ON ? HIT 'y' OR 'n'.\n");
		GET_RESPONSE2:		
		switch(WaitkbHit())
			{
			case 'y':
				Micro_Print 	= TRUE;
				Micro_Print_F 	= TRUE;			
				break;
			case 'n':
				Micro_Print 	= FALSE;
				Micro_Print_F 	= FALSE;
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE2;
			}
		}	
}

int Find_Level_Transformations_Of_The_Initial_Presentation()
{
	int		FLTRV;
	
    if(Length == 0)
        {
        if(Batch != 7) printf("\n");
        printf("\nThis is an empty presentation!!");
        return(1);
        }
    ReadPres = NumFilled;    
    if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG)
        {
        if(Batch != 7) printf("\n");
        printf("\nPresentation %d is too long!",ReadPres + 1);
        return(2);
        }
    if(Automorphisms)
        {
        if(Batch != 7) printf("\n");
        printf("\nPresentation %d does not have minimal length.",ReadPres + 1);
        printf("\nHeegaard will only find level transformations for presentations ");
        printf("that have minimal length.");
        return(3);
        }
    if(Batch != 7) 
    	{     
		printf("\nPRINT STABILIZER INFO ? HIT 'y' OR 'n'.");
		GET_RESPONSE1:   
		switch(WaitkbHit())
			{
			case 'y':
				Compute_Stabilizers = TRUE;
				break;
			case 'n':
				Compute_Stabilizers = FALSE;
				break;
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE1;        
			}
		printf("\n\nListing members of the orbit of presentation 1 under level transformations. . .\n");
    	}
    if(Batch == 7) Compute_Stabilizers = FALSE;	
    Num_Level_Transformations = 0L;
    Left[0] = Right[0] = INFINITE;    
    if(Save_Pres(ReadPres,0,Length,1,2,1,0,0))
        {
        Compute_Stabilizers = FALSE;
        return(0);
        }
    if(Init_Find_Level_Transformations(TRUE))
        {
        Compute_Stabilizers = FALSE;
        return(0);
        }
    do
        {
        OnStack --;
        FLTRV = Find_Level_Transformations(FALSE,TRUE);
        if(FLTRV == 5 || FLTRV == FULL_HOUSE)
            {
            if(Batch != 7) printf("\n\nWe have run out of memory set aside for presentations!!!!");
            Compute_Stabilizers = FALSE;
            break;
            }
        ReadPres ++;
        }
    while(ReadPres < NumFilled);
    
    if(Batch != 7) SysBeep(5); 
      
    Compute_Stabilizers = FALSE;
    if(ReadPres < NumFilled && Batch != 7) printf("\n\nThe search for level transformations was interrupted.");
    if(Batch != 7)
    	{
		printf("\nHeegaard performed %lu level-transformation",Num_Level_Transformations);
		if(Num_Level_Transformations == 1)
			printf(".");
		else
			printf("s.");
		printf("\n\nThe orbit of presentation 1 under level transformations has ");
		switch(NumFilled - 1)
			{
			case 0:
				printf("1 member.\n");
				break;
			default:
				if(NumFilled >= MAX_SAVED_PRES -3)
					printf("at least %u members.\n",NumFilled);
				else
					printf("%u members.\n",NumFilled);
				break;
			}
        }   
    
    if(Batch == 7) 
    	{
    	printf("\n|LTs| %lu, ",Num_Level_Transformations);
    	if(H_Results != NULL) 
    		fprintf(H_Results,"\n\n%s\n|LTs| %lu, ",PresName,Num_Level_Transformations);
    	if(NumFilled >= MAX_SAVED_PRES -3)
    		{
    		printf("|Orbit| >= %u, ",NumFilled);
    		if(H_Results != NULL) fprintf(H_Results,"|Orbit| >= %u, ",NumFilled);
    		}
    	else
    		{
    		printf("|Orbit| %u, ",NumFilled);
    		if(H_Results != NULL) fprintf(H_Results,"|Orbit| %u, ",NumFilled);
    		}
    	return(0);
    	}
        
    if(Batch != 7)
    	{
    	printf("\n\nHIT 'v' TO REVIEW THE PRESENTATIONS IN THE ORBIT.");
    	if(NumFilled > 1)
			{
			printf("\nHIT 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY BY SUMMAND NUMBER,");
			printf("\n        NUMGENERATORS, NUMRELATORS, LENGTH AND 'LEXICOGRAPHIC' ORDER.");
			}
    	printf("\nHIT ANY OTHER KEY TO CONTINUE.");
    	}   
    switch(WaitkbHit())
        {
        case 'v':
            REVIEW:
            Report(0,0,0,0,0,0,0,0,1,NULL);
            printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
            GET_RESPONSE4:            
            switch(WaitkbHit())
                {
                case 'y':
                    goto REVIEW;
                case 'n':
                    break;
                default:
                    if(Batch == FALSE) SysBeep(5);
                    goto GET_RESPONSE4;
                }
            break;
        case 'w':
            printf("\n\n     Sorting presentations. . . .");
            Sort_Presentations_In_Memory(FALSE);
        	break;
        default:
            break;    
        }
        
    printf("\n\nDISCARD THE PRESENTATIONS IN THE ORBIT ?  HIT 'y' OR 'n'.");        
    GET_RESPONSE5:   
    switch(WaitkbHit())
        {
        case 'y':
            return(0);
        case 'n':
            printf("\n\nTEST THESE PRESENTATIONS FOR PSEUDO-MINIMALITY, SEPARATING VERTICES AND REALIZABILITY ?");
            printf("\nHIT 'y' OR 'n'.");
            GET_RESPONSE6:            
            switch(WaitkbHit())
                {
                case 'y':
                    Test_LT_For_Pseudo_Min();
                case 'n':
                    return(1);
                default:
                    if(Batch == FALSE) SysBeep(5);
                    goto GET_RESPONSE6;
                }
        default:
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE5;
        }    
}                        

int Reduce_The_Initial_Presentation_To_Minimal_Length()
{
    int   	i;
    
    long    Scratch;
	if(Batch != 16 && Batch != 53) Turn_Micro_Print_On();
    Scratch = Length;
    i = NumGenerators;
    switch(Find_Flow_A(NORMAL,FALSE))
        {
        case 1:
            if(Batch != 53) Rewrite_Input();
            break;
        case TOO_LONG:
            printf("\n\n     This presentation may be too long for Heegaard to handle. Sorry!");
            return(1);    
        }        
    if(Length == Scratch)
    	{
        printf("\n\nThis presentation has minimal length.");
        if(Batch == 16 && H_Results != NULL)
        	{				
			fprintf(H_Results,"\n\n%s",PresName);
			Print_Relators2(Relators,NumRelators);        	
        	}
        if(Batch == 53) return(0);	
        }
    else
        {
        printf("\n\n%lu automorphism(s) reduced the length from %ld to %lu.\n",Automorphisms,Scratch,Length);
        if(i > NumGenerators)
            printf("\nand reduced the number of generators from %d to %d.\n",i,NumGenerators);    
        Canonical_Rewrite(Relators,FALSE,FALSE);
		if(Batch == 16)
			{
			Print_Relators(Relators,NumRelators);
			if(H_Results != NULL) 
				{
				fprintf(H_Results,"\n\n%s",PresName);
				Print_Relators2(Relators,NumRelators);
				}
			return(0); 
			} 
		if(Batch == 53 && H_Results != NULL)
			{
			fprintf(H_Results,"\n\n%s <-- Not a HS Rep! (IP does not have minimal length.)",PresName);
			return(1);
			}	    
        if(Save_Pres(ReadPres,0,Length,1,2,1,0,0)) return(0);
        UDV[NumFilled - 1] = 0;
        Report(0,0,0,0,0,0,0,0,1,NULL);
        printf("\n\nDISCARD THIS MINIMAL LENGTH PRESENTATION ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:        
        switch(WaitkbHit())
            {
            case 'y':
                return(0);
            case 'n':
                return(1);
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE1;
            }
        }
    return(0);        
}

int Initial_Realizability_Check()
{
    int     i;

    /*****************************************************************************************
        If the Whitehead graph of a presentation P is connected, and has no cut-vertices,
        then the Whitehead graph of P must be planar if P is realizable. This routine checks
        whether the Whitehead graph of P is connected, has no cut-vertices, and is planar.
        The routine returns TRUE if it determines that P is not realizable, and otherwise
        returns FALSE.
    *****************************************************************************************/
            
    Fill_A(NumRelators);
    Saved_Vertices = 0;
    Get_Matrix();
    for(i = 0; i < Vertices; i++) ZZ[i] = 0;
    if(Connected_(0,0) == FALSE) return(FALSE);
    if(Find_Cut_Vertices()) return(FALSE);
    if(Planar(TRUE,FALSE))
        {        
        if(Save_Pres(ReadPres,0,Length,1,2,1,0,0)) return(TRUE);
        ReadPres = NumFilled - 1;
        Fatal_Error();
        printf("\n\n     The Whitehead graph is connected, has no cut-vertices and is non-planar.");
        printf("\n     This is impossible if the presentation is realizable.");        
        return(TRUE);
        }
    return(FALSE);    
}

void Print_Realizability(int Del_Only_Triv_Rel, unsigned int WhichPres)
{
    printf("\n\n                    The data appears consistent.");
    if(Del_Only_Triv_Rel)
    	{
        printf("\n\n                    The initial presentation is realizable.");
        if(Num_Saved_LPres && (Batch == 4 || Batch == 10 || Batch == 11 || Batch == 53))
			{
			printf("\n\nAnd is now:");
			Print_Relators(Relators,NumRelators);
			}
        }
    else
    	{
        printf("\n\n                    Presentation %d is realizable.",WhichPres);
        if(Batch == 4 || Batch == 10 || Batch == 11 || Batch == 53)
			{
			printf("\n\nPresentation %d is:",WhichPres);
			Print_Relators(Relators,NumRelators);
			}
        }
    if(Batch == 3) 
    	{
    	NumRealizable ++;
    	if(H_Results != NULL) fprintf(H_Results,"\n\n%s\nRealizable.",PresName);
    	}   
}

void Realization_Warning()
{
    printf("\n\n    NOTE: Heegaard has not directly verified that the initial presentation is realizable.");
    printf("\n          This means it is possible that the initial presentation is not realizable,");
    printf("\n          even though 'derived' presentations are realizable.");  
}

char mykbhit()
{
#ifdef MAC
    char            	charhit;
    EventRecord        theEvent;
    Boolean            gotEvent;
    
    gotEvent = WaitNextEvent(keyDownMask + autoKeyMask,&theEvent,180,NULL);
    if(gotEvent)
        {
        charhit = theEvent.message & charCodeMask;
        return(charhit);
        }
    else
        return(EOS);
#else
    int c = 0;
    struct termios oldtermio, newtermio;
    tcgetattr(0, &oldtermio);
    newtermio = oldtermio;
    newtermio.c_cc[VTIME]=0;
    newtermio.c_cc[VMIN]=0;
    tcsetattr(0, 0, &newtermio);
    c = getchar();
    tcsetattr(0, 0, &oldtermio);
    if (c == EOF){
      c = 0;
      clearerr(stdin);
    }
    return (char)c;
#endif
}

char WaitkbHit()
{
#ifdef MAC
    char            charhit;
    EventRecord        theEvent;
    
    FlushEvents(everyEvent,0);
    while(!WaitNextEvent(keyDownMask + autoKeyMask,&theEvent,10,NULL))        ;
    charhit = theEvent.message & charCodeMask;
    return(charhit);
#else
    int c;
    c = getchar();
    if (c == EOF){
      c = 0;
      clearerr(stdin);
    }
    return (char)c;
#endif
}
