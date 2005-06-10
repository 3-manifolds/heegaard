#include "Heegaard.h"
#include <ctype.h>
#include <string.h>

#ifndef MAC
#include <termios.h>
#include <stdio.h>
struct termios normal_termio;
#endif

FILE 
    *fptr,            
    *myout,
    *input_relators;

char
    *ER;
    
unsigned char
    *CBC[MAXNUMCOMPONENTS],
    *BCF,
    *BCWG,
    *CO[VERTICES],
    **Copy_Of_Input[MAXNUMRELATORS + 1],
    **Copy_Of_Rel_1[MAXNUMRELATORS + 1],
    **Copy_Of_Rel_2[MAXNUMRELATORS + 1],    
    *CS,
    **CD_Surgery_Rel[MAXNUMRELATORS + 1],
    *DeletedEdgePtr,
    *DeletedEdges,
    **DelRelators[MAXNUMRELATORS + 1],
    **DualRelators[MAXNUMRELATORS + 1],
    **Exp_Surgery_Rel[MAXNUMRELATORS + 1],    
    *Face[2*VERTICES],    
    *FV,
    *GBC,
    *Inst,
    **KorLRelators[MAXNUMRELATORS + 1],
    *N1H,
    *NCS,
    *NS1XD2,
    *NS1XS2,
    **OutRelators[MAXNUMRELATORS + 1],
    *PresName,
    *QPM,            
    **Relators[MAXNUMRELATORS + 1],
    **TopOfChain[MAXNUMRELATORS + 1],
    **RWR[MAXNUMDUPS],
    ***SLR[MAX_SAVED_LEVELS],
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
    *TP;

int 
    *BDY,
    BdryData,
    Boundary,
    Compute_Stabilizers,
    Connected,
    CopyNumGenerators,
    CopyNumRelators,
    Count,
    CurrentComp,
    CycleDiagrams,
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
    GoingUp,
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
    NumRelators,
    NumSepComps,
    NumTimes,
    OnlyReducingBandsums,
    *PRIM,
    ReadPres,
    *SaveBdry,
    Save_Init_Pres,
    SaveMinima,
    Saved_Vertices,
    SepPairs,
    SRError,
    SReadPres,
    TestRealizability1,
    TestRealizability2,
    TestRealizability3,    
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
    *BSV1,
    *BSV2,
    *ComponentNum,
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
    *InDisk,
    *InQueue,
    *IV,
    *Left,
    LensSpace,    
    *Lowpt,
    MaxLength,
    *NEBC,    
    *NEX[(VERTICES)/2],
    *NFBC,    
    NotNewPres,
    *NRBC,    
    *Number,
    NumCalled,
    NumFilled,
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
    *SV,
    This_Pres,
    *TV,
    *UDV,
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
    Length,
    *LR,
    *LSP,
    *LSQ,
    *MLC[MAXNUMCOMPONENTS],
    NumDualized,
    Num_Level_Transformations,
    OrigLength,
    P,
    Q,
    SLength,
    *SURL,
    TOCLength,
    TotalAuts;

main()
{
    char            c;
        
    unsigned char    *p,
                     *q;
                            
    unsigned int     i;
    
    long             Scratch;

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
    
    if(i = Do_Initialization())
        {
        printf("\n\nThe program was unable to allocate memory for all of its data structures.");
        printf("\nTry increasing the amount of memory which the program can use.");
        printf("\nFailure at allocation number %u.",i);
#ifndef MAC
	tcsetattr(0, 0, &normal_termio);
#endif
        return(0);
        }
    
    Input = INITIAL_PRES;
    
_BEGIN:
    if(Input == BEGIN)
        {
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;
        printf("\n");
        #ifdef DEBUGGING
            printf("\nHIT 'b' TO PRINT DEBUGGING INFO INTO THE FILE SIMPLIFY_H_ERRORS.");
        #endif
        printf("\nHIT 'd' TO SEE THE HEEGAARD DIAGRAMS.");
        if(NoReport)
            printf("\nHIT 'F' TO SAVE THE PRESENTATIONS NOW IN MEMORY TO THE FILE 'Heegaard_Results'.");
        printf("\nHIT 'i' TO SEE THE FILE 'Input_Presentations'.");
        printf("\nHIT 'n' TO ENTER A NEW PRESENTATION.");
        printf("\nHIT 'o' TO SEE THE FILE 'Heegaard_Results'.");
        printf("\nHIT 'p' TO FIND CANCELLATION PATHS.");        
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
	printf("\n");
        if(Knot_Or_Link || Did_Exponent_Surgery || Did_Cutting_Disk_Surgery)
            printf("\nHIT 'r' TO TRY ANOTHER SURGERY ON A PREVIOUS PRESENTATION, OR TO RERUN A PRESENTATION.");        
        else
            printf("\nHIT 'r' TO RERUN A PRESENTATION.");
        printf("\nHIT 's' TO LOOK FOR SYMMETRIES.");
        printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS NOW IN MEMORY.");
        if(NumFilled > 1)
            {
            printf("\nHIT 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY BY SUMMAND NUMBER,");
            printf("\n        NUMGENERATORS, NUMRELATORS, LENGTH AND 'LEXICOGRAPHIC' ORDER.");
            }
        printf("\nHIT 'x' TO CLOSE 'Heegaard_Results' FOR 'EXTERIOR' EDITING.");    
        printf("\nHIT '?' FOR HELP.");
	printf("\n");
        GET_RESPONSE1:
        c = WaitkbHit();
        switch(c)
            {
            case 'b':
                #ifdef DEBUGGING
                    DrawingDiagrams = FALSE;
                    fptr = myout;
                    Debug();
                    Input = BEGIN;
                    goto _BEGIN;
                #else
                    SysBeep(5);
                    DrawingDiagrams = FALSE;
                    Input = BEGIN;
                    goto GET_RESPONSE1;                    
                #endif
                
            case 'd':
                Display_Diagrams();
                Input = BEGIN;
                goto _BEGIN;
                        
            case 'F':
                if(NoReport)
                    Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,1,1,1,0,1,0);
                Input = BEGIN;
                goto _BEGIN;
                        
            case 'i':
                Display_File_Input_Presentations();
                Input = BEGIN;
                goto _BEGIN;
                    
            case 'n': 
                Input = INITIAL_PRES;
                DrawingDiagrams = FALSE;
                Did_Exponent_Surgery = FALSE;
                Did_Cutting_Disk_Surgery = FALSE;
                Delete_Old_Presentations();
                break;
                
            case 'o':
                Display_File_Simplify_Heegaard();
                Input = BEGIN;
                goto _BEGIN;

            case 'p':
                Find_Cancellation_Paths();
                goto _BEGIN;
                            
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
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,0);
                Input = BEGIN;
                goto _BEGIN;
                
            case 'w':
                printf("\n\n     Sorting presentations. . . .");
                Sort_Presentations_In_Memory();
                Input = BEGIN;
                goto _BEGIN;
                
            case 'x':
                Edit_MyOut();
                Input = BEGIN;
                goto _BEGIN;
                        
            case '?':
                Display_Help_File();
                Input = BEGIN;
                goto _BEGIN;
                                            
            default:
                SysBeep(5);
                DrawingDiagrams = FALSE;
                Input = BEGIN;
                goto GET_RESPONSE1;                                
            }            
        }    
    if(Input == INITIAL_PRES)
        {
        PRINT_INPUT_OPTIONS:            
        printf("\n\nHIT 'f' IF THE PRESENTATION WILL COME FROM THE FILE 'Input_Presentations'.");
        printf("\nHIT 'r' TO REVIEW THE FILE 'Input_Presentations'.");
        printf("\nHIT 'k' IF THE PRESENTATION WILL BE ENTERED FROM THE KEYBOARD.");
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
        printf("\nHIT '?' FOR HELP.");
	printf("\n");
        GET_RESPONSE2:
        c = WaitkbHit();
        switch(c)
            {
            case 'f':
                if(Get_Presentation_From_File()) goto PRINT_INPUT_OPTIONS;
                break;
            case 'r':
                Display_File_Input_Presentations();
                goto PRINT_INPUT_OPTIONS;        
            case 'k':
                if(Get_Presentation_From_KeyBoard()) goto PRINT_INPUT_OPTIONS;
                break;
            case 'q':
                goto _QUIT;
            case '?':
                Display_Help_File();
                goto PRINT_INPUT_OPTIONS;    
            default:
                SysBeep(5);
                goto GET_RESPONSE2;
            }
            
        /**************************************************************************************
            Call Wirtinger() and check whether this presentation seems to meet enough criteria
            to be the Wirtinger presentation of a knot or link and, if so, give the user the
            option of performing Dehn-fillings.
        **************************************************************************************/
                
        if(Wirtinger()) Knot_Or_Link = FALSE;
        if(Knot_Or_Link)
            {
            printf("\n\nThe surgered presentation is:\n");
            Print_Relators(Relators,NumRelators,stdout);            
            }                    
        }
    if(Input == RERUN) Input = INITIAL_PRES;        
    if(Input == INITIAL_PRES)
        {
        /**************************************************************************************
            Echo the initial relators to the output so we will have a copy of them. Then call
            Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
            which serves as the initial presentation for the program. 
        **************************************************************************************/    
        
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;
        ObscureCursor();
        fprintf(myout,"\n\nThe initial presentation was:");
        fprintf(myout," %s\n",PresName);
        Print_Relators(Relators,NumRelators,myout);    
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            Scratch += GetHandleSize((char **) Relators[i]);
        Scratch -= NumRelators;
        printf("\n\nThis presentation has length %ld ",Scratch);
        fprintf(myout,"\n\nThis presentation has length %ld ",Scratch);
        if(Freely_Reduce() == TOO_LONG)
            {
            printf("\n\nThis presentation is too long!!");
            SysBeep(5);
            Input = BEGIN;
            goto _BEGIN;
            }
        if(Scratch > OrigLength)
            {
            printf("\nand freely reduces to the following presentation of length %lu.\n",
                OrigLength);
            fprintf(myout,"\nand freely reduces to the following presentation of length %lu.\n",
                OrigLength);
            Print_Relators(Relators,NumRelators,stdout);
            Print_Relators(Relators,NumRelators,myout);    
            Scratch = OrigLength;
            }
        else
            {
            printf("and is freely reduced.");
            fprintf(myout,"and is freely reduced.");
            }
        Micro_Print = TRUE;
        Micro_Print_F = TRUE;    
        if(Rewrite_Input())
            {
            printf("\n\nThere must be at least one non-empty relator!!");
            SysBeep(5);
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
        CopyNumGenerators     = NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            ReallocateHandle((char **) Copy_Of_Input[i],GetHandleSize((char **) Relators[i]));                
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while(*p++ = *q++) ;
            }
                                        
        /**************************************************************************************
                Call Init_G_Variables() to initialize some global variables. Then call
                Canonical_Rewrite() to rewrite the presentation in canonical form.
        **************************************************************************************/
        
        Init_G_Variables();
        Length = Scratch;
        printf("\n\nNumRelators = %d, NumGenerators = %d\n",NumRelators, NumGenerators);
        Micro_Print = TRUE;
        Micro_Print_F = TRUE;    
        if(Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG)
            {
            printf("\n\nThis presentation has too many symmetries!!");
            SysBeep(5);
            Micro_Print = FALSE;
            Micro_Print_F = FALSE;
        /*    Input = BEGIN;
            goto _BEGIN;    */    /* Overrode the error message from Canonical_Rewrite(). */
            }
        Micro_Print = FALSE;
        Micro_Print_F = FALSE;                            
        if(Compare_Input_Pres() == FALSE)
            {
            printf("\n\nThe rewritten initial presentation is:\n");
            Print_Relators(Relators,NumRelators,stdout);
            fprintf(myout,"\n\nThe rewritten initial presentation is:\n");
            Print_Relators(Relators,NumRelators,myout);
            }
        
        /**************************************************************************************
                        Update the copy of the initial set of relators.
        **************************************************************************************/
        
        CopyNumRelators     = NumRelators;
        CopyNumGenerators     = NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            ReallocateHandle((char **) Copy_Of_Input[i],GetHandleSize((char **) Relators[i]));                
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while(*p++ = *q++) ;
            }
    
        /**************************************************************************************
                                Present the user with some options.
        **************************************************************************************/
                                        
    _OPTIONS:
        printf("\n");
        printf("\nHIT 'c' TO CHECK REALIZABILITY OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'd' TO SEE THE HEEGAARD DIAGRAM OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'e' TO TRY EXPONENT SURGERY ON THE INITIAL PRESENTATION.");
        printf("\nHIT 'g' TO GENERATE GENUS TWO DIAGRAMS.");            
        printf("\nHIT 'h' TO FIND THE INTEGRAL FIRST HOMOLOGY OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'i' TO SEE THE FILE 'Input_Presentations'.");
        printf("\nHIT 'k' TO TRY CUTTING DISK SURGERY ON THE INITIAL PRESENTATION.");
        printf("\nHIT 'l' TO FIND LEVEL TRANSFORMATIONS OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'm' TO LOOK FOR ALL MINIMAL PRESENTATIONS.");
        printf("\nHIT 'n' TO ENTER A NEW PRESENTATION.");
        printf("\nHIT 'q' TO QUIT RUNNING THE PROGRAM.");
        printf("\nHIT 'r' TO REDUCE AND SIMPLIFY THE INITIAL PRESENTATION.");
        printf("\nHIT 's' TO FIND SYMMETRIES OF THE INITIAL PRESENTATION.");
        printf("\nHIT 'v' TO REVIEW THE INITIAL PRESENTATION.");
        printf("\nHIT 'x' TO CLOSE 'Heegaard_Results' FOR 'EXTERIOR' EDITING.");
        printf("\nHIT 'z' TO REDUCE THE INITIAL PRESENTATION TO MINIMAL LENGTH.");
        printf("\nHIT '?' FOR HELP.\n");
        GET_RESPONSE3:
        c = WaitkbHit();
        switch(c)
            {
            case 'c':
                if(Check_Realizability_Of_The_Initial_Presentation())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                break;
                
            case 'd':
                Display_Diagram_Of_The_Initial_Presentation();
                break;
                
            case 'e':
                if(Try_Exponent_Surgery()) break;
                Did_Exponent_Surgery = TRUE;
                goto _OPTIONS;
            
            case 'g':
                Generate_Orbits_Under_Auts();
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                    
            case 'h':
                printf("\n\nComputing the integral first homology of the initial presentation . . .");
                Compute_Homology();
                break;
            
            case 'i':
                Display_File_Input_Presentations();
                goto _OPTIONS;
            
            case 'k':
                if(Try_Cutting_Disk_Surgery()) break;
                goto _OPTIONS;
                
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
                Find_All_Min_Pres = TRUE;
                goto _GET_INITIAL_DIAGRAM;    
            
            case 'n':
                Input = INITIAL_PRES;
                goto _BEGIN;
                            
            case 'q':
                goto _QUIT;    
            
            case 'r':
                Get_Simplification_Parameters_From_User(FALSE,TRUE);
                Find_All_Min_Pres = FALSE;                
                goto _GET_INITIAL_DIAGRAM;
            
            case 's':
                printf("\n\n    This is the initial presentation:\n");
                Print_Relators(Relators,NumRelators,stdout);
                fprintf(myout,"\n\n   This was the initial presentation:\n");
                Print_Relators(Relators,NumRelators,myout);
                i = NoReport;
                NoReport = FALSE;
                Find_Symmetries(TRUE);
                NoReport = i;
                break;
            
            case 'v':
                printf("\n\n    This is the initial presentation:\n");
                Print_Relators(Relators,NumRelators,stdout);    
                goto _OPTIONS;
            
            case 'x':
                Edit_MyOut();
                goto _OPTIONS;
                    
            case 'z':
                if(Reduce_The_Initial_Presentation_To_Minimal_Length())
                    {
                    Input = BEGIN;
                    goto _BEGIN;
                    }
                break;    
            
            case '?':
                Display_Help_File();
                goto _OPTIONS;
                
            default:
                SysBeep(5);
                goto GET_RESPONSE3;
            }
        
        NumRelators     = CopyNumRelators;
        NumGenerators     = CopyNumGenerators;
        Vertices         = 2*NumGenerators;
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            {
            ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));
            Scratch += GetHandleSize((char **) Copy_Of_Input[i]);                
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while(*q++ = *p++) ;
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
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);
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
    fclose(myout);
#ifndef MAC
    tcsetattr(0, 0, &normal_termio);
#endif
    return(0);            
}        

Rewrite_Input()
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
                            
    register unsigned int     i,
                            j,
                            k;
                            
    int                        h,
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
        If k = 0, the presentation the program was given has reduced to the empty presentation.
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
        while(z = *p++) TT[z]++;
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
                {
                printf("\n\nMade the following replacements of generators in the presentation:\n");
                if(Micro_Print_F)
                    fprintf(myout,"\n\nMade the following replacements of generators in the presentation:\n");
                }
            if(m > 80)
                {
                printf("\n");
                if(Micro_Print_F)                
                    fprintf(myout,"\n");
                m = 0;
                }    
            m += printf("%d) %c -> %c ",h,s,x);
            if(Micro_Print_F)            
            fprintf(myout,"%d) %c -> %c ",h,s,x);    
            }
        for(k = 1; k <= NumRelators; k++)
            {
            p = *Relators[k];
            while(z = *p)
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
    
    register unsigned char     **Temp;
        
    register int             i,
                            j;
    
    int                     SNumRelators;
    
    unsigned long            HS;
    
    SNumRelators = NumRelators;
    NumEmptyRels = 0;
    for(i = 1; i <= NumRelators; i++) HLock((char **) Relators[i]);
    
    /*****************************************************************************************
                        First, find and remove any empty relators.
    *****************************************************************************************/
    
    for(i = 1; i < SNumRelators; i++)
        {
        HS = GetHandleSize((char **) Relators[i]);
        if(HS == 1L)                        
            {
            SNumRelators--;
            NumEmptyRels ++;
                
            /*********************************************************************************
            Relators[i] is empty, find the last nonempty relator and swap it with Relators[i].
            *********************************************************************************/
                
            for(j = SNumRelators + 1; j > i; j--)
                {
                if(GetHandleSize((char **) Relators[j]) > 1L)
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
    for(i = 1; i <= NumRelators; i++) HUnlock((char **) Relators[i]);             
    return(SNumRelators);                
}

Compare_Input_Pres()
{
    /******************************************************************************************
                This routine compares the presentation in Copy_Of_Input[] and the
            presentation given by Relators[]. It returns TRUE if they are identical and
            FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char     *p,
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

Compare_Pres(k)
register int k;
{
    /******************************************************************************************
                This routine compares the presentation in SUR[k][] and the
            presentation given by Relators[]. It returns TRUE if they are identical and
            FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char     *p,
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

Compare_Dual_Pres(k)
register int k;
{
    /******************************************************************************************
            This routine compares the presentation in SUR[k][] and the presentation given by
            DualRelators[]. It returns TRUE if they are identical and FALSE otherwise.
    ******************************************************************************************/
    
    register unsigned char     *p,
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

Get_Initial_Diagram(PrintFlag)
int        PrintFlag;
{
    /******************************************************************************************
                        Get_Initial_Diagram() is supposed to determine 
            if the initial presentation is realizable by a Heegaard diagaram.
    ******************************************************************************************/
        
    register unsigned char     *p,
                            *q;
    
    register int             i,
                            j;
                            
    unsigned int            Flag1,
                            Flag2,
                            Flag3,
                            Flag4,
                            SNumFilled;                                                
    
    int                        Del_Only_Triv_Rel,
                            DistinctNonEmpty,
                            FirstPass,
                            NumReTries1,
                            NumReTries2,
                            SMicro_Print,
                            SMicro_Print_F,
                            SNumRelators;
                            
    unsigned long            HS;                        
    
    unsigned int Whitehead_Graph();
    unsigned int Reduce_Genus();        
    
    Input = NORMAL;
    TestRealizability1 = TRUE;
    TestRealizability2 = FALSE;
    Del_Only_Triv_Rel = TRUE;
    NumReTries1 = 0;
    NumReTries2 = 0;
    FirstPass = TRUE;
    BdryData = TRUE;

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
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,0);
                goto LIST_OPTIONS2;
            case 'x':
                printf("\n\n    HIT 'b' TO MICRO_PRINT TO BOTH THE SCREEN AND 'Heegaard_Results'.");
                printf("\n    HIT 's' TO MICRO_PRINT ONLY TO THE SCREEN.");
                if(SMicro_Print)
                    printf("\n    HIT 'o' TO TURN MICRO_PRINTING COMPLETELY OFF.");        
                GET_RESPONSE1:
                switch(WaitkbHit())
                    {
                    case 'b':
                        SMicro_Print = TRUE;
                        SMicro_Print_F = TRUE;
                        break;
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
                        goto GET_RESPONSE1;
                    }    
                break;
            default:
                SysBeep(5);
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
                    {
                    printf("\n\n                    Deleted a trivial generator.\n");
                    if(Micro_Print_F)                    
                        fprintf(myout,"\n\n                    Deleted a trivial generator.\n");
                    }
                else
                    {    
                    printf("\n\n                    Deleted %d trivial generators.\n",
                        SNumRelators - NumRelators);
                    if(Micro_Print_F)
                        fprintf(myout,"\n\n                    Deleted %d trivial generators.\n",
                            SNumRelators - NumRelators);
                    }
                }
            break;
        case TOO_LONG:
            printf("\n\n                    This presentation may be too long!");
            fprintf(myout,"\n\n                    This presentation may be too long!");
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
        
    if(Length == 0L)
        {
        printf("\n\nThe presentation has reduced to a presentation of length zero.");
        fprintf(myout,"\n\nThe presentation has reduced to a presentation of length zero.");
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
            Print_Relators(Relators,NumRelators,stdout);
            if(Micro_Print_F)
                {    
                fprintf(myout,"\n\n%lu automorphism(s) reduced the length to %lu.",
                    Automorphisms,Length);
                fprintf(myout,"\n\nThe presentation is currently:\n");
                Print_Relators(Relators,NumRelators,myout);
                }
            }
        else
            {
            printf("\n\nThe current set of relators has minimal length of %lu.",Length);
            if(Micro_Print_F)
                fprintf(myout,"\n\nThe current set of relators has minimal length of %lu.",Length);
            }        
        }
            
    if(Flag4 == TOO_LONG)
        {
        printf("\n\n                    This presentation may be too long!");
        fprintf(myout,"\n\n                    This presentation may be too long!");
        return(1);
        }                
    if(Flag4 == 1)
        {
        Modified_Init_Pres = TRUE;
        if(NumRelators == 1)
            {
            printf("\n\n                    This relator is not of full rank.");
            printf("\n\n                    Reducing the number of generators. . .\n");
            fprintf(myout,"\n\n                    This relator is not of full rank.");
            fprintf(myout,"\n\n                    Reducing the number of generators. . .\n");
            }
        else
            {    
            printf("\n\n                    The initial relators are not of full rank.");
            printf("\n\n                    Reducing the number of generators. . .\n");
            fprintf(myout,"\n\n                    The initial relators are not of full rank.");
            fprintf(myout,"\n\n                    Reducing the number of generators. . .\n");
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
                fprintf(myout,"\n\n                    Unable to determine whether the presentation is realizable.");                        
                if(Do_Not_Reduce_Genus || Delete_Only_Short_Primitives)
                    {
                    printf("\n\n          Suggest rerunning the example and allowing the program to delete all primitives.");    
                    fprintf(myout,"\n\n          Suggest rerunning the example and allowing the program to delete all primitives.");    
                    }    
                return(1);                
                }
            printf("\n\n                    Retrying the initial presentation.");
            NumRelators         = CopyNumRelators;
            NumGenerators         = CopyNumGenerators;
            Vertices             = 2*NumGenerators;
            Del_Only_Triv_Rel     = TRUE;
            FirstPass            = TRUE;
            for(i = 1,Length = 0L; i <= NumRelators; i++)
                {
                ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));
                Length += GetHandleSize((char **) Copy_Of_Input[i]);                
                p = *Copy_Of_Input[i];
                if((q = *Relators[i]) == NULL)
                    {
                    printf("\n\n    Memory Error. Sorry!");
                    return(1);
                    }
                while(*q++ = *p++) ;
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
        Print_Relators(Relators,NumRelators,stdout);        
        if(Micro_Print_F)
            {        
            fprintf(myout,"\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");        
            Print_Relators(Relators,NumRelators,myout);
            }
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
    
    if(Flag1 = Whitehead_Graph()) switch(Flag1)
        {
        case NON_PLANAR:
            printf("\n\n                    The Whitehead graph is nonplanar.");
            fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");
        case FATAL_ERROR:
            Fatal_Error();        
            return(1);
        case TOO_LONG:
            printf("\n\n                    This presentation may be too long!");
            fprintf(myout,"\n\n                    This presentation may be too long!");
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
                return(1);
            return(REDUCE_GENUS);
        case REDUCE_GENUS:
            goto _RESTART;    
        case NOT_CONNECTED:
            /**********************************************************************************
                The Whitehead graph corresponding to the initial presentation is not
                connected. Hence the corresponding Heegaard diagram is reducible.
            **********************************************************************************/
            
            if(!Del_Only_Triv_Rel)
                {
                Realization_Warning(stdout);
                Realization_Warning(myout);
                }
            printf("\n\n                    The Whitehead graph is not connected.");
            printf("\n\n                    Please check each summand separately.");
            fprintf(myout,"\n\n                    The Whitehead graph is not connected.");
            BdryData = FALSE;
            return(1);
                        
        /**************************************************************************************
                If the program gets to this point, then the program cannot find the diagram
            corresponding to the initial presentation. Before we give up, we let the program
            try to find some level transformations which might change the presentation into
            one for which the program can find the corresponding diagram. If the program can't
            find any level transformations that work,then we look for ways to reduce the genus.
        **************************************************************************************/    
                
        case SEP_PAIRS:
            printf("\n\n                    The Whitehead graph has a pair of separating vertices.\n");
            fprintf(myout,"\n\n                    The Whitehead graph has a pair of separating vertices.");
            UDV[This_Pres] = SEP_PAIRS;
            if(V1 & 1)
                LSP[This_Pres] = V1/2 + 97;
            else
                LSP[This_Pres] = V1/2 + 65;
            if(V2 & 1)
                LSQ[This_Pres] = V2/2 + 97;
            else
                LSQ[This_Pres] = V2/2 + 65;                
            NumCalled = 0;
            NotNewPres = 0;
            ReadPres = This_Pres;    
            switch(Flag3 = Level_Transformations(TRUE,TRUE,FALSE,0,FALSE))
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
                            fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");
                        case FATAL_ERROR:
                            Fatal_Error();        
                            return(1);
                        case TOO_LONG:
                            printf("\n\n                    This presentation may be too long!");
                            fprintf(myout,"\n\n                    This presentation may be too long!");
                            return(1);
                        case NON_UNIQUE_1:
                        case NON_UNIQUE_2:
                        case NON_UNIQUE_3:
                        case NON_UNIQUE_4:
                            UDV[This_Pres] = Flag2;
                        case V2_ANNULUS_EXISTS:         
                            for(i = 0; i < NumCalled; i++)
                            for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                            if(Non_Unique_Initial_Diagram())
                                return(1);
                            return(REDUCE_GENUS);    
                        case NO_ERROR:
                            printf("\n\n                    After some level transformations:");
                            fprintf(myout,"\n\n                    After some level transformations:");
                            Print_Realizability(Del_Only_Triv_Rel,ReadPres + 1);
                            if(UDV[This_Pres] == V2_ANNULUS_EXISTS)
                                {
                                printf("\n\n                    However, the realization is not unique because an annulus exists.");    
                                fprintf(myout,"\n\n                    However, the realization is not unique because an annulus exists.");
                                }
                            for(i = 0; i < NumCalled; i++)
                            for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                            if(Micro_Print)
                                {
                                printf("\n\nStarted with Presentation %d:\n",ReadPres + 1);
                                Print_Relators(Relators,NumRelators,stdout);        
                                if(Micro_Print_F)
                                    {        
                                    fprintf(myout,"\n\nStarted with Presentation %d:\n",ReadPres + 1);        
                                    Print_Relators(Relators,NumRelators,myout);
                                    }
                                }                                
                            goto DELETE_REDUNDANT;
                        break;    
                        }        
                    break;
                case 3:
                    for(i = 0; i < NumCalled; i++)
                    for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                    printf("\n\n                    The Whitehead graph is nonplanar.");
                    fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");
                    Fatal_Error();        
                    return(1);
                case 4:
                    break;
                case 5:
                    printf("\n\n          After some level transformations, ");
                    printf("the presentation contains a trivial generator.");
                    fprintf(myout,"\n\n          After some level transformations, ");
                    fprintf(myout,"the presentation contains a trivial generator.");        
                    if(!Do_Not_Reduce_Genus)
                        {
                        for(i = 0; i < NumCalled; i++)
                        for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
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
                }
            if(NumCalled) for(i = 1; i <= NumRelators; i++)
                {
                ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SLR[0][i]));
                if((p = *Relators[i]) == NULL)
                    {
                    printf("\n\nAn out of memory condition has arisen. Sorry!");
                    for(i = 0; i < NumCalled; i++)
                    for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                    return(1);
                    }
                q = *SLR[0][i];
                while(*p++ = *q++) ;
                }        
            for(i = 0; i < NumCalled; i++)
            for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
            if(Micro_Print)
                {
                printf("\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
                Print_Relators(Relators,NumRelators,stdout);                
                if(Micro_Print_F)
                    {
                    fprintf(myout,"\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");                
                    Print_Relators(Relators,NumRelators,myout);
                    }
                }    
            if(Flag3 != 5)
                {                
                printf("\n\n                    Unable to remove the separating pair of vertices.");
                fprintf(myout,"\n\n                    Unable to remove the separating pair of vertices.");        
                if(NumReTries1 <= 3 && !Do_Not_Reduce_Genus)
                    {
                    printf("\n\n                    Trying to reduce the genus of the presentation.");
                    fprintf(myout,"\n\n                    Trying to reduce the genus of the presentation.");
                    }
                }
            
        default:
        if(NumRelators > 1 && NumReTries1 <= 3 && !Do_Not_Reduce_Genus)
            {
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
                    printf("\n\n                    Unable to delete a primitive.");
                    if(Micro_Print_F)
                        fprintf(myout,"\n\n                    Unable to delete a primitive.");                    
                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                    NumReTries1 ++;
                    printf("\n\n                    Retrying the initial presentation.");
                    NumRelators         = CopyNumRelators;
                    NumGenerators         = CopyNumGenerators;
                    Vertices             = 2*NumGenerators;
                    Del_Only_Triv_Rel     = TRUE;
                    FirstPass            = TRUE;
                    for(i = 1,Length = 0L; i <= NumRelators; i++)
                        {
                        ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));
                        Length += GetHandleSize((char **) Copy_Of_Input[i]);                
                        p = *Copy_Of_Input[i];
                        if((q = *Relators[i]) == NULL)
                            {
                            printf("\n\n    Memory Error. Sorry!");
                            return(1);
                            }
                        while(*q++ = *p++) ;
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
                        {
                        Realization_Warning(stdout);
                        Realization_Warning(myout);
                        }
                    printf("\n\n                    This manifold is a lens space.");
                    fprintf(myout,"\n\n                    This manifold is a lens space.");
                    return(REDUCE_GENUS);
                    }
                if(FoundPrimitive)
                    {
                    FoundPrimitive = FALSE;
                    printf("\n\n                    The program found relator(s) which are primitive");
                    printf("\n                    and deleted their consequences from the presentation.");
                    fprintf(myout,"\n\n                    The program found relator(s) which are primitive");
                    fprintf(myout,"\n                    and deleted their consequences from the presentation.");
                    }
                if(FoundPower)
                    {
                    FoundPower = FALSE;
                    printf("\n\n                    The program found relator(s) which are proper powers");
                    printf("\n                    and deleted their consequences from the presentation.");
                    fprintf(myout,"\n\n                    The program found relator(s) which are proper powers");
                    fprintf(myout,"\n                    and deleted their consequences from the presentation.");
                    }
                if(EmtyRel)
                    {
                    if(SNumRelators > NumRelators)
                        {
                        Realization_Warning(stdout);
                        Realization_Warning(myout);
                        }
                    EmtyRel = FALSE;
                    printf("\n\n                    This manifold is a connected sum.");
                    printf("\n                    Please check each summand separately.\n");
                    fprintf(myout,"\n\n                    This manifold is a connected sum.");
                    fprintf(myout,"\n                    Please check each summand separately.\n");
                    return(1);
                    }
                printf("\n");
                fprintf(myout,"\n");        
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
            NumCalled = 0;
            NotNewPres = 0;
            if(Micro_Print)
                {
                printf("\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
                Print_Relators(Relators,NumRelators,stdout);
                if(Micro_Print_F)
                    {                                
                    fprintf(myout,"\n\nAt Get_Initial_Diagram(), the presentation is currently:\n");
                    Print_Relators(Relators,NumRelators,myout);
                    }
                }
            Fill_A(NumRelators);
            Get_Matrix();
            switch(Flag3 = Level_Transformations(TRUE,FALSE,FALSE,0,Del_Only_Triv_Rel))
                {
                case 3:
                    for(i = 0; i < NumCalled; i++)
                    for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                    printf("\n\n                    The Whitehead graph is nonplanar.");
                    fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");
                    Fatal_Error();        
                    return(1);                                        
                case 5:
                    for(i = 0; i < NumCalled; i++)
                    for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                    if(Micro_Print)
                        {
                        printf("\n\nThe presentation is currently:\n");
                        Print_Relators(Relators,NumRelators,stdout);                        
                        if(Micro_Print_F)
                            {                        
                            fprintf(myout,"\n\nThe presentation is currently:\n");
                            Print_Relators(Relators,NumRelators,myout);
                            }
                        }
                    printf("\n\n          After some level transformations, ");
                    printf("the presentation contains a trivial generator.");
                    fprintf(myout,"\n\n          After some level transformations, ");
                    fprintf(myout,"the presentation contains a trivial generator.");                    
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
            for(i = 0; i < NumCalled; i++)
            for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);        
            }
        
        if(!Delete_Only_Short_Primitives && NumReTries2 < 3
            && Flag1 == SEP_PAIRS && Flag3 == 4)
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
                fprintf(myout,"\n\n                    Unable to remove the annulus which is present.");
                if(j == 1)
                    {
                    SysBeep(5);
                    printf("\n\n    Modifying the initial presentation by deleting a duplicated or empty relator.\n",j);
                    fprintf(myout,"\n\n    Modifying the initial presentation by deleting a duplicated or empty relator.\n",j);
                    printf("\n    NOTE: This modification may change the homeomorphism type of M.\n");
                    fprintf(myout,"\n    NOTE: This modification may change the homeomorphism type of M.\n");
                    }                
                else
                    {
                    SysBeep(5);
                    printf("\n\n    Modifying the initial presentation by deleting %d duplicated or empty relators.\n",j);
                    fprintf(myout,"\n\n    Modifying the initial presentation by deleting %d duplicated or empty relators.\n",j);
                    printf("\n    NOTE: This modification may change the homeomorphism type of M.\n");
                    fprintf(myout,"\n    NOTE: This modification may change the homeomorphism type of M.\n");
                    }
                NumRelators = DistinctNonEmpty;
                for(i = 1,Length = 0L; i <= NumRelators; i++)
                    Length += GetHandleSize((char **) Relators[i]);
                Length -= NumRelators;        
                printf("\n\nThe modified presentation is:\n");
                fprintf(myout,"\n\nThe modified presentation is:\n");
                Print_Relators(Relators,NumRelators,stdout);
                Print_Relators(Relators,NumRelators,myout);    
                goto _RESTART;
                }
            }
                    
        printf("\n\n                    Unable to determine whether the presentation is realizable.");
        fprintf(myout,"\n\n                    Unable to determine whether the presentation is realizable.");                        
        if(Do_Not_Reduce_Genus || Delete_Only_Short_Primitives)
            {
            printf("\n\n          Suggest rerunning the example and allowing the program to delete all primitives.");    
            fprintf(myout,"\n\n          Suggest rerunning the example and allowing the program to delete all primitives.");    
            }    
        return(1);    
        }
        
    /******************************************************************************************
            Otherwise, if the program gets to here, then the program found the diagram
        corresponding to the initial presentation. So we are in business.
    ******************************************************************************************/    
    
    Print_Realizability(Del_Only_Triv_Rel,Save_Init_Pres + 1);
    if(UDV[This_Pres] == V2_ANNULUS_EXISTS)
        {
        printf("\n\n                    However, the realization is not unique because an annulus exists.");    
        fprintf(myout,"\n\n                    However, the realization is not unique because an annulus exists.");
        }
        
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
        Save_Pres() is called when the program has determined that a presentation should be
        put on file for further processing. It saves a copy of the presentation in the array
        SUR[][] and saves some other data about the presentation in the arrays listed below.
    ******************************************************************************************/
    
    register unsigned char     *p,
                            *q;
                            
    unsigned char             ***MyRelators;
                    
    register int             i,
                            j;
                    
    int                        Error = FALSE;
    
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
    
    for(i = 1; i <= NumRelators; i++)
        {
        SUR[NumFilled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) MyRelators[i]));            
        if((q = *SUR[NumFilled][i]) == NULL)
            {
            Error = TRUE;
            goto END;
            }
        p = *MyRelators[i];
        while(*q++ = *p++) ;                                    
        }    

    ComponentNum[NumFilled] = CurrentComp;
    Daughters[NumFilled]    = Daut;
    ER[NumFilled]            = F5;
    FR[NumFilled]            = From;
    NG[NumFilled]             = NumGenerators;
    NR[NumFilled]             = NumRelators;
    PRIM[NumFilled]         = F2;
    SURL[NumFilled]         = Len;
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
    if(Num_Level_Transformations && ! Find_All_Min_Pres)
        OnStack = NumFilled - ReadPres;
    else    
        OnStack += 2*NumGenerators;
        
    if(Micro_Print)
        {
        printf("\n\nSaved the current presentation as: Presentation %u.\n",NumFilled);
        if(Micro_Print_F)        
            fprintf(myout,"\n\nSaved the current presentation as: Presentation %u.\n",NumFilled);
        }
    printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,CurrentComp);
    printf("Gen%3d  Rel%3d  Length%6lu  From%6u  ",NumGenerators,NumRelators,Len,From + 1);
    
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
    END:
    if(Error)
        {
        for(j = 1; j < i; j++) DisposeHandle((char **) SUR[NumFilled][j]);
        printf("\n    Memory Error. Sorry!");
        return(TOO_LONG);
        }
    return(NO_ERROR);        
}

int Non_Unique_Initial_Diagram(void)
{
    /******************************************************************************************
        This routine is called when the program has determined that if the initial presentation
        is realizable by a Heegaard diagram, then the realization is not unique. The routine
        determines whether we have a presentation that is realizable, but not uniquely so, or
        we have a presentation that is not realizable at all. It returns 0 if the presentation
        is realizable and returns 1 if it is not realizable.
    ******************************************************************************************/
    
    int                i;
                    
    unsigned int    AnnulusExists,
                    Flag1;                
    
    unsigned int Whitehead_Graph();
    unsigned int Reduce_Genus();
    
    /******************************************************************************************
                            Check whether a valence-two annulus exists.
    ******************************************************************************************/
    
    AnnulusExists         = FALSE;    
    TestRealizability2     = TRUE;
    TestRealizability3     = TRUE;
    DrawingDiagrams     = FALSE;
    
    Flag1 = Whitehead_Graph();
    
    TestRealizability3     = FALSE;
    
    if(Flag1 == V2_ANNULUS_EXISTS)
        {
        Flag1 = Whitehead_Graph();
        AnnulusExists = TRUE;
        }
        
    switch(Flag1)
        {
        case NO_ERROR:
            printf("\n\n                    Presentation %d is realizable.",This_Pres + 1);
            fprintf(myout,"\n\n                    Presentation %d is realizable.",This_Pres + 1);
            printf("\n\n                    But the realization is not unique.");
            fprintf(myout,"\n\n                    But the realization is not unique.");
            
            if(AnnulusExists)
                {
                printf("\n\n    NOTE: The program is investigating one diagram that realizes the initial presentation.");
                printf("\n          Subsequent results may not apply to all realizations.");
                fprintf(myout,"\n\n    NOTE: The program is investigating one diagram that realizes the initial presentation.");
                fprintf(myout,"\n          Subsequent results may not apply to all realizations.");
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
                    printf("\n\n                    Unable to delete a primitive.");
                    if(Micro_Print_F)
                        fprintf(myout,"\n\n                    Unable to delete a primitive.");                    
                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                    return(1);        
                }
            if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
                {
                if(!AnnulusExists)
                    {
                    printf("\n\nNOTE: The program is investigating one diagram that realizes the initial presentation.");
                    printf("\n      Subsequent results may not apply to all realizations.");
                    fprintf(myout,"\n\nNOTE: The program is investigating one diagram that realizes the initial presentation.");
                    fprintf(myout,"\n      Subsequent results may not apply to all realizations.");
                    }
                if(FoundPrimitive)
                    {
                    printf("\n\n                    The program found a relator which is a primitive");
                    printf("\n                    and deleted its consequences from the presentation.");
                    fprintf(myout,"\n\n                    The program found a relator which is a primitive");
                    fprintf(myout,"\n                    and deleted its consequences from the presentation.");    
                    }
                if(FoundPower)
                    {
                    printf("\n\n                    The program found a relator which is a proper power");
                    printf("\n                    and deleted its consequences from the presentation.");
                    fprintf(myout,"\n\n                    The program found a relator which is a proper power");
                    fprintf(myout,"\n                    and deleted its consequences from the presentation.");                    
                    }
                if(LensSpace)
                    {
                    printf("\n\n                    This manifold is a lens space.");
                    fprintf(myout,"\n\n                    This manifold is a lens space.");
                    }
                if(EmtyRel)
                    {
                    printf("\n\n                    This manifold is a connected sum.");
                    fprintf(myout,"\n\n                    This manifold is a connected sum.");
                    }
                printf("\n");
                fprintf(myout,"\n");    
                return(0);
                }    
            return(1);
        case FATAL_ERROR:
            printf("\n\n                    The initial presentation is not realizable.\n");
            Print_Relators(Relators,NumRelators,stdout);
            fprintf(myout,"\n\n                    The initial presentation is not realizable.\n");
            Print_Relators(Relators,NumRelators,myout);
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
        
    register int     i,
                    j;
    
    if(Length > 0L)    Canonical_Rewrite(Relators,FALSE,FALSE);        
     for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
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
         }
     
     if(Micro_Print)
         {
         if(i < NumFilled)
             {
             printf("\n\nA copy of the current presentation is already on file as presentation %d.",
                 i+1);
             if(Micro_Print_F)    
                 fprintf(myout,"\n\nA copy of the current presentation is already on file as presentation %d.",
                     i+1);
             }
         else if(Dup_On_File < INFINITE)
             {
             printf("\n\nThe current presentation is a duplicate of presentation %d of summand %d.",
                 i+1,ComponentNum[Dup_On_File]);
             if(Micro_Print_F)    
                 fprintf(myout,"\n\nThe current presentation is a duplicate of presentation %d of summand %d.",
                     i+1,ComponentNum[Dup_On_File]);
             }                         
         }    
     return(i);                    
}

void Init_G_Variables()
{
    /******************************************************************************************
        This routine initializes some of the global variables used by the program.
    ******************************************************************************************/
        
    register int     i,
                    j;
                    
    Band_Sums                         = 0L;
    BDY[0]                            = 2;
    BdryData                        = FALSE;
    BytesAvailable                    = 2097152L;
    BytesUsed                        = 0L;
    Compute_Stabilizers                = FALSE;
    Count                             = 0;
    CurrentComp                     = 1;
    Delete_Only_Short_Primitives     = FALSE;
    Do_Not_Reduce_Genus                = TRUE;
    DrawingDiagrams                 = FALSE;
    EmtyRel                            = FALSE;
    Find_All_Min_Pres                = FALSE;
    FormBandsums                    = FALSE;
    FoundPower                         = FALSE;
    FoundPrimitive                     = FALSE;
    From_BANDSUM                    = 0;
    From_DUALIZE                     = 0;
    Length                             = 0L;
    LensSpace                         = FALSE;
    Level_Interrupt                    = FALSE;
    Micro_Print                        = FALSE;
    Micro_Print_F                    = FALSE;
    Minimum                         = 0L;
    Modified_Init_Pres                 = FALSE;
    NoReport                         = TRUE;
    NumDiagrams                     = 0L;
    NumDualized                        = 0L;
    NumFilled                         = 0;
    Num_Level_Transformations         = 0L;
    NumTimes                        = 0;
    OnStack                            = 0;
    ReadPres                         = 0;
    SaveMinima                         = FALSE;
    Saved_Vertices                    = 0;
    Start_Level_Search                = 0;
    Starting_Pres                    = 0;
    TotalAuts                        = 0L;
    TotalComp                         = 1;
    UserSaidQuit                    = FALSE;
    Vertices                         = 2*NumGenerators;
    for(i = 0; i < MAX_SAVED_PRES; i++)
        {
        SURL[i] = 0L;
        QPM[i] = EOS;
        }
    for(i = 0; i < MAXNUMCOMPONENTS; i++)
    for(j = 1; j <= MAXNUMGENERATORS; j++) MLC[i][j] = BIG_NUMBER;
    for(i = 0; i < MAXNUMCOMPONENTS; i++) CBC[i][0] = BDRY_UNKNOWN;
    for(i = 0; i <= MAXNUMCOMPONENTS; i++) CS[i] = EOS;
}

void Delete_Old_Presentations(void)
{
    unsigned int    i,
                    j;
                            
    for(i = 0; i < NumFilled; i++)
        {
        if(UDV[i] == ANNULUS_EXISTS || UDV[i] == V2_ANNULUS_EXISTS)
            DisposeHandle((char **) SUR[i][0]);
        UDV[i] = PRIM[i] = 0;
        for(j = 1; j <= NR[i]; j++) DisposeHandle((char **) SUR[i][j]);
        }
    NumFilled         = 0;
    BytesAvailable    = 2097152L;
    BytesUsed         = 0L;
    UserSaidQuit     = FALSE;                    
}

ReRun_A_Presentation()
{
    unsigned char   *p,
                    *q,
                    *r,
                    revnum[8];
    
    unsigned int    PresNum;
    
    unsigned int    h,
                    i,
                    j,
                    k;

    DrawingDiagrams = FALSE;

    if(NoReport)
        {
        SysBeep(5);
        printf("\n\n                                    Caution!!!");
        printf("\n\n     All the presentations in memory, except the one you rerun, ");
        printf("\n     will be vaporized unless you save them now!");
        printf("\n\n     HIT 'v' TO REVIEW THESE PRESENTATIONS.");
        if(NoReport)
            printf("\n     HIT 's' TO SAVE THESE PRESENTATIONS TO THE FILE 'Heegaard_Results'.");
        printf("\n     OR HIT ANY OTHER KEY TO CONTINUE.");
        switch(WaitkbHit())
            {
            case 'v':
                REVIEW:
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);
                printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
                GET_RESPONSE2:
                switch(WaitkbHit())
                    {
                    case 'y':
                        goto REVIEW;
                    case 'n':
                        break;
                    default:
                        SysBeep(5);
                        goto GET_RESPONSE2;
                    }
                break;
            case 's':
                if(NoReport)
                    Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,1,1,1,0,1,0);
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
                    ReallocateHandle((char **) Relators[i],GetHandleSize((char **) CD_Surgery_Rel[i]));                
                    p = *CD_Surgery_Rel[i];
                    if((q = *Relators[i]) == NULL)
                        {
                        printf("\n\n    Memory Error. Sorry!");
                        return(1);
                        }
                    while(*q++ = *p++) ;
                    }
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators,stdout);
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                fptr = myout;
                Input = RERUN;
                Did_Cutting_Disk_Surgery = FALSE;
                return(0);
            case 'n':
                Did_Cutting_Disk_Surgery = FALSE;
                break;
            default:
                SysBeep(5);
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
                    ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Exp_Surgery_Rel[i]));                
                    p = *Exp_Surgery_Rel[i];
                    if((q = *Relators[i]) == NULL)
                        {
                        printf("\n\n    Memory Error. Sorry!");
                        Did_Exponent_Surgery = FALSE;
                        return(1);
                        }
                    while(*q++ = *p++) ;
                    }
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators,stdout);
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                fptr = myout;
                Input = RERUN;
                Did_Exponent_Surgery = FALSE;
                return(0);
            case 'n':
                Did_Exponent_Surgery = FALSE;
                break;
            default:
                SysBeep(5);
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
                NumRelators = NumKnot_Or_Link_Rel;
                for(i = 1; i <= NumRelators; i++)
                    {
                    ReallocateHandle((char **) Relators[i],GetHandleSize((char **) KorLRelators[i]));                
                    p = *KorLRelators[i];
                    if((q = *Relators[i]) == NULL)
                        {
                        printf("\n\n    Memory Error. Sorry!");
                        return(1);
                        }
                    while(*q++ = *p++) ;
                    }
                printf("\n\nThe initial presentation is: %s\n",PresName);                    
                Print_Relators(Relators,NumRelators,stdout);                                    
                if(Wirtinger())
                    {
                    printf("\nSurgery is not possible. Sorry!");
                    Knot_Or_Link = FALSE;
                    break;
                    }        
                Delete_Old_Presentations();
                NoReport = TRUE;
                WhichInput = 0;
                printf("\nThe surgered presentation is:\n");
                Print_Relators(Relators,NumRelators,stdout);
                fptr = myout;
                Input = RERUN;
                return(0);        
            case 'n':
                Knot_Or_Link = FALSE;
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE3;
            }
        }

    r = (unsigned char *) NewPtr(100L);        
    printf("\n\nENTER A PRESENTATION FROM 0 TO %u THAT YOU WANT RERUN AND HIT 'return'.     ",NumFilled);
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
        printf("\n\nRerunning the original presentation.\n");
        NumRelators = CopyNumRelators;
        for(i = 1; i <= NumRelators; i++)
            {
            ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Input[i]));                
            p = *Copy_Of_Input[i];
            if((q = *Relators[i]) == NULL)
                {
                DisposePtr((char *) r);
                printf("\n\n    Memory Error. Sorry!");
                return(1);
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
            ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[WhichInput][i]));                
            p = *SUR[WhichInput][i];
            if((q = *Relators[i]) == NULL)
                {
                DisposePtr((char *) r);
                printf("\n\n    Memory Error. Sorry!");
                return(1);
                }            
            while(*q++ = *p++) ;
            }
        }
    DisposePtr((char *) r);        
    Delete_Old_Presentations();
    NoReport = TRUE;
    
    Print_Relators(Relators,NumRelators,stdout);
    fptr = myout;
    Input = RERUN;
    return(0);
}

int Display_File_Input_Presentations(void)
{
    char            c;
    
    register int     x;                                        
    
    if((input_relators = fopen("Input_Presentations","r+")) == NULL)
        {
        SysBeep(5);
        printf("\nUnable to open the file 'Input_Presentations'.\n");
        return(1);
        }
    
    ObscureCursor();
    
    printf("\n\n                    This is the file 'Input_Presentations'.");
    printf("\n\nNote: Hit 'return' to stop reviewing this file.");
    printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.\n\n");
        
    while((x = getc(input_relators)) != EOF)
        {
        if(x == '\n' && (c = mykbhit())) switch(c)
            {
            case ' ':
                WAIT:
                switch(WaitkbHit())
                    {
                    case ' ':
                        break;
                    case '\n':
                        fclose(input_relators);
                        return(NO_ERROR);
                    case '\r':
                        fclose(input_relators);
                        return(NO_ERROR);    
                    default:
                        SysBeep(5);
                        printf("\n\nNote: Hit 'return' to stop reviewing this file.");
                        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                        goto WAIT;        
                    }
                break;    
            case '\n':
                fclose(input_relators);
                return(NO_ERROR);
            case '\r':
                fclose(input_relators);
                return(NO_ERROR);    
            default:
                break;        
            }
        putc((char)x,stdout);
        }
    fclose(input_relators);
    return(NO_ERROR);
}

int Display_File_Simplify_Heegaard(void)
{
    char            c;

    register int     x;                                        
    
    if(myout == NULL)
        {
        SysBeep(5);
        printf("\nUnable to display the file 'Heegaard_Results'.\n");
        return(1);
        }
    
    ObscureCursor();
    
    printf("\n\n                    This is the file 'Heegaard_Results'.\n\n");
    printf("\n\nNote: Hit 'return' to stop reviewing this file.");
    printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.\n\n");
    rewind(myout);    
    while((x = getc(myout)) != EOF)
        {
        if(x == '\n' && (c = mykbhit())) switch(c)
            {
            case ' ':
                WAIT:
                switch(WaitkbHit())
                    {
                    case ' ':
                        break;
                    case '\n':
                        return(NO_ERROR);
                    case '\r':
                        return(NO_ERROR);    
                    default:
                        SysBeep(5);
                        printf("\n\nNote: Hit 'return' to stop reviewing this file.");
                        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                        goto WAIT;        
                    }
                break;    
            case '\n':
                return(NO_ERROR);
            case '\r':
                return(NO_ERROR);    
            default:
                break;        
            }
        putc((char)x,stdout);
        }
    return(NO_ERROR);    
}

int Display_Help_File(void)
{
    char            c;
    
    register int     x;
    
    FILE            *fptr;

    if((fptr = fopen("Heegaard_Help","r")) == NULL)
        {
        SysBeep(5);
        printf("\nUnable to open the Help file 'Heegaard_Help'.");
        printf("\nPlease locate the file 'Heegaard_Help' and put it in the same folder with");
        printf("\nthe program.");
        return(1);
        }
    
    ObscureCursor();
    while((x = getc(fptr)) != EOF)
        {
        if(x == '\n' && (c = mykbhit())) switch(c)
            {
            case ' ':
                WAIT:
                switch(WaitkbHit())
                    {
                    case ' ':
                        break;
                    case '\n':
                        fclose(fptr);
                        return(NO_ERROR);
                    case '\r':
                        fclose(fptr);
                        return(NO_ERROR);    
                    default:
                        SysBeep(5);
                        printf("\n\nNote: Hit 'return' to stop reviewing this file.");
                        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                        goto WAIT;        
                    }
                break;    
            case '\n':
                fclose(fptr);
                return(NO_ERROR);
            case '\r':
                fclose(fptr);
                return(NO_ERROR);    
            default:
                SysBeep(5);
                printf("\n\nNote: Hit 'return' to stop reviewing this file.");
                printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                break;        
            }        
        putc((char)x,stdout);
        }
    fclose(fptr);
    return(NO_ERROR);    
}

Get_Presentation_From_File()
{
    register unsigned char    *p,
                            *q,
                            *r,
                            t;
    
    unsigned int            h,
                            i;
                            
    long                    StrLength;                                                
    
    if((input_relators = fopen("Input_Presentations","r+")) == NULL)
        {
        SysBeep(5);
        printf("\nUnable to open the file 'Input_Presentations'.");
        printf("\nPlease locate the file 'Input_Presentations', make sure it is closed,");
        printf("\nand place it in the same folder that the program is in.");
        return(1);
        }
    r = (unsigned char *) NewPtr((Size)(MAXLENGTH + 1));
    printf("\n\nPLEASE INDICATE WHICH PRESENTATION YOU WANT USED BY ENTERING ITS IDENTIFIER:\n\n    ");        
    ReadString((char *)PresName, GetPtrSize(PresName));

    /******************************************************************************************
                    Look for a presentation with the given identifier.
    ******************************************************************************************/    
    
    rewind(input_relators);
    do
        {
        if(fgets((char *) r,MAXLENGTH,input_relators) == NULL)
            {
            printf("\n\nUnable to find this presentation in the file 'Input_Presentations'.");
            SysBeep(5);
            fclose(input_relators);
            DisposePtr((char *) r);
            printf("\n\nREVIEW THE FILE 'Input_Presentations' ?  HIT 'y' OR 'n'.");
            GET_RESPONSE1:
            switch(WaitkbHit())
                {
                case 'y':
                    Display_File_Input_Presentations();
                    break;
                case 'n':
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE1;
                }
            return(1);
            }
        p = r;
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
        if(fgets((char *) r,MAXLENGTH,input_relators) == NULL)
            {
            printf("\n\nUnable to find any non-trivial relators in this presentation!");
            SysBeep(5);
            fclose(input_relators);
            DisposePtr((char *) r);
            return(1);
            }
    while(*r == '\n' || *r == '\r');

    /******************************************************************************************
        Read in the relators, at one relator to each non-empty line, stripping off leading
        spaces and tabs.
    ******************************************************************************************/
    
    ObscureCursor();
                        
    for(i = 1,NumRelators = 0; i <= MAXNUMRELATORS; i++)
        {
        p = r;
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
        while(t = *p)
            {
            if(t == '\n' || t == ' ' || t == '\t')
                {
                *p = EOS;
                break;
                }     
            putc(t,stdout);
            if(!isalpha(t))
                {
                SysBeep(5);
                printf("\nPlease check relator %2u.\n",i);
                printf("\nA relator can only contain upper and lower case letters!\n");
                printf("\nPlease make sure each relator is on its own line!");
                printf("\nAnd the presentation is terminated by an empty line!");
                fclose(input_relators);
                DisposePtr((char *) r);
                return(1);
                }
            p++;        
            }                                
        if(*r == '\n') break;
        StrLength = strlen((char *) r);
        if(StrLength >= MAXLENGTH)
            {
            SysBeep(5);
            printf("\nRelator %d is too long!!",NumRelators + 1);
            fclose(input_relators);
            DisposePtr((char *) r);
            return(1);
            }
        if(StrLength == h) break;
        NumRelators ++;
        ReallocateHandle((char **) Relators[i],StrLength + 1 - h);            
        if((p = *Relators[i]) == NULL)
            {
            SysBeep(5);
            fclose(input_relators);
            DisposePtr((char *) r);
            printf("\n\n    Memory Error. Sorry!");
            return(1);
            }
        q = r + h;
        while(*p++ = *q++) ;
        if(fgets((char *) r,MAXLENGTH,input_relators) == NULL) break;        
        }
    if(NumRelators == MAXNUMRELATORS)
        printf("\n\nThe program will not accept any more than %d relators.",NumRelators);            
    fclose(input_relators);
    DisposePtr((char *) r);
    return(0);
}

Get_Presentation_From_KeyBoard()
{
    unsigned char   *p,
                    *q,
                    *r,
                    t;
    
    int             Flag1;
    
    unsigned int    i;
                            
    unsigned long   StrLength;                                                
    
/*    r = (unsigned char *) NewPtr((Size)(MAXLENGTH + 1));
    if(r == NULL)
        {
        SysBeep(5);
        printf("\n\n    Memory Error. Sorry!");
        return(1);
        }                */
    printf("\n\nPLEASE ENTER THE RELATORS.");
    printf("        (Note: Hit 'return' twice to terminate entry of the presentation.)\n\n");
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
        /*    printf("R%3u)   ",i);
            ungetc((char)t,stdin);
            gets((char *) r);    */
            Flag1 = FALSE;
            while(t = *p++)
                {
                if(!isalpha(t))
                    {
                    Flag1 = TRUE;
                    SysBeep(5);
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
        if(StrLength == 0L) break;
        NumRelators ++;
        ReallocateHandle((char **) Relators[i],StrLength + 1);            
        if((p = *Relators[i]) == NULL)
            {
            printf("\n\n    Memory Error. Sorry!");
            return(1);
            }
        q = r;
        while(*p++ = *q++) ;            
        }
    if(NumRelators == MAXNUMRELATORS)
        printf("\n\nThe program will not accept more than %d relators.",NumRelators);        
    printf("\n\nPlease enter a name by which the program can refer to this presentation,");
    printf("\nand then hit 'return'.");
    printf("\n\nSAVE THIS PRESENTATION AS: ");
    ReadString((char *)PresName, GetPtrSize(PresName));
    return(0);
}

Check_Realizability_Of_The_Initial_Presentation()
{
    Get_Simplification_Parameters_From_User(FALSE,FALSE);
    if(Get_Initial_Diagram(FALSE) == 2)
        NoReport = FALSE;
    else    
        Report(0,NumDiagrams,0,0,0,0,1,0,1,0);
    printf("\n\nDISCARD ALL BUT THE INITIAL PRESENTATION ?  HIT 'y' OR 'n'.");
    GET_RESPONSE1:
    switch(WaitkbHit())
        {
        case 'y':
            return(0);
        case 'n':
            return(1);
        default:
            SysBeep(5);
            goto GET_RESPONSE1;
        }        
}

void Display_Diagram_Of_The_Initial_Presentation(void)
{
    int             Flag1,
                    i,
                    SNumGenerators;
    
    unsigned long    Scratch;
    
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
    if(Length == 0L)
        {
        printf("\n\nThe presentation has reduced to the empty presentation.");
        printf("\nThere is no diagram to display.");
        DrawingDiagrams = FALSE;
        return;
        }    
    if(Automorphisms)
        {
        SysBeep(5);
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
        printf("\nA copy of this presentation will be saved in the file 'Heegaard_Results'.");
        Canonical_Rewrite(Relators,FALSE,FALSE);
        Fill_A(NumRelators);
        Saved_Vertices = 0;
        fprintf(stdout,"\n\nThe minimal length presentation is:\n");
        Print_Relators(Relators,NumRelators,stdout);
        fprintf(myout,"\n\nThe minimal length presentation is:\n");
        Print_Relators(Relators,NumRelators,myout);    
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            Scratch += GetHandleSize((char **) Relators[i]);
        Scratch -= NumRelators;
        printf("\nThis presentation has length %ld ",Scratch);
        fprintf(myout,"\nThis presentation has length %ld ",Scratch);
        printf("\n\n    HIT ANY KEY TO SEE THE HEEGAARD DIAGRAM.");
        WaitkbHit();
        }
    Get_Matrix();
    Check_Connected();
    SepPairs = Sep_Pairs(0,0);
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
    NonPlanar = Planar(FALSE,FALSE);
    Print_Graph(FALSE);
    DrawingDiagrams = FALSE;
}

void Get_Simplification_Parameters_From_User(int Flag1,int Flag2)
{
    if(Flag2 == FALSE) OnlyReducingBandsums = FormBandsums = FALSE;
    if(Flag1)
        {
        if(FormBandsums)
            {
            if(OnlyReducingBandsums)
                printf("\nThe program is currently forming only length reducing bandsums.\n");            
            else
                printf("\nThe program is currently forming all possible bandsums.\n");
            }
        else
            printf("\nThe program is not currently forming any bandsums.\n");
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
                        FormBandsums = TRUE;
                        OnlyReducingBandsums = FALSE;
                        break;
                    case 'r':
                        FormBandsums = TRUE;
                        OnlyReducingBandsums = TRUE;
                        break;
                    case 'n':
                        OnlyReducingBandsums = FormBandsums = FALSE;
                        break;    
                    default:
                        SysBeep(5);
                        goto GET_RESPONSE2;
                    }
                printf("\n");    
                break;        
            case 'n':
                OnlyReducingBandsums = FormBandsums = FALSE;                
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE1;    
            }    
        }
    if(Flag1)
        {
        if(!Do_Not_Reduce_Genus && !Delete_Only_Short_Primitives)
            {
            printf("\nThe program is currently trying to delete all primitive relators.\n");
            printf("\nCONTINUE TO DELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
            }
        else
            {
            if(Delete_Only_Short_Primitives && !Do_Not_Reduce_Genus)
            printf("\nThe program is currently deleting only primitive relators of length 1 and 2.\n");
            else
            printf("\nThe program is not currently attempting to delete any primitive relators.\n");
            printf("\nCHANGE TO DELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
            }
        }
    else        
    printf("\nDELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.\n");
    GET_RESPONSE3:
    switch(WaitkbHit())
        {
        case 'y':
            Do_Not_Reduce_Genus = FALSE;
            Delete_Only_Short_Primitives = FALSE;
            break;
        case 'n':    
            printf("\nDELETE PRIMITIVE RELATORS OF LENGTH 1 AND LENGTH 2 ? HIT 'y' OR 'n'.\n");
            GET_RESPONSE4:
            switch(WaitkbHit())
                {
                case 'y':
                    Delete_Only_Short_Primitives = TRUE;
                    Do_Not_Reduce_Genus = FALSE;
                    break;
                case 'n':
                    Delete_Only_Short_Primitives = FALSE;
                    Do_Not_Reduce_Genus = TRUE;
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE4;
                }
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE3;
        }
        
    if(!Flag1) Turn_Micro_Print_On();
    
    if(Flag1)
        {
        if(Find_All_Min_Pres)
            {
            printf("\nThe program is currently trying to make broad use of level-transformations.\n");
            printf("\nCONTINUE TO MAKE BROAD USE OF LEVEL-TRANSFORMATIONS ? HIT 'y' OR 'n'.\n");
            }
        else
            {
            printf("\nThe program is currently trying to make only limited use of level-transformations.\n");
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
                SysBeep(5);
                goto GET_RESPONSE5;
            }
        }                        
}

void Turn_Micro_Print_On(void)
{
    printf("\nTURN Micro_Printing ON ? HIT 'y' OR 'n'.\n");
    GET_RESPONSE1:
    switch(WaitkbHit())
        {
        case 'y':
            printf("\n    HIT 'b' TO MICRO_PRINT TO BOTH THE SCREEN AND 'Heegaard_Results'.");
            printf("\n    HIT 's' TO MICRO_PRINT ONLY TO THE SCREEN.    ");
            if(Micro_Print)
                printf("\n    HIT 'o' TO TURN MICRO_PRINTING COMPLETELY OFF.    ");
            printf("\n");    
            GET_RESPONSE2:
            switch(WaitkbHit())
                {
                case 'b':
                    Micro_Print = TRUE;
                    Micro_Print_F = TRUE;
                    break;
                case 's':
                    Micro_Print = TRUE;
                    Micro_Print_F = FALSE;
                    break;                            
                case 'o':
                    if(!Micro_Print)
                        {
                        SysBeep(5);
                        goto GET_RESPONSE2;                            
                        }
                    Micro_Print = FALSE;
                    Micro_Print_F = FALSE;
                    break;                        
                default:
                    SysBeep(5);
                    goto GET_RESPONSE2;
                }    
            break;
        case 'n':
	  //	    printf("\n");
            Micro_Print = FALSE;
            Micro_Print_F = FALSE;
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE1;
        }
}

Find_Level_Transformations_Of_The_Initial_Presentation()
{
    if(Length == 0L)
        {
        printf("\n\nThis is an empty presentation!!");
        return(0);
        }
    ReadPres = NumFilled;    
    if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG)
        {
        printf("\n\nPresentation %d is too long!",ReadPres + 1);
        return(0);
        }
    if(Automorphisms)
        {
        printf("\n\nPresentation %d does not have minimal length.",ReadPres + 1);
        printf("\n\nThe program will only find level transformations for presentations ");
        printf("that have minimal length.");
        return(0);
        }    
    printf("\nPRINT STABILIZER INFO INTO 'Heegaard_Results' ? HIT 'y' OR 'n'.");
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
            SysBeep(5);
            goto GET_RESPONSE1;        
        }
    printf("\n\nListing members of the orbit of presentation 1 under level transformations. . .\n");
    Num_Level_Transformations = 1L;
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
        if(Find_Level_Transformations(FALSE,TRUE) == 5)
            {
            printf("\n\nWe have run out of memory set aside for presentations!!!!");
            Compute_Stabilizers = FALSE;
            break;
            }
        if(mykbhit() == ' ')
            {
            printf("\n\n   There are currently %u unexamined presentations.",OnStack);            
            printf("\n      HIT 't' TO TERMINATE THIS RUN.");
            printf("\n         HIT 'r' TO RESUME LOOKING FOR LEVEL TRANSFORMATIONS.");
            GET_RESPONSE2:
            switch(WaitkbHit())
                {
                case 't':
                    goto _STOP_LT;
                case 'r':
                    if(Compute_Stabilizers)
                        {
                        printf("\nCONTINUE TO PRINT STABILIZER INFO ?  HIT 'y' OR 'n'.\n");
                        GET_RESPONSE3:
                        switch(WaitkbHit())
                            {
                            case 'y':
                                break;
                            case 'n':
                                Compute_Stabilizers = FALSE;
                                break;
                            default:
                                SysBeep(5);
                                goto GET_RESPONSE3;    
                            }
                        }
                    break;                        
                default:
                    SysBeep(5);
                    goto GET_RESPONSE2;
                }
            }        
        ReadPres ++;
        }
    while(ReadPres < NumFilled);
    SysBeep(5);        
_STOP_LT:        
    Free_Memory_For_Find_Level_Transformations(TRUE,1000);
    Compute_Stabilizers = FALSE;
    if(ReadPres < NumFilled)
        {
        printf("\n\nThe search for level transformations was interrupted.");
        fprintf(myout,"\n\nThe search for level transformations was interrupted.");
        }
    switch(NumFilled - 1)
        {
        case 0:
            printf("\n\nThe orbit of presentation 1 under level transformations has 1 member.\n");
            fprintf(myout,"\n\nThe orbit of presentation 1 under level transformations has 1 member.\n");
            break;
        default:
            printf("\n\nThe orbit of presentation 1 under level transformations has %u members.\n",
                NumFilled);
            fprintf(myout,"\n\nThe orbit of presentation 1 under level transformations has %u members.\n",
                NumFilled);
            break;
        }
    
    printf("\n\nHIT 'v' TO REVIEW THE PRESENTATIONS IN THE ORBIT.");    
    printf("\nHIT 'f' TO FILE THESE PRESENTATIONS IN THE FILE 'Heegaard_Results'.");
    printf("\nHIT ANY OTHER KEY TO CONTINUE.");
    switch(WaitkbHit())
        {
        case 'v':
            REVIEW:
            Report(0,0,0,0,0,0,0,0,1,0);
            printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
            GET_RESPONSE4:
            switch(WaitkbHit())
                {
                case 'y':
                    goto REVIEW;
                case 'n':
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE4;
                }
            break;
        case 'f':
            Report(0,0,0,0,1,1,0,0,1,0);
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
                    SysBeep(5);
                    goto GET_RESPONSE6;
                }
        default:
            SysBeep(5);
            goto GET_RESPONSE5;
        }    
}                        

Reduce_The_Initial_Presentation_To_Minimal_Length()
{
    int        i;
    
    long    Scratch;

    Turn_Micro_Print_On();
    Scratch = Length;
    i = NumGenerators;
    switch(Find_Flow_A(NORMAL,FALSE))
        {
        case 1:
            Rewrite_Input();
            break;
        case TOO_LONG:
            printf("\n\n     This presentation may be too long for the program to handle. Sorry!");
            return(0);    
        }        
    if(Length == Scratch)
        printf("\n\nThis presentation has minimal length.");
    else
        {
        printf("\n\n%lu automorphism(s) reduced the length from %ld to %lu.\n",
            Automorphisms,Scratch,Length);
        if(i > NumGenerators)
            printf("\nand reduced the number of generators from %d to %d.\n",
            i,NumGenerators);    
        Canonical_Rewrite(Relators,FALSE,FALSE);
        if(Save_Pres(ReadPres,0,Length,1,2,1,0,0)) return(0);
        UDV[NumFilled - 1] = 0;
        Report(0,0,0,0,0,0,0,0,1,0);
        printf("\n\nDISCARD THIS MINIMAL LENGTH PRESENTATION ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:
        switch(WaitkbHit())
            {
            case 'y':
                return(0);
            case 'n':
                return(1);
            default:
                SysBeep(5);
                goto GET_RESPONSE1;
            }
        }
    return(0);        
}

Initial_Realizability_Check()
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
        fprintf(myout,"\n\n     The Whitehead graph is connected, has no cut-vertices and is non-planar.");
        fprintf(myout,"\n     This is impossible if the presentation is realizable.");
        return(TRUE);
        }
    return(FALSE);    
}
        
void Edit_MyOut(void)
{
    fflush(myout);
    fclose(myout);
    printf("\n\nThe file 'Heegaard_Results' is now closed.");
    printf("\nYou may now transfer to an external program and edit this file.");
    printf("\nWHEN DONE EDITING, HIT ANY KEY TO RESUME RUNNING HEEGAARD.");
    WaitkbHit();
    if((myout = fopen("Heegaard_Results","a+")) == NULL)
        printf("\nUnable to open the file 'Heegaard_Results'.");
}

void Print_Realizability(int Del_Only_Triv_Rel, unsigned int WhichPres)
{
    printf("\n\n                    The data appears consistent.");
    fprintf(myout,"\n\n                    The data appears consistent.");            
    if(Del_Only_Triv_Rel)
        {
        printf("\n\n                    The initial presentation is realizable.");
        fprintf(myout,"\n\n                    The initial presentation is realizable.");
        }
    else
        {    
        printf("\n\n                    Presentation %d is realizable.",WhichPres);
        fprintf(myout,"\n\n                    Presentation %d is realizable.",WhichPres);
        }
}

void Realization_Warning(FILE * fptr)
{
    fprintf(fptr,"\n\n    NOTE: The program has not directly verified that the initial presentation is realizable.");
    fprintf(fptr,"\n          This means it is possible that the initial presentation is not realizable,");
    fprintf(fptr,"\n          even though 'derived' presentations are realizable.");  
}

char mykbhit()
{
#ifdef MAC
    char            charhit;
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
    newtermio.c_cc[VTIME]=2;
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
