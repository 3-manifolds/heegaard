#include "Heegaard.h"
#include "Heegaard_Dec.h"

int     *Table;

void Sort_Presentations_In_Memory(void)
{    
    int         i;    
    int         qkst_compare();
    void        qkst_swap();
    
    Table = (int*) NewPtr((sizeof(int)*NumFilled));
    if(Table == NULL) return;
    
    for(i = 0; i < NumFilled; i++) Table[i] = i;
    
    qksort(NumFilled);
    
    printf("\n\nHIT 'v' TO REVIEW THESE SORTED PRESENTATIONS.");
    if(NoReport)
        printf("\nHIT 's' TO SAVE THESE SORTED PRESENTATIONS TO THE FILE 'Heegaard_Results'.");
    printf("\nOR HIT ANY OTHER KEY TO CONTINUE.");
    switch(WaitkbHit())
        {
        case 'v':
            qksort_Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,NULL);
            break;
        case 's':
            if(NoReport)
            qksort_Report(Band_Sums,NumDiagrams,OnStack,0,1,1,1,0,1,NULL);
            break;    
        default:
            break;        
        }
    
    DisposePtr((char *) Table);
}

int    qkst_compare(int i,int j)
{
    register unsigned char    *p,
                            *q;
                            
    int                        Ti,
                            Tj,
                            k;
                            
                
    Ti = Table[i];
    Tj = Table[j];
    if(Ti == Tj) return(0);
    if(ComponentNum[Ti] > ComponentNum[Tj]) return(-1);
    if(ComponentNum[Ti] < ComponentNum[Tj]) return(1);
    if(NG[Ti] < NG[Tj]) return(-1);
    if(NG[Ti] > NG[Tj]) return(1);
    if(NR[Ti] < NR[Tj]) return(-1);
    if(NR[Ti] > NR[Tj]) return(1);
    if(SURL[Ti] < SURL[Tj]) return(-1);
    if(SURL[Ti] > SURL[Tj]) return(1);
    for(k = 1; k <= NR[Ti]; k++)
        {
        if(GetHandleSize((char **) SUR[Ti][k]) > GetHandleSize((char **) SUR[Tj][k])) return(-1);
        if(GetHandleSize((char **) SUR[Ti][k]) < GetHandleSize((char **) SUR[Tj][k])) return(1);
        }
    for(k = 1; k <= NR[Ti]; k++)
        {
        p = *SUR[Ti][k];
        q = *SUR[Tj][k];
        while(*p && *p == *q)
            {
            p++;
            q++;
            }
        if(*p < *q) return(-1);
        if(*p > *q) return(1);
        }    
    if(Ti < Tj) return(-1);
    if(Ti > Tj) return(1);
    return(0);
}

void qkst_swap(i,j)
int            i,
            j;
{
    int            Temp;
    
    Temp        = Table[i];
    Table[i]     = Table[j];
    Table[j]     = Temp;
}            

int qksort_Report(long Band_Sums,long NumDiagrams,unsigned int OnStack,unsigned int Starting_Pres,
                unsigned int Flag1,unsigned int Flag2,unsigned int Flag3,unsigned int Flag4,
                unsigned int Flag5,unsigned char * Ptr1)
{
    /******************************************************************************************
        This variant of Report() is an output routine called on termination of the sorting
        routine which sorts the presentations presently in memory.
    ******************************************************************************************/    

    char            c;
            
    unsigned char     *p,
                    x,
                    y;
    
    int                NumRelators,
                    One_By_One,
                    SNumFilled,
                    SStarting_Pres;
                    
    unsigned int     h,
                    i,
                    j,
                    k,
                    m,
                    n;                                
    
    unsigned long    Length;
    
    if(Flag5) Update_Bdry_Data();

LIST_PRESENTATIONS:
    
    One_By_One = FALSE;
    ObscureCursor();
    for(k = Flag1; k <= Flag2; k++)
    {
    if(k == 0)
        fptr = stdout;
    else
        fptr = myout;
    
    if(Flag5)            
    fprintf(fptr,"\n\n<------------------------------------   REPORT   ---------------------------------------->");
    
    if(k == 0 && NumFilled > 1 && Flag5)
        {
        printf("\n\n      Hit 'return' to stop reviewing this report.");
        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
        printf("\n      Use the 'd','i','n','p' and 's' keys to: view diagrams, save a presentation in");
        printf("\n      'Input_Presentations', see the next presentation, see the previous presentation,");
        printf("\n      or save a presentation in 'Heegaard_Results'.");
        }
        
    if(BdryData && Flag5)
        {
        fprintf(fptr,"\n");
        Print_Bdry_Data(Starting_Pres);
        }

    fprintf(fptr,"\n\nThe initial presentation was: %s",PresName);    
    
    for(n = Starting_Pres; n < NumFilled; n++)
        {
        i = Table[n];
        if(Flag4 == 1 && !Ptr1[i]) continue;
        if(Flag4 && Ptr1[i] != Flag4) continue;
        NumRelators = NR[i];
        Length = SURL[i];
        fprintf(fptr,"\n\nPresentation %d  of Summand %u:  %d Generator(s)  Length  %lu  From Presentation %u  ",
        i+1,ComponentNum[i],NG[i],Length,FR[i]+1);
        
        switch(PRIM[i])
            {
            case 1:
            case 101:
                fprintf(fptr,"DR");
                break;
            case 2:
            case 102:
                fprintf(fptr,"IP");
                break;
            case 3:
            case 103:
                fprintf(fptr,"LS");
                break;
            case 4:
            case 104:
                fprintf(fptr,"1G");
                break;
            case 5:
            case 105:
                fprintf(fptr,"S3");
                break;
            case 6:
            case 106:
                fprintf(fptr,"FP");
                break;
            case 7:
            case 107:
                fprintf(fptr,"BC");
                break;
            case 8:
            case 108:
                fprintf(fptr,"PM");
                break;
            case 9:
            case 109:
                fprintf(fptr,"PM");
                break;
            case 10:
            case 110:
                fprintf(fptr,"BC");
                break;
            case 11:
            case 111:
                fprintf(fptr,"CF");
                break;
            case 12:
            case 112:
                fprintf(fptr,"ER");
                break;
            case 13:
            case 113:
                fprintf(fptr,"Er");
                break;                    
            case 20:
            case 120:
                fprintf(fptr,"NC");
                break;
            case 30:
            case 130:
            case 40:
            case 140:
                fprintf(fptr,"MG");
                break;
            case 60:
            case 160:
                fprintf(fptr,"PP");
                break;    
            case 70:
            case 170:
                fprintf(fptr,"Lt");
                break;
            case 75:
            case 175:
                fprintf(fptr,"LT");
                break;    
            case 80:
            case 180:
                fprintf(fptr,"A2");
                break;    
            default:
                break;
            }
                
        switch(UDV[i])
            {
            case SPLIT:
                j = Daughters[i];
                m = j + NCS[i] - 1;
                switch(PRIM[j])
                    {
                    case 40:
                    case 140:
                        fprintf(fptr,"\nThe presentation dual to presentation %d 'split' into presentations %u through %u of M",
                                    i + 1,j + 1,m + 1);
                        fprintf(fptr,"\nwhich correspond to summands %u through %u of M.",ComponentNum[j],ComponentNum[m]);                    
                        break;
                    default:
                        fprintf(fptr,"\n'split' into presentations %u through %u corresponding to summands %u through %u of M.",
                        j + 1,m + 1,ComponentNum[j],ComponentNum[m]);
                        break;                    
                    }
                if(BdryData) for(m = 0; m < NCS[i]; m++) Print_Bdry_Data(j + m);
                break;
                
            case GENERIC_LENS_SPACE:
                if(LSP[i] > 4L)
                    fprintf(fptr,"\npresents a Lens space of the form L(%lu,Q).",LSP[i]);
                else
                if(LSP[i] == 1L)
                    fprintf(fptr,"\npresents the 3-Sphere.");
                else
                    fprintf(fptr,"\npresents a Lens space of the form L(%lu,1).",LSP[i]);
                break;        
            
            case THREE_SPHERE:
                fprintf(fptr,"\npresents the 3-Sphere.");
                break;
            
            case NOT_CONNECTED:
                fprintf(fptr,"\nthe diagram is not connected.");
                break;
                
            case S1_X_S2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X S2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X S2.",NG[i]);
                break;
            
            case S1_X_D2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X D2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X D2.",NG[i]);
                break;
            
            case S1_X_X2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X ?2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X ?2.",NG[i]);
                break;

            case MISSING_GEN_DONE2:
                fprintf(fptr,"\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of I X D2, ");
                else
                    fprintf(fptr,"%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of S1 X S2, ");
                else
                    fprintf(fptr,"%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    fprintf(fptr,"and 1 copy of S1 X D2.");
                else
                    fprintf(fptr,"and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                                        
                fprintf(fptr,"\nThe program was not able to unambiguously determine the boundary components.");
                break;
                
            case MISSING_GEN_DONE1:
                fprintf(fptr,"\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of I X D2, ");
                else
                    fprintf(fptr,"%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of S1 X S2, ");
                else
                    fprintf(fptr,"%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    fprintf(fptr,"and 1 copy of S1 X D2.");
                else
                    fprintf(fptr,"and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                            
                break;
                
            case KNOWN_LENS_SPACE:
                switch(LSP[i])
                    {
                    case 0L:
                        fprintf(fptr,"\npresents S1 X S2.");
                        break; 
                    case 1L:
                        fprintf(fptr,"\npresents the 3-Sphere.");
                        break;
                    default:
                        fprintf(fptr,"\npresents the Lens space L(%lu,%lu).",LSP[i],LSQ[i]);
                        break;
                    }
                break;
            
            case SEP_PAIRS:
                if(PRIM[i] >= 100)
                    {
                    fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                    fprintf(fptr," Vertex '%c' and vertex '%c' separate.",    
                        x = LSP[i],y = LSQ[i]);                    
                    }
                else                        
                    fprintf(fptr,"\nvertex '%c' and vertex '%c' separate the diagram.",
                    x = LSP[i],y = LSQ[i]);
                break;
            
            case ANNULUS_EXISTS:
                p = *SUR[i][0];
                x = *p++;
                y = *p++;                
                fprintf(fptr,"\nVertices '%c' and '%c' separate the diagram.",x,y);
                fprintf(fptr,"\nThe component consisting of vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    fprintf(fptr,",%c",x);
                    p++;
                    }
                fprintf(fptr,"}");    
                p++;        
                fprintf(fptr,"\nlies in an annulus which swallows the component and otherwise follows the curve:\n");
                h = 0;
                while(*p)
                    {
		    fputc(*p++,fptr);
                    if(h++ >= 90 && *p)
                        {
                        fprintf(fptr,"\n");
                        h = 0;
                        }
                    }
                break;
            
            case V2_ANNULUS_EXISTS:
                p = *SUR[i][0];
                fprintf(fptr,"\nThere exists an annulus which swallows vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    fprintf(fptr,",%c",x);
                    p++;
                    }
                fprintf(fptr,"}");    
                p++;        
                fprintf(fptr,"\nand otherwise follows the curve:\n");
                h = 0;
                while(*p)
                    {
		    fputc(*p++,fptr);
                    if(h++ >= 90 && *p)
                        {
                        fprintf(fptr,"\n");
                        h = 0;
                        }
                    }
                break;                
            
            case DELETED_RELATOR:
                break;
            
            case NON_UNIQUE_4:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nwith only one exponent and that exponent is greater than 6.");
                break;
            
            case NON_UNIQUE_3:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nonly with exponent 5.");
                break;
            
            case NON_UNIQUE_2:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nonly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
                break;
            
            case NON_UNIQUE_1:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
                break;
            
            case DUPLICATE:
                fprintf(fptr,"\nis a duplicate of presentation %d of summand %u.",
                Daughters[i] + 1,ComponentNum[Daughters[i]]);
                break;                                                                                                                break;        
            
            default:
                {
                j = PRIM[i];
                switch(j)
                    {
                    case 8:
                    case 108:
                        fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 70:
                        if(QPM[i])
                            {
                            fprintf(fptr,"\nvia level transformations of sep_pairs,");
                            fprintf(fptr," is quasi-pseudo-minimal.");
                            }
                        else        
                            fprintf(fptr,"\nvia level transformations of sep_pairs.");
                        break;
                    case 75:
                        if(QPM[i])
                            {
                            fprintf(fptr,"\nvia general level transformations,");
                            fprintf(fptr," is quasi-pseudo-minimal.");
                            }
                        else
                            fprintf(fptr,"\nvia general level transformations.");
                        break;    
                    case 170:
                        fprintf(fptr,"\nvia level transformations of sep_pairs and ");
                        fprintf(fptr,"dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 175:
                        fprintf(fptr,"\nvia general level transformations and ");
                        fprintf(fptr,"dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;    
                    default:
                        if(j >= 100)
                            fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                        if(QPM[i])
                            fprintf(fptr,"\nis quasi-pseudo-minimal.");    
                        break;
                    }
                break;
                }                                                                                
            }
                
        fprintf(fptr,"\n");
        for(j = 1; j <= NumRelators; j++)
            {    
            HLock((char **) SUR[i][j]);
            fprintf(fptr,"\nR%3u)    ",j);
            p = *SUR[i][j];
            h = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(h++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    h = 9;
                    }
                }
            HUnlock((char **) SUR[i][j]);
            }
        
        GET_RESPONSE1:
        if(k == 0 && Flag5 && One_By_One && (c = mykbhit())) switch(c)
            {
            case ' ':
                One_By_One = FALSE;
                ObscureCursor();
                break;
            case 'd':
                WhichInput = i;
                if(Display_A_Diagram())
                    {
                    if(n > Starting_Pres) n -= 2;
                    else n --;                    
                    }
                fptr = stdout;
                ObscureCursor();
                break;
            case 'i':
                Save_Pres_To_Input_Presentations(Table[n]);
                ObscureCursor();
                goto GET_RESPONSE1;                    
            case 'n':
                ObscureCursor();
                break;
            case 'p':
                if(n > Starting_Pres) n -= 2;
                else
                    {
                    n --;
                    SysBeep(5);
                    }
                ObscureCursor();
                break;        
            case '\n':
                return(NO_ERROR);
            case '\r':
                return(NO_ERROR);    
            case 's':
                SNumFilled = NumFilled;
                SStarting_Pres = Starting_Pres;
                Starting_Pres = n;
                NumFilled = Starting_Pres + 1;
                qksort_Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,
                    1,1,Flag3,Flag4,0,Ptr1);                        
                NumFilled = SNumFilled;
                Starting_Pres = SStarting_Pres;
                printf("\n\n    Saved Presentation %d in 'Heegaard_Results'.",i + 1);
                fptr = stdout;
                ObscureCursor();
                goto GET_RESPONSE1;                    
            default:
                SysBeep(5);
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);            
                printf("\n      Hit 'n' to see the next presentation.");
                if(n > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to start scrolling.");
                goto GET_RESPONSE1;                    
            }        
            
        if(k == 0 && Flag5 && (c = mykbhit())) switch(c)
            {
            case ' ':
                WAIT:
                switch(WaitkbHit())
                    {
                    case ' ':
                        ObscureCursor();
                        break;
                    case 'd':
                        WhichInput = i;
                        if(Display_A_Diagram())
                            {
                            if(n > Starting_Pres) n -= 2;
                            else n --;                            
                            }
                        fptr = stdout;
                        ObscureCursor();
                        One_By_One = TRUE;
                        break;
                    case 'i':
                        Save_Pres_To_Input_Presentations(Table[n]);
                        ObscureCursor();
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;                                                    
                    case 'n':
                        One_By_One = TRUE;
                        ObscureCursor();
                        break;
                    case 'p':
                        if(n > Starting_Pres) n -= 2;
                        else
                            {
                            n --;
                            SysBeep(5);
                            }
                        One_By_One = TRUE;
                        ObscureCursor();
                        break;                                    
                    case '\n':
                        return(NO_ERROR);
                    case '\r':
                        return(NO_ERROR);    
                    case 's':
                        SNumFilled = NumFilled;
                        SStarting_Pres = Starting_Pres;
                        Starting_Pres = n;
                        NumFilled = Starting_Pres + 1;
                        qksort_Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,
                            1,1,Flag3,Flag4,0,Ptr1);                        
                        NumFilled = SNumFilled;
                        Starting_Pres = SStarting_Pres;
                        printf("\n\n    Saved Presentation %d in 'Heegaard_Results'.",i + 1);
                        fptr = stdout;
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;                        
                    default:
                        SysBeep(5);
                        printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                        printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);            
                        printf("\n      Hit 'n' to see the next presentation.");
                        if(n > Starting_Pres)
                        printf("\n      Hit 'p' to see the previous presentation.");
                        printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                        printf("\n      Hit 'return' to stop reviewing this report.");
                        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;    
                    }
                break;
            case 'n':
            case 'p':
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);            
                printf("\n      Hit 'n' to see the next presentation.");
                if(n > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");            
                One_By_One = TRUE;
                goto GET_RESPONSE1;                
            case '\n':
                if(Flag2 == 0)
                    return(NO_ERROR);
                else
                    {
                    k = 1;
                    goto LIST_PRESENTATIONS;
                    }
            case '\r':
                if(Flag2 == 0)
                    return(NO_ERROR);
                else
                    {
                    k = 1;
                    goto LIST_PRESENTATIONS;
                    }        
            default:
                SysBeep(5);
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);            
                printf("\n      Hit 'n' to see the next presentation.");
                if(n > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                One_By_One = TRUE;
                goto GET_RESPONSE1;        
            }            
        }                    
    
    if(Flag3 && Flag5)
        {
        fprintf(fptr,"\n\nThe program performed %lu automorphism(s), examined %ld bandsum(s),",
            TotalAuts,Band_Sums);
        fprintf(fptr,"\nexamined %ld diagram(s), and dualized %lu diagram(s).",
            NumDiagrams,NumDualized);
        fprintf(fptr," Left to do %u.\n",OnStack);
        }
    }
    if(Flag2 && Flag5)
        {
        printf("\n\nThese results are printed in the file 'Heegaard_Results'.");    
        NoReport = FALSE;
        }
    
    if(Flag5)
        {    
        printf("\n\nHIT 'v' TO REVIEW THESE PRESENTATIONS.");
        if(NoReport)
            printf("\nHIT 's' TO SAVE THESE PRESENTATIONS TO THE FILE 'Heegaard_Results'.");
        printf("\nOR HIT ANY OTHER KEY TO CONTINUE.");
        switch(WaitkbHit())
            {
            case 'v':
                Flag1 = 0;
                Flag2 = 0;
                Flag3 = 1;
                goto LIST_PRESENTATIONS;    
            case 's':
                if(NoReport)
                    {
                    Flag1 = 1;
                    Flag2 = 1;
                    Flag3 = 1;
                    goto LIST_PRESENTATIONS;
                    }
                break;    
            default:
                break;        
            }
        }
    return(NO_ERROR);    
}

int Report(long Band_Sums,long NumDiagrams,unsigned int OnStack,unsigned int Starting_Pres,unsigned int Flag1,
            unsigned int Flag2,unsigned int Flag3,unsigned int Flag4,unsigned int Flag5,unsigned char * Ptr1)
{
    /******************************************************************************************
        Report() is an output routine called automatically on termination of the program. It
        can also be called by the user when the program has been interrupted.
    ******************************************************************************************/    

    char            c;
            
    unsigned char     *p,
                    x,
                    y;
    
    int                NumRelators,
                    One_By_One,
                    SNumFilled,
                    SStarting_Pres;
                    
    unsigned int     h,
                    i,
                    j,
                    k,
                    m;                                
    
    unsigned long     Length;
    
    if(Flag5) Update_Bdry_Data();

LIST_PRESENTATIONS:
    
    One_By_One = FALSE;
    ObscureCursor();
    for(k = Flag1; k <= Flag2; k++)
    {
    if(k == 0)
        fptr = stdout;
    else
        fptr = myout;    
    
    if(Flag5)            
    fprintf(fptr,"\n\n<------------------------------------   REPORT   ---------------------------------------->");
    
    if(k == 0 && NumFilled > 1 && Flag5)
        {
        printf("\n\n      Hit 'return' to stop reviewing this report.");
        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
        printf("\n      Use the 'd','i','n','p' and 's' keys to: view diagrams, save a presentation in");
        printf("\n      'Input_Presentations', see the next presentation, see the previous presentation,");
        printf("\n      or save a presentation in 'Heegaard_Results'."); 
        }
        
    if(BdryData && Flag5)
        {
        fprintf(fptr,"\n");
        Print_Bdry_Data(Starting_Pres);
        }

    fprintf(fptr,"\n\nThe initial presentation was: %s",PresName);    
            
    for(i = Starting_Pres; i < NumFilled; i++)
        {
        if(Flag4 == 1 && !Ptr1[i]) continue;
        if(Flag4 && Ptr1[i] != Flag4) continue;
        NumRelators = NR[i];
        Length = SURL[i];
        fprintf(fptr,"\n\nPresentation %d  of Summand %u:  %d Generator(s)  Length  %lu  From Presentation %u  ",
        i+1,ComponentNum[i],NG[i],Length,FR[i]+1);
        
        switch(PRIM[i])
            {
            case 1:
            case 101:
                fprintf(fptr,"DR");
                break;
            case 2:
            case 102:
                fprintf(fptr,"IP");
                break;
            case 3:
            case 103:
                fprintf(fptr,"LS");
                break;
            case 4:
            case 104:
                fprintf(fptr,"1G");
                break;
            case 5:
            case 105:
                fprintf(fptr,"S3");
                break;
            case 6:
            case 106:
                fprintf(fptr,"FP");
                break;
            case 7:
            case 107:
                fprintf(fptr,"BC");
                break;
            case 8:
            case 108:
                fprintf(fptr,"PM");
                break;
            case 9:
            case 109:
                fprintf(fptr,"PM");
                break;
            case 10:
            case 110:
                fprintf(fptr,"BC");
                break;
            case 11:
            case 111:
                fprintf(fptr,"CF");
                break;
            case 12:
            case 112:
                fprintf(fptr,"ER");
                break;
            case 13:
            case 113:
                fprintf(fptr,"Er");
                break;                
            case 20:
            case 120:
                fprintf(fptr,"NC");
                break;
            case 30:
            case 130:
            case 40:
            case 140:
                fprintf(fptr,"MG");
                break;
            case 60:
            case 160:
                fprintf(fptr,"PP");
                break;    
            case 70:
            case 170:
                fprintf(fptr,"Lt");
                break;
            case 75:
            case 175:
                fprintf(fptr,"LT");
                break;    
            case 80:
            case 180:
                fprintf(fptr,"A2");
                break;    
            default:
                break;
            }
                
        switch(UDV[i])
            {
            case SPLIT:
                j = Daughters[i];
                m = j + NCS[i] - 1;
                switch(PRIM[j])
                    {
                    case 40:
                    case 140:
                        fprintf(fptr,"\nThe presentation dual to presentation %d 'split' into presentations %u through %u of M",
                                    i + 1,j + 1,m + 1);
                        fprintf(fptr,"\nwhich correspond to summands %u through %u of M.",ComponentNum[j],ComponentNum[m]);                    
                        break;
                    default:
                        fprintf(fptr,"\n'split' into presentations %u through %u corresponding to summands %u through %u of M.",
                        j + 1,m + 1,ComponentNum[j],ComponentNum[m]);
                        break;                    
                    }
                if(BdryData) for(m = 0; m < NCS[i]; m++) Print_Bdry_Data(j + m);
                break;
            
            case GENERIC_LENS_SPACE:
                if(LSP[i] > 4L)
                    fprintf(fptr,"\npresents a Lens space of the form L(%lu,Q).",LSP[i]);
                else
                if(LSP[i] == 1L)
                    fprintf(fptr,"\npresents the 3-Sphere.");
                else
                    fprintf(fptr,"\npresents a Lens space of the form L(%lu,1).",LSP[i]);
                break;        
            
            case THREE_SPHERE:
                fprintf(fptr,"\npresents the 3-Sphere.");
                break;

            case NOT_CONNECTED:
                fprintf(fptr,"\nthe diagram is not connected.");
                break;
                    
            case S1_X_S2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X S2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X S2.",NG[i]);
                break;
            
            case S1_X_D2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X D2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X D2.",NG[i]);
                break;
            
            case S1_X_X2:
                if(NG[i] == 1)
                    fprintf(fptr,"\npresents S1 X ?2.");
                else
                    fprintf(fptr,"\npresents %d copies of S1 X ?2.",NG[i]);
                break;
            
            case MISSING_GEN_DONE2:
                fprintf(fptr,"\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of I X D2, ");
                else
                    fprintf(fptr,"%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of S1 X S2, ");
                else
                    fprintf(fptr,"%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    fprintf(fptr,"and 1 copy of S1 X D2.");
                else
                    fprintf(fptr,"and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                                        
                fprintf(fptr,"\nThe program was not able to unambiguously determine the boundary components.");
                break;
                
            case MISSING_GEN_DONE1:
                fprintf(fptr,"\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of I X D2, ");
                else
                    fprintf(fptr,"%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    fprintf(fptr,"1 copy of S1 X S2, ");
                else
                    fprintf(fptr,"%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    fprintf(fptr,"and 1 copy of S1 X D2.");
                else
                    fprintf(fptr,"and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                            
                break;
                                
            case KNOWN_LENS_SPACE:
                switch(LSP[i])
                    {
                    case 0L:
                        fprintf(fptr,"\npresents S1 X S2.");
                        break; 
                    case 1L:
                        fprintf(fptr,"\npresents the 3-Sphere.");
                        break;
                    default:
                        fprintf(fptr,"\npresents the Lens space L(%lu,%lu).",LSP[i],LSQ[i]);
                        break;
                    }
                break;
            
            case SEP_PAIRS:
                if(PRIM[i] >= 100)
                    {
                    fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                    fprintf(fptr," Vertex '%c' and vertex '%c' separate.",
                        x = LSP[i],y = LSQ[i]);                    
                    }
                else                        
                    fprintf(fptr,"\nvertex '%c' and vertex '%c' separate the diagram.",
                    x = LSP[i],y = LSQ[i]);            
                break;
            
            case ANNULUS_EXISTS:
                p = *SUR[i][0];
                x = *p++;
                y = *p++;                
                fprintf(fptr,"\nVertices '%c' and '%c' separate the diagram.",x,y);
                fprintf(fptr,"\nThe component consisting of vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    fprintf(fptr,",%c",x);
                    p++;
                    }
                fprintf(fptr,"}");    
                p++;        
                fprintf(fptr,"\nlies in an annulus which swallows the component and otherwise follows the curve:\n");
                h = 0;
                while(*p)
                    {
		    fputc(*p++,fptr);
                    if(h++ >= 90 && *p)
                        {
                        fprintf(fptr,"\n");
                        h = 0;
                        }
                    }
                break;
            
            case V2_ANNULUS_EXISTS:
                p = *SUR[i][0];
                fprintf(fptr,"\nThere exists an annulus which swallows vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    fprintf(fptr,",%c",x);
                    p++;
                    }
                fprintf(fptr,"}");    
                p++;        
                fprintf(fptr,"\nand otherwise follows the curve:\n");
                h = 0;
                while(*p)
                    {
		    fputc(*p++,fptr);
                    if(h++ >= 90 && *p)
                        {
                        fprintf(fptr,"\n");
                        h = 0;
                        }
                    }
                break;                
            
            case DELETED_RELATOR:
                break;
            
            case NON_UNIQUE_4:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nwith only one exponent and that exponent is greater than 6.");
                break;
            
            case NON_UNIQUE_3:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nonly with exponent 5.");
                break;
            
            case NON_UNIQUE_2:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nonly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
                break;
            
            case NON_UNIQUE_1:
                fprintf(fptr,"\nthe diagram is not unique because there is a generator which appears");
                fprintf(fptr,"\nwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
                break;
            
            case DUPLICATE:
                fprintf(fptr,"\nis a duplicate of presentation %d of summand %u.",
                Daughters[i] + 1,ComponentNum[Daughters[i]]);
                break;                                                                                                                break;        
            
            default:
                {
                j = PRIM[i];
                switch(j)
                    {
                    case 8:
                    case 108:
                        fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 70:
                        if(QPM[i])
                            {
                            fprintf(fptr,"\nvia level transformations of sep_pairs,");
                            fprintf(fptr," is quasi-pseudo-minimal.");
                            }
                        else        
                            fprintf(fptr,"\nvia level transformations of sep_pairs.");
                        break;
                    case 75:
                        if(QPM[i])
                            {
                            fprintf(fptr,"\nvia general level transformations,");
                            fprintf(fptr," is quasi-pseudo-minimal.");
                            }
                        else
                            fprintf(fptr,"\nvia general level transformations.");
                        break;    
                    case 170:
                        fprintf(fptr,"\nvia level transformations of sep_pairs and ");
                        fprintf(fptr,"dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 175:
                        fprintf(fptr,"\nvia general level transformations and ");
                        fprintf(fptr,"dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;    
                    default:
                        if(j >= 100)
                            fprintf(fptr,"\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                        if(QPM[i])
                            fprintf(fptr,"\nis quasi-pseudo-minimal.");                                
                        break;
                    }
                break;
                }                                                                                
            }
                
        fprintf(fptr,"\n");
        for(j = 1; j <= NumRelators; j++)
            {    
            HLock((char **) SUR[i][j]);
            fprintf(fptr,"\nR%3u)    ",j);
            p = *SUR[i][j];
            h = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(h++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    h = 9;
                    }
                }
            HUnlock((char **) SUR[i][j]);
            }
        
        GET_RESPONSE1:
        if(k == 0 && Flag5 && One_By_One && (c = mykbhit())) switch(c)
            {
            case ' ':
                One_By_One = FALSE;
                ObscureCursor();
                break;
            case 'd':
                WhichInput = i;
                if(Display_A_Diagram())
                    {
                    if(i > Starting_Pres) i -= 2;
                    else i --;                            
                    }                
                fptr = stdout;
                ObscureCursor();
                break;
            case 'i':
                Save_Pres_To_Input_Presentations(i);
                ObscureCursor();
                goto GET_RESPONSE1;                        
            case 'n':
                ObscureCursor();
                break;
            case 'p':
                if(i > Starting_Pres) i -= 2;
                else
                    {
                    i--;
                    SysBeep(5);
                    }
                ObscureCursor();
                break;    
            case '\n':
                return(NO_ERROR);
            case 0:
                return(NO_ERROR);        
            case 's':
                SNumFilled = NumFilled;
                SStarting_Pres = Starting_Pres;
                Starting_Pres = i;
                NumFilled = Starting_Pres + 1;
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,
                    1,1,Flag3,Flag4,0,Ptr1);                        
                NumFilled = SNumFilled;
                Starting_Pres = SStarting_Pres;
                printf("\n\n    Saved Presentation %d in 'Heegaard_Results'.",i + 1);
                fptr = stdout;
                ObscureCursor();
                goto GET_RESPONSE1;                    
            default:
                SysBeep(5);
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);
                printf("\n      Hit 'n' to see the next presentation.");
                if(i > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to start scrolling.");
                goto GET_RESPONSE1;                    
            }
                
        if(k == 0 && Flag5 && (c = mykbhit())) switch(c)
            {
            case ' ':
                WAIT:
                switch(WaitkbHit())
                    {
                    case ' ':
                        ObscureCursor();
                        break;
                    case 'd':
                        WhichInput = i;
                        if(Display_A_Diagram())
                            {
                            if(i > Starting_Pres) i -= 2;
                            else i --;                            
                            }    
                        fptr = stdout;
                        ObscureCursor();
                        One_By_One = TRUE;
                        break;
                    case 'i':
                        Save_Pres_To_Input_Presentations(i);
                        ObscureCursor();
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;                                                    
                    case 'n':
                        One_By_One = TRUE;
                        ObscureCursor();
                        break;
                    case 'p':
                        if(i > Starting_Pres) i -= 2;
                        else
                            {
                            i--;
                            SysBeep(5);
                            }
                        ObscureCursor();
                        One_By_One = TRUE;
                        break;                                
                    case '\n':
                        return(NO_ERROR);
                    case '\r':
                        return(NO_ERROR);    
                    case 's':
                        SNumFilled = NumFilled;
                        SStarting_Pres = Starting_Pres;
                        Starting_Pres = i;
                        NumFilled = Starting_Pres + 1;
                        Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,
                            1,1,Flag3,Flag4,0,Ptr1);                        
                        NumFilled = SNumFilled;
                        Starting_Pres = SStarting_Pres;
                        printf("\n\n    Saved Presentation %d in 'Heegaard_Results'.",i + 1);
                        fptr = stdout;
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;    
                    default:
                        SysBeep(5);
                        printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                        printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);
                        printf("\n      Hit 'n' to see the next presentation.");
                        if(i > Starting_Pres)
                        printf("\n      Hit 'p' to see the previous presentation.");
                        printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                        printf("\n      Hit 'return' to stop reviewing this report.");
                        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                        One_By_One = TRUE;
                        goto GET_RESPONSE1;        
                    }
                break;
            case 'n':
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);
                printf("\n      Hit 'n' to see the next presentation.");
                if(i > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");                
                One_By_One = TRUE;
                goto GET_RESPONSE1;                        
            case '\n':
                if(Flag2 == 0)
                    return(NO_ERROR);
                else
                    {
                    k = 1;
                    goto LIST_PRESENTATIONS;
                    }
            case '\r':
                if(Flag2 == 0)
                    return(NO_ERROR);
                else
                    {
                    k = 1;
                    goto LIST_PRESENTATIONS;
                    }        
            default:
                SysBeep(5);
                printf("\n\n      Hit 'd' to see the diagram of presentation %d.",i + 1);
                printf("\n      Hit 'i' to save presentation %d to 'Input_Presentations'.",i + 1);
                printf("\n      Hit 'n' to see the next presentation.");
                if(i > Starting_Pres)
                printf("\n      Hit 'p' to see the previous presentation.");
                printf("\n      Hit 's' to save presentation %d to 'Heegaard_Results'.",i + 1);
                printf("\n      Hit 'return' to stop reviewing this report.");
                printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
                One_By_One = TRUE;
                goto GET_RESPONSE1;            
            }            
        }                    
    
    if(Flag3 && Flag5)
        {
        fprintf(fptr,"\n\nThe program performed %lu automorphism(s), examined %ld bandsum(s),",
            TotalAuts,Band_Sums);
        fprintf(fptr,"\nexamined %ld diagram(s), and dualized %lu diagram(s).",
            NumDiagrams,NumDualized);
        fprintf(fptr," Left to do %u.\n",OnStack);
        }
    }
    if(Flag2 && Flag5)
        {
        printf("\n\nThese results are printed in the file 'Heegaard_Results'.");    
        NoReport = FALSE;
        }
    
    if(Flag5)
        {    
        printf("\n\nHIT 'v' TO REVIEW THESE PRESENTATIONS.");
        if(NoReport)
            printf("\nHIT 's' TO SAVE THESE PRESENTATIONS TO THE FILE 'Heegaard_Results'.");
        printf("\nOR HIT ANY OTHER KEY TO CONTINUE.");
        switch(WaitkbHit())
            {
            case 'v':
                Flag1 = 0;
                Flag2 = 0;
                Flag3 = 1;
                goto LIST_PRESENTATIONS;
            case 's':
                if(NoReport)
                    {
                    Flag1 = 1;
                    Flag2 = 1;
                    Flag3 = 1;
                    goto LIST_PRESENTATIONS;
                    }
                break;    
            default:
                break;        
            }
        }
    return(NO_ERROR);    
}        

void Update_Bdry_Data(void)
{
    int            h,
                i,
                j,
                k,
                m,
                n,
                NotClosed,
                NumNotClosed,
                NumUpdates;
    
    if(TotalComp == 1 || NoReport == FALSE) return;
                
    do
        {
        for(i = 0, NumUpdates = 0; i < NumFilled; i++)
            {
            if(UDV[i] == DUPLICATE)
                {
                m = ComponentNum[i];
                h = ComponentNum[Daughters[i]];
                if(CBC[m][0] == BDRY_UNKNOWN && CBC[h][0] < BDRY_UNKNOWN
                    && CBC[h][1] == BDRY_UNKNOWN)
                    {
                    CBC[m][0] = 1;
                    CBC[m][1] = BDRY_UNKNOWN;
                    NumUpdates ++;
                    }
                }
            if(UDV[i] == SPLIT)
                {
                m = ComponentNum[i];
                if(CBC[m][0] < BDRY_UNKNOWN && CBC[m][1] == BDRY_UNKNOWN)
                    {
                    n = ComponentNum[Daughters[i]] + NCS[i];
                    for(j = ComponentNum[Daughters[i]]; j < n; j++)
                    if(CBC[j][0] == BDRY_UNKNOWN)
                        {
                        NumUpdates ++;
                        CBC[j][0] = 1;
                        CBC[j][1] = BDRY_UNKNOWN;
                        }
                    }
                if(CBC[ComponentNum[i]][0] == BDRY_UNKNOWN)
                    {
                    /************************************************************************
                        Check whether all but at most one of the summands of M is closed.
                    ************************************************************************/
                        
                    h = Daughters[i];
                    n = ComponentNum[h] + NCS[i];
                    for(j = ComponentNum[h]; j < n; j++)
                        if(CBC[j][0] == BDRY_UNKNOWN) break;
                    if(j < n) continue;
                    for(j = ComponentNum[h],NumNotClosed = 0; j < n; j++)
                        if(CBC[j][1] != BDRY_UNKNOWN)
                            {
                            NotClosed = j;
                            if(++NumNotClosed > 1)
                            break;
                            }
                    if(j < n) continue;                
                    NumUpdates ++;
                    CBC[m][0] = 1;
                    if(NumNotClosed == 1)
                        for(k = 1; (CBC[m][k] = CBC[NotClosed][k]) < BDRY_UNKNOWN; k++) ;
                    else
                        CBC[m][1] = BDRY_UNKNOWN;    
                    for(j = ComponentNum[h]; j < n; j++) if(CBC[j][0])
                        CBC[m][0] += CBC[j][0] - 1;
                    }    
                }
            }    
        }
    while(NumUpdates);        
}

void Print_Bdry_Data(unsigned int WhichPres)
{
    unsigned int    i,
                    j,
                    k,
                    kk;
    
    k = kk = ComponentNum[WhichPres];
    if(UDV[WhichPres] == DUPLICATE)
        k = ComponentNum[Daughters[WhichPres]];
        
    if(CBC[k][0] == BDRY_UNKNOWN)
        {
        if(TotalComp == 1)
            fprintf(fptr,"\n    The program was unable to determine the boundary components of this manifold.");
        else    
            fprintf(fptr,"\n    The program was unable to determine the boundary components of 'summand' %u of M.",
                kk);
        return;
        }
        
    for(i = 1,j = 0; CBC[k][i] < BDRY_UNKNOWN; i++) j += CBC[k][i];
    switch(j)
        {
        case 0:
            if(TotalComp == 1)
                fprintf(fptr,"\n    This manifold is closed.");
            else
                fprintf(fptr,"\n    'Summand' %u of M is closed.",kk);
            break;
        case 1:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i])
                {
                if(TotalComp == 1)
                    fprintf(fptr,"\n    This manifold has one boundary component of genus %u.",i);                
                else
                    fprintf(fptr,"\n    'Summand' %u of M has one boundary component of genus %u.",
                        kk,i);
                break;
                }
            break;        
        default:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i] == j)
                {
                if(TotalComp == 1)
                    fprintf(fptr,"\n    This manifold has %u boundary components of genus %u.",
                        j,i);                
                else
                    fprintf(fptr,"\n    'Summand' %u of M has %u boundary components of genus %u.",
                        kk,j,i);
                return;
                }
            if(TotalComp == 1)
                fprintf(fptr,"\n    This manifold has the following boundary components:");
            else    
                fprintf(fptr,"\n    'Summand' %u of M has the following boundary components:",kk);
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(j = CBC[k][i])
                {
                if(j == 1)
                    fprintf(fptr,"\n        %2u component  of genus %2u.",j,i);
                else
                    fprintf(fptr,"\n        %2u components of genus %2u.",j,i);
                }
            break;    
        }
}

void Fatal_Error(void)
{
    /******************************************************************************************
        Fatal_Error() is called when the program has discovered that a presentation is not
        realizable by a Heegaard diagram. It prints a message to this effect, and also prints
        the offending presentation.
    ******************************************************************************************/
    
    SysBeep(5);
    BdryData = FALSE;
    if(NR[ReadPres] == NumRelators && Compare_Pres(ReadPres))
        {
        fprintf(stdout,"\n\nPresentation %d is not realizable.\n",
            ReadPres + 1);
        Print_Relators(Relators,NumRelators,stdout);        
        fprintf(myout,"\n\nPresentation %d is not realizable.\n",
            ReadPres + 1);
        }
    else
        {    
        fprintf(stdout,"\n\nThis presentation, obtained from presentation %d, is not realizable.\n",
            ReadPres + 1);
        Print_Relators(Relators,NumRelators,stdout);        
        fprintf(myout,"\n\nThis presentation obtained, from presentation %d, is not realizable.\n",
            ReadPres + 1);
        }    
    Print_Relators(Relators,NumRelators,myout);    

    UDV[ReadPres] = DONE;
}

void Print_Relators(unsigned char ***MyRelators,int MyNumRelators,FILE *fptr)
{
    register int    i,
                    j;
        
    register unsigned char *p;
    
    for(i = 1; i <= MyNumRelators; i++)    
        {
        HLock((char **) MyRelators[i]);        
        fprintf(fptr,"\n R%3u:   ",i);
            p = *MyRelators[i];
            j = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(j++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    j = 9;
                    }
                }
        HUnlock((char **) MyRelators[i]);                                                        
        }        
    fprintf(fptr,"\n");
}

void Micro_Print_Reset(void)
{
    register int     i,
                    j;
        
    register unsigned char *p;
    
    printf("\n\nStarted with Presentation %d of Summand %d, Length %lu:\n",
        ReadPres + 1,CurrentComp,SURL[ReadPres]);
    if(Micro_Print_F)    
        fprintf(myout,"\n\nStarted with Presentation %d of Summand %d, Length %lu:\n",
            ReadPres + 1,CurrentComp,SURL[ReadPres]);
    
    fptr = stdout;
    for(i = 1; i <= NumRelators; i++)    
        {
        HLock((char **) SUR[ReadPres][i]);        
        fprintf(fptr,"\n R%3u:   ",i);
            p = *SUR[ReadPres][i];
            j = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(j++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    j = 9;
                    }
                }
        HUnlock((char **) SUR[ReadPres][i]);                                                        
        }        
    fprintf(fptr,"\n");
    
    if(Micro_Print_F)
        {
        fptr = myout;
        for(i = 1; i <= NumRelators; i++)    
            {
            HLock((char **) SUR[ReadPres][i]);        
            fprintf(fptr,"\n R%3u:   ",i);
                p = *SUR[ReadPres][i];
                j = 9;
                while(*p)
                    {
		    fputc(*p++,fptr);
                    if(j++ >= 90 && *p)
                        {
                        fprintf(fptr,"\n         ");
                        j = 9;
                        }
                    }
            HUnlock((char **) SUR[ReadPres][i]);                                                        
            }        
        fprintf(fptr,"\n");
        }
}

void Micro_Print_Freely_Reduce(unsigned long length, unsigned long origlength)
{
    printf("\n\nFree reductions reduced the length of the current presentation from %lu to %lu.",
        length,origlength);
    printf("\nThe reduced presentation is:\n");
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        {
        fprintf(myout,"\n\nFree reductions reduced the length of the current presentation from %lu to %lu.",
            length,origlength);
        fprintf(myout,"\nThe reduced presentation is:\n");
        Print_Relators(Relators,NumRelators,myout);
        }    
}

void Micro_Print_Dualize(void)
{
    printf("\n\nDualized the current relators to get the following dual relators:\n");
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        {
        fprintf(myout,"\n\nDualized the current relators to get the following dual relators:\n");
        Print_Relators(Relators,NumRelators,myout);
        }
}

void Micro_Print_Bandsum(void)
{
    int        SNumRelators;
    
    SNumRelators = NumRelators;    
    printf("\n\nReplaced Relator %u with the following bandsum of Relator %u and Relator %u.",
        Word1,Word1,Word2);
    printf(" Delta Length = %ld.\n",Length - SLength);    
    if(Micro_Print_F)
        {
        fprintf(myout,"\n\nReplaced Relator %u with the following bandsum of Relator %u and Relator %u.",
            Word1,Word1,Word2);
        fprintf(myout," Delta Length = %ld.\n",Length - SLength);
        }        
    NumRelators = 1;
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        Print_Relators(Relators,NumRelators,myout);
    NumRelators = SNumRelators;
    if(Word1 != 1)
        {
        printf("\n\nAnd then swapped Relator %u and Relator 1.",Word1);
        if(Micro_Print_F)
            fprintf(myout,"\n\nAnd then swapped Relator %u and Relator 1.",Word1);
        }
    printf("\n\nThe current presentation is:\n");
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        {
        fprintf(myout,"\n\nThe current presentation is:\n");
        Print_Relators(Relators,NumRelators,myout);
        }        
}    

void Micro_Print_Level_Transformations_Reset(void)
{
    printf("\n\nLooking for level-transformations of:\n");
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        {
        fprintf(myout,"\n\nLooking for level-transformations of:\n");
        Print_Relators(Relators,NumRelators,myout);
        }
}

void Micro_Print_Level_Transformations(unsigned int TheComp,
				       unsigned int V1,
				       unsigned int V2,
				       unsigned int Type,
				       unsigned int NumReps)
{
    char            x,
                    y;
                                        
    unsigned char     *p;
    
    int             i;
    
    unsigned int    j;
    
    if(V1 & 1)
        x = V1/2 + 97;
    else
        x = V1/2 + 65;
    if(V2 & 1)
        y = V2/2 + 97;
    else
        y = V2/2 + 65;                
    printf("\n\nVertices %c and %c form a Type %u separating pair.",x,y,Type);
    if(Micro_Print_F)
        fprintf(myout,"\n\nVertices %c and %c form a Type %u separating pair.",x,y,Type);
    
    ReallocateHandle((char **) Temp9,2*VERTICES);
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
    HLock((char **) Temp8);
    HLock((char **) Temp9);    
    printf("\nPerformed a level-transformation by sliding vertice(s):");
    printf("\n{%s}",*Temp9);
    printf("\nalong a path represented by: ");
    for(j = 0; j < NumReps; j++) printf("%s",*Temp8);
    printf("\nto obtain the presentation:\n");
    if(Micro_Print_F)
        {
        fprintf(myout,"\nPerformed a level-transformation by sliding vertice(s):");
        fprintf(myout,"\n{%s}",*Temp9);
        fprintf(myout,"\nalong a path represented by: ");
        for(j = 0; j < NumReps; j++) fprintf(myout,"%s",*Temp8);
        fprintf(myout,"\nto obtain the presentation:\n");
        }
    HUnlock((char **) Temp8);
    HUnlock((char **) Temp9);
    Print_Relators(Relators,NumRelators,stdout);
    if(Micro_Print_F)
        Print_Relators(Relators,NumRelators,myout);
}

void Micro_Print_Do_Aut(unsigned int Source, unsigned int NumReps)
{
    unsigned char    A,
                    a,
                    x;
                    
    int                i,
                    j,
                    k;
                                    
    A = ((Source >> 1) + 65);
    a = A + 32;
    
    if(Micro_Print)
        {
        fprintf(stdout,"\nDo Aut %u time(s): ",NumReps);
        if(Micro_Print_F)
            fprintf(myout,"\nDo Aut %u time(s): ",NumReps);
        }
    else
        {
        fprintf(stdout,"\n%6lu) ",Num_Level_Transformations + 1);
        fprintf(myout,"\n%6lu) ",Num_Level_Transformations + 1);
        }    
        
    for(i = j = k = 0; i < Vertices; i+= 2)
        {
        if(VA[i >> 1] == 0) continue;
        if(!ZZ[i] && !ZZ[i+1])
            j++;
        else
        if(ZZ[i] && ZZ[i+1])
            k++;
        }
    if(j <= k) for(i = j = 0; i < Vertices; i += 2)
        {
        if(i == Source) continue;
        if(VA[i >> 1] == 0) continue;
        x = (i >> 1) + 65;
        if(!ZZ[i])
            {
            if(ZZ[i+1])
                {
                if(++j > 12)
                    {
                    j = 1;
                    fprintf(stdout,"\n        ");
                    if(Compute_Stabilizers || Micro_Print_F)                    
                        fprintf(myout,"\n        ");
                    }
                fprintf(stdout,"%c->%c%c  ",x,A,x);
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"%c->%c%c  ",x,A,x);
                }
            else
                {
                if(++j > 12)
                    {
                    j = 1;
                    fprintf(stdout,"\n        ");
                    if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"\n        ");
                    }
                fprintf(stdout,"%c->%c%c%c ",x,A,x,a);
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"%c->%c%c%c ",x,A,x,a);
                }
            }
        else
        if(!ZZ[i+1])
            {
            if(++j > 12)
                {
                j = 1;
                fprintf(stdout,"\n        ");
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"\n        ");
                }
            fprintf(stdout,"%c->%c%c  ",x,x,a);
            if(Compute_Stabilizers || Micro_Print_F)
                fprintf(myout,"%c->%c%c  ",x,x,a);
            }
        }
    else  for(i = j = 0; i < Vertices; i += 2)
        {
        if(i == Source) continue;
        if(VA[i >> 1] == 0) continue;
        x = (i >> 1) + 65;
        if(!ZZ[i])
            {
            if(ZZ[i+1])
                {
                if(++j > 12)
                    {
                    j = 1;
                    fprintf(stdout,"\n        ");
                    if(Compute_Stabilizers || Micro_Print_F)
                        fprintf(myout,"\n        ");
                    }
                fprintf(stdout,"%c->%c%c  ",x,x,A);
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"%c->%c%c  ",x,x,A);
                }
            }
        else
            {
            if(!ZZ[i+1])
                {
                if(++j > 12)
                    {
                    j = 1;
                    fprintf(stdout,"\n        ");
                    if(Compute_Stabilizers || Micro_Print_F)
                        fprintf(myout,"\n        ");
                    }
                fprintf(stdout,"%c->%c%c  ",x,a,x);
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"%c->%c%c  ",x,a,x);
                }
            else
                {
                if(++j > 12)
                    {
                    j = 1;
                    fprintf(stdout,"\n        ");
                    if(Compute_Stabilizers || Micro_Print_F)
                        fprintf(myout,"\n        ");
                    }
                fprintf(stdout,"%c->%c%c%c ",x,a,x,A);
                if(Compute_Stabilizers || Micro_Print_F)
                    fprintf(myout,"%c->%c%c%c ",x,a,x,A);
                }    
            }    
        }                
}

        
#ifdef PRINT
void Print_DelRelators(void)
{
    int i;
        
    fprintf(stdout,"\n");    
    for(i = 1; i <= NumRelators; i++)
        {
        HLock((char **) DelRelators[i]);
        fprintf(stdout,"\nR %2u:    %s",i,*DelRelators[i]);
        HUnlock((char **) DelRelators[i]);
        }
}

void Print_DualRelators(int F1)
{
    register int             i,
                            j;
        
    register unsigned char *p;
    
    FILE                    *fptr;
    
    fptr = stdout;
    
    RERUN_AND_SAVE:
    
    fprintf(fptr,"\n\nThe 'Dual' Relators of Diagram %d are:\n",WhichInput + 1);
    
    for(i = 1; i <= NumGenerators; i++)    
        {
        HLock((char **) DualRelators[i]);        
        fprintf(fptr,"\n R%3u:   ",i);
            p = *DualRelators[i];
            j = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(j++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    j = 9;
                    }
                }
        HUnlock((char **) DualRelators[i]);                                                        
        }        
    fprintf(fptr,"\n");
    
    if(fptr == stdout)
        {
        printf("\n    SAVE A COPY OF THIS DATA IN 'Heegaard_Results' ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:
        switch(WaitkbHit())
            {
            case 'y':
                fptr = myout;
                goto RERUN_AND_SAVE;
            case 'n':
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE1;    
            }
        }
        
    printf("\n\nHIT 'P' TO PRINT THE DIAGRAM.");
    printf("\n   HIT 'm' TO RETURN TO DIAGRAM %d.",WhichInput + 1);
    printf("\n      HIT 'v' TO REVIEW PRESENTATION %d.",WhichInput + 1);
    if(F1)
        {
        printf("\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
        printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
        printf("\n               HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");
        printf("\n                  HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
        printf("\n                     HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");                                            
        }
    else
        {    
        printf("\n         HIT 'n' TO SEE THE NEXT DIAGRAM.");
        printf("\n            HIT 'p' TO SEE THE PREVIOUS DIAGRAM.");
        printf("\n               HIT 'q' TO QUIT VIEWING DIAGRAMS.");
        printf("\n                  HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");
        printf("\n                     HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
        printf("\n                        HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");        
        }                
}

void Print_OutRelators(F1)
{
    register int             i,
                            j;
        
    register unsigned char *p;
    
    FILE                    *fptr;
    
    fptr = stdout;
    
    RERUN_AND_SAVE:
    
    fprintf(fptr,"\n\nThe 'Out' Relators of Diagram %d are:\n",WhichInput + 1);
    
    for(i = 1; i <= NumRelators; i++)    
        {
        HLock((char **) OutRelators[i]);        
        fprintf(fptr,"\n R%3u:   ",i);
            p = *OutRelators[i];
            j = 9;
            while(*p)
                {
		fputc(*p++,fptr);
                if(j++ >= 90 && *p)
                    {
                    fprintf(fptr,"\n         ");
                    j = 9;
                    }
                }
        HUnlock((char **) OutRelators[i]);                                                        
        }        
    fprintf(fptr,"\n");
    
    if(fptr == stdout)
        {
        printf("\n    SAVE A COPY OF THIS DATA IN 'Heegaard_Results' ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:
        switch(WaitkbHit())
            {
            case 'y':
                fptr = myout;
                goto RERUN_AND_SAVE;
            case 'n':
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE1;    
            }
        }
        
    printf("\n\nHIT 'P' TO PRINT THE DIAGRAM.");
    printf("\n   HIT 'm' TO RETURN TO DIAGRAM %d.",WhichInput + 1);
    printf("\n      HIT 'v' TO REVIEW PRESENTATION %d.",WhichInput + 1);
    if(F1)
        {
        printf("\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
        printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
        printf("\n               HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");                                    
        printf("\n                  HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
        printf("\n                     HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");        
        }
    else
        {    
        printf("\n         HIT 'n' TO SEE THE NEXT DIAGRAM.");
        printf("\n            HIT 'p' TO SEE THE PREVIOUS DIAGRAM.");
        printf("\n               HIT 'q' TO QUIT VIEWING DIAGRAMS.");
        printf("\n                  HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");
        printf("\n                     HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
        printf("\n                        HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");        
        }        
}

void Print_SLR(int i)
{
    int j;
    
    fprintf(stdout,"\n");    
    for(j = 1; j <= NumRelators; j++)
        {
        HLock((char **) SLR[i][j]);
        fprintf(stdout,"\nSLR [%2u][%2u]:    %s",i,j,*SLR[i][j]);
        HUnlock((char **) SLR[i][j]);
        }
}
#endif        

Display_A_Diagram()
{    
    int                Reply;
    
    unsigned int    SaveUDV;            
    
    unsigned long    SLSP,
                    SLSQ;

    DrawingDiagrams = TRUE;
    CycleDiagrams     = FALSE;

    if(Get_Relators_From_SUR(WhichInput))
        {
        printf("\n\n    Memory Error. Sorry!");
        goto _ERROR;        
        }
                
    if(Length == 0L)
        {
        printf("\n\nThis is an empty presentation. There is nothing to display.");
        goto _ERROR;
        }

    if(WhichInput == 0 && UDV[0] == SPLIT)
        {
        switch(Find_Flow_A(NORMAL,FALSE))
            {
            case 1:
                break;
            case TOO_LONG:
                SysBeep(5);
                printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                goto _ERROR;
            }
        if(Automorphisms)
            {
            SysBeep(5);
            printf("\n\n                    NOTE!");
            printf("\n\n    Presentation 1 does not have minimal length.");
            printf("\n    The program will only display diagrams of minimal length presentations.");
            printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
            goto _ERROR;
            }    
        }
    
    Fill_A(NumRelators);
    if(ComputeValences_A())
        {
        SysBeep(5);
        printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
        goto _ERROR;    
        }    
    Get_Matrix();
    Check_Connected();
    SepPairs = Sep_Pairs(0,0);
    if(SepPairs == TOO_LONG)
        {
        SysBeep(5);
        printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
        goto _ERROR;
        }
    SaveUDV = UDV[WhichInput];
    SLSP = LSP[WhichInput];
    SLSQ = LSQ[WhichInput];    
    if(SepPairs)
        {
        if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
        if(V1 & 1)
            LSP[WhichInput] = V1/2 + 97;
        else
            LSP[WhichInput] = V1/2 + 65;
        if(V2 & 1)
            LSQ[WhichInput] = V2/2 + 97;
        else
            LSQ[WhichInput] = V2/2 + 65;        
        }
    NonPlanar = Planar(FALSE,FALSE);
    Reply = Print_Graph(TRUE);
    if(UDV[WhichInput] != ANNULUS_EXISTS && UDV[WhichInput] != V2_ANNULUS_EXISTS)
        UDV[WhichInput] = SaveUDV;
    LSP[WhichInput] = SLSP;
    LSQ[WhichInput] = SLSQ;    
    DrawingDiagrams = FALSE;            
    printf("\f");
    return(Reply);

_ERROR:
    printf("\n\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
    printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
    GET_RESPONSE:
    switch(WaitkbHit())
        {
        case 'n':
            Reply = 0;
            break;
        case 'p':
            Reply = 1;
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE;
        }
    DrawingDiagrams = FALSE;            
    printf("\f");
    return(Reply);                        
}

void Display_Diagrams(void)
{
    unsigned char   DispList[MAX_SAVED_PRES],
                    *r;                    
                            
    int             NoSepPairs,
                    NumConnected,
                    PP,
                    Response,
                    SWhichInput;
    
    unsigned int    h,
                    i,
                    j,
                    k,
                    SaveUDV;
                    
    unsigned long    SLSP,
                    SLSQ;                                

    printf("\n\n                    Displaying Diagrams. . .");
REDRAW:
    if(NumFilled > 1)
        {
        printf("\n\n    REVIEW ALL PRESENTATIONS AVAILABLE ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:
        switch(WaitkbHit())
            {
            case 'y':
                REVIEW:
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,0);
                printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
                GET_RESPONSE5:
                switch(WaitkbHit())
                    {
                    case 'y':
                        goto REVIEW;
                    case 'n':
                        break;
                    default:
                        SysBeep(5);
                        goto GET_RESPONSE5;
                    }
                break;    
            case 'n':
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE1;
            }
        }        
    DrawingDiagrams = TRUE;
    
    printf("\n\n    a) SHOW ALL DIAGRAMS,");
    printf("\n    b) SHOW ALL DIAGRAMS THAT ARE CONNECTED AND HAVE NO PAIRS OF SEPARATING VERTICES,");
    printf("\n    c) OR SHOW YOUR CHOICE OF A PARTICULAR DIAGRAM ?");
    printf("\n\n    HIT 'a','b', OR 'c'");
    GET_RESPONSE2:
    switch(Response = WaitkbHit())
        {
        case 'a':
        case 'b':
        case 'c':
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE2;
        }
    if(Response == 'c')
        {
        CycleDiagrams = FALSE;
        r = (unsigned char *) NewPtr(100L);
        printf("\n\nENTER A DIAGRAM FROM 1 TO %u THAT YOU WANT TO SEE AND HIT 'return'.      ",NumFilled);
        for(i = j = 0; j < NumFilled; j++) if(SURL[j] == 0L) i ++;
        if(i)
            {
            if(i == 1)
                {
                for(i = 0; i < NumFilled && SURL[i]; i++) ;
                i++;
                printf("\n\nExcept for diagram %d which is a diagram of S1 X S2 or S1 X D2 and is empty.     ",i);
                }
            else
                {
                j = 0;
                j += printf("\n\nExcept for diagrams: ");
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
                printf("\nThese are diagrams of S1 X S2 (s) or S1 X D2 (s) and are empty.      ");
                }
            }
        GET_RESPONSE3:
        WhichInput = 0;                
        ReadString((char *)r, GetPtrSize(r));
        sscanf((char *) r,"%d",&WhichInput);        
        if(WhichInput < 1 || WhichInput > NumFilled || SURL[WhichInput-1] == 0L)
            {
            SysBeep(5);
            goto GET_RESPONSE3;
            }    
        DisposePtr((char *) r);    
        WhichInput --;
                
        if(Get_Relators_From_SUR(WhichInput))
            {
            printf("\n\n    Memory Error. Sorry!");
            goto REDRAW;        
            }
                
        if(Length == 0L)
            {
            printf("\n\nThis is an empty presentation. There is nothing to display.");
            goto REDRAW;
            }
        
        if(WhichInput == 0 && UDV[0] == SPLIT)
            {
            switch(Find_Flow_A(NORMAL,FALSE))
                {
                case 1:
                    break;
                case TOO_LONG:
                    SysBeep(5);
                    printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                    goto REDRAW;
                }
            if(Automorphisms)
                {
                SysBeep(5);
                printf("\n\n                    NOTE!");
                printf("\n\n    Presentation 1 does not have minimal length.");
                printf("\n    The program will only display diagrams of minimal length presentations.");
                printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
                printf("\n\n    HIT ANY KEY TO CONTINUE.");
                WaitkbHit();
                goto REDRAW;
                }    
            }
        
        Fill_A(NumRelators);
        if(ComputeValences_A())
            {
            SysBeep(5);
            printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
            goto REDRAW;
            }    
        Get_Matrix();
        Check_Connected();
        SepPairs = Sep_Pairs(0,0);
        if(SepPairs == TOO_LONG)
            {
            SysBeep(5);
            printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
            goto REDRAW;
            }
        SaveUDV = UDV[WhichInput];
        SLSP = LSP[WhichInput];
        SLSQ = LSQ[WhichInput];
        SWhichInput = WhichInput;        
        if(SepPairs)
            {
            if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
            if(V1 & 1)
                LSP[WhichInput] = V1/2 + 97;
            else
                LSP[WhichInput] = V1/2 + 65;
            if(V2 & 1)
                LSQ[WhichInput] = V2/2 + 97;
            else
                LSQ[WhichInput] = V2/2 + 65;
            }
        NonPlanar = Planar(FALSE,FALSE);
        Print_Graph(FALSE);
        if(UDV[SWhichInput] != ANNULUS_EXISTS && UDV[SWhichInput] != V2_ANNULUS_EXISTS)
            UDV[SWhichInput] = SaveUDV;
        LSP[SWhichInput] = SLSP;
        LSQ[SWhichInput] = SLSQ;        
        }
    if(Response == 'a' || Response == 'b')
        {
        CycleDiagrams = TRUE;
        NoSepPairs = FALSE;
        NumConnected = 0;
        for(i = 0; i < NumFilled; i++) DispList[i] = EOS;
        for(WhichInput = 0; WhichInput < NumFilled; WhichInput ++)
        if(SURL[WhichInput] != 0L)
            {
            h = UDV[WhichInput];
            if(Response == 'b' && (h == SEP_PAIRS || h == ANNULUS_EXISTS)) continue;        

            if(Get_Relators_From_SUR(WhichInput))
                {
                printf("\n\n    Memory Error. Sorry!");
                goto REDRAW;        
                }
            
            if(WhichInput == 0 && h == SPLIT)
                {
                switch(Find_Flow_A(NORMAL,FALSE))
                    {
                    case 1:
                        break;
                    case TOO_LONG:
                        SysBeep(5);
                        printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                        continue;
                    }
                if(Automorphisms)
                    {
                    SysBeep(5);
                    printf("\n\n                    NOTE!");
                    printf("\n\n    Presentation 1 does not have minimal length.");
                    printf("\n    The program will only display diagrams of minimal length presentations.");
                    printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
                    printf("\n\nHIT ANY KEY TO SEE THE NEXT DIAGRAM.");
                    WaitkbHit();
                    continue;
                    }    
                }
            
            Fill_A(NumRelators);
            if(ComputeValences_A())
                {
                SysBeep(5);
                printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                continue;
                }            
            Get_Matrix();
            Check_Connected();
            if(!Connected && Response == 'b') continue;
            NumConnected ++;
            SepPairs = Sep_Pairs(0,0);
            if(SepPairs == TOO_LONG)
                {
                SysBeep(5);
                printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                continue;
                }    
            SaveUDV = UDV[WhichInput];
            SLSP = LSP[WhichInput];
            SLSQ = LSQ[WhichInput];
            if(SepPairs)
                {
                if(Response == 'b') continue;                
                if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
                if(V1 & 1)
                    LSP[WhichInput] = V1/2 + 97;
                else
                    LSP[WhichInput] = V1/2 + 65;
                if(V2 & 1)
                    LSQ[WhichInput] = V2/2 + 97;
                else
                    LSQ[WhichInput] = V2/2 + 65;        
                }
            NoSepPairs = TRUE;    
            NonPlanar = Planar(FALSE,FALSE);
            DispList[WhichInput] = TRUE;
            SWhichInput = WhichInput;
            if(Print_Graph(FALSE))
                {
                if(Response == 'a')
                    {
                    if(WhichInput) WhichInput -= 2;
                    else WhichInput --;    
                    }
                if(Response == 'b')
                    {
                    for(PP = WhichInput - 1; PP >= 0; PP--)
                    if(DispList[PP])
                        {
                        WhichInput = PP - 1;
                        break;
                        }
                    if(PP < 0) WhichInput --;                            
                    }
                }
            if(UDV[SWhichInput] != ANNULUS_EXISTS && UDV[SWhichInput] != V2_ANNULUS_EXISTS)
                UDV[SWhichInput] = SaveUDV;
            LSP[SWhichInput] = SLSP;
            LSQ[SWhichInput] = SLSQ;                    
            }
        }

    if(Response == 'b' && (!NoSepPairs || !NumConnected))
        {
        SysBeep(5);
        printf("\n\n    There are no diagrams which are connected and without pairs of separating vertices!");
        }
    if(Response == 'c' || WhichInput <= NumFilled)
        {
        printf("\n\n    CONTINUE DISPLAYING DIAGRAMS  ?  HIT 'y' OR 'n'.");    
        GET_RESPONSE4:
        switch(WaitkbHit())
            {
            case 'y':
                goto REDRAW;
            case 'n':
                break;
            default:
                SysBeep(5);
                goto GET_RESPONSE4;
            }
        }    
    DrawingDiagrams = FALSE;            
    printf("\f");                    
}

int Save_Pres_To_Input_Presentations(WhichPres)
{
    register unsigned char    *p,
                            *q,
                            *r;
    
    unsigned char            *NewPresName;
    
    int                     i,
                            NumRelators;                                    
    
    if((input_relators = fopen("Input_Presentations","r+")) == NULL)
        {
        SysBeep(5);
        printf("\nUnable to open the file 'Input_Presentations'.\n");
        return(1);
        }
    
    if((r = (unsigned char *) NewPtr((Size)(MAXLENGTH + 1))) == NULL)
        {
        SysBeep(5);
        printf("\nMemory error. Sorry!\n");
        fclose(input_relators);        
        return(1);        
        }
    if((NewPresName = (unsigned char *) NewPtr(1000L)) == NULL)
        {
        SysBeep(5);
        printf("\nMemory error. Sorry!\n");
        DisposePtr((char *) r);
        fclose(input_relators);        
        return(1);        
        }

    GET_ID:
    printf("\n\nPlease enter a name by which the program can refer to this presentation,");
    printf("\nand then hit 'return'. Or just hit 'return' to skip saving this presentation.");
    printf("\n\nSAVE THIS PRESENTATION AS: ");        
    ReadString((char *)NewPresName, GetPtrSize(NewPresName));
    
    if(*NewPresName == EOS || *NewPresName == '\n')
        {
        DisposePtr((char *) r);
        DisposePtr((char *) NewPresName);
        fclose(input_relators);
        printf("\n\n      Hit 'return' to stop reviewing this report.");
        printf("\n      Hit the 'space-bar' to alternately stop and start scrolling.");
        printf("\n      Use the 'd','i','n','p' and 's' keys to: view diagrams, save a presentation in");
        printf("\n      Input_Presentations, see the next presentation, see the previous presentation,");
        printf("\n      or save a presentation in Heegaard_Results.");        
        return(1);
        }
        
    rewind(input_relators);
    do
        {
        if(fgets((char *) r,MAXLENGTH,input_relators) == NULL) goto ID_IS_UNIQUE;
        p = r;
        q = NewPresName;
        while(*p && *q && *p == *q)
            {
            p++;
            q++;
            }
        if(*p == '\n' || *p == ' ' || *p == '\t') *p = EOS;        
        }
    while(*p != *q);    

    SysBeep(5);
    printf("\nInput_Presentations already contains a presentation with this identifier!");
    goto GET_ID;
        
    ID_IS_UNIQUE:
    DisposePtr((char *) r);
    fseek(input_relators,0L,2);        
    fprintf(input_relators,"\n\n%s",NewPresName);
    NumRelators = NR[WhichPres];
    for(i = 1; i <= NumRelators; i++)    
        {
        HLock((char **) SUR[WhichPres][i]);        
        fprintf(input_relators,"\n    %s",*SUR[WhichPres][i]);
        HUnlock((char **) SUR[WhichPres][i]);                                                        
        }        
    fprintf(input_relators,"\n");

    fflush(input_relators);
    fclose(input_relators);
    printf("\n    Saved Presentation %d in 'Input_Presentations' as: %s",
        WhichPres + 1,NewPresName);
    DisposePtr((char *) NewPresName);
    return(0);
}
