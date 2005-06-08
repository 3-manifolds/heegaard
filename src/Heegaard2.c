#include "Heegaard.h"
#include "Heegaard_Dec.h"

#define MAX_COUNT    5                    
#define MTLK        100

Get_Diagrams()
{
    char                    c;
    
    register unsigned char  *p,
                            *q,
                            **Temp;
                            
    register unsigned int   h,
                            i,
                            j;
    
    int                     DistinctNonEmpty,
                            hh,
                            jj,
                            kk,
                            SMicro_Print,
                            SMicro_Print_F,
                            SRNumGenerators,
                            SRNumRelators,
                            SRReadPres;
                             
    unsigned int            DoingDup,
                            Flag2,
                            ii,
                            MaxUDV,
                            MinNR,
                            NumBandSum,
                            SNumFilled,
                            TheComp;
                            
    unsigned long           HS,
                            MinLength,
                            SRLength,
                            SRSLength;                            

    unsigned int Whitehead_Graph();
    unsigned int Reduce_Genus();

/**********************************************************************************************
                            THIS IS THE BEGINNING OF THE MAIN PROGRAM.
**********************************************************************************************/

    OnStack = 0;                                                                                    
    if(Input == RESET) goto _RESET;
    if(Input == REDUCE_GENUS)
        {
        SReadPres = ReadPres;
        Input = NORMAL;
        SRError = TRUE;
        goto _REDUCE_GENUS;
        }    
_DUALIZE:                                                                                                                    
    if(Input == DUALIZE)            
        {
        if(NumGenerators != NumRelators)
            {
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }
            
        /**************************************************************************************
            If the program is here, then the manifold is closed, and we are going to swap
            relators and dual relators to see if reductions in length are possible.
        **************************************************************************************/
        
        for(i = 1; i <= NumRelators; i++)
            { 
            Temp             = Relators[i];
            Relators[i]     = DualRelators[i];
            DualRelators[i] = Temp;    
            }
        NumDualized ++;        
        SLength = Length;
        if(Micro_Print) Micro_Print_Dualize();
        i = Freely_Reduce();
        if(i == TOO_LONG)
            {
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }
        if(i)
            {
            Length = OrigLength;
            Empty_Relator_D();
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }    
        Length = OrigLength;
          }            
_BANDSUM:        
    if(Input == BANDSUM)
        {    
        /**************************************************************************************
            Call New_Relator() to produce a new relator which is a bandsum of two previous
            relators.
        **************************************************************************************/
        
        SLength = Length;
        if(OnlyReducingBandsums) GoingUp = FALSE;
        if(FormBandsums)
            New_Relator(TRUE);
        else
            {
            Minimum = BIG_NUMBER;
            SaveMinima = TRUE;
            }
            
        /**************************************************************************************
            If FormBandsums is false, the program is put into a mode in which it will attempt
            to simplify a presentation only by dualizing (when that is possible) and by
            looking for primitives, proper powers and lens spaces.
        **************************************************************************************/    
        
        if(Minimum < BIG_NUMBER)
            {
            if(OnlyReducingBandsums && Minimum >= 0L)
                {
                if(Micro_Print)
                    {
                    printf("\n\nThere is no reducing bandsum for the current presentation.");
                    if(Micro_Print_F)
                        fprintf(myout,"\n\nThere is no reducing bandsum for the current presentation.");
                    }
                if(Boundary)
                    {
                    Length += Minimum;
                    goto QUASI_PSEUDO_MINIMAL;
                    }
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;    
                }
            SRError = TRUE;    
            LR[0] = GetHandleSize((char **) OutRelators[Word1]) - 1;
            i = Word1;
            HLock((char **) OutRelators[i]);
            Word1 = abs(Compare(*OutRelators[i]));
            HUnlock((char **) OutRelators[i]);
            if(Word1 == 0)
                {
                SysBeep(5);
                printf("\n\nError in forming a bandsum. Relators are probably too long.");
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }
            if(Micro_Print)
                {
                LR[0] = GetHandleSize((char **) OutRelators[Word2]) - 1;
                i = Word2;
                HLock((char **) OutRelators[i]);
                Word2 = abs(Compare(*OutRelators[i]));
                HUnlock((char **) OutRelators[i]);
                if(Word2 == 0)
                    {
                    SysBeep(5);
                    printf("\n\nError in forming a bandsum. Relators are probably too long.");
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }
                }        
            Temp             = Temp3;
            Temp3             = Relators[Word1];    
            Relators[Word1] = Relators[1];
            Relators[1]     = Temp2;
            Temp2             = Temp;
            LR[Word1]        = LR[1];
            LR[1]            = GetHandleSize((char **) Relators[1]) - 1;
            OrigLength         = Length - GetHandleSize((char **) Temp3) + GetHandleSize((char **) Relators[1]);
            Length            = OrigLength;
            if(Micro_Print) Micro_Print_Bandsum();
            }
        else
            {
            Input = RESET;
            if(OnlyReducingBandsums) SaveMinima = TRUE;
            goto _RESET;
            }
        
        From_BANDSUM ++;
        NumBandSum ++;    
        
        if(EmtyRel)
            {
            EmtyRel = FALSE;
            Empty_Relator_BS();
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }
        }
        
    /******************************************************************************************
        If there is more than one relator, we call Reduce_Genus() which looks for primitives,
        proper powers of free generators and lens spaces. If Reduce_Genus() finds such relators
        it deletes other relators which are consequences and, in the case of primitives,
        it deletes the primitive.
    ******************************************************************************************/    
                
    if((Input != DUALIZE && NumRelators > 1 && UDV[ReadPres] > 2)
        || (Input == DUALIZE && NumGenerators == 2))
        {
        if(Input == BANDSUM)
            {
            /**********************************************************************************
                                Save a copy of the current relators.
            **********************************************************************************/
            
            SRLength = Length;
            SRSLength = SLength;
            SRNumRelators = NumRelators;
            SRNumGenerators = NumGenerators;
            SRError = FALSE;
            SRReadPres = ReadPres;
            for(i = 1; i <= NumRelators; i++)
                {
                ReallocateHandle((char **) Copy_Of_Rel_2[i],GetHandleSize((char **) Relators[i]));
                if((p = *Copy_Of_Rel_2[i]) == NULL)
                    {
                    SRError = TRUE;
                    break;
                    }
                q = *Relators[i];
                while(*p++ = *q++) ;                    
                }    
            }
        else
            SRError = TRUE;    
        switch(Reduce_Genus(Input,FALSE,FALSE))
            {
            case NO_ERROR:
                break;
            case FATAL_ERROR:
                Fatal_Error();
                return(1);
            case TOO_LONG:
            case CAN_NOT_DELETE:
                if(Micro_Print)
                    {
                    printf("\n\nUnable to delete a primitive.");
                    if(Micro_Print_F)
                        fprintf(myout,"\n\nUnable to delete a primitive.");
                    }
                if(TP[ReadPres]) TP[ReadPres] --;    
                FoundPower = FoundPrimitive = LensSpace = EmtyRel = FALSE;
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;    
            }
_REDUCE_GENUS:
        if(LensSpace)
            {
            /********************************************************************************** 
                If LensSpace = NOT_CONNECTED, the program found that a genus two diagram of a
            lens space is not connected and split it into two genus one diagrams. Otherwise,
            LensSpace is true only when the routine Lens_Space() was passed a presentation of
            a lens space and succesfully discovered which lens space it was.
            **********************************************************************************/
            
            FoundPower = FoundPrimitive = EmtyRel = FALSE;
            if(LensSpace == NOT_CONNECTED)                             
                {
                LensSpace = FALSE;
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }
            LensSpace = FALSE;
            Length = GetHandleSize((char **) Relators[1]) + GetHandleSize((char **) Relators[2]) - 2;
            if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                {
                CBC[ComponentNum[ReadPres]][0] = 1;
                CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                }
            h = On_File();
            if(h == NumFilled)    
                {
                if(Dup_On_File < INFINITE)
                    {
                     if(Save_Pres(ReadPres,Dup_On_File,Length,1,3,0,0,0)) return(1);
                    UDV[NumFilled - 1] = DUPLICATE;                    
                    if(UDV[Dup_On_File] <= UNKNOWN)
                        {
                        if(UDV[Dup_On_File] == ANNULUS_EXISTS
                            || UDV[Dup_On_File] == V2_ANNULUS_EXISTS)
                            DisposeHandle((char **) SUR[Dup_On_File][0]);   /* Can this happen? */
                        BDY[Dup_On_File] = 0;        
                        UDV[Dup_On_File] = KNOWN_LENS_SPACE;
                        LSP[Dup_On_File] = P;
                        LSQ[Dup_On_File] = Q;
                        Mark_As_Found_Elsewhere(ComponentNum[Dup_On_File]);
                        }    
                     }
                 else
                     {
                    if(Save_Pres(ReadPres,0,Length,1,3,0,0,0)) return(1);                                
                    BDY[NumFilled - 1] = 0;        
                    UDV[NumFilled - 1] = KNOWN_LENS_SPACE;
                    LSP[NumFilled - 1] = P;
                    LSQ[NumFilled - 1] = Q;                    
                    Mark_As_Found_Elsewhere(CurrentComp);
                    }                                
                }            
            else if(UDV[h] <= UNKNOWN)
                {
                if(UDV[h] == ANNULUS_EXISTS
                    || UDV[h] == V2_ANNULUS_EXISTS)
                        DisposeHandle((char **) SUR[h][0]);   /* Can this happen? */                    
                BDY[h]  = 0;    
                UDV[h]  = KNOWN_LENS_SPACE;
                LSP[h]  = P;
                LSQ[h]  = Q;
                Mark_As_Found_Elsewhere(CurrentComp);
                }            
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;        
            }
        if(EmtyRel)
            {
            FoundPower = FoundPrimitive = EmtyRel = FALSE;                
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }                    
        if(FoundPower || FoundPrimitive)
            Input = NORMAL;
        }
                    
    /******************************************************************************************
        At this point, we call Find_Flow_A(), which reduces the length of the relators in
        Relators[] via automorphisms if possible. It may happen that a generator which
        oppeared in the relators no longer appears after some automorphisms have been
        performed. In this case, we have a "missing" generator. When this happens,
        Find_Flow_A() returns TRUE, and we call Missing_Gen().
    ******************************************************************************************/    

_FIND_FLOW:    
    if(i = Find_Flow_A(Input,FALSE))
        {
        FoundPrimitive = FoundPower = FALSE;        
        if(i == 1 && Missing_Gen() == TOO_MANY_COMPONENTS) return(1);
        SaveMinima = TRUE;
        Input = RESET;
        goto _RESET;
        }
    
    if(SRError == 2) Input = BANDSUM;
        
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
        
    /******************************************************************************************
        After Find_Flow_A() has reduced the presentation to minimal length, we check whether
        any  "termination" condition has occured. If so, we want to stop processing 
        some presentations. We also check whether the current presentation is shorter or
        simpler than previous presentations, in which case, we want to save a copy.
    ******************************************************************************************/    
                
    if(NumGenerators == 1)
        {
        /**************************************************************************************
            If the present presentation has only one generator, then it corresponds to a 
            "lens-space" of some sort. Try to determine which "lens-space", and save a copy
            of the presentation.
        **************************************************************************************/    
        
        FoundPrimitive = FoundPower = FALSE;
        TP[ReadPres] = FALSE;
        if((DistinctNonEmpty = Delete_Dups()) > 1)
            {
            Fatal_Error();        
            return(1);
            }
        HS = GetHandleSize((char **) Relators[1]) - 1;        
        if(On_File() == NumFilled)
            {    
            switch(HS)
                {
                case 0:
                    if(BDY[ReadPres] == TRUE)
                        {
                        if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                            {
                            CBC[ComponentNum[ReadPres]][0] = EOS;
                            CBC[ComponentNum[ReadPres]][1] = 1;
                            CBC[ComponentNum[ReadPres]][2] = BDRY_UNKNOWN;
                            }
                        UDV[NumFilled] = S1_X_D2;
                        BDY[NumFilled] = TRUE;
                        }
                    if(BDY[ReadPres] == FALSE)
                        {
                        if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                            {
                            CBC[ComponentNum[ReadPres]][0] = NumRelators - NumEmptyRels;
                            CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                            }
                        UDV[NumFilled] = S1_X_S2;
                        BDY[NumFilled] = FALSE;
                        }
                    if(BDY[ReadPres] > 1)
                        {
                        UDV[NumFilled] = S1_X_X2;
                        BDY[NumFilled] = BDY[ReadPres];
                        }    
                    break;
                case 1:
                    if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                        {
                        CBC[ComponentNum[ReadPres]][0] = NumRelators - NumEmptyRels;
                        CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                        }
                    UDV[NumFilled] = THREE_SPHERE;
                    BDY[NumFilled] = FALSE;    
                    break;
                case 2:
                case 3:
                case 4:
                    if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                        {
                        CBC[ComponentNum[ReadPres]][0] = NumRelators - NumEmptyRels;
                        CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                        }
                    LSP[NumFilled] = HS;
                    LSQ[NumFilled] = 1;    
                    UDV[NumFilled] = KNOWN_LENS_SPACE;
                    BDY[NumFilled] = FALSE;
                    break;
                default:
                    if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                        {
                        CBC[ComponentNum[ReadPres]][0] = NumRelators - NumEmptyRels;
                        CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                        }
                    LSP[NumFilled] = HS;
                    LSQ[NumFilled] = 1;                                                        
                    UDV[NumFilled] = GENERIC_LENS_SPACE;
                    BDY[NumFilled] = FALSE;
                    if(Micro_Print)
                        {
                        printf("\n\nThe current presentation presents a manifold of Heegaard genus one.");
                        if(Micro_Print_F)                        
                            fprintf(myout,"\n\nThe current presentation presents a manifold of Heegaard genus one.");
                        }                    
                    if(Dup_On_File < INFINITE)
                        {
                        if(Save_Pres(ReadPres,Dup_On_File,Length,1,4,0,0,0)) return(1);
                        Mark_As_Duplicate(Dup_On_File);                        
                        Input = RESET;
                        goto _RESET;
                        }
                    if(Save_Pres(ReadPres,0,Length,1,4,0,0,0)) return(1);                                
                    Input = RESET;
                    goto _RESET;                
                }
            if(Micro_Print)
                {
                printf("\n\nThe current presentation presents a manifold of Heegaard genus one.");
                if(Micro_Print_F)                
                    fprintf(myout,"\n\nThe current presentation presents a manifold of Heegaard genus one.");
                }                
            if(Dup_On_File < INFINITE)
                {
                if(Save_Pres(ReadPres,Dup_On_File,Length,1,4,0,0,0)) return(1);
                Mark_As_Duplicate(Dup_On_File);                
                Input = RESET;
                goto _RESET;
                }
            if(Save_Pres(ReadPres,0,Length,1,4,0,0,0)) return(1);
            Mark_As_Found_Elsewhere(CurrentComp);                                        
            Input = RESET;
            goto _RESET;    
            }
        }
        
    /******************************************************************************************    
        If the total length of the relators is equal to the number of generators, then this
        may be a canonical presentation of the 3-Sphere. Check whether this is the case, and
        if so terminate further processing of the appropriate associated presentations.
    ******************************************************************************************/    
    
    if(Length == NumGenerators && NumRelators == NumGenerators && !Boundary)
        {
        FoundPrimitive = FoundPower = FALSE;        
        for(i = 1; i <= NumRelators; i++)
            if(GetHandleSize((char **) Relators[i]) != 2L) break;
        if(i > NumRelators && Delete_Dups() == NumRelators)
            {
            if(Micro_Print)
                {
                printf("\n\nThe current presentation presents the 3-Sphere.");
                if(Micro_Print_F)                
                    fprintf(myout,"\n\nThe current presentation presents the 3-Sphere.");
                }            
            if(On_File() == NumFilled)    
                {
                if(Dup_On_File < INFINITE)
                    {
                     if(Save_Pres(ReadPres,Dup_On_File,Length,1,5,0,0,0)) return(1);
                     if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                        {
                        CBC[CurrentComp][0] = 1;
                        CBC[CurrentComp][1] = BDRY_UNKNOWN;
                        }
                    Mark_As_Duplicate(Dup_On_File);                    
                    if(UDV[Dup_On_File] <= UNKNOWN)
                        {
                        UDV[Dup_On_File] = THREE_SPHERE;
                        Mark_As_Found_Elsewhere(ComponentNum[Dup_On_File]);
                        }
                     }
                 else
                     {
                    if(Save_Pres(ReadPres,0,Length,1,5,0,0,0)) return(1);                                                     
                    UDV[NumFilled - 1] = THREE_SPHERE;
                    BDY[NumFilled - 1] = FALSE;
                     if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                        {
                        CBC[CurrentComp][0] = 1;
                        CBC[CurrentComp][1] = BDRY_UNKNOWN;
                        }                    
                    }
                Mark_As_Found_Elsewhere(CurrentComp);                    
                }
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }
        SaveMinima = TRUE;    
        Input = RESET;
        goto _RESET;    
        }
        
    /******************************************************************************************
        If FoundPrimitive is TRUE, then the subroutine Reduce_Genus() has found and deleted a
        primitive relator. If FoundPower is TRUE, then Reduce_Genus() found a relator which is
        a proper power of a free generator. Check whether the resulting presentation is new
        and save it for further investigation if it is new.
    ******************************************************************************************/    
        
    if(FoundPrimitive == TRUE || FoundPower == TRUE)
        {            
        if(On_File() == NumFilled)
            {
            if(Dup_On_File < INFINITE)
                {                 
                 if(FoundPrimitive && Save_Pres(ReadPres,Dup_On_File,Length,1,6,1,0,0))
                     {
                     FoundPrimitive = FoundPower = FALSE;
                     return(1);
                     }
                 else if(FoundPower && Save_Pres(ReadPres,Dup_On_File,Length,1,60,1,0,0))
                     {
                     FoundPrimitive = FoundPower = FALSE;
                     return(1);
                     }
                 Mark_As_Duplicate(Dup_On_File);                                                                      
                 }
             else
                 {
                if(FoundPrimitive && Save_Pres(ReadPres,0,Length,1,6,1,0,0))
                    {
                    FoundPrimitive = FoundPower = FALSE;
                    return(1);
                    }
                else if(FoundPower && Save_Pres(ReadPres,0,Length,1,60,1,0,0))
                    {
                    FoundPrimitive = FoundPower = FALSE;
                    return(1);
                    }
                BDY[NumFilled - 1] = BDY[ReadPres];
                UDV[NumFilled - 1] = 0;                    
                }    
            }
        else
            {
            if(TP[ReadPres]) TP[ReadPres] --;
            if(SRError == FALSE)
                {
                FoundPrimitive = FoundPower = FALSE;
                Length = SRLength;
                SLength = SRSLength;
                NumRelators = SRNumRelators;
                NumGenerators = SRNumGenerators;
                Vertices = 2*NumGenerators;
                ReadPres = SRReadPres;
                for(i = 1; i <= NumRelators; i++)
                    { 
                    Temp                 = Relators[i];
                    Relators[i]         = Copy_Of_Rel_2[i];
                    Copy_Of_Rel_2[i]     = Temp;    
                    }
                SRError = 2;
                if(Micro_Print)
                    {
                    printf("\n\nReverting to the presentation formed by the last bandsum:\n");
                    Print_Relators(Relators,NumRelators,stdout);
                    if(Micro_Print_F)
                        {
                        fprintf(myout,"\nReverting to the presentation formed by the last bandsum:\n");
                        Print_Relators(Relators,NumRelators,myout);
                        }
                    }                                                
                goto _FIND_FLOW;
                }
            }                
        FoundPrimitive = FoundPower = FALSE;                
        Input = RESET;
        goto _RESET;
        }
        
    /******************************************************************************************
        Here we test for the analog of pseudo-minimality, in the case where the manifold is
        not closed. We want to save a copy of the previous presentation P1, which has minimal
        length equal to SLength, provided the current presentation P2, which was obtained from
        P1 by forming a bandsum and then reducing via automorphisms, has minimal length greater
        than or equal to SLength. (Note that the program automatically saves a copy of P1 in
        the array OutRelators[]. Note this copy is valid provided SRError != 3)
    ******************************************************************************************/    
        
    if(Boundary && Input == BANDSUM && !GoingUp && SLength <= Length && SRError != 3)
        {
        QUASI_PSEUDO_MINIMAL:
        Length = SLength;
        for(i = 1; i <= NumRelators; i++)
            {
            Temp             = Relators[i];
            Relators[i]     = OutRelators[i];
            OutRelators[i]     = Temp;
            }
        if(Micro_Print)
            {
            printf("\n\nRecalled the previous set of relators and made them current.\n");
            if(Micro_Print_F)            
                fprintf(myout,"\n\nRecalled the previous set of relators and made them current.\n ");
            }            
        j = On_File();        
        if(j == NumFilled)
            {
            if(Micro_Print)
                {
                printf("\n\nSaving the current presentation, which has length %lu and is quasi-pseudo-minimal:\n",Length);
                Print_Relators(Relators,NumRelators,stdout);
                if(Micro_Print_F)
                    {
                    fprintf(myout,"\n\nSaving the current presentation, which has length %lu and is quasi-pseudo-minimal:\n",Length);    
                    Print_Relators(Relators,NumRelators,myout);
                    }
                }                
            if(Dup_On_File < INFINITE)
                {
                 if(Save_Pres(ReadPres,Dup_On_File,Length,1,7,1,0,0)) return(1);
                 QPM[Dup_On_File] = TRUE;                                
                 Mark_As_Duplicate(Dup_On_File);                 
                 }
             else
                 {
                if(Save_Pres(ReadPres,0,Length,1,7,1,0,0)) return(1);                             
                BDY[NumFilled - 1] = TRUE;
                UDV[NumFilled - 1] = 0;
                QPM[NumFilled - 1] = TRUE;
                }                    
            }
        QPM[j] = TRUE;    
        if(UDV[ReadPres] <= 1 || OnlyReducingBandsums) SaveMinima = TRUE;        
        Input = RESET;
        goto _RESET;
        }
        
    /******************************************************************************************
        If the manifold is closed, and this is a balanced presentation, then this presentation
        may be pseudo-minimal. This is the case if SLength = Length. If this has happened and
        this is a new presentation we want to save copies of both this presentation and its
        dual presentation.
    ******************************************************************************************/    
            
    if(Input == DUALIZE)        
        {            
        if(SLength == Length)
            {
            /**********************************************************************************
                            Check whether the dual-presentation is on file.
            **********************************************************************************/    
            
            h = On_File();    
            if(h == NumFilled)
                {
                if(Dup_On_File < INFINITE)
                    {
                     if(Save_Pres(ReadPres,Dup_On_File,Length,1,108,1,0,0)) return(1);
                     Mark_As_Duplicate(Dup_On_File);
                     Input = RESET;
                     goto _RESET;                            
                     }
                 else
                     {
                     /**************************************************************************
                         The dual presentation is not on file. So we will save it. Next, check
                         whether the presentation given by the input relators is on file.
                     **************************************************************************/    
                     
                     if(Micro_Print)
                        {
                        printf("\n\nNow rewriting the original presentation.");
                        if(Micro_Print_F)
                            fprintf(myout,"\n\nNow rewriting the original presentation.");
                        }    
                             
                    BDY[NumFilled] = BDY[ReadPres];
                    UDV[NumFilled] = 0;
                    Canonical_Rewrite(DualRelators,FALSE,FALSE);
                    for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
                    if(SURL[i] == Length 
                        && NG[i] == NumGenerators
                         && NR[i] == NumRelators)
                         {
                         for(j = 1; j <= NumRelators; j++)
                             if(GetHandleSize((char **) DualRelators[j]) != GetHandleSize((char **) SUR[i][j])) break;
                         if(j > NumRelators && Compare_Dual_Pres(i))
                             {
                             if(ComponentNum[i] != CurrentComp)
                                 {
                                 if(Dup_On_File == INFINITE) Dup_On_File = i;
                                 }
                             else    
                                 break;
                             }    
                         }            
                    if(i == NumFilled)
                        {
                        /**********************************************************************
                            The presentation corresponding to the relators is not on file.
                        **********************************************************************/
                            
                        if(Dup_On_File < INFINITE)
                            {
                            if(Micro_Print)
                                {
                                printf("\n\nReverted to the undualized presentation.");
                                if(Micro_Print_F)
                                    fprintf(myout,"\n\nReverted to the undualized presentation.");
                                }                            
                             if(Save_Pres(ReadPres,Dup_On_File,Length,0,108,1,0,0)) return(1);                                
                             Mark_As_Duplicate(Dup_On_File);                             
                             }
                         else
                             {                                     
                            if(Save_Pres(ReadPres,NumFilled + 1,Length,1,108,1,0,0)) return(1);
                            
                            /******************************************************************
                                        Check whether the presentation is self-dual.
                            ******************************************************************/
                            
                            if(Compare_Dual_Pres(NumFilled - 1))
                                Daughters[NumFilled - 1] = NumFilled - 1;
                            else
                                {
                                if(Micro_Print)
                                    {
                                    printf("\n\nReverted to the undualized presentation.");
                                    if(Micro_Print_F)
                                        fprintf(myout,"\n\nReverted to the undualized presentation.");
                                    }                
                                if(Save_Pres(ReadPres,NumFilled - 1,Length,0,108,1,0,0)) return(1);         
                                BDY[NumFilled - 1] = BDY[ReadPres];
                                UDV[NumFilled - 1] = 0;
                                }
                            }
                        }
                    else
                        {
                        /**********************************************************************
                            The presentation corresponding to the relators is on file, but the
                            dual presentation is not. Save the dual presentation. Update
                            PRIM[i] to indicate that presentation i is pseudo-minimal and
                            update Daughters[i] to point to its new found dual. 
                        **********************************************************************/

                        if(Micro_Print)
                            {
                            printf("\n\nThe current presentation is on file, but the dual presentation is not.");
                            printf("\nSwapping the current presentation with the dual presentation.");
                            if(Micro_Print_F)
                                {
                                fprintf(myout,"\n\nThe current presentation is on file, but the dual presentation is not.");
                                fprintf(myout,"\nSwapping the current presentation with the dual presentation.");
                                }
                            }
                        if(Save_Pres(ReadPres,i,Length,1,108,1,0,0)) return(1);                        
                        if(PRIM[i] < 100) PRIM[i] += 100;
                        if(UDV[i] != SPLIT && UDV[i] != DUPLICATE) Daughters[i] = NumFilled - 1;    
                        }
                    }                
                }
            else
            if(Dup_On_File == INFINITE)
                {
                /******************************************************************************
                    The dual-presentation was previously found and is already on file.
                    Check whether the original presentation is on file. Update PRIM[h] to
                    reflect the fact that we now know this presentation is psuedo-minimal.
                    Update Daughters[h] to make it point to its new found dual.
                ******************************************************************************/
        
                if(Micro_Print)
                    {
                    printf("\n\nReverted to the undualized presentation.");
                    if(Micro_Print_F)
                        fprintf(myout,"\n\nReverted to the undualized presentation.");
                    }    
                
                if(PRIM[h] < 100) PRIM[h] += 100;        
                Canonical_Rewrite(DualRelators,FALSE,FALSE);
                for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
                if(SURL[i] == Length 
                    && NG[i] == NumGenerators
                     && NR[i] == NumRelators)
                     {
                     for(j = 1; j <= NumRelators; j++)
                         if(GetHandleSize((char **) DualRelators[j]) != GetHandleSize((char **) SUR[i][j])) break;
                     if(j > NumRelators && Compare_Dual_Pres(i))
                         {
                         if(ComponentNum[i] != CurrentComp)
                             {
                             if(Dup_On_File == INFINITE) Dup_On_File = i;
                             }
                         else    
                             break;
                         }    
                     }            
                if(i == NumFilled)
                    {
                    /**************************************************************************
                                The original presentation is not on file. So save it.
                    **************************************************************************/    
                    
                    if(Dup_On_File < INFINITE)
                        {
                         if(Save_Pres(ReadPres,Dup_On_File,Length,0,108,1,0,0)) return(1);                                
                        BDY[NumFilled - 1] = BDY[ReadPres];
                         Mark_As_Duplicate(Dup_On_File);
                         }
                     else
                         {
                        if(Save_Pres(ReadPres,h,Length,0,108,1,0,0)) return(1);             
                        BDY[NumFilled - 1] = BDY[ReadPres];
                        UDV[NumFilled - 1] = 0;
                        if(UDV[h] != SPLIT && UDV[h] != DUPLICATE) Daughters[h] = NumFilled - 1;
                        }
                    }
                else
                    {    
                    /**************************************************************************
                        Update Daughters[h] and Daughters[i]. Note that if i = h the
                                            presentation is self-dual. 
                    **************************************************************************/    

                    if(Micro_Print)
                        {
                        printf("\n\nBoth the current presentation and its dual presentation are already on file.");
                        if(Micro_Print_F)
                            fprintf(myout,"\n\nBoth the current presentation and its dual presentation are already on file.");
                        }
                                
                    if(UDV[h] != SPLIT && UDV[h] != DUPLICATE) Daughters[h] = i;
                    if(PRIM[i] < 100) PRIM[i] += 100;
                    if(UDV[i] != SPLIT && UDV[i] != DUPLICATE) Daughters[i] = h;
                    }
                }        
            Input = RESET;
            goto _RESET;                            
            }
            
        /**************************************************************************************
            If this presentation is not pseudo-minimal, but its length is shorter than
            that of any previous presentation we have for this summand, then we want to save
            a copy of this presentation.
        **************************************************************************************/    
                    
        if(Length < MLC[CurrentComp][NumGenerators])    
            {
            if(On_File() == NumFilled)
                {
                if(Dup_On_File < INFINITE)
                    {
                     if(Save_Pres(ReadPres,Dup_On_File,Length,1,10,1,0,0)) return(1);                                
                     Mark_As_Duplicate(Dup_On_File);
                     }
                 else
                     {
                    if(Save_Pres(ReadPres,0,Length,1,10,1,0,0)) return(1);                                 
                    BDY[NumFilled - 1]     = BDY[ReadPres];
                    UDV[NumFilled - 1]     = 0;
                    }    
                }
            Input = RESET;
            goto _RESET;
            }
        Input = NORMAL;        
        From_DUALIZE ++;                        
        }
        
    /******************************************************************************************
        Next, we do some crude checking designed to keep the program from getting "hung up" or
        otherwise spending too much time looking for something that may not be there.
    ******************************************************************************************/
                    
    if(Input == BANDSUM) 
        {
        if(Boundary)
            {
            if(GoingUp)
                {
                if(From_BANDSUM > (Vertices << 3))
                    {
                    Input = RESET;
                    goto _RESET;
                    }
                }
            else
                {
                if(Length >= SLength || (NumBandSum > (Vertices << 4) && Minimum >= 0L))
                    {
                    Input = RESET;
                    goto _RESET;
                    }                
                }
            }
        else if(NumBandSum > (Vertices << 2) || (!GoingUp && SLength == Length))
            {
            SaveMinima = TRUE;
            Input = RESET;
            goto _RESET;
            }
        }        
    if(Length > MaxTotalLength)
        {
        Input = RESET;
        goto _RESET;
        }    

_RESET:                        
    /******************************************************************************************
        When the program is finished processing one presentation and is ready for another it
        jumps to or returns to this point. Here we update some of the running data kept by the
        program, select the next presentation to be investigated and set some of the flags and
        parameters used by the rest of the program in processing a presentation.
    ******************************************************************************************/
    
    if((c = mykbhit()) || Level_Interrupt)
        {
        Level_Interrupt = FALSE;
        switch(c)
            {
            #ifdef DEBUGGING
            case 'd':
                fptr = stdout;
                Debug();
                break;
            case 'r':
                Print_Relators(Relators,NumRelators,stdout);
                break;
            #endif    
            case 's':
                printf("\n  Status: TotalAuts %lu  BandSums %ld  NumDiagrams %ld  NumDualized %lu  ToDo %u",
                    TotalAuts,Band_Sums,NumDiagrams,NumDualized,OnStack);
                break;                        
            case ' ':
                SMicro_Print = Micro_Print;
                SMicro_Print_F = Micro_Print_F;
                Micro_Print = FALSE;
                Micro_Print_F = FALSE;
                LIST_OPTIONS2:
                printf("\n\nHIT 'c' TO CHANGE SIMPLIFICATION PARAMETERS.");
                printf("\nHIT 'd' TO VIEW HEEGAARD DIAGRAMS OF THESE PRESENTATIONS.");
                printf("\nHIT 'p' TO PRUNE THE SEARCH TREE.");
                printf("\nHIT 't' TO TERMINATE THIS RUN.");
                printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS.");
                if(NumFilled > 1)
                    {
                    printf("\nHIT 'w' TO SORT THE PRESENTATIONS NOW IN MEMORY BY SUMMAND NUMBER,");
                    printf("\n        NUMGENERATORS, NUMRELATORS, LENGTH AND 'LEXICOGRAPHIC' ORDER.");
                    }
                printf("\nHIT 'x' TO TOGGLE MICRO_PRINTING ON AND OFF.");
                printf("\nHIT 'r' TO RESUME RUNNING THIS EXAMPLE.");
                GET_RESPONSE:
                switch(WaitkbHit())
                    {
                    case 'c':
                        printf("\n");
                        Get_Simplification_Parameters_From_User(TRUE,TRUE);
                        break;
                    case 'd':
                        Display_Diagrams();
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto LIST_OPTIONS2;                        
                    case 'p':
                        Prune_Search_Tree();
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto LIST_OPTIONS2;
                    case 'r':
                        break;    
                    case 't':
                        return(1);
                    case 'v':
                        Report(Band_Sums,NumDiagrams,OnStack,0,0,0,1,0,1,0);
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto LIST_OPTIONS2;
                    case 'w':
                        printf("\n\n     Sorting presentations. . . .");
                        Sort_Presentations_In_Memory();
                        SaveMinima = TRUE;
                        Input = RESET;
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
                        goto GET_RESPONSE;    
                    }
                printf("\n\n     Resumed running. . . .\n");
                NoReport = TRUE;
                Micro_Print = SMicro_Print;
                Micro_Print_F = SMicro_Print_F;
                break;            
            }
        }    
                    
    if(Input == RESET)
        {
        if(SaveMinima || ++Count >= MAX_COUNT)                                                
            {
            if(NumFilled >= MAX_SAVED_PRES - 3)
                {
                SysBeep(5);
                printf("\n\nThe program has saved the maximum number of presentations presently allowed.");
                return(STOP);
                }
            if(BytesUsed > BytesAvailable)
                {
                if(UserSaidQuit)
                    {
                    UserSaidQuit = FALSE;
                    return(STOP);
                    }
                if(UserSaidQuit = User_Says_Quit())
                    {
                    UserSaidQuit = FALSE;
                    return(STOP);
                    }
                }
            RESTART:            
            Count         = 0;
            SaveMinima     = FALSE;
            NumBandSum     = 0;
            TheComp     = 0;
            DoingDup     = FALSE;
            for(h = i = 0,j = VERTICES,MinLength = BIG_NUMBER; i < NumFilled; i++)
                {
                if(UDV[i] == SPLIT && CS[ComponentNum[i]] < 2)
                    {
                    jj = ComponentNum[Daughters[i]];
                    hh = jj + NCS[i];
                    for( ; jj < hh && CS[jj] == 2; jj++) ;
                    if(jj >= hh)
                        {
                        Mark_As_Found_Elsewhere(ComponentNum[i]);
                        goto RESTART;
                        }
                    }    
                if(UDV[i] == DUPLICATE)
                    {
                    ii = Daughters[i];
                    if(CS[ComponentNum[i] + 1] == 3 &&
                        (CBC[ComponentNum[ii]][0] < BDRY_UNKNOWN ||
                        CBC[ComponentNum[i]][0] < BDRY_UNKNOWN))
                        {
                        hh = ComponentNum[i];
                        jj = ComponentNum[ii];
                        if(CBC[jj][0] != BDRY_UNKNOWN && CBC[hh][0] == BDRY_UNKNOWN)
                        for(kk = 0; (CBC[hh][kk] = CBC[jj][kk]) < BDRY_UNKNOWN; kk++) ;    
                        if(MG_Bdry_Comp_Data(i)) goto RESTART;
                        }                    
                    if(UDV[ii] > UNKNOWN && CS[ComponentNum[i]] < 2)
                        {
                        Mark_As_Found_Elsewhere(ComponentNum[i]);
                        goto RESTART;
                        }
                    if(UDV[ii] < 2*NG[ii])
                        {
                        if(ComponentNum[i] < TheComp) continue;
                        if(ComponentNum[i] > TheComp)
                            {
                            TheComp = ComponentNum[i];
                            j = VERTICES;
                            DoingDup = TRUE;
                            }
                        if(NG[ii] < j)
                            {
                            j             = NG[ii];
                            MinNR        = NR[ii];
                            MinLength     = SURL[ii];
                            MaxUDV         = UDV[ii];
                            ReadPres     = ii;
                            }
                        else
                        if(NG[ii] == j)
                            {
                            if(NR[ii] < MinNR)
                                {
                                MinNR         = NR[ii];
                                MinLength     = SURL[ii];
                                MaxUDV         = UDV[ii];
                                ReadPres     = ii;
                                }
                            else
                            if(NR[ii] == MinNR)
                                {
                                if(SURL[ii] < MinLength)
                                    {
                                    MinLength     = SURL[ii];
                                    MaxUDV         = UDV[ii];
                                    ReadPres     = ii;
                                    }
                                else
                                if(SURL[ii] == MinLength && UDV[ii] < MaxUDV)
                                    {
                                    MaxUDV         = UDV[ii];
                                    ReadPres     = ii;
                                    }
                                }    
                            }    
                        }
                    else
                        if(UDV[ii] == 2*NG[ii]) UDV[ii] = DONE;
                    continue;
                    }    
                if(UDV[i] < 2*NG[i])
                    {
                    h += 2*NG[i] - UDV[i];
                    if(ComponentNum[i] < TheComp && !DoingDup) continue;
                    if(ComponentNum[i] > TheComp)
                        {
                        TheComp = ComponentNum[i];
                        j = VERTICES;
                        DoingDup = FALSE;
                        }
                    if(NG[i] < j)
                        {
                        j             = NG[i];
                        MinNR        = NR[i];
                        MinLength     = SURL[i];
                        MaxUDV         = UDV[i];
                        ReadPres     = i;
                        }
                    else
                    if(NG[i] == j)
                        {
                        if(NR[i] < MinNR)
                            {
                            MinNR         = NR[i];
                            MinLength     = SURL[i];
                            MaxUDV         = UDV[i];
                            ReadPres     = i;
                            }
                        else
                        if(NR[i] == MinNR)
                            {
                            if(SURL[i] < MinLength)
                                {
                                MinLength     = SURL[i];
                                MaxUDV         = UDV[i];
                                ReadPres     = i;
                                }
                            else
                            if(SURL[i] == MinLength && UDV[i] < MaxUDV)
                                {
                                MaxUDV         = UDV[i];
                                ReadPres     = i;
                                }
                            }    
                        }    
                    }
                else
                    if(UDV[i] == 2*NG[i]) UDV[i] = DONE;
                }    
            OnStack = h - 1;
            SReadPres = ReadPres;
                                        
            if(j == VERTICES && MinLength == BIG_NUMBER)
                {
                OnStack = 0;
                if(Find_All_Min_Pres) switch(Check_Level_Transformations())
                    {
                    case FALSE:
                        break;
                    case TRUE:
                        Input = RESET;
                        SaveMinima = TRUE;
                        goto _RESET;
                    case 2:
                        Input = NORMAL;
                        SRError = TRUE;
                        From_BANDSUM = 0;    
                        From_DUALIZE = 0;
                        goto _REDUCE_GENUS;
                    }        
                SysBeep(5);
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);
                return(1);
                }
                                
            CurrentComp         = ComponentNum[ReadPres];
            NumGenerators         = NG[ReadPres];
            NumRelators         = NR[ReadPres];
            Vertices             = NumGenerators << 1;
            if(SURL[ReadPres] < MTLK)
                MaxTotalLength     = SURL[ReadPres] << 4;
            else
                MaxTotalLength  = SURL[ReadPres] << 3;    
                        
            if(CS[CurrentComp] >= 2) goto RESTART;
            
            if(CS[TheComp] == 1)
                {
                for(i = 0; i < NumFilled; i++) if(UDV[i] == SPLIT && ComponentNum[i] == TheComp)
                    {
                    jj = ComponentNum[Daughters[i]];
                    hh = jj + NCS[i];
                    for( ; jj < hh && CS[jj] == 2; jj++) ;
                    if(jj >= hh) break;
                    }
                if(i < NumFilled)    
                    {
                    Mark_As_Found_Elsewhere(TheComp);
                    goto RESTART;    
                    }
                else
                    {
                    printf("\n    Rechecking presentations of Summand %u. This Summand previously 'split'.",
                        TheComp);
                    CS[TheComp] = EOS;    
                    }    
                }
                
            From_BANDSUM = 0;    
            From_DUALIZE = 0;    
            
            /**********************************************************************************
                If this is a new presentation, as indicated by the fact that UDV[ReadPres] =
                0, then call Test_New_Pres() to run some initializing routines on this
                presentation.
            **********************************************************************************/    
            
            if(UDV[ReadPres] == 0)
                {
                if(Micro_Print) Micro_Print_Reset();                
                switch(Test_New_Pres())
                    {
                    case STOP:
                        Fatal_Error();
                    case TOO_MANY_COMPONENTS:    
                        return(STOP);
                    case DUALIZE:
                        Input = DUALIZE;
                        goto _DUALIZE;
                    case BANDSUM:
                        Input = BANDSUM;
                        From_BANDSUM = 0;
                        NumBandSum = 0;
                        goto _BANDSUM;    
                    case REDUCE_GENUS:
                        Input = NORMAL;
                        SRError = TRUE;
                        goto _REDUCE_GENUS;
                    case RESET:
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto _RESET;    
                    }
                }            
            Boundary = BDY[ReadPres];
            if(Boundary)
                {
                if(UDV[ReadPres])
                    GoingUp = TRUE;
                else
                    GoingUp = FALSE;
                }
            else
                GoingUp = TRUE;    
            UDV[ReadPres] ++;
            
            if(ER[ReadPres] < 0 && NumGenerators > 1)
                {
                /******************************************************************************
                    Since ER[ReadPres] < 0, this presentation contains redundant relator(s).
                    Check whether the program can find a new way to delete one or more of these
                    relators.
                ******************************************************************************/
                
                if(Micro_Print) Micro_Print_Reset();
            
                if(Get_Relators_From_SUR(ReadPres))
                    {
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;                    
                    }
    
                Fill_A(NumRelators);
                NumDiagrams ++;
                DrawingDiagrams = TRUE;
                WhichInput = ReadPres;
                if(Whitehead_Graph())
                    {
                    DrawingDiagrams = FALSE;
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }
                DrawingDiagrams = FALSE;                                
                Get_Bdry_Comps(FALSE,FALSE,ReadPres);
                if(BCWG[0] > 1 || (BCWG[0] && NumBdryComps > BCWG[0]))                            
                    Delete_Redundant_Relators();
                else
                    ER[ReadPres] = 0;    
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }
                
            if(TP[ReadPres] && NumRelators > 1 && !Do_Not_Reduce_Genus && UDV[ReadPres] > 2)                    
                {
                /******************************************************************************
                    Since TP[ReadPres] is TRUE, check whether there are primitives, which the
                        program has not previously discovered, in this presentation.
                ******************************************************************************/    
                
                if(Micro_Print) Micro_Print_Reset();
                
                if(Get_Relators_From_SUR(ReadPres))
                    {
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }
                
                if(SetUp_TopOfChain())
                    {
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }
                    
                switch(Reduce_Genus(NORMAL,FALSE,TRUE))
                    {
                    case NO_ERROR:
                        break;
                    case FATAL_ERROR:
                        Fatal_Error();    
                        Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);    
                        return(STOP);
                    case TOO_LONG:
                    case CAN_NOT_DELETE:
                        if(Micro_Print)
                            {
                            printf("\n\nUnable to delete a primitive.");
                            if(Micro_Print_F)
                                fprintf(myout,"\n\nUnable to delete a primitive.");
                            }
                        TP[ReadPres] --;                            
                        FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto _RESET;    
                    }    
                if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
                    {
                    Input = NORMAL;
                    SRError = TRUE;
                    goto _REDUCE_GENUS;
                    }                    
                TP[ReadPres] = FALSE;
                UDV[ReadPres] --;
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;            
                }
                
            if((!FormBandsums && UDV[ReadPres] > 2)
                || (OnlyReducingBandsums && (UDV[ReadPres] > 2 || QPM[ReadPres])))                
                {
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;                
                }
            if(Micro_Print) Micro_Print_Reset();        
            }    
        From_BANDSUM = 0;    
        From_DUALIZE = 0;        
        if(Count == 0)
            {
            if(Get_Relators_From_SUR(ReadPres))
                {
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }            

            if(SetUp_TopOfChain())
                {
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }            
            }                    
        else    
            {
            NumGenerators     = NG_TOC;
            NumRelators     = NR_TOC;
            Vertices         = 2*NG_TOC;
            Length             = TOCLength;
            for(i = 1; i <= NumRelators; i++)
                {
                ReallocateHandle((char **) Relators[i],GetHandleSize((char **) TopOfChain[i]));
                if((q = *Relators[i]) == NULL)
                    {
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }                
                p = *TopOfChain[i];
                while(*q++ = *p++) ;        
                }
            GoingUp     = TRUE;    
            NumBandSum     = 0;
            if(Micro_Print)
                {
                printf("\n\nSet Relators[] = TopOfChain[], Length = %lu:\n",Length);
                Print_Relators(Relators,NumRelators,stdout);
                if(Micro_Print_F)
                    {
                    fprintf(myout,"\n\nSet Relators[] = TopOfChain[], Length = %lu:\n",Length);
                    Print_Relators(Relators,NumRelators,myout);
                    }
                }    
            }        
        Fill_A(NumRelators);
        Automorphisms = 0L;
        }        
    NumDiagrams ++;
    
    /******************************************************************************************
        At this point, we call Whitehead_Graph() which determines whether the presentation
        given in Relators[] is realizable and if the presentation is realizable finds the
        Heegaard diagram.
    ******************************************************************************************/    
    
    if(Flag2 = Whitehead_Graph())
        {
        /**************************************************************************************
            If Whitehead_Graph() returns a non-zero value, it indicates either a fatal error
            of some sort has occured, or the program had some difficulty finding the
            diagram. The actual return value indicates what the problem was.
        **************************************************************************************/    
        
        switch(Flag2)
            {
            case NON_PLANAR:
            case FATAL_ERROR:
                if(Flag2 == NON_PLANAR)
                    {
                    fprintf(stdout,"\n\n                    The Whitehead graph is nonplanar.");
                    fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");
                    }
                Fatal_Error();        
                return(1);
            case TOO_LONG:
            case V2_ANNULUS_EXISTS:
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
            case TOO_MANY_COMPONENTS:
                return(1);    
            default:
                break;
            }                                    
        if((From_DUALIZE || (Boundary && From_BANDSUM)) && Length < SURL[ReadPres]) switch(Flag2)
            {
            case NON_UNIQUE_1:
            case NON_UNIQUE_2:
            case NON_UNIQUE_3:
            case NON_UNIQUE_4:
            case SEP_PAIRS:
            
            /**********************************************************************************
                We are in a descending chain of diagrams, and the program cannot determine
                what the diagram associated with the current presentation should be.
                So this descending chain must be terminated. However we want to save this
                presentation if it is not already on file and if the problem is that the 
                diagram has separating pairs of vertices, then we give Level_Transformations()
                a try.
            **********************************************************************************/
                                        
            if(On_File() == NumFilled)
                {
                if(NumFilled >= MAX_SAVED_PRES - 3)
                    {
                    SaveMinima = TRUE;
                    Input = RESET;
                    goto _RESET;
                    }
                if(Dup_On_File < INFINITE)
                    {
                     if(Save_Pres(ReadPres,Dup_On_File,Length,1,11,1,0,0)) return(1);
                     Mark_As_Duplicate(Dup_On_File);                 
                     SaveMinima = TRUE;                                
                     }
                 else
                     {
                    if(Save_Pres(ReadPres,0,Length,1,11,1,0,0)) return(1);        
                    BDY[NumFilled - 1] = BDY[ReadPres];
                    switch(Flag2)
                        {
                        case SEP_PAIRS:
                            Fill_A(NumRelators);
                            Get_Matrix();
                            Sep_Pairs(0,0);
                            UDV[NumFilled - 1] = SEP_PAIRS;
                            if(V1 & 1)
                                LSP[NumFilled - 1] = V1/2 + 97;
                            else
                                LSP[NumFilled - 1] = V1/2 + 65;
                            if(V2 & 1)
                                LSQ[NumFilled - 1] = V2/2 + 97;
                            else
                                LSQ[NumFilled - 1] = V2/2 + 65;    
                            break;    
                        case NON_UNIQUE_1:    
                        case NON_UNIQUE_2:    
                        case NON_UNIQUE_3:    
                        case NON_UNIQUE_4:
                            UDV[NumFilled - 1] = Flag2;
                            break;    
                        }    
                    if(Flag2 == SEP_PAIRS)
                        {
                        NumCalled     = 0;
                        NotNewPres     = 0;
                        SNumFilled     = NumFilled;
                        ReadPres     = NumFilled - 1;
                        h = Level_Transformations(TRUE,!Find_All_Min_Pres,FALSE,0,0);                                                                                
                        for(i = 0; i < NumCalled; i++)
                        for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
                        if(Level_Interrupt == 1)
                            {
                            Input = RESET;
                            goto _RESET;    
                            }
                        if(h == FATAL_ERROR)
                            {
                            Fatal_Error();
                            return(1);
                            }
                        if(Micro_Print && h != 2)
                            {
                            printf("\n\nThe current Presentation is:\n");
                            Print_Relators(Relators,NumRelators,stdout);
                            if(Micro_Print_F)
                                {
                                fprintf(myout,"\n\nThe current Presentation is:\n");
                                Print_Relators(Relators,NumRelators,myout);
                                }
                            }
                        if(h == 5 && !Do_Not_Reduce_Genus)
                        switch(Delete_Trivial_Generators(FALSE))
                            {
                            case 0:
                                break;
                            case 1:
                                FoundPrimitive = TRUE;
                                SRError = TRUE;                                    
                                goto _REDUCE_GENUS;
                            case TOO_LONG:
                                SaveMinima = TRUE;
                                Input = RESET;
                                goto _RESET;    
                            }
                        if(Find_All_Min_Pres && SNumFilled == NumFilled
                            && !Init_Find_Level_Transformations(FALSE))
                            {
                            switch(Find_Level_Transformations(Delete_Only_Short_Primitives,0))
                                {
                                case 0:
                                case 1:
                                case 2:
                                    break;
                                case 3:
                                    Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                                    ReadPres = NumFilled - 1;
                                    switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                                        {
                                        case NO_ERROR:
                                            break;
                                        case FATAL_ERROR:
                                            Fatal_Error();
                                            return(1);
                                        case TOO_LONG:
                                        case CAN_NOT_DELETE:
                                            if(Micro_Print)
                                                {
                                                printf("\n\nUnable to delete a primitive.");
                                                if(Micro_Print_F)
                                                    fprintf(myout,"\n\nUnable to delete a primitive.");
                                                }
                                            FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;                                        
                                            SaveMinima = TRUE;
                                            Input = RESET;
                                            goto _RESET;    
                                        }
                                    SRError = TRUE;                
                                    goto _REDUCE_GENUS;
                                case 5:
                                    Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                                    SaveMinima = TRUE;
                                    Input = RESET;
                                    goto _RESET;
                                default:
                                    break;    
                                }
                            Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                            }            
                        if(SNumFilled < NumFilled) SaveMinima = TRUE;        
                        }        
                    Input = RESET;
                    goto _RESET;
                    }
                }
            Input = RESET;
            goto _RESET;    
            }
        if(Flag2 == NOT_CONNECTED)
            {
            /**********************************************************************************
                The Whitehead graph corresponding to the current presentation is not connected.
                Hence the corresponding Heegaard diagram is reducible and the presentation has
                been "split" by Whitehead_Graph().
            **********************************************************************************/
            
            if(NumFilled >= MAX_SAVED_PRES - 3)
                {
                fprintf(stdout,"\n\nWe have run out of memory set aside for presentations!!!!");
                SysBeep(5);
                Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);
                return(1);
                }    
            Input = RESET;
            goto _RESET;
            }
        if(Flag2 == REDUCE_GENUS)
            {
            /**********************************************************************************
                The Whitehead graph is not connected, but there are primitives present which
                means that the genus of the diagram can be reduced.
            **********************************************************************************/    
            
            Input = NORMAL;
            if(TestRealizability1 || TestRealizability2)
                Delete_Trivial_Generators(FALSE);
            FoundPrimitive = TRUE;
            SRError = TRUE;
            goto _REDUCE_GENUS;
            }        
        Input = RESET;
        goto _RESET;
        }    
        
        /**************************************************************************************
            After returning from Whitehead_Graph(), if the program was able to find the
            Heegaard diagram without problems, execution returns to this point.
        **************************************************************************************/    
                                                    
        if(Automorphisms)
            {
            if(Boundary || NumGenerators != NumRelators)
                {
                if(GoingUp)
                    {
                    GoingUp = FALSE;
                    if(SetUp_TopOfChain())
                        {
                        SaveMinima = TRUE;
                        Input = RESET;
                        goto _RESET;
                        }    
                    }    
                Input = BANDSUM;
                From_BANDSUM = 0;
                goto _BANDSUM;
                }
            if(!From_DUALIZE && SetUp_TopOfChain())
                {
                SaveMinima = TRUE;
                Input = RESET;
                goto _RESET;
                }
            Input = DUALIZE;
            goto _DUALIZE;
            }    
        Input = BANDSUM;
        goto _BANDSUM;
            
/**********************************************************************************************
                     THIS IS THE END OF THE MAIN LOOP OF THE PROGRAM.
**********************************************************************************************/                         
}                

void Fill_A(int NumRelators)
{
    /******************************************************************************************
        This routine takes the cyclic relator pointed to by "p" and examines each pair of 
        consecutive letters that appear in that relator, including the pair comprised of the
        last letter and the first letter. For each such pair of letters, it adds the 
        appropriate edge to the array A[][].
    ******************************************************************************************/
    
    register unsigned char     i,
                            j,
                            *p,
                            x;
                            
    register unsigned int    *q;                            
    
    for(i = 0; i < Vertices; i++)
    for(j = 0,q = A[i]; j < Vertices; j++,q++) *q = 0;
    for(j = 1; j <= NumRelators; j++)
        {
        p = *Relators[j];
        if(*p == EOS) continue;
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
            A[i][x]++;
            if(x & 1) i = x - 1;
            else i = x + 1;
            p++;
            }
        x = **Relators[j];
        x = x << 1;            
        if(x < 194) x -= 130;
        else x -= 193;
        A[i][x]++;
        }
    for(i = 0; i < Vertices - 1; i++)
    for(j = i + 1; j < Vertices; j++)
        {
        A[i][j] += A[j][i];
        A[j][i] = A[i][j];
        }                        
}

void Fill_AA(int NumRelators)
{
    /******************************************************************************************
        This routine takes the cyclic relator pointed to by "p" and examines each pair of 
        consecutive letters that appear in that relator, including the pair comprised of the
        last letter and the first letter. For each such pair of letters, it adds the 
        appropriate edge to the array AA[][].
    ******************************************************************************************/
    
    register unsigned char     i,
                            j,
                            *p,
                            x;
                            
    register unsigned int     *q;                        
    
    for(i = 0; i < Vertices; i++)
    for(j = 0,q = AA[i]; j < Vertices; j++,q++) *q = 0;
    for(j = 1; j <= NumRelators; j++)
        {
        p = *Relators[j];
        if(*p == EOS) continue;
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
        x = **Relators[j];
        x = x << 1;            
        if(x < 194) x -= 130;
        else x -= 193;
        AA[i][x]++;
        }
    for(i = 0; i < Vertices - 1; i++)
    for(j = i + 1; j < Vertices; j++)
        {
        AA[i][j] += AA[j][i];
        AA[j][i] = AA[i][j];
        }                    
}

void UpDate_Fill_A(void)
{
    /******************************************************************************************
        This routine takes the cyclic relator pointed to by "p" and examines each pair of 
        consecutive letters that appear in that relator, including the pair comprised of the
        last letter and the first letter. For each such pair of letters, it adds the 
        appropriate edge to the array A[][]. This version of FillA() is used upon returning
        from New_Relator().
    ******************************************************************************************/
    
    register unsigned char     i,
                            *p,
                            x;    
    
    p = *Temp3;
    if(*p)
        {
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
            A[i][x] --;
            A[x][i] --;
            if(x & 1) i = x - 1;
            else i = x + 1;
            p++;
            }
        x = **Temp3;
        x = x << 1;            
        if(x < 194) x -= 130;
        else x -= 193;
        A[i][x] --;
        A[x][i] --;
        }    
    p = *Relators[1];
    if(*p)
        {
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
            A[i][x] ++;
            A[x][i] ++;
            if(x & 1) i = x - 1;
            else i = x + 1;
            p++;
            }
        x = **Relators[1];
        x  = x << 1;            
        if(x < 194) x -= 130;
        else x -= 193;
        A[i][x] ++;
        A[x][i] ++;
        }                
}

Find_Flow_A(Input,F1)
int     Input,
        F1;
{
    /******************************************************************************************
        If there is a Whitehead T-transformation that reduces the length of the relators in
        Relators[], this routine will find one. This is a reasonably efficient polynomial time
        procedure, which finds a minimal "cut set" of edges of a network by finding a 
        maximal "flow" in the network. The procedure works even if the presentation is not 
        realizable. While if the presentation is realizable, then the T-transformation
        produced by this procedure is a geometric T-transformation--provided the Whitehead
        graph is connected, and there are no edges of the Heegaard diagram which form "waves".
    ******************************************************************************************/
            
    register unsigned int     i,
                            j,
                            max,
                            min,
                            *p,
                            *q,
                            *r;
                            
    unsigned int             Flow,
                            GL[VERTICES/2],
                            h,
                            k,
                            MaxFlow,
                            MinExp,
                            S[VERTICES],
                            Sink,
                            Source;
                            
    unsigned int            Get_MinExp();                        
    
    if(NumGenerators == 0) return(TOO_LONG);
    
    if(Input == BANDSUM)
        UpDate_Fill_A();
    else
        Fill_A(NumRelators);        
    if(ComputeValences_A()) return(TOO_LONG);                    
    for(i = 0; i < Vertices; i++)
        {
        for(j = k = 0; j < Vertices; j++) if(A[i][j])
            {
            AJ3[i][k] = j;
            k++;
            }
        AJ3[i][k] = VERTICES;
        }    
    for(i = 0; i < NumGenerators; i++) GL[i] = i;
    k = NumGenerators;
    Automorphisms = 0L;        
    while(1)
        {
        do
            {
            i = abs(rand()) % k;
            k--;
            h = GL[i];
            GL[i] = GL[k];
            GL[k] = h;
            }
        while(k && VA[h] == 0);        /* Choose a random generator with non-zero valance to work on. */
        if(k == 0 && VA[h] == 0) return(TRUE);
        Source             = h << 1;
        Sink             = Source + 1;
        MaxFlow         = VA[h];
        Flow             = A[Source][Sink];
        A[Source][Sink] = 0;
        A[Sink][Source] *= 2;
        if(Flow && (MaxFlow > Flow))
            MinExp = Flow/(MaxFlow - Flow) + 1;
        else
            MinExp = 1;                
        while(MaxFlow)
            {        
            for(i = 0,p = ZZ,q = InQueue; i < Vertices; i++,p++,q++) *p = *q = 0;
            ZZ[Sink] = INFINITE;
            InQueue[Source] = TRUE;        /*  Setting InQueue[Source] = TRUE prevents unnecessary 
                                        attempts to change paths from vertices to the Sink to
                                        paths which pass through the Source.                 */
                                        
            for(r = UpDate,*r = Sink,p = r + 1; r < p; r++)
                {
                /******************************************************************************
                    ZZ[j] holds the current value of a maximal width path of edges from vertex
                    j to the Sink. S[j] is the neighbor of vertex j to which this path proceeds
                    from vertex j. UpDate is a queue of vertices whose width of path to the 
                    sink has been increased and whose neighbors need to be examined to see if 
                    there is a path to the sink of greater width than their current path to the
                    Sink. InQueue[] is used to prevent redundant entries of a vertex in UpDate.
                        Note that there is no point in updating paths from vertices to the
                    Sink unless the updated path has width greater than the current width of
                    the path from Source to Sink. Hence no updating is done unless the path has
                    width greater than the current value of ZZ[Source].
                ******************************************************************************/
                
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
                        if(min > ZZ[j])    
                            {
                            /******************************************************************
                                (Note : A[j][i] gives the current number of edges available for 
                            a path from vertex j to vertex i.)
                                There is a new path of 'min' parallel edges from vertex j to 
                            the Sink which leaves vertex j and goes to vertex i. And 'min' is
                            greater than the current flow ZZ[j] from vertex j to the Sink.
                                Increase the flow from vertex j to the sink to the value 'min', 
                            and set S[j] equal to i. If 'min' is greater than the current flow 
                            from Source to Sink, and the flows from neighbors of vertex j to j
                            are not already slated to be updated, set InQueue[j] 'TRUE' so we 
                            will check to see if the flow from neighboring vertices of vertex 
                            j can be increased. Also enter vertex j in the UpDate queue.
                                (The function of InQueue[] is to keep a vertex from being
                            entered more than once in the array UpDate[] between times it is
                            reexamined.)
                            ******************************************************************/
                            
                            ZZ[j] = min;
                            S[j] = i;
                            if(!InQueue[j] && min > ZZ[Source])
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
            if(Flow >= MaxFlow) break;
            if(max == 0) break;
            i = Source;
            while(i != Sink)
                {
                /******************************************************************************
                    We have found a path of 'max' parallel edges from Source to Sink. Update
                    the entries of A[][] to reflect this.
                ******************************************************************************/
                j = S[i];
                A[i][j] -= max;                                                                
                A[j][i] += max;                                    
                i = j;
                }                                                
            }
                                                                    
        /**************************************************************************************
                If MaxFlow is greater than Flow,then the flow from source to sink is less 
                than the valence of the source. Hence there is a T-transformation which 
                reduces the length of the relators.Increment the number of automorphisms. 
                Set k = NumGenerators.
        **************************************************************************************/
        
        if(MaxFlow > Flow)    
            {
            if(F1) return(2);    
            k = NumGenerators;
            if(MinExp > 1) MinExp = Get_MinExp(Source,NumRelators);
            Automorphisms += MinExp;
            if(MinExp > 2)
                {
                if(Do_Auts(Source,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                }
            else
                if(Do_Aut(Source,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
            }    
        if(k == NumGenerators)
            {                                                
            Fill_A(NumRelators);        
            Length -= MinExp*(MaxFlow - Flow);
            VA[h] -= MinExp*(MaxFlow - Flow);
            for(i = 0; i < Vertices; i++)
                {
                for(j = 0,p = AJ3[i]; j < Vertices; j++) if(A[i][j])
                    {
                    *p = j;
                    p++;
                    }
                *p = VERTICES;
                }
            }
        else
            {
            for(i = 1; i < Vertices; i++)
            for(p = AJ3[i]; (j = *p) < i; p++)
            A[i][j] = A[j][i] = (A[i][j] + A[j][i]) >> 1;
            if(k == 0)
                {
                /******************************************************************************
                    If k is equal to zero, then we have checked each generator in turn without
                    finding any (new) automorphism that reduces the length. Hence the current
                    set of relators has minimal length. Next we check whether there are any
                    generators which no longer appear in the relators. This condition may have
                    occured through the action of an automorphism or by the replacement of a
                    relator with a new relator via New_Relator() or Reduce_Genus(). If this has
                    happened, Find_Flow_A returns TRUE and we call Missing_Gen() in order to
                    guarantee that the Whitehead graph of the presentation Relators[] has no
                    isolated vertices.
                ******************************************************************************/
                
                for(i = 0 ; i < NumGenerators; i++) if(VA[i] == 0) return(TRUE);
                return(FALSE);        
                }
            }                                                        
        }
}

Find_Primitives(flag)
int flag;
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
        know that Relators[1] is a proper power. As in the case when the parameter flag = 1, 
        we want to perform enough automorphisms to reduce Relators[1] to a proper power 
        and then make the appropriate substitutions in the remaining relators.
    ******************************************************************************************/
    
    register unsigned char    *p;
        
    register unsigned int     i,
                            j,
                            k,
                            *q;
                            
    unsigned int             C[125],
                            MinExp,
                            VG[VERTICES],
                            NumIV,
                            root;
                            
    unsigned long            SLength1,
                            SLength2,
                            SLength3;
                            
    unsigned int            Get_MinExp();                                                    
    
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
        Fill_AA(1);
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
            Count the number of isolated vertices, i.e. the number of AJ3[][0] with AJ3[][0] =
            VERTICES. Set NumIV equal to the number of isolated vertices, and set up the array 
            IV[] to hold a list of these.
        **************************************************************************************/
        
        for(i = k = 0; i < Vertices; i++)
            {
            if(AJ3[i][0] == VERTICES)
                {
                IV[i] = 2;
                k++;
                }
            else
                IV[i] = 0;
            }    
        NumIV = k;
        if(NumIV >= Vertices - 2) return(2);
            
        /**************************************************************************************
             Check whether the graph obtained by ignoring the isolated vertices is connected.
         *************************************************************************************/
         
        for(i = 0; i < Vertices; i++) ZZ[i] = IV[i];
        for(i = 0; i < Vertices && ZZ[i]; i++) ;
        
        if(Connected_AJ3(i,NumIV) == FALSE)
            {
            /**********************************************************************************
                If the graph is not connected, then there are two vertices corresponding to the
                same generator which are in different components. After we have found such a
                pair of vertices, we perform an automorphism which gets rid of all appearances
                of the corresponding generator in Relators[1].
            **********************************************************************************/
            
            for(i = 0; i < Vertices; i += 2)
                if((!ZZ[i] && ZZ[i+1]) || (ZZ[i] && !ZZ[i+1])) break;
            if(ZZ[i])
                {
                for(j = 0; j < Vertices; j++)
                    {
                    if(ZZ[j] == 1)
                        ZZ[j] = 0;
                    else
                    if(ZZ[j] == 0)
                        ZZ[j] = 1;
                    }
                }
            switch(flag)
                {
                case 1:
                    if(Do_Aut(i,1,NumRelators) == TOO_LONG) return(TOO_LONG);
                    break;
                case 2:
                    if(Do_Aut(i,1,1) == TOO_LONG) return(TOO_LONG);
                    break;
                case 3:
                    if(Do_Aut(i,1,NumRelators) == TOO_LONG) return(TOO_LONG);
                    break;
                }
            if(flag == 1 || flag == 2)
                {
                SLength2 = GetHandleSize((char **) Relators[1]);
                C[i] -= SLength1 - SLength2;            
                SLength1 = SLength2;
                }                    
            }    
        else
            {
            /**********************************************************************************
                When the program gets here, the graph is connected, but may have isolated
                vertices. Find the first non-isolated vertex and use this vertex as the root
                of a depth-first-search of the graph.
            **********************************************************************************/    
            
            for(i = 0; i < Vertices && IV[i]; i += 2) ;
            for(j = 0; j < Vertices; j++) Number[j] = 0;
            NumVert     = 1;
            Number[i]     = 1;
            Lowpt[i]     = 1;
            Father[i]     = i;
            VG[i]         = 0;
            root        = i;
            for(q = UpDate,*q = root; q >= UpDate; q--)
                {
                NEW_VERT:
                i = *q;
                for(k = VG[i]; (j = AJ3[i][k]) < VERTICES; k++)
                    {
                    if(Number[j] == 0)
                        {
                        NumVert     ++;
                        Number[j]     = NumVert;
                        Lowpt[j]     = NumVert;
                        Father[j]     = i;
                        VG[j]         = 0;
                        VG[i]         = k + 1;
                        q             ++;
                        *q             = j;
                        goto         NEW_VERT;        
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
                Otherwise, we look for other cut vertices.
            **********************************************************************************/
                        
            if(k <= 1)
                {        
                for(j = 0; j < Vertices; j++)
                    {
                    k = Father[j];
                    if(IV[j] == 0 && Lowpt[j] >= Number[k] && k != root)
                        break;
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
                
            for(i = 0; i < Vertices; i++) ZZ[i] = IV[i];
            i = Father[j];
            ZZ[i] = 3;
            if(i & 1)
                {
                Connected_AJ3(i-1,NumIV+1);
                ZZ[i] = 0;
                for(j = 0; j < Vertices; j++)
                    {
                    if(ZZ[j] == 1)
                        ZZ[j] = 0;
                    else
                    if(ZZ[j] == 0)
                        ZZ[j] = 1;
                    }
                SLength3 = AA[i-1][i];
                if(SLength3 == 0L)
                    MinExp = 1;
                else
                    {
                    for(j = 0,SLength2 = 0L; (k = AJ3[i][j]) < VERTICES; j++)
                    SLength2 += AA[i][k];
                    SLength2 -= SLength3;
                    if(SLength2)
                        MinExp = SLength3/SLength2 + 1;
                    else
                        MinExp = 1;    
                    }        
                switch(flag)
                    {
                    case 1:
                        if(MinExp > 1) MinExp = Get_MinExp(i-1,NumRelators);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i-1,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i-1,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                        break;
                    case 2:
                        if(MinExp > 1) MinExp = Get_MinExp(i-1,1);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i-1,MinExp,1) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i-1,MinExp,1) == TOO_LONG) return(TOO_LONG);
                        break;
                    case 3:
                        if(MinExp > 1) MinExp = Get_MinExp(i-1,NumRelators);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i-1,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i-1,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                        break;
                    }
                if(flag == 1 || flag == 2)
                    {
                    SLength2 = GetHandleSize((char **) Relators[1]);
                    C[i-1] -= SLength1 - SLength2;            
                    if(C[i-1] == 1) return(TRUE);
                    SLength1 = SLength2;
                    }            
                }
            else
                {
                Connected_AJ3(i+1,NumIV+1);
                ZZ[i] = 0;
                SLength3 = AA[i][i+1];
                if(SLength3 == 0L)
                    MinExp = 1;
                else
                    {
                    for(j = 0,SLength2 = 0L; (k = AJ3[i][j]) < VERTICES; j++)
                    SLength2 += AA[i][k];
                    SLength2 -= SLength3;
                    if(SLength2)
                        MinExp = SLength3/SLength2 + 1;
                    else
                        MinExp = 1;    
                    }    
                switch(flag)
                    {
                    case 1:
                        if(MinExp > 1) MinExp = Get_MinExp(i,NumRelators);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);            
                        break;
                    case 2:
                        if(MinExp > 1) MinExp = Get_MinExp(i,1);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i,MinExp,1) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i,MinExp,1) == TOO_LONG) return(TOO_LONG);
                        break;
                    case 3:
                        if(MinExp > 1) MinExp = Get_MinExp(i,NumRelators);
                        if(MinExp > 2)
                            {
                            if(Do_Auts(i,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                            }
                        else
                            if(Do_Aut(i,MinExp,NumRelators) == TOO_LONG) return(TOO_LONG);
                        break;
                    }
                if(flag == 1 || flag == 2)
                    {
                    SLength2 = GetHandleSize((char **) Relators[1]);
                    C[i] -= SLength1 - SLength2;            
                    if(C[i] == 1) return(TRUE);
                    SLength1 = SLength2;
                    }            
                }    
            }    
        }    
}
    
Connected_AJ3(i,k)
register unsigned int     i,
                        k;
{
    /******************************************************************************************
        This routine finds those vertices in the component of vertex i in the graph specified
        in the adjacency lists AJ3[]. The array ZZ[] is initialized by the calling routine
        which sets the entries of vertices which should be deleted from the adjacency lists to
        a non-zero value and passes the number of deleted vertices as the parameter k. 
        The routine returns FALSE if the graph is not connected and TRUE if it is connected.
    ******************************************************************************************/    
     
    register unsigned int     h,
                            j,
                            *p,
                            *r;
     
    ZZ[i] = 1;
    k ++;
    for(r = UpDate,*r = i,p = r + 1; r < p; r++)
        {
        i = *r;
        for(h = 0; (j = AJ3[i][h]) < VERTICES; h++)
            {
            if(i == j) continue;
            if(ZZ[j] == 0)
                {
                ZZ[j] = 1;
                *p++ = j;
                if(++k >= Vertices) return(TRUE);
                }
            }
        }                
    return(FALSE);        
}
                        
ComputeValences_A()
{
    /******************************************************************************************
        This routine is called by Find_Flow_A(). It computes the valences of the vertices of
        the Whitehead graph of the presentation given by the relators in Relators[]. The
        valence of vertex j is left in the array VA[] at VA[j]. The routine also checks the
        total length of the relators and returns TOO_LONG in case of an error.
    ******************************************************************************************/
        
    register unsigned int     i,
                            j;
    
    register unsigned long     CheckLength,
                            Valence;                        
    
    for(i = 0; i < Vertices; i++) A[i][i] = 0;
    for(i = 0,CheckLength = 0L; i < Vertices; i += 2)
        {
        Valence = 0L;
        for(j = 0; j < Vertices; j++) Valence += A[i][j];
        CheckLength += Valence;    
        j = (i >> 1);
        if(Valence > MAXLENGTH) return(TOO_LONG);
        VA[j] = Valence;
        }
    if(CheckLength != Length) return(TOO_LONG);
    return(0);    
}                

Do_Aut(Source,NumReps,NumRelators)
unsigned int     Source,
                NumReps;
int             NumRelators;
{        
    /******************************************************************************************
        This routine performs the automorphism(s) determined by the routines Find_Flow_A()
        Find_Primitives() and Level_Transformations(). The exact automorphism to be performed
        is specified by the entries of the global array ZZ[]. The routine performs a T-
        transformation corresponding to a splitting of the vertices of the Whitehead graph
        into two subsets: those vertices accessible from the "sink" and the remaining vertices.
        If 'A' is the source vertex and 'a' is the sink vertex. and X is a generator, (X != A),
        then the T-transformation acts on X according to the following table.
            1)    If 'X' not accessible and 'x' not accessible, then X: --> AXa.
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
    
    register unsigned char     A,
                            a,
                            *p,
                            *q,
                            *r,
                            x,
                            y;
                            
    register                 int i;
    
    unsigned char             TX[125],
                            TY[125],
                            **Temp;
                            
    unsigned int            j;
    
    unsigned long            HS;                        
            
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
    
    if(Micro_Print || Compute_Stabilizers) Micro_Print_Do_Aut(Source,NumReps);
    
    for(j = 0; j < NumReps; j++)
        {        
        for(i = 1; i <= NumRelators; i++)
            {
            HS = GetHandleSize((char **) Relators[i]);
            if(HS > MAXLENGTH) return(TOO_LONG);
            HS += HS;
            ReallocateHandle((char **) Temp5,HS);
            if((p = *Temp5) == NULL) return(TOO_LONG);
            q = *Relators[i];
            r = q;
            if(*r == EOS) continue;
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
            if(TX[x] && !TY[*r])
                *p++ = a;
            else
            if(!TX[x] && TY[*r])
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
    TotalAuts += NumReps;    
    return(0);                
}    

Do_Auts(Source, NumReps, NumRelators)
  unsigned int Source, NumReps,  NumRelators;
{        
    /******************************************************************************************
            Do_Auts() is used when the program has determined that it can perform the same
        T-transformation several times in succession. This can be done when the smallest
        absolute value with which the generator, corresponding to the "Source", appears in
        the relators is greater than one. Before calling Do_Auts(), the program called
        Get_MinExp() which determined the value of this minimal exponent. If the minimal
        exponent has value MinExp, then the T-transformation can be composed with itself
        MinExp times.
            Since using Do_Auts() requires calling Get_MinExp() and also requires calling
        Do_Aut() once in order to determine how the length of each individual relator is
        changed by a single application of Do_Aut(), it is more efficient to use Do_Aut()
        unless MinExp is greater than two. 
    ******************************************************************************************/
    
    register unsigned char     A,
                            a,
                            *p,
                            *q,
                            *r,
                            x,
                            y;
                            
    register long int        h;
    
    unsigned char             TX[125],
                            TY[125],
                            **Temp;
                            
    unsigned int            i;
    
    long                    Diff[MAXNUMRELATORS + 1],
                            HS,
                            sexp;
            
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
    
    if(Micro_Print || Compute_Stabilizers) Micro_Print_Do_Aut(Source,NumReps);
            
    for(i = 1; i <= NumRelators; i++)
        {
        h = 0;    
        q = *Relators[i];
        r = q;
        if(*r == EOS)
            continue;
        x = *q++;
        while(y = *q++)
            {
            if(x != A && x != a) h++;
            if(TX[x] && !TY[y])
                h++;
            else
            if(!TX[x] && TY[y])
                h++;
            x = y;
            }
        if(x != A && x != a) h++;
        if(TX[x] && !TY[*r])
            h++;
        else
        if(!TX[x] && TY[*r])
            h++;
        Diff[i] = h + 1 - GetHandleSize((char **) Relators[i]);                    
        }
    
    for(i = 1; i <= NumRelators; i++) if(**Relators[i])
        {
        HS = GetHandleSize((char **) Relators[i]) + NumReps*Diff[i];
        if(HS > MAXLENGTH) return(TOO_LONG);
        ReallocateHandle((char **) Temp5,HS);
        if((p = *Temp5) == NULL ) return(TOO_LONG);
        q = *Relators[i];
        r = q;
        x = *q;
        h = 0;
        if(x == A)
            {
            while(*q == x)
                {
                h++;
                q++;
                }
            if(*q == EOS) continue;
            r += GetHandleSize((char **) Relators[i]) - 2;
            while(*r == x)
                {
                h++;
                r--;
                }
            h -= NumReps;    
            r++;
            *r = EOS;
            r = *Relators[i];    
            if(TX[x] && !TY[*q])
                h -= NumReps;
            else
            if(!TX[x] && TY[*q])
                h += NumReps;    
            }    
        if(x == a)
            {
            while(*q == x)
                {
                h--;
                q++;
                }
            if(*q == EOS) continue;
            r += GetHandleSize((char **) Relators[i]) - 2;
            while(*r == x)
                {
                h--;
                r--;
                }
            h += NumReps;    
            r++;
            *r = EOS;
            r = *Relators[i];    
            if(TX[x] && !TY[*q])
                h -= NumReps;
            else
            if(!TX[x] && TY[*q])
                h += NumReps;
            }    
        sexp = h;
        h = 0;    
        x = *q++;
        while(y = *q++)
            {
            if(x != A && x != a)
                *p++ = x;
            if(TX[x] && !TY[y])
                h -= NumReps;
            else
            if(!TX[x] && TY[y])
                h += NumReps;    
            if(y != A && y != a)
                {
                while(h > 0)
                    {
                    *p++ = A;
                    h--;
                    }
                while(h < 0)
                    {
                    *p++ = a;
                    h++;
                    }    
                }                
            x = y;
            if(y == A)
                {
                while(*q == A)
                    {
                    q++;
                    h++;
                    }
                h -= NumReps - 1;    
                }
            if(y == a)
                {
                while(*q == a)
                    {
                    q++;
                    h--;
                    }
                h += NumReps - 1;    
                }            
            }
        if(x != A && x != a)
            *p++ = x;
        if(TX[x] && !TY[*r])
            h -= NumReps;
        else
        if(!TX[x] && TY[*r])
            h += NumReps;    
        h += sexp;
        while(h > 0)
            {
            *p++ = A;
            h--;
            }
        while(h < 0)
            {
            *p++ = a;
            h++;
            }
        *p = EOS;
        q = *Temp5;
        HS = p + 1 - q;
        if(HS > MAXLENGTH) return(TOO_LONG);
        SetHandleSize((char **) Temp5,HS);
        Temp = Relators[i];
        Relators[i] = Temp5;
        Temp5 = Temp;                    
        }
    TotalAuts += NumReps;    
    return(0);                
}    

Freely_Reduce()
{
    /******************************************************************************************
                        This routine freely reduces the input relators.
    ******************************************************************************************/
    
    register unsigned char     *p,
                            *q;
    
    register char            x;
                            
    register int             i;
            
    int                        NumEmptyRelators = 0;
    
    unsigned long            HS,
                            length;
    
    if(Micro_Print) for(i = 1,length = 0L; i <= NumRelators; i++)
        length += GetHandleSize((char **) Relators[i]) - 1;
            
    for(i = 1,OrigLength = 0L; i <= NumRelators; i++)
        {
        HS = GetHandleSize((char **) Relators[i]);        
        ReallocateHandle((char **) Temp5,HS + 2);
        if((p = *Temp5) == NULL) return(TOO_LONG);
        q = *Relators[i];
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
        q = *Temp5;
        q++;
        if(p > q) while((*p - *q) == 32 || (*p - *q) == -32)
            {
            p--;
            q++;
            }
        p++;
        *p = EOS;
        if(p == q) NumEmptyRelators ++;                            
        OrigLength += p - q;
        if(HS > p + 1 - q)
            {
            SetHandleSize((char **) Relators[i],p + 1 - q);
            if((p = *Relators[i]) == NULL) return(TOO_LONG);
            while(*p++ = *q++) ;
            }                    
        }
    if(Micro_Print && OrigLength < length) Micro_Print_Freely_Reduce(length,OrigLength);    
    return(NumEmptyRelators);            
}

int CheckPrimitivity(void)
{
    /******************************************************************************************
        This routine checks whether Relators[1] is a primitive or a proper power of a generator.
        It returns 1 if Relators[1] is a primitive, returns -1 if Relators[1] is a proper
        power of a generator and otherwise returns 0.
    ******************************************************************************************/            

    register unsigned char     *p,
                            x;
    
    if(GetHandleSize((char **) Relators[1]) <= 1) return(0);        
    switch(Find_Primitives(2))
        {
        case 0:
            return(0);
        case 1:
            return(1);
        case 2:
            p = *Relators[1];
            x = *p;
            while(x == *p) p++;
            if(!*p) return(-1);
            return(0);
        case TOO_LONG:
            return(TOO_LONG);    
        }
    return(0);                                                                        
}

Test_New_Pres()
{
    /******************************************************************************************
        Test_New_Pres() is called when RESET has chosen a new presentation to investigate.
        That is: a presentation for which UDV[ReadPres] = 0. This routine then runs certain
        tests on the new presentation as detailed below.
    ******************************************************************************************/
        
    register unsigned char     *p,
                            *q;
                            
    register unsigned int     i,
                            j;
    
    unsigned char            **Temp;
    
    int                        DistinctNonEmpty,
                            k,
                            SRBoundary,
                            SRDelete_Only_Short_Primitives,
                            SRNumGenerators,
                            SRNumRelators,
                            SRReadPres;
    
    unsigned int            SNumFilled,
                            RVWG;
                            
    unsigned long            HS,
                            SRLength;                        

    unsigned int Whitehead_Graph();
    unsigned int Reduce_Genus();
    
    NumDiagrams ++;
    Input = NORMAL;
    UDV[ReadPres] = 1;

    if(Get_Relators_From_SUR(ReadPres)) return(RESET);
    
    if(NumGenerators == 1)
        {
        /**************************************************************************************
            If this presentation is a presentation on only one generator, then we are
            essentially done.
        **************************************************************************************/    
        
        TP[ReadPres] = FALSE;
        DistinctNonEmpty = Delete_Dups();
        if(DistinctNonEmpty > 1)
            {
            UDV[ReadPres] = DONE;
            return(STOP);    
            }
        HS = GetHandleSize((char **) Relators[1]) - 1;    
        switch(HS)
            {
            case 0:
                if(BDY[ReadPres] == TRUE)
                    {
                    if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                        {
                        CBC[CurrentComp][0] = EOS;
                        CBC[CurrentComp][1] = 1;
                        CBC[CurrentComp][2] = BDRY_UNKNOWN;
                        }
                    UDV[NumFilled] = S1_X_D2;
                    Mark_As_Found_Elsewhere(CurrentComp);
                    }
                if(BDY[ReadPres] == FALSE)
                    {
                    if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                        {
                        CBC[CurrentComp][0] = NumRelators - NumEmptyRels;
                        CBC[CurrentComp][1] = BDRY_UNKNOWN;
                        }
                    UDV[NumFilled] = S1_X_S2;
                    Mark_As_Found_Elsewhere(CurrentComp);
                    }
                if(BDY[ReadPres] > 1)
                    UDV[NumFilled] = S1_X_X2;        
                break;
            case 1:
                if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                    {
                    CBC[CurrentComp][0] = NumRelators - NumEmptyRels;
                    CBC[CurrentComp][1] = BDRY_UNKNOWN;
                    }
                UDV[ReadPres] = THREE_SPHERE;
                Mark_As_Found_Elsewhere(CurrentComp);    
                break;
            case 2:
            case 3:
            case 4:
                if(CBC[ComponentNum[ReadPres]][0] == BDRY_UNKNOWN)
                    {
                    CBC[ComponentNum[ReadPres]][0] = NumRelators - NumEmptyRels;
                    CBC[ComponentNum[ReadPres]][1] = BDRY_UNKNOWN;
                    }
                LSP[ReadPres] = HS;
                LSQ[ReadPres] = 1;    
                UDV[ReadPres] = KNOWN_LENS_SPACE;
                BDY[ReadPres] = FALSE;
                Mark_As_Found_Elsewhere(CurrentComp);
                break;    
            default:
                if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                    {
                    CBC[CurrentComp][0] = NumRelators - NumEmptyRels;
                    CBC[CurrentComp][1] = BDRY_UNKNOWN;
                    }
                LSP[ReadPres] = HS;
                LSQ[ReadPres] = 1;                    
                UDV[ReadPres] = GENERIC_LENS_SPACE;
                break;
            }
        if(Micro_Print)
            {
            printf("\n\nThe current presentation presents a manifold of Heegaard genus one.");
            if(Micro_Print_F)
                fprintf(myout,"\n\nThe current presentation presents a manifold of Heegaard genus one.");
            }                
        return(RESET);            
        }

    if(Length == NumGenerators && NumRelators == NumGenerators)
        {
        /**************************************************************************************
                If Length == NumGenerators and NumRelators == NumGenerators then we check
                whether this is a canonical presentation of the 3-Sphere.
        **************************************************************************************/
            
        for(i = 1; i <= NumRelators; i++) if(GetHandleSize((char **) Relators[i]) != 2L) break;
        if(i > NumRelators && Delete_Dups() == NumRelators)
            {
            UDV[ReadPres] = THREE_SPHERE;
            if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                {
                CBC[CurrentComp][0] = 1;
                CBC[CurrentComp][1] = BDRY_UNKNOWN;
                }
            if(Micro_Print)
                {
                printf("\n\nThe current presentation presents the 3-Sphere.");
                if(Micro_Print_F)
                    fprintf(myout,"\n\nThe current presentation presents the 3-Sphere.");
                }    
            Mark_As_Found_Elsewhere(CurrentComp);
            return(RESET);
            }
        return(RESET);        
        }
        
    From_BANDSUM = 0;
    From_DUALIZE = 0;
    
    /******************************************************************************************
                Just for good measure, call Freely_Reduce(). Then call Find_Flow_A().
    ******************************************************************************************/
            
    if(Freely_Reduce() == TOO_LONG) return(RESET);
    Length = OrigLength;
    switch(Find_Flow_A(NORMAL,FALSE))
        {
        case 1:
            switch(Missing_Gen())
                {
                case TOO_LONG:
                    return(RESET);
                case TOO_MANY_COMPONENTS:
                    return(TOO_MANY_COMPONENTS);
                case NO_ERROR:
                    return(RESET);
                }    
        case TOO_LONG:    
            return(RESET);
        }
        
    if(SetUp_TopOfChain()) return(RESET);    

    Saved_Vertices = 0;
    
    /******************************************************************************************
        Set Saved_Vertices = 0. (This acts as a flag for Whitehead_Graph() and causes
        WhiteHead_Graph() to update the adjacency matrices it uses.)
        Call Whitehead_Graph() to look for the Heegaard diagram corresponding to the
        presentation in Relators[].
    ******************************************************************************************/
        
    switch(RVWG = Whitehead_Graph())
        {
        case NO_ERROR:
            if((CBC[CurrentComp][0] == BDRY_UNKNOWN || ER[ReadPres] < 0) 
                && NumGenerators > 1)
                {
                Get_Bdry_Comps(FALSE,FALSE,ReadPres);
                if(NumBdryComps == BCWG[0])
                    {
                    BDY[ReadPres] = FALSE;
                    Boundary = FALSE;        
                    }
                else
                    {
                    BDY[ReadPres] = TRUE;
                    Boundary = TRUE;
                    }
                for(i = 0; (CBC[CurrentComp][i] = BCWG[i]) < BDRY_UNKNOWN; i++) ;
                if(CS[CurrentComp + 1] == 3) MG_Bdry_Comp_Data(ReadPres);        
                if((BCWG[0] > 1 || (BCWG[0] && NumBdryComps > BCWG[0]))
                    && Delete_Redundant_Relators()) return(RESET);
                ER[ReadPres] = 0;        
                }
            else
                Boundary = BDY[ReadPres];
            if(NumRelators == 1)
                {
                /******************************************************************************
                    If there is only one relator, once we have reduced it and found its
                    Heegaard diagram, there is nothing more to do. So flag this presentation
                    as done and return RESET.
                ******************************************************************************/
                
                UDV[ReadPres] = DONE;
                return(RESET);
                }

            if(Automorphisms)
                {
                if(Boundary || NumGenerators != NumRelators) 
                    {
                    GoingUp = FALSE;
                    return(BANDSUM);
                    }
                return(DUALIZE);
                }    

            /**********************************************************************************
                    Save a copy of the current relators so we can restore them later.
            **********************************************************************************/
            
            SRLength = Length;
            SRNumRelators = NumRelators;
            SRNumGenerators = NumGenerators;
            SRReadPres = ReadPres;
            for(i = 1; i <= NumRelators; i++)
                {
                ReallocateHandle((char **) Copy_Of_Rel_2[i],GetHandleSize((char **) Relators[i]));
                if((p = *Copy_Of_Rel_2[i]) == NULL) return(RESET);
                q = *Relators[i];
                while(*p++ = *q++) ;                    
                }                
            
            SRBoundary = Boundary;
            SRDelete_Only_Short_Primitives = Delete_Only_Short_Primitives;
                
            if(TP[ReadPres] && NumRelators > 1 && !Do_Not_Reduce_Genus)
                {
                /******************************************************************************
                    If TP[ReadPres] is TRUE, then we want to check at this point whether:
                        1) There are primitives of length one or two which can be deleted.
                        2) There is a primitive relator or proper power relator whose deletion
                            creates some empty relators.
                        3) This is a presentation of a lens space.
                ******************************************************************************/    
                
                Delete_Only_Short_Primitives = TRUE;
                switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                    {
                    case NO_ERROR:
                        break;
                    case FATAL_ERROR:
                        Delete_Only_Short_Primitives = SRDelete_Only_Short_Primitives;
                        return(STOP);
                    case TOO_LONG:
                    case CAN_NOT_DELETE:
                        FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                        ReadPres = SRReadPres;
                        NumGenerators = SRNumGenerators;
                        NumRelators = SRNumRelators;
                        Vertices = 2*NumGenerators;
                        Length = SRLength;
                        Boundary = SRBoundary;                    
                        for(i = 1; i <= NumRelators; i++)
                            {
                            ReallocateHandle((char **) Relators[i],GetHandleSize((char **) Copy_Of_Rel_2[i]));
                            if((q = *Relators[i]) == NULL) return(RESET);        
                            p = *Copy_Of_Rel_2[i];
                            while(*q++ = *p++) ;                                    
                            }                        
                        if(Micro_Print)
                            {
                            printf("\n\nUnable to delete a primitive. Reverting to Presentation %d: Length %lu",
                                ReadPres + 1,Length);
                            if(Micro_Print_F)
                                fprintf(myout,"\n\nUnable to delete a primitive. Reverting to Presentation %d: Length %lu",
                                    ReadPres + 1,Length);
                            }
                    }
                Delete_Only_Short_Primitives = SRDelete_Only_Short_Primitives;                    
                if(FoundPrimitive || FoundPower || LensSpace || EmtyRel) return(REDUCE_GENUS);        
                }


            if(TP[ReadPres] && NumRelators > 1 && !Do_Not_Reduce_Genus
                && !Delete_Only_Short_Primitives)
                {
                /******************************************************************************
                    If TP[ReadPres] is TRUE, we make another call to Reduce_Genus(). This time
                    checking whether:
                        1) There are any primitives or proper powers which can be deleted.
                        2) There is a primitive relator or proper power relator whose deletion
                            creates some empty relators.
                        3) This is a presentation of a lens space.
                ******************************************************************************/    
                
                switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                    {
                    case NO_ERROR:
                        break;
                    case FATAL_ERROR:
                        return(STOP);
                    case TOO_LONG:
                    case CAN_NOT_DELETE:
                        FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                        TP[ReadPres] --;            
                        if(Micro_Print)
                            {
                            printf("\n\nUnable to delete a primitive. Reverting to Presentation %d: Length %lu.",
                                ReadPres + 1,SRLength);
                            if(Micro_Print_F)
                                fprintf(myout,"\n\nUnable to delete a primitive. Reverting to Presentation %d: Length %lu.",
                                    ReadPres + 1,SRLength);
                            }
                    }
                    
                if(LensSpace || EmtyRel)
                    {
                    FoundPrimitive = FoundPower = FALSE;
                    return(REDUCE_GENUS);
                    }                        
                if(FoundPrimitive || FoundPower)
                    {
                    if(SRNumGenerators == SRNumRelators && BDY[SRReadPres] == FALSE
                        && PRIM[SRReadPres] == 108)
                        {
                        /**********************************************************************
                                The initial presentation was part of a pseudo-minimal pair.
                        **********************************************************************/    
                        
                        return(REDUCE_GENUS);
                        }
                    if(SRBoundary && QPM[SRReadPres])
                        {
                        /**********************************************************************
                                The initial presentation was "quasi-pseudo-minimal".
                        **********************************************************************/
                        
                        return(REDUCE_GENUS);
                        }        
                    }        
                }

            /**********************************************************************************
                            Restore the initial presentation we saved.
            **********************************************************************************/

            ReadPres = SRReadPres;
            NumGenerators = SRNumGenerators;
            NumRelators = SRNumRelators;
            Vertices = 2*NumGenerators;
            Length = SRLength;
            Boundary = SRBoundary;                    
            for(i = 1; i <= NumRelators; i++)
                {
                Temp = Relators[i];
                Relators[i] = Copy_Of_Rel_2[i];
                Copy_Of_Rel_2[i] = Temp;                                    
                }            

            if(Micro_Print && (FoundPrimitive || FoundPower || LensSpace || EmtyRel) )
                {
                printf("\n\nReverting to Presentation %d: Length %lu.",
                    ReadPres + 1,Length);
                if(Micro_Print_F)
                    fprintf(myout,"\n\nReverting to Presentation %d: Length %lu.",
                        ReadPres + 1,Length);
                }
        
            FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
            Fill_A(NumRelators);
            if(Whitehead_Graph()) return(RESET);
            if(NumGenerators == NumRelators && BDY[ReadPres] == FALSE)
                {
                if(PRIM[ReadPres] != 108) return(DUALIZE);
                GoingUp = TRUE;
                return(BANDSUM);    
                }
            GoingUp = FALSE;            
            return(BANDSUM);             
        case NON_UNIQUE_1:
        case NON_UNIQUE_2:
        case NON_UNIQUE_3:
        case NON_UNIQUE_4:
            UDV[ReadPres] = RVWG;
        case V2_ANNULUS_EXISTS:                
            if(TP[ReadPres] && NumRelators > 1 && !Do_Not_Reduce_Genus)
                {
                switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                    {
                    case NO_ERROR:
                        break;
                    case FATAL_ERROR:
                        return(STOP);
                    case TOO_LONG:
                    case CAN_NOT_DELETE:
                        if(Micro_Print)
                            {
                            printf("\n\nUnable to delete a primitive.");
                            if(Micro_Print_F)
                                fprintf(myout,"\n\nUnable to delete a primitive.");
                            }
                        TP[ReadPres] --;                        
                        FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
                        return(RESET);    
                    }    
                if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
                    return(REDUCE_GENUS);
                TP[ReadPres] = FALSE;    
                }        
            return(RESET);            
        case NOT_CONNECTED:
            if(UDV[ReadPres] < DONE) UDV[ReadPres] = NOT_CONNECTED;
            return(RESET);
        case REDUCE_GENUS:
            if(UDV[ReadPres] < DONE) UDV[ReadPres] = NOT_CONNECTED;
            if(TestRealizability1 || TestRealizability2)
                Delete_Trivial_Generators(FALSE);
            FoundPrimitive = TRUE;
            return(REDUCE_GENUS);            
        case SEP_PAIRS:
        
            /**********************************************************************************
                The diagram corresponding to Relators[] has a pair of separating vertices.
                Record the two vertices that separate in LSP and LSQ.
                    Then, if this presentation hasn't been checked for primitives etc. call
                Reduce_Genus(). If Reduce_Genus turns up something return REDUCE_GENUS. 
                Otherwise, call Level_Transformations to see if it can find a modified
                presentation that doesn't have any separating vertices. Return RESET.
            **********************************************************************************/
            
            UDV[ReadPres] = SEP_PAIRS;
            if(V1 & 1)
                LSP[ReadPres] = V1/2 + 97;
            else
                LSP[ReadPres] = V1/2 + 65;
            if(V2 & 1)
                LSQ[ReadPres] = V2/2 + 97;
            else
                LSQ[ReadPres] = V2/2 + 65;

            NumCalled = 0;
            NotNewPres = 0;
            SNumFilled = NumFilled;
            k = Level_Transformations(TRUE,!Find_All_Min_Pres,FALSE,0,0);                            
            for(i = 0; i < NumCalled; i++)
            for(j = 1; j <= NumRelators; j++) DisposeHandle((char **) SLR[i][j]);
            if(Level_Interrupt == 1) return(RESET);
            if(k == FATAL_ERROR)
                {
                Fatal_Error();
                return(STOP);
                }        
            if(Micro_Print && k != 2)
                {
                printf("\n\nThe current Presentation is:\n");
                Print_Relators(Relators,NumRelators,stdout);
                if(Micro_Print_F)
                    {
                    fprintf(myout,"\n\nThe current Presentation is:\n");
                    Print_Relators(Relators,NumRelators,myout);
                    }
                }            
            if(k == 5 && !Do_Not_Reduce_Genus) switch(Delete_Trivial_Generators(FALSE))
                {
                case 0:
                    break;
                case 1:
                    FoundPrimitive = TRUE;                                            
                    return(REDUCE_GENUS);
                case TOO_LONG:
                    return(RESET);    
                }
            if(Find_All_Min_Pres && SNumFilled == NumFilled)
                {
                if(Get_Relators_From_SUR(ReadPres)) return(RESET);
                if(Micro_Print) Micro_Print_Reset();                    
                if(Freely_Reduce() == TOO_LONG) return(RESET);
                Length = OrigLength;
                if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG) return(RESET);
                if(!Init_Find_Level_Transformations(FALSE))
                    {
                    switch(Find_Level_Transformations(Delete_Only_Short_Primitives,0))
                        {
                        case 0:
                        case 1:
                        case 2:
                            break;
                        case 3:
                            Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                            ReadPres = NumFilled - 1;
                            switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                                {
                                case NO_ERROR:
                                    break;
                                case FATAL_ERROR:
                                    Fatal_Error();
                                    return(STOP);
                                case TOO_LONG:
                                case CAN_NOT_DELETE:
                                    if(Micro_Print)
                                        {
                                        printf("\n\nUnable to delete a primitive.");
                                        if(Micro_Print_F)
                                            fprintf(myout,"\n\nUnable to delete a primitive.");
                                        }                                    
                                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;                            
                                    return(RESET);    
                                }            
                            return(REDUCE_GENUS);
                        case 5:
                            Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                            return(RESET);
                        default:
                            break;    
                        }
                    Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                    }
                }    
                
            if(TP[ReadPres] && SNumFilled == NumFilled && NumRelators > 1
                && !Do_Not_Reduce_Genus)
                {
                if(Get_Relators_From_SUR(ReadPres)) return(RESET);    
                if(Micro_Print) Micro_Print_Reset();                    
                if(Freely_Reduce() == TOO_LONG) return(RESET);
                Length = OrigLength;
                if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG) return(RESET);
                switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                    {
                    case NO_ERROR:
                        break;
                    case FATAL_ERROR:
                        return(STOP);
                    case TOO_LONG:
                    case CAN_NOT_DELETE:
                        TP[ReadPres] --;
                        if(Micro_Print)
                            {
                            printf("\n\nUnable to delete a primitive.");
                            if(Micro_Print_F)
                                fprintf(myout,"\n\nUnable to delete a primitive.");
                            }                        
                        FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;                        
                        return(RESET);    
                    }
                if(FoundPrimitive || FoundPower || LensSpace || EmtyRel)
                    return(REDUCE_GENUS);        
                }
                                                            
            return(RESET);        
        case TOO_MANY_COMPONENTS:
            return(TOO_MANY_COMPONENTS);
        case NON_PLANAR:
            fprintf(stdout,"\n\n                    The Whitehead graph is nonplanar.");
            fprintf(myout,"\n\n                    The Whitehead graph is nonplanar.");    
        case FATAL_ERROR:
            return(STOP);        
        case TOO_LONG:
            if(UDV[ReadPres] < DONE) UDV[ReadPres] = DONE;
            return(RESET);        
        }
    return(NO_ERROR);                        
}

Check_Level_Transformations()
{
    /******************************************************************************************
        This routine looks for new presentations obtainable by performing level
        T-transformations. In order to reduce the proliferation of presentations, the routine
        is currently configured to look for level transformations only for presentations
        which have a minimal number of generators among the presentations of its summand and,
        among these presentations, have shortest length.If new presentations, which are not on
        file, are found these are saved and we return TRUE. Otherwise we return FALSE.
    ******************************************************************************************/
        
    unsigned int            i,
                            SNumFilled;

    for(ReadPres = Start_Level_Search; ReadPres < NumFilled; ReadPres++)
        {
        /**************************************************************************************
            This "for" loop filters out presentations for which we do not want to find level-
            transformations.
        **************************************************************************************/    
        if(NG[ReadPres] < 2) continue;
        if(SURL[ReadPres] == 0L) continue;
        for(i = 1; i < NG[ReadPres]; i++) if(MLC[ComponentNum[ReadPres]][i] < BIG_NUMBER)
            break;
        if(i < NG[ReadPres]) continue;    
        if(MLC[ComponentNum[ReadPres]][NG[ReadPres]] < SURL[ReadPres]) continue;
        switch(UDV[ReadPres])
            {
            case SPLIT:
            case GENERIC_LENS_SPACE:
            case KNOWN_LENS_SPACE:
            case S1_X_S2:
            case S1_X_D2:
            case MISSING_GEN_DONE1:
            case MISSING_GEN_DONE2:
            case THREE_SPHERE:
            case FOUND_ELSEWHERE:
                break;
            default:
                Start_Level_Search    = ReadPres + 1;
                NumGenerators         = NG[ReadPres];
                NumRelators         = NR[ReadPres];
                Vertices             = 2*NumGenerators;
                CurrentComp            = ComponentNum[ReadPres];    
                if(Init_Find_Level_Transformations(FALSE) == FALSE)
                    {
                    SNumFilled = NumFilled;
                    if(Micro_Print) Micro_Print_Reset();
                    switch(Find_Level_Transformations(Delete_Only_Short_Primitives,0))
                        {
                        case 0:
                        case 1:
                        case 2:
                            break;
                        case 3:
                            ReadPres = NumFilled - 1;
                            Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                            if(Micro_Print) Micro_Print_Reset();
                            switch(Reduce_Genus(NORMAL,FALSE,FALSE))
                                {
                                case NO_ERROR:
                                    break;
                                case FATAL_ERROR:
                                    Fatal_Error();
                                    return(FALSE);
                                case TOO_LONG:
                                case CAN_NOT_DELETE:
                                    if(Micro_Print)
                                        {
                                        printf("\n\nUnable to delete a primitive.");
                                        if(Micro_Print_F)
                                            fprintf(myout,"\n\nUnable to delete a primitive.");
                                        }                                    
                                    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;                                
                                    return(FALSE);    
                                }            
                            return(2);        
                        case 5:
                            Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                            return(TRUE);
                        default:
                            break;    
                        }        
                    Free_Memory_For_Find_Level_Transformations(FALSE,1000);
                    if(NumFilled > SNumFilled) return(TRUE);
                    }
                break;
            }
        }
    return(FALSE);        
}

unsigned int Get_MinExp(unsigned int Source,int MyNumRelators)
{
    register unsigned char     *p,
                            x,
                            y,
                            z;
                            
    register unsigned int    j,
                            MinExp;                        
                            
    unsigned int             h,
                            i;
                            
    /******************************************************************************************
        Let G be the generator which corresponds to the vertex "Source" of the Whitehead graph.
        This routine finds and returns the smallest absolute value of the exponents with which
        G appears in the relators.
    ******************************************************************************************/    
    
    x = Source/2 + 65;
    y = x + 32;
    
    for(i = 1,MinExp = INFINITE; i <= MyNumRelators; i++)
        {
        p = *Relators[i];
        h = 0;
        z = *p;
        if(z == x || z == y) while(*p == z)
            {
            p++;
            h++;
            }
        if(*p == EOS && (h > 0) && (h < MinExp))
            {
            MinExp = h;
            if(MinExp < 2) return(1);    
            continue;
            }        
        while(z = *p)
            {
            p++;
            if(z != x && z != y)
                {
                if(*p == EOS && (h > 0) && (h < MinExp))
                    {
                    MinExp = h;
                    if(MinExp < 2) return(1);
                    }            
                if(*p == x || *p == y)
                    {
                    for(j = 0, z = *p; z == *p; j++, p++) ;
                    if(*p)
                        {
                        if(j < MinExp)
                            {
                            MinExp = j;
                            if(MinExp < 2) return(1);
                            }
                        }    
                    else
                        {
                        if(j + h < MinExp)
                            {
                            MinExp = j + h;
                            if(MinExp < 2) return(1);
                            }
                        }        
                    }
                }    
            }
        }
    return(MinExp);    
}

void Mark_As_Duplicate(Dup_On_File)
     unsigned int Dup_On_File;
{
    register int        i,
                        j;
    
    UDV[NumFilled - 1] = DUPLICATE;
    if(UDV[Dup_On_File] > UNKNOWN)
        {
        for(i = 0; i < NumFilled; i++) if(ComponentNum[i] == CurrentComp && UDV[i] <= DONE)
            UDV[i] = FOUND_ELSEWHERE;
        CS[CurrentComp] = 2;
        }
    j = ComponentNum[Dup_On_File];    
    if(CBC[CurrentComp][0] == BDRY_UNKNOWN && CBC[j][0] != BDRY_UNKNOWN)
        for(i = 0; (CBC[CurrentComp][i] = CBC[j][i]) < BDRY_UNKNOWN; i++) ;
    if(CBC[CurrentComp][0] != BDRY_UNKNOWN && CBC[j][0] == BDRY_UNKNOWN)
        for(i = 0; (CBC[j][i] = CBC[CurrentComp][i]) < BDRY_UNKNOWN; i++) ;    
    if(CS[CurrentComp + 1] == 3) MG_Bdry_Comp_Data(NumFilled - 1);        
}

void Mark_As_Found_Elsewhere(TheComp)
{
    register int        i;
    
    for(i = 0; i < NumFilled; i++)
    if(ComponentNum[i] == TheComp && UDV[i] <= DONE)
        UDV[i] = FOUND_ELSEWHERE;
    CS[TheComp] = 2;
    if(CS[TheComp + 1] == 3)
        {
        for(i = 0; i < NumFilled; i++) if(ComponentNum[i] == TheComp) break;
        if(i < NumFilled) MG_Bdry_Comp_Data(i);
        }    
}

int SetUp_TopOfChain(void)
{
    register unsigned char    *p,
                            *q;
                            
    register int            i;                            
    
    TOCLength     = Length;
    NG_TOC         = NumGenerators;
    NR_TOC         = NumRelators;
    for(i = 1; i <= NumRelators; i++)
        { 
        ReallocateHandle((char **) TopOfChain[i],GetHandleSize((char **) Relators[i]));
        if((q = *TopOfChain[i]) == NULL) return(TRUE);            
        p = *Relators[i];
        while(*q++ = *p++) ;    
        }
    if(Micro_Print)
        {
        printf("\n\nSaved a copy of the current Relators[] as TopOfChain[].");
        if(Micro_Print_F)
            fprintf(myout,"\n\nSaved a copy of the current Relators[] as TopOfChain[].");
        }    
    return(FALSE);
}

Get_Relators_From_SUR(MyReadPres)
register int    MyReadPres;
{
    register unsigned char    *p,
                            *q;
                            
    register int            i;
                            
    NumRelators     = NR[MyReadPres];
    NumGenerators     = NG[MyReadPres];
    Length             = SURL[MyReadPres];
    Vertices         = 2*NumGenerators;
    
    for(i = 1; i <= NumRelators; i++)
        {
        ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[MyReadPres][i]));
        if((q = *Relators[i]) == NULL) return(TRUE);        
        p = *SUR[MyReadPres][i];
        while(*q++ = *p++) ;        
        }
    return(FALSE);
}
