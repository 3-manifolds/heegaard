#include "Heegaard.h"
#include "Heegaard_Dec.h"

/*  #define MANUALLY_RENUMBER_PRESENTATIONS    */
/*  #define PRINT_CYCLES    */
/*  #define PRINT_TRIAL_CYCLES    */
                            
unsigned int Whitehead_Graph()
{
    int     i;
    
    /******************************************************************************************
            Call Get_Matrix() to determine the valences of the vertices of the "reduced"
            Whitehead graph, and then delete extraneous components if the "reduced"
            Whitehead graph is not connected.
    ******************************************************************************************/    
    
    i = Get_Matrix();
    if(i == FALSE || (i == TRUE && !Connected) || (i == TRUE && MajorVert == 2))    
        {
        Saved_Vertices = Vertices;
        if(Check_Connected() == FALSE)
            {
            if(NotConnectedError == TOO_LONG)
                return(TOO_LONG);
            if(NotConnectedError == TOO_MANY_COMPONENTS)
                return(TOO_MANY_COMPONENTS);
            if(NotConnectedError)    
                return(REDUCE_GENUS);
            else
                return(NOT_CONNECTED);
            }
            
        /**************************************************************************************
            Next, call Sep_Pairs() to find out if there are any pairs of vertices whose removal
            gives a non-trivial separation of the "reduced" Whitehead graph.
        **************************************************************************************/
        
        SepPairs = Sep_Pairs(0,0);
        switch(SepPairs)
            {
            case 0:
                break;
            case 1:
                return(SEP_PAIRS);
            case TOO_LONG:
                return(TOO_LONG);
            }
                    
        /**************************************************************************************
            We now have a graph which is 3-connected. Call Planar(TRUE,FALSE) to determine
            whether this graph is planar.
        **************************************************************************************/    
                         
        if(NonPlanar = Planar(TRUE,FALSE)) return(NON_PLANAR);
        }
        else
        {
        if(SepPairs)     return(SEP_PAIRS);
        if(NonPlanar)     return(NON_PLANAR);
        }                     
    return(Diagram_Main());                                                            
}    

Get_Matrix()
{
    register int     i,
                    j,
                    k;
                    
    unsigned int     *Temp;
    
    /****************************************************************************************** 
            Set VWG[i] equal to the valence of vertex i in the "reduced" Whitehead graph. 
            Set NumEdges equal to the number of edges in the "reduced" Whitehead graph.
    ******************************************************************************************/
    
    for(i = NumEdges = 0;i < Vertices; i ++)
        {
        A[i][i] = 0;
        for(j = k = 0; j < Vertices; j ++) if(A[i][j])
            {
            AJ2[i][k] = j;
            k++;
            }
        VWG[i] = k;
        NumEdges += k;
        AJ2[i][k] = VERTICES;
        }                            
    NumEdges /= 2;
    
    /******************************************************************************************
        Compare the adjacency matrices AJ1 and AJ2. If they are identical and Saved_Vertices
        equals Vertices return TRUE.
        Otherwise swap the arrays AJ1 and AJ2 and return FALSE.
    ******************************************************************************************/
    
    for(i = 0; i < Vertices; i++)
        {
        for(j = 0; ; j++)
            {
            if(AJ1[i][j] != AJ2[i][j])
                {
                for(i = 0; i < Vertices; i ++)
                    {
                    Temp = AJ1[i];
                    AJ1[i] = AJ2[i];
                    AJ2[i] = Temp;
                    }
                return(FALSE);
                }
            if(AJ2[i][j] == VERTICES) break;
            }
        }
    if(Saved_Vertices != Vertices) return(FALSE);            
    return(TRUE);        
}    
        
Check_Connected()
{    
    register unsigned char     *p,
                            *q;
                            
    register int             i,
                            j;
                            
    unsigned char             **Temp,
                            x;                        
    
    int                        h,
                            SaveCS,
                            SaveNumRelators;
    
    unsigned int            SaveUDV;
        
    /****************************************************************************************** 
                    Check whether the "reduced" Whitehead graph is connected. 
    ******************************************************************************************/    
    
    for(i = 0; i < Vertices; i++) ZZ[i] = 0;
    Connected = Connected_(0,0);    
    if(DrawingDiagrams == TRUE) return(TRUE);
    
    /****************************************************************************************** 
        If Connected is FALSE, then the "reduced" Whitehead graph is not connected, hence the 
        Heegaard diagram is reducible. So we call Split_Relators() to "split" the relators into
        two subsets corresponding to this splitting of the Heegaard diagram. 
    ******************************************************************************************/    
        
    if(Connected == FALSE)
        {
        NotConnectedError = FALSE;    

        /**************************************************************************************
                If there are more than two generators and more than one relator, call
            Delete_Trivial_Generators(TRUE) to see if there are trivial generators which can
            be removed. If we are testing realizability, then we just check whether there are
            any trivial generators which can be deleted.
        **************************************************************************************/
        
        if(NumGenerators > 2 && NumRelators > 1 && !Do_Not_Reduce_Genus)
            {
            if(TestRealizability1 || TestRealizability2)
                h = Delete_Trivial_Generators(TRUE);
            else
                h = Delete_Trivial_Generators(FALSE);    
            switch(h)
                {
                case 0:
                    break;
                case 1:
                    NotConnectedError = TRUE;                                    
                    return(FALSE);
                case TOO_LONG:
                    NotConnectedError = TOO_LONG;
                    return(FALSE);    
                }
            }
            
        /**************************************************************************************
            Check whether the presentation that is "splitting" is already on file, if it is,
            we don't want to save another copy of it.
        **************************************************************************************/    
        
        SaveUDV = UDV[ReadPres];
        Canonical_Rewrite(Relators,FALSE,FALSE);        
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
                 else if(i != ReadPres)
                     {
                     NotConnectedError = FALSE;
                     return(FALSE);
                     }    
                 break;        
                 }    
             }
        if(i == NumFilled && Dup_On_File < INFINITE)
            {
             if((NumFilled >= MAX_SAVED_PRES - 3) || 
                 Save_Pres(ReadPres,Dup_On_File,Length,1,20,1,0,0))
                 {
                 NotConnectedError = TOO_LONG;
                 return(FALSE);
                 }        
             Mark_As_Duplicate(Dup_On_File);
             NotConnectedError = FALSE;        
            return(FALSE);                        
             }
 
         if(i < NumFilled) switch(SaveUDV)
            {
            case SPLIT:
            case ANNULUS_EXISTS:
            case V2_ANNULUS_EXISTS:
                NotConnectedError = FALSE;
                return(FALSE);
            }
                                         
        /**************************************************************************************
                Since we have rewritten the relators, we need to update the relevant arrays.
        **************************************************************************************/    
        
        Fill_A(NumRelators);
        Get_Matrix();
        for(j = 0; j < Vertices; j++) ZZ[j] = 0;
        Connected_(0,0);

        if(Micro_Print)
            {
            printf("\n\nThe diagram of the following presentation from Presentation %d is not connected.\n",
                ReadPres + 1);
            Print_Relators(Relators,NumRelators,stdout);
            if(Micro_Print_F)
                {
                fprintf(myout,"\n\nThe diagram of the following presentation from Presentation %d is not connected.\n",
                    ReadPres + 1);
                Print_Relators(Relators,NumRelators,myout);
                }
            }
            
        if(i == NumFilled)
            {
            /**********************************************************************************
                The presentation that is "splitting" is new, so we want to save a copy.
            **********************************************************************************/
            
            /**********************************************************************************
                If Whitehead_Graph() was called by Lens_Space(), we want to check whether the
                presentation of the "lens-space" is the standard presentation of the 3-Sphere
                at this point. If we don't do this, the program will eventually discover that
                it has the 3-Sphere, but only after producing some redundant presentations and
                diagrams.
            **********************************************************************************/
            
            if(Length == NumGenerators && NumRelators == NumGenerators && !Boundary)
                {
                for(i = 1; i <= NumRelators; i++) if(GetHandleSize((char **) Relators[i]) != 2L) break;
                if(i > NumRelators && Delete_Dups() == NumRelators)
                    {
                    if(Micro_Print)
                        {
                        printf("\n\nThe current presentation presents the 3-Sphere.");
                        if(Micro_Print_F)
                            fprintf(myout,"\n\nThe current presentation presents the 3-Sphere.");
                        }
                    if((NumFilled >= MAX_SAVED_PRES - 3) || 
                        Save_Pres(ReadPres,0,Length,1,5,0,0,0))
                        {
                        NotConnectedError = TOO_LONG;
                        return(FALSE);
                        }        
                    BDY[NumFilled - 1] = BDY[ReadPres];
                    UDV[NumFilled - 1] = THREE_SPHERE;
                    if(CBC[CurrentComp][0] == BDRY_UNKNOWN)
                        {
                        CBC[CurrentComp][0] = 1;
                        CBC[CurrentComp][1] = BDRY_UNKNOWN;
                        }
                    Mark_As_Found_Elsewhere(CurrentComp);
                    }
                return(FALSE);    
                }
                    
            if((NumFilled >= MAX_SAVED_PRES - 3) ||Save_Pres(ReadPres,0,Length,1,20,1,0,0))
                {
                NotConnectedError = TOO_LONG;
                return(FALSE);
                }
            BDY[NumFilled - 1] = BDY[ReadPres];
            UDV[NumFilled - 1] = 0;            
            SaveUDV = 0;
            ReadPres = NumFilled - 1;
            }
            
        /**************************************************************************************
             If the program already has as many summands as it can handle, flag any other
             presentations corresponding to this summand so that we will quit processing them.
         **************************************************************************************/    
         
         if(TotalComp > MAXNUMCOMPONENTS - 3)
             {
            Mark_As_Found_Elsewhere(CurrentComp);    
            NotConnectedError = TOO_MANY_COMPONENTS;
            SysBeep(5);
            printf("\n\nStopping because the program cannot deal with any more summands. Sorry!");    
             return(FALSE);
            }
                         
        /**************************************************************************************
            Call Split_Relators() to partition the relators into two subsets. One subset will
            consist of those relators which contain generators corresponding to vertices
            occuring in the component of the Whitehead graph containing vertex number zero, 
            while the other subset will consist of the remaining relators.
        **************************************************************************************/
                
        for(i = NumDelRelators = 0; i < Vertices; i += 2) if(ZZ[i])
            {
            x = i/2 + 65;
            Split_Relators(x);                                        
            }
                                
        /**************************************************************************************
            Call Rewrite_Input() to rewrite the undeleted relators. Save the rewritten
            relators in SUR[]. Increase TotalComp by 1 and update Daughters[]. Then increase
            TotalComp by 1 more, and call Rewrite_Input() to rewrite the deleted relators.
            Save the rewritten deleted relators in SUR[].
        **************************************************************************************/
                
        Rewrite_Input();
        for(i = 1, Length = 0L; i <= NumRelators; i++)
            Length += GetHandleSize((char **) Relators[i]);
        Length -= NumRelators;
        
        SaveCS = CS[CurrentComp];
        if(!CS[CurrentComp]) CS[CurrentComp] = TRUE;
        
        TotalComp                         ++;
        UDV[ReadPres]                     = SPLIT;
        NCS[ReadPres]                    = 2;                
        Daughters[ReadPres]             = NumFilled;
        ComponentNum[NumFilled]         = TotalComp;        
        ER[NumFilled]                    = -3;
        FR[NumFilled]                    = ReadPres;        
        MLC[TotalComp][NumGenerators]     = Length;
        NG[NumFilled]                     = NumGenerators;
        NR[NumFilled]                     = NumRelators;
        PRIM[NumFilled]                    = 20;
        SURL[NumFilled]                 = Length;
        UDV[NumFilled]                     = 0;
        TP[NumFilled]                    = NumRelators;
        BDY[NumFilled]                     = BDY[ReadPres];
        OnStack                            += 2*NumGenerators;
        
        Canonical_Rewrite(Relators,FALSE,FALSE);
        
        for(i = 1; i <= NumRelators; i++)
            {
            SUR[NumFilled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));            
            if((q = *SUR[NumFilled][i]) == NULL)
                {
                for(j = 1; j < i; j++) DisposeHandle((char **) SUR[NumFilled][j]);
                MLC[TotalComp][NumGenerators]     = BIG_NUMBER;
                UDV[ReadPres]                     = SaveUDV;
                CS[CurrentComp]                 = SaveCS;
                OnStack                         -= 2*NumGenerators;
                TotalComp                         --;
                NotConnectedError                 = TOO_LONG;
                return(FALSE);
                }
            p = *Relators[i];
            while(*q++ = *p++) ;                                    
            }
            
        BytesUsed += Length;        

        for(i = 0; i < NumFilled; i++)
            if(SURL[i] == Length  
                && NG[i] == NumGenerators
                && NR[i] == NumRelators)
            {
             for(j = 1; j <= NumRelators; j++)
                 if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;    
             if(j > NumRelators && Compare_Pres(i))
                 {
                 UDV[NumFilled] = DUPLICATE;
                 Daughters[NumFilled] = i;
                 if(CBC[ComponentNum[i]][0] != BDRY_UNKNOWN)
                     {
                     h = ComponentNum[i];
                     for(j = 0; (CBC[TotalComp][j] = CBC[h][j]) < BDRY_UNKNOWN; j++) ;
                     }
                 break;
                 }    
             }
                                                                 
        NumFilled ++;
        SaveMinima = TRUE;
        
        if(Micro_Print)
            {
            printf("\nThe presentation of the first summand is:\n");
            Print_Relators(Relators,NumRelators,stdout);
            printf("\n\nSaved this presentation as: Presentation %u\n",NumFilled);
            if(Micro_Print_F)
                {
                fprintf(myout,"\nThe presentation of the first summand is:\n");
                Print_Relators(Relators,NumRelators,myout);
                fprintf(myout,"\n\nSaved this presentation as: Presentation %u\n",NumFilled);
                }
            }
                
        printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
        printf("Gen%3d  Rel%3d  Length%6lu  From%6d  NC",
            NumGenerators,NumRelators,Length,ReadPres + 1);
        
        SaveNumRelators = NumRelators;
        NumRelators = NumDelRelators;
        
        for(i = 1; i <= NumRelators; i++)
            {
            Temp = Relators[i];
            Relators[i] = DelRelators[i];
            DelRelators[i] = Temp;
            }
        
        Rewrite_Input();
        
        for(i = 1, Length = 0L; i <= NumRelators; i++)
            Length += GetHandleSize((char **) Relators[i]);
        Length -= NumRelators;
        
        TotalComp                         ++;
        ComponentNum[NumFilled]         = TotalComp;    
        ER[NumFilled]                    = -3;
        FR[NumFilled]                    = ReadPres;    
        MLC[TotalComp][NumGenerators]     = Length;
        NG[NumFilled]                     = NumGenerators;
        NR[NumFilled]                     = NumRelators;
        PRIM[NumFilled]                    = 20;    
        SURL[NumFilled]                 = Length;
        UDV[NumFilled]                     = 0;
        TP[NumFilled]                    = NumRelators;
        BDY[NumFilled]                    = BDY[ReadPres];
        OnStack                            += 2*NumGenerators;
        
        Canonical_Rewrite(Relators,FALSE,FALSE);
        
        for(i = 1; i <= NumRelators; i++)
            {
            SUR[NumFilled][i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));            
            if((q = *SUR[NumFilled][i]) == NULL)
                {
                for(j = 1; j < i; j++) DisposeHandle((char **) SUR[NumFilled][j]);
                NumFilled --;
                for(j = 1; j <= SaveNumRelators; j++) DisposeHandle((char **) SUR[NumFilled][j]);
                MLC[TotalComp][NumGenerators]     = BIG_NUMBER;
                UDV[ReadPres]                     = SaveUDV;
                CS[CurrentComp]                 = SaveCS;
                OnStack                         -= 2*NumGenerators;
                TotalComp                         --;
                NotConnectedError                 = TOO_LONG;
                return(FALSE);
                }
            p = *Relators[i];
            while(*q++ = *p++) ;                                    
            }
        
        BytesUsed += Length;
            
        for(i = 0; i < NumFilled; i++)
            if(SURL[i] == Length  
                && NG[i] == NumGenerators
                && NR[i] == NumRelators)
            {
             for(j = 1; j <= NumRelators; j++)
                 if(GetHandleSize((char **) Relators[j]) != GetHandleSize((char **) SUR[i][j])) break;    
             if(j > NumRelators && Compare_Pres(i))
                 {
                 UDV[NumFilled] = DUPLICATE;
                 Daughters[NumFilled] = i;
                 if(CBC[ComponentNum[i]][0] != BDRY_UNKNOWN)
                     {
                     h = ComponentNum[i];
                     for(j = 0; (CBC[TotalComp][j] = CBC[h][j]) < BDRY_UNKNOWN; j++) ;
                     }    
                 break;
                 }    
             }
                                                         
        NumFilled ++;
        SaveMinima = TRUE;
        
        if(Micro_Print)
            {
            printf("\nThe presentation of the second summand is:\n");
            Print_Relators(Relators,NumRelators,stdout);
            printf("\n\nSaved this presentation as: Presentation %u\n",NumFilled);
            if(Micro_Print_F)
                {
                fprintf(myout,"\nThe presentation of the second summand is:\n");
                Print_Relators(Relators,NumRelators,myout);
                fprintf(myout,"\n\nSaved this presentation as: Presentation %u\n",NumFilled);
                }
            }        
        
        if(BDY[ReadPres] == FALSE)
            {
            if(CBC[TotalComp - 1][0] == BDRY_UNKNOWN)
                {
                CBC[TotalComp - 1][0] == 1;
                CBC[TotalComp - 1][1] == BDRY_UNKNOWN;
                }
            if(CBC[TotalComp][0] == BDRY_UNKNOWN)
                {
                CBC[TotalComp][0] == 1;
                CBC[TotalComp][1] == BDRY_UNKNOWN;
                }    
            }
            
        for(i = 1; i <= NumDelRelators; i++) EmptyHandle((char **) DelRelators[i]);
        NumDelRelators = 0;
                        
        printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
        printf("Gen%3d  Rel%3d  Length%6lu  From%6d  NC",
            NumGenerators,NumRelators,Length,ReadPres + 1);            
        }
    return(Connected);                
}    

Connected_(i,k)
register unsigned int     i,
                        k;
{    
    /******************************************************************************************
        This routine finds those vertices in the component of vertex i in the graph specified
        in the adjacency lists AJ1[]. The array ZZ[] is initialized by the calling routine
        which sets the entries of vertices which should be deleted from the adjacency lists to
        a non-zero value and passes the number of deleted vertices as the parameter k. 
        The routine returns FALSE if the graph is not connected and TRUE if it is connected.
    ******************************************************************************************/    
     
    register unsigned int     h,
                            j,
                            *p,
                            *r,
                            *zz;
                            
    zz = ZZ;                         
    zz[i] = 1;
    k ++;
    for(r = UpDate,*r = i,p = r + 1; r < p; r++)
        {
        i = *r;
        for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
            {
            if(zz[j] == 0)
                {
                zz[j] = 1;
                *p++ = j;
                if(++k == Vertices) return(TRUE);
                }
            }
        }            
    if(k < Vertices) return(FALSE);        
    return(TRUE);        
}                                

void Split_Relators(register unsigned char x)
{
    /******************************************************************************************
        This routine searches through the relators in Relators[] looking for relators that
        contain the specified character x or its inverse. Any relator that contains x or its
        inverse is removed from the set of relators and saved in the array DelRelators[].
    ******************************************************************************************/
        
    register unsigned char     *p,
                            y,
                            z;
                            
    register int             i;
    
    unsigned char             **Temp;
    
    y = x + 32;
    for(i = 1; i <= NumRelators; i++)
        { 
        if(GetHandleSize((char **) Relators[i]) > 1L)
            {
            p = *Relators[i];
            while(z = *p++)
                {
                if((z == x) || (z == y))
                    {
                    Temp = DelRelators[NumDelRelators + 1];
                    DelRelators[NumDelRelators + 1] = Relators[i];
                    Relators[i] = Temp;
                    ReallocateHandle((char **) Relators[i],1L);
                    p = *Relators[i];
                    *p = EOS;        
                    NumDelRelators++;
                    break;
                    }    
                }    
            }
        }    
}

Sep_Pairs(VI,VJ)
int     VI,
        VJ;
{    
    /******************************************************************************************
        This routine determines whether the graph specified by the adjacency lists in AJ1[],
        has a pair of separating vertices. It deletes, in turn, each vertex, of valence greater
        than two, from the graph and then uses a stack-based version of depth-first-search to
        determine if the resulting graph has a separating vertex. The routine returns TRUE if
        it finds a pair of vertices which "essentially" separates the graph and otherwise
        returns FALSE. If the original graph has more than two major vertices, the routine
        calls Sep_Pairs_Sub() which returns the first pair of separating vertices (V1,V2),
        which follows the ordered pair (VI,VJ) in lexicographic order. Otherwise, it returns
        the ordered pair of major vertices.
    ******************************************************************************************/
    
    register unsigned int     i,
                            j,
                            k,
                            K,
                            m,
                            *p;
    
    unsigned int            VG[VERTICES],
                            root;                
    
    /******************************************************************************************
        First, count the number of vertices of the graph which have valence greater than two,
        and deal with the special case where there are exactly two of these.
    ******************************************************************************************/
    
    for(m = MajorVert = 0; m < Vertices; m++)
        if(VWG[m] > 2)
            {
            MajorVert++;
            if(MajorVert == 1)
                i = m;
            else
                {    
                if(MajorVert == 2)
                    j = m;
                else
                    break;
                }
            }    
    if(MajorVert == 2)
        {
        /**************************************************************************************    
                The graph has exactly two major vertices given by the values of i and j.
                If the valence of vertex i is greater than 3 or there is more than one edge 
                joining vertex i and j then the graph does not have a unique planar embedding.
        **************************************************************************************/    
                        
        if((VWG[i] > 3) || (A[i][j] > 1))
            { 
            V1 = i;
            V2 = j;
            NumComps = VWG[i];
            return(TRUE);
            }
        else
            {
            V1 = i;
            V2 = j;
            return(FALSE);                                                    
            }                
        }
        
    /******************************************************************************************
        If we are here, the graph has more than two major vertices, and we go into the main
        loop of the routine.
    ******************************************************************************************/
                            
    for(i = VI; i < Vertices - 1; i ++)
        {
        if(VWG[i] < 3) continue;
        
        for(j = 0; j < Vertices; j++) Number[j] = 0;
        if(i == 0)
            root = 1;
        else
            root = 0;
        Lowpt[i]          = 0;        
        Father[i]         = i;
        NumVert         = 1;
        Number[root]     = 1;
        Lowpt[root]      = 1;
        Father[root]     = root;
        VG[root]        = 0;
        for(p = UpDate,*p = root; p >= UpDate; p--)
            {
            NEW_VERT:
            m = *p;
            for(k = VG[m]; (j = AJ1[m][k]) < VERTICES; k++)
                {
                if(j == i) continue;
                if(Number[j] == 0)
                    {
                    NumVert     ++;
                    Number[j]     = NumVert;
                    Lowpt[j]      = NumVert;
                    Father[j]     = m;
                    VG[j]          = 0;
                    VG[m]          = k + 1;
                    p             ++;
                    *p             = j;
                    goto         NEW_VERT;        
                    }
                if(j != Father[m] && Number[j] < Lowpt[m])
                    Lowpt[m] = Number[j];        
                }
            if(Lowpt[m] < Lowpt[Father[m]]) Lowpt[Father[m]] = Lowpt[m];    
            }

        if(i == 0 && VJ == 0 && root == 1 && VWG[root] >= 3)
            {
            /**********************************************************************************
                Check whether the root of the depth-first-search tree is a cut vertex. This 
                will be the case iff the root has more than one son.
            **********************************************************************************/
        
            for(m = k = 0; (j = AJ1[root][m]) < VERTICES; m++)
                if(Father[j] == root && ++k > 1) break;
            
            if(k > 1 && Sep_Pairs_Sub(i,root)) return(TRUE);
            
            /**********************************************************************************
                If k > 1, the root of the depth-first-search tree is a cut vertex.
                Otherwise, we look for other cut vertices.
            **********************************************************************************/
            
            }

        for(j = 0,K = VERTICES; j < Vertices; j ++)
            {
            k = Father[j];
            if(VWG[k] < 3) continue;
            if(i == VI && k <= VJ) continue;
            if(i != VI && k <= i) continue;
            if(k >= K) continue;                        
            if(Lowpt[j] >= Number[k] && k != root && Sep_Pairs_Sub(i,k))
                K = k;
            }
        if(K < VERTICES) return(TRUE);        
        }    
        
    /******************************************************************************************
        We have checked each pair of distinct major vertices of the graph without finding a
                            pair which produced a proper separation.
    ******************************************************************************************/
                                                        
    return(FALSE);    
}

int Sep_Pairs_Sub(v1,v2)
{
    register int            i,
                            j,
                            k,
                            m,
                            n;

    register unsigned int    *p,
                            *r,
                            *zz;
    
    /******************************************************************************************
        We have found a pair of distinct vertices v1 and v2 which separate the graph and both
        v1 and v2 have valence greater than two. Determine how many components arise when v1
        and v2 are deleted from the graph.
    ******************************************************************************************/
    
    zz = ZZ;
    for(k = 0; k < Vertices; k++) zz[k] = 0;
    zz[v1] = zz[v2] = VERTICES;
    for(k = 0; (k < Vertices) && zz[k]; k++) ;
    
    p             = UpDate;
    *p++         = k;
    zz[k]         = 1;
    k             = 3;
    NumComps     = 1;
    while(1)
        {                
        for(r = UpDate; r < p; r++)
            {
            i = *r;
            for(m = 0; (j = AJ1[i][m]) < VERTICES; m++)
                {
                if(zz[j] == 0)
                    {
                    zz[j] = NumComps;
                    *p++ = j;
                    if(++k >= Vertices) goto FOUND_COMPS;
                    }
                }
            }
        if(k >= Vertices) break;
        if(++NumComps > 2) break;
        for(m = 0, p = UpDate; m < Vertices; m++) if(!zz[m])
            {
            zz[m] = NumComps;
            k++;
            *p++ = m;
            break;
            }    
        }
        
    FOUND_COMPS:
    
    i = v1;
    j = v2;
    
    /******************************************************************************************
        If removing the vertices i and j produced a graph with more than two components, or a
        graph with two components and vertices i and j are joined by an edge, then the graph
        is not 3-connected. (In the 2-component case, at least one component contains a
        major vertex.)
    ******************************************************************************************/    
            
    if((NumComps >= 3) || ((NumComps == 2) && A[i][j]))
        {
        V1 = i; 
        V2 = j;
        return(TRUE);
        }
        
    /******************************************************************************************
        If execution gets to this point, the graph has two components and vertices i and j are
        not joined by an edge. Check whether each component contains a vertex with valence
        greater than two. If so, the graph is not 3-connected.
    ******************************************************************************************/
        
    for(m = 0,k = FALSE,n = -1; m < Vertices; m++)
        {
        if(m == i || m == j) continue; 
        if(VWG[m] > 2)
            {
            if(n < 0) n = zz[m];
            if(n != zz[m])
                {
                k = TRUE;
                break;
                }
            }    
        }
        
    /******************************************************************************************
                If k is TRUE, each component of the graph contains a vertex with valence
                greater than two and the graph is not 3-connected.
    ******************************************************************************************/
            
    if(k)
        {
        V1 = i;
        V2 = j;
        return(TRUE);
        }
        
    return(FALSE);            
}

Planar(Flag,SaveFaces)
int     Flag,
        SaveFaces;
{
    /******************************************************************************************
        This routine determines whether the graph specified by the set of adjacency lists
        AJ1[i][] is planar. It also locates a cycle of vertices which will be the boundary of
        the "infinite" face of the graph. It returns FALSE if the graph is planar and TRUE if
        the graph is nonplanar.
    ******************************************************************************************/
    
    register int            h,
                            i,
                            j,
                            k;
                    
    register unsigned int   *MyBdry,
                            *q,
                            length;
                    
    int                     Bdry_Major_Vert,
                            Max_Bdry_Major_Vert,
                            MaxNumFaces,
                            NumInPS,
                            s1,
                            s2;
                            
    /******************************************************************************************
        First take care of the special case where the number of vertices is less than or equal
        to two. 
        Next, if the graph has more than 3(Vertices - 2) edges, return TRUE, since the graph
        is nonplanar.
    ******************************************************************************************/
                            
    if(Vertices <= 2)
        {
        FV[0] = 1;
        FV[1] = 0;
        CO[0][1] = 1;
        CO[1][0] = 0;
        NumFaces = 1;
        MaxLength = 1;
        return(FALSE);
        }        
    if(NumEdges > 3*(Vertices - 2)) return(TRUE);
                                                                                                                                    
    MaxLength             = VERTICES + 1;
    NumFaces             = 0;
    Max_Bdry_Major_Vert = -1;
    MaxNumFaces         = NumEdges - Vertices + 2;
    for(i = 0; i < Vertices; i++)
    for(q = AJ1[i]; (j = *q) < VERTICES; q++) GB[i][j] = TRUE;
    for(i = 0; i < Vertices; i++)
        {
        InDisk[i]     = 0;
        InPS[i]     = FALSE;
        ZZ[i]         = 0;
        }
        
    /******************************************************************************************
        Find an oriented edge of the graph to serve as the first edge of the boundary of the
        initial face of the graph. Then delete the oppositely oriented edge from the array
        GB[][] so that the path-finding routine, Find_Minimal_Path(), will not return with a
        minimal path which yields a degenerate cycle. Planar() now works essentially by 
        adjoining faces one-by-one to a planar surface yielding a new planar surface which
        eventually contains all of the faces of the graph -- provided the graph is planar.
        Note that Planar() works whether or not the graph has pairs of separating vertices.
    ******************************************************************************************/
            
    for(i = 0; (i < Vertices) && (VWG[i] < 3); i++) ;
    if(i == Vertices) for(i = 0; (i < Vertices) && !VWG[i]; i++) ;
    j                     = AJ1[i][0];
    MyBdry                = Bdry;
    MyBdry[0]             = i;
    MyBdry[1]             = j;
    GB[j][i]            = FALSE;
    DeletedEdgePtr         = DeletedEdges;
    *DeletedEdgePtr++     = j;
    *DeletedEdgePtr++    = i;
    InPS[i]                = 1;
    InPS[j]                = 1;
    NumInPS                = 2;
    s1                     = 0;
    s2                    = 0;
    
FIND_BDRY:

    length = Find_Minimal_Path();
    
    if(length == 0)
        {
        /**************************************************************************************
            There is no directed path joining Bdry[1] to Bdry[0]. This can only occur when the
            graph is not planar. So set NumFaces = MaxNumFaces + 1, as a flag, and goto OUTPUT.
        **************************************************************************************/    
        
        NumFaces = MaxNumFaces + 1;
        goto OUTPUT;
        }

    /******************************************************************************************
        The program has found a cycle of length at least three. See if removing the vertices
        of this cycle disconnects the graph.
    ******************************************************************************************/    
    
    #ifdef PRINT_TRIAL_CYCLES
        printf("\n");
        for(k = 0; k <= length; k++)
            {
            h = MyBdry[k];
            if(h & 1)
                c = h/2 + 97;
            else
                c = h/2 + 65;    
            printf("%c",c);
            }
    #endif
        
    if(NumInPS < Vertices && (length + 1 < Vertices))
        {
        for(h = 0; h <= length; h++) ZZ[MyBdry[h]] = VERTICES;
        k = Planar_Connected_(length + 1);
        }    
    else
        k = 1;
                                            
    if(k == 1)            
        {
        /**************************************************************************************
            Removing the vertices of the cycle does not disconnect the graph. Hence, it is a
            boundary of a face. Restore any edges which where deleted from GB[][] while looking
            for the boundary of this face. Then increment the number of faces.      
        **************************************************************************************/
        
        while(DeletedEdgePtr > DeletedEdges)
            {
            DeletedEdgePtr --;
            j = *DeletedEdgePtr--;
            GB[*DeletedEdgePtr][j] = TRUE;
            }
                    
        length      ++;                            
        NumFaces ++;
        
        /**************************************************************************************
            If Flag is FALSE, we intend to draw the graph. In this case,if the number of major
            vertices in this cycle is greater than the maximal number of major vertices in any
            previous cycle, we save this cycle in the array SaveBdry[]. Then, if the graph is
            planar, the boundary of the infinite face will contain the maximal possible number
            of major vertices and, subject to this constraint, as many vertices as possible.
        **************************************************************************************/
            
        if(Flag == FALSE)
            {
            for(Bdry_Major_Vert = h = 0; h < length; h++)
            if(VWG[MyBdry[h]] > 2) Bdry_Major_Vert ++;         
            if(Bdry_Major_Vert > Max_Bdry_Major_Vert)
                { 
                Max_Bdry_Major_Vert = Bdry_Major_Vert;
                MaxLength = length;
                for(k = 0;k < MaxLength; k++) SaveBdry[k] = MyBdry[k];    
                }
            else
            if(Bdry_Major_Vert == Max_Bdry_Major_Vert && length > MaxLength)
                { 
                MaxLength = length;
                for(k = 0;k < MaxLength; k++) SaveBdry[k] = MyBdry[k];    
                }                
            }
        
        if(SaveFaces)
            {
            for(k = 0; k < length; k++) Face[NumFaces][k] = MyBdry[k];
            Face[NumFaces][k] = VERTICES;    
            }
                        
        #ifdef PRINT_CYCLES    
            printf("\n %d) ",NumFaces);
            for(k = 0; k < length; k++)
                {
                h = MyBdry[k];
                if(h & 1)
                    c = h/2 + 97;
                else
                    c = h/2 + 65;    
                printf("%c",c);
                }
            for( ; k < 26; k++)
                printf(".");    
            printf(" bounds.");
        #endif
                                                
        /**************************************************************************************
            For each oriented edge of this cycle, set the corresponding entry of GB[][] equal
            to zero. (If the graph is planar, each oriented edge will appear in the boundary of
            exactly one oriented face.) Then update the entries of the array CO[][]. This array
            gives the cyclic order in which vertices appear around a vertex in the embedding
            of the graph.
        **************************************************************************************/

        for(k = 0,h = MyBdry[length - 1]; k < length; k++)
            {
            j = MyBdry[k];
            GB[h][j] = 0;
            h = j;
            ZZ[j] = 1;
            if(!InPS[j])
                {
                InPS[j] = TRUE;
                NumInPS ++;
                }
            }
                
        for(k = 0; k < length - 2; k++) CO[MyBdry[k+1]][MyBdry[k]] = MyBdry[k+2];
        CO[MyBdry[k+1]][MyBdry[k]] = MyBdry[0];
        CO[MyBdry[0]][MyBdry[k+1]] = MyBdry[1];
        
        /**************************************************************************************
            Look for an oriented edge, in the boundary of the current planar surface, which
            can serve as the first edge of the boundary of the next face. If no such edge
            exists, goto OUTPUT.
        **************************************************************************************/
            
        for(h = s1,i = s2; (j = AJ1[h][i]) < VERTICES; i++) if(GB[h][j])
            {
            s1 = h;
            s2 = i;
            goto FIND_FIRST_EDGE;
            }
        for(h = s1 + 1; h < Vertices; h++)
        for(i = 0; (j = AJ1[h][i]) < VERTICES; i++) if(GB[h][j])
            {
            s1 = h;
            s2 = i;
            goto FIND_FIRST_EDGE;
            }    
        FIND_FIRST_EDGE:
        
        for(h = s1,i = s2; (j = AJ1[h][i]) < VERTICES; i++) if(GB[h][j] && !GB[j][h])
            {
            MyBdry[0] = h;
            MyBdry[1] = j;
            goto FIND_BDRY;
            }
        for(h = s1 + 1; h < Vertices; h++)
        for(i = 0; (j = AJ1[h][i]) < VERTICES; i++) if(GB[h][j] && !GB[j][h])
            {
            MyBdry[0] = h;
            MyBdry[1] = j;
            goto FIND_BDRY;
            }
                    
        goto OUTPUT;                                                        
        }
    
    goto FIND_BDRY;    
    
OUTPUT:
    
    if(MaxLength == Vertices || NumFaces == MaxNumFaces)
        {
        for(i = 0; i < Vertices; i++) if(VWG[i])
            FV[i] = AJ1[i][0];
            
        /**************************************************************************************
                FV[i] is the first vertex in the counterclockwise cyclic ordering of the
                                vertices joined to vertex i.
        **************************************************************************************/
        
        #ifdef PRINT_LINKS
            for(i = 0; i < Vertices; i++) if(VWG[i])
                {
                if(i & 1)
                    c = i/2 + 97;
                else
                    c = i/2 + 65;
                printf("\nThe vertices in cyclic order about vertex %c are: ",c);
                for(j = 0; j < Vertices; j++) if(A[i][j] && i != j) break;
                if(j & 1)
                    c = j/2 + 97;
                else
                    c = j/2 + 65;
                printf("%c",c);
                k = CO[i][j];
                while(k != j) 
                    {
                    if(k & 1)
                        c = k /2 + 97;
                    else
                        c = k/2 + 65; 
                    printf("%c",c);
                    k = CO[i][k];
                    }
                }
        #endif
        return(FALSE);
        }            
    else
        return(TRUE);                                                                
}                

unsigned int Find_Minimal_Path(void)        
{
    register unsigned int    h,
                            i,
                            j,
                            *p,
                            *r;
                            
    unsigned int            *Beg1,
                            *Beg2,
                            *End1,
                            *End2,
                            FreeEdges1,
                            FreeEdges2,
                            From[VERTICES],
                            Disk1[VERTICES],
                            Disk2[VERTICES],
                            length,
                            Radius1,
                            Radius2;
    
    /******************************************************************************************
        Find_Minimal_Path() is called by Planar(). Find_Minimal_Path() finds a minimal length
        directed path, in the graph specified by GB[][], which joins the vertices Bdry[1] and
        Bdry[0]. The path is returned in the array Bdry[], while the value returned by the
        routine is equal to the length of this path. If no path can be found, the routine
        returns 0. This routine is called frequently, and it uses a couple of ideas which
        attempt to to make it efficient. First, the routine uses breadth-first-search
        starting from both vertices which are terminal vertices of the path, instead of
        searching from just one of these vertices. Secondly, the routine keeps track of the
        total number of edges it would have to search in order to increase the "radius" of
        each "disk" by one unit, and then the routine always chooses to increase the radius
        of the disk whose radius can be increased by searching the fewest edges. Finally,
        before it returns, the routine clears any entries of the array InDisk[] which it set.
        This allows the routine to find subsequent paths, in a graph with V vertices, without
        having to do V work at each call. 
    ******************************************************************************************/

    End1 = Beg1         = Disk1;
    End2 = Beg2         = Disk2;
    *End1++             = Bdry[0];
    *End2++             = Bdry[1];
    InDisk[Bdry[0]]     = 1;
    InDisk[Bdry[1]]     = 2;
    FreeEdges1             = VWG[Bdry[0]];
    FreeEdges2             = VWG[Bdry[1]];
    Radius1 = Radius2     = 0;
    
    while(1)
        {
        if(FreeEdges1 <= FreeEdges2)
            {
            Radius1 ++;
            FreeEdges1 = 0;
            for(r = Beg1,p = End1; r < p; r++)
                {
                i = *r;
                for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
                    {
                    if(InDisk[j] == 2 && GB[j][i])
                        {
                        Bdry[Radius2 + 1] = j;
                        length = Radius1 + Radius2;
                        for(h = Radius2; h > 1; h--) Bdry[h] = From[Bdry[h+1]];
                        From[j] = i;
                        for(h = Radius2 + 1; h < length; h++)
                            Bdry[h+1] = From[Bdry[h]];
                        goto FOUND_PATH;
                        }
                    if(InDisk[j] == 0 && GB[j][i])
                        {
                        From[j]     = i;
                        InDisk[j]     = 1;
                        *End1++     = j;
                        FreeEdges1     += VWG[j];
                        }
                    }
                }
            Beg1 = r;
            if(r == End1) return(0);
            }
        else
            {
            Radius2 ++;
            FreeEdges2 = 0;    
            for(r = Beg2,p = End2; r < p; r++)
                {
                i = *r;
                for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
                    {
                    if(InDisk[j] == 1 && GB[i][j])
                        {
                        Bdry[Radius2 + 1] = j;
                        length = Radius1 + Radius2;
                        for(h = Radius2 + 1; h < length; h++)
                            Bdry[h+1] = From[Bdry[h]];
                        From[j] = i;
                        for(h = Radius2; h > 1; h--) Bdry[h] = From[Bdry[h+1]];
                        goto FOUND_PATH;
                        }
                    if(InDisk[j] == 0 && GB[i][j])
                        {
                        From[j]     = i;
                        InDisk[j]     = 2;
                        *End2++     = j;
                        FreeEdges2     += VWG[j];
                        }
                    }
                }
            Beg2 = r;
            if(r == End2) return(0);
            }        
        }
    
    FOUND_PATH:
    for(r = Disk1; r < End1; r++) InDisk[*r] = 0;
    for(r = Disk2; r < End2; r++) InDisk[*r] = 0;
    return(length);                                                                        
}

int Planar_Connected_(unsigned int length)                                            
{    
    /******************************************************************************************
        This routine is used to find the components produced when the cycle specified in
        Bdry[] is deleted from the graph. It is sufficient to know which component each vertex
        neighboring the deleted cycle belongs to. Given such a vertex, the routine uses
        breadth-first search to find the component in which the vertex lies. Note that
        vertices which are known to belong to the planar-surface, and are not in the cycle,
        have their entries in ZZ[] set to 1.
    ******************************************************************************************/    
     
    register unsigned int     h,
                            i,
                            j,
                            k,
                            *p,
                            *r,
                            *zz;
    
    unsigned int            NumComps,
                            *q,
                            V,
                            W;                                            
                            
    NumComps = 1;
    p          = UpDate;                        
    zz          = ZZ;
    
    for(V = 0; V < length; V++)
    for(W = 0; (j = AJ1[Bdry[V]][W]) < VERTICES; W++) if(zz[j] == 0)
        {
        h         = NumComps + 1;
        q         = p;
        zz[j]     = h;
        *p++     = j;
        for(r = q; r < p; r++)
            {
            for(k = 0; (j = AJ1[*r][k]) < VERTICES; k++) if((i = zz[j]) < h)
                {
                if(i == 0)
                    {
                    zz[j] = h;
                    *p++  = j;
                    }
                else 
                    {
                    for(r = q; r < p; r++) zz[*r] = i;
                    goto NEXT_BDRY_VERTEX;
                    }
                }
            }
        NumComps ++;    
        NEXT_BDRY_VERTEX:
        continue;
        }
        
    if(NumComps > 1 && !Check_Bridge_Interlacing(NumComps,length))
        {
        for(r = UpDate; r < p; r++)
            zz[*r] = 0;
        for(h = 0; h < length; h++)
            {
            j = Bdry[h];
            if(InPS[j])
                zz[j] = 1;
            else
                zz[j] = 0;
            }    
        return(FALSE);
        }
    for(r = UpDate; r < p; r++) zz[*r] = 0;    
    return(TRUE);
}

int Check_Bridge_Interlacing(unsigned int NumComps, unsigned int length)            
{
    register unsigned int    h,
                            i,
                            j,
                            k,
                            kk,
                            *MyBdry,
                            *q,
                            *zz;
    
    char                    SideC[VERTICES];
    
    unsigned char            Comp[VERTICES],
                            MaxC[VERTICES],
                            MinC[VERTICES];
                            
    unsigned int            CompsFound,
                            Found_h,
                            Found_i,
                            MaxC_h,
                            MaxC_i,
                            MinC_h,
                            MinC_i,
                            V,
                            W;
    
    MyBdry             = Bdry;
    zz                 = ZZ;
    MyBdry[length]     = MyBdry[0];
                                        
    /******************************************************************************************
        For each "bridge" or component formed when the cycle specified by the path in Bdry[]
        is deleted from the graph, find the first and last points at which the bridge meets
        the path.
    ******************************************************************************************/
    
    MinC[1] = 1;
    for(h = 2; h <= NumComps; h++) MinC[h] = 0;
    for(h = i = 1; h <= length; h++)
    for(q = AJ1[MyBdry[h]]; (j = *q) < VERTICES; q++)
        {
        if((k = zz[j]) != VERTICES && MinC[k] == 0)
            {
            MinC[k] = h;
            if(++i == NumComps) goto FOUND_MINS;
            }                
        }
    FOUND_MINS:
    
    MaxC[1] = length;
    for(h = 2; h <= NumComps; h++) MaxC[h] = 0;
    for(h = length,i = 1; h > 0; h--)
    for(q = AJ1[MyBdry[h]]; (j = *q) < VERTICES; q++)
        {
        if((k = zz[j]) != VERTICES && MaxC[k] == 0)
            {
            MaxC[k] = h;
            if(++i == NumComps) goto FOUND_MAXS;
            }                
        }
    FOUND_MAXS:
    
    /******************************************************************************************
            Determine whether components are forced to be embedded on the "inside" or on the
        "outside" of this cycle.
            Start by putting component 1 on the "outside". Continue by putting other components
        on the side forced by the interlacing. Then, if there are components whose location is
        still indeterminate, put the first such component on the "outside". Continue until the
        locations of all components have been determined.
    ******************************************************************************************/
    
    for(i = 1; i <= NumComps; i++) Comp[i] = i;
    CompsFound = 0;
    
    do
        {
        CompsFound ++;
        SideC[Comp[CompsFound]] = 1;
        for(V = CompsFound; V <= CompsFound; V++)
        for(W = CompsFound + 1; W <= NumComps; W++)
            {
            h         = Comp[V];
            i         = Comp[W];
            MinC_h     = MinC[h];
            MaxC_h     = MaxC[h];
            MinC_i     = MinC[i];
            MaxC_i     = MaxC[i];
            if(MaxC_h <= MinC_i || MaxC_i <= MinC_h)
                continue;
            if(MinC_h < MinC_i && MinC_i < MaxC_h && MaxC_h < MaxC_i)
                goto FOUND_COMP;
            if(MinC_i < MinC_h && MinC_h < MaxC_i && MaxC_i < MaxC_h)
                goto FOUND_COMP;
            if(MinC_h == MinC_i && MaxC_h == MaxC_i)
                {
                Found_h = Found_i = FALSE;
                if(h == 1) for(j = MinC_h + 1; j < MaxC_h; j++)
                    {
                    if(InPS[MyBdry[j]])
                        goto FOUND_COMP;
                    for(q = AJ1[MyBdry[j]]; (kk = *q) < VERTICES; q++) if(zz[kk] == h)
                        goto FOUND_COMP;
                    }
                else for(j = MinC_h + 1; j < MaxC_h; j++)
                    {
                    for(q = AJ1[MyBdry[j]]; (kk = *q) < VERTICES; q++)
                        {
                        if(zz[kk] == h)
                            Found_h = TRUE;
                        if(zz[kk] == i)
                            Found_i = TRUE;    
                        }
                    if(Found_h && Found_i) goto FOUND_COMP;
                    }
                continue;        
                }
            if(MinC_i <= MinC_h && MaxC_h <= MaxC_i)
                {
                for(j = MinC_h + 1; j < MaxC_h; j++)
                for(q = AJ1[MyBdry[j]]; (kk = *q) < VERTICES; q++)
                    if(zz[kk] == i) goto FOUND_COMP;
                continue;
                }
            if(MinC_h <= MinC_i && MaxC_i <= MaxC_h)
                {
                if(h == 1)
                    {
                    if((MinC_i + 1) == MaxC_i
                        && !GB[MyBdry[MaxC_i]][MyBdry[MinC_i]]) goto FOUND_COMP;
                    for(j = MinC_i + 1; j < MaxC_i; j++)
                        {
                        if(InPS[MyBdry[j]]) goto FOUND_COMP;
                        for(q = AJ1[MyBdry[j]]; (kk = *q) < VERTICES; q++)
                            if(zz[kk] == h) goto FOUND_COMP;
                        }
                    }
                else
                    {
                    for(j = MinC_i + 1; j < MaxC_i; j++)
                    for(q = AJ1[MyBdry[j]]; (kk = *q) < VERTICES; q++)
                        if(zz[kk] == h) goto FOUND_COMP;
                    }
                continue;
                }
            FOUND_COMP:
            SideC[i]             = -SideC[h];
            CompsFound             ++;
            Comp[W]             = Comp[CompsFound];
            Comp[CompsFound]     = i;                        
            }
        }
    while(CompsFound < NumComps);
    
    /******************************************************************************************
        Check whether it was possible to put all of the components on the "outside" of this
        cycle.    If it was, return TRUE.
    ******************************************************************************************/

    for(h = 2; h <= NumComps && SideC[h] == 1; h++) ;
        if(h > NumComps) return(TRUE);
    
    /******************************************************************************************
        Temporarily delete the oriented edges leading from this cycle to components lying on
        the "outside". This insures that the next cycle will lie "inside" of this one.
    ******************************************************************************************/
        
    for(h = 0; h < length; h++)
        {
        i = MyBdry[h];
        for(q = AJ1[i]; (k = *q) < VERTICES; q++)
            {
            if((j = zz[k]) == VERTICES) continue;
            if(SideC[j] == 1 && GB[i][k])
                {
                *DeletedEdgePtr++    = i;
                *DeletedEdgePtr++     = k;
                GB[i][k]             = FALSE;    
                }
            }
        }
        
    /******************************************************************************************
        Temporarily delete the oriented edge of this cycle which leaves the first point at
        which an interior bridge is connected. Similarly, temporarily delete the oriented edge
        of this cycle which enters the last point at which an interior bridge is connected.
        Doing this insures that the next cycle will lie properly "inside" this one.
    ******************************************************************************************/

    for(h = 2,i = INFINITE,j = 0; h <= NumComps; h++) if(SideC[h] == -1)
        { 
        if(MinC[h] < i) i = MinC[h];
        if(MaxC[h] > j) j = MaxC[h];
        }
        
    *DeletedEdgePtr++             = MyBdry[i];
    *DeletedEdgePtr++             = MyBdry[i+1];
    GB[MyBdry[i]][MyBdry[i+1]]    = FALSE;    
    
    *DeletedEdgePtr++             = MyBdry[j-1];
    *DeletedEdgePtr++             = MyBdry[j];
    GB[MyBdry[j-1]][MyBdry[j]]     = FALSE;
    
    return(FALSE);        
}

Find_Cut_Vertices()
{
    register unsigned int     i,
                            j,
                            k,
                            *q;
                            
    unsigned int            VG[VERTICES],
                            root;
    
    /*****************************************************************************************
        This routine uses depth-first-search to determine whether the graph specified by the
        adjacency-matrix AJ1[][] has a cut-vertex. The routine returns TRUE if the graph has
        a cut-vertex, and otherwise returns FALSE.
    *****************************************************************************************/
                        
    for(j = 0; j < Vertices; j++) Number[j] = 0;
    NumVert         = 1;
    root            = 0;
    Number[root]     = 1;
    Lowpt[root]     = 1;
    Father[root]     = 0;
    VG[root]         = 0;

    for(q = UpDate,*q = root; q >= UpDate; q--)
        {
        NEW_VERT:
        i = *q;
        for(k = VG[i]; (j = AJ1[i][k]) < VERTICES; k++)
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
    
    for(i = k = 0; (j = AJ1[root][i]) < VERTICES; i++)
        if(Father[j] == root && ++k > 1) return(TRUE);
    
    /**********************************************************************************
        If k > 1, the root of the depth-first-search tree is a cut vertex.
        Otherwise, we look for other cut vertices.
    **********************************************************************************/
                        
    for(j = 0; j < Vertices; j++)
        {
        k = Father[j];
        if(Lowpt[j] >= Number[k] && k != root)
            return(TRUE);
        }    
    return(FALSE);
    
    /**********************************************************************************
        If j = Vertices, the Whitehead graph has no "cut vertex". So we return FALSE.
        Otherwise, the vertex Father[j] is a "cut vertex" and we return TRUE.
    **********************************************************************************/
}
    
Plot_Graph(k,m,F1)
int     k,
        m,
        F1;
{
#ifndef MAC
  printf("\n************ Function Plot_Graph is disabled. ************\n");
  return 0;
#else
    /****************************************************************************************** 
            This routine plots the vertices of the graph using the locations in X[] and Y[]
                determined by Gauss-Seidel() and draws the edges using the info in A[][]. 
    ******************************************************************************************/
    
    register int     i,
                    j,
                    x,
                    y;
                    
    char             c,
                    Legend[200];
                    
    Rect             *r,
                    rect;
    
    r = &rect;
    if(WhichInput == 0 && Length < SURL[WhichInput])
        {
        if(k == 1)
            sprintf(Legend,"DIAGRAM OF PRESENTATION %d':",m);
        else    
            sprintf(Legend,"DIAGRAM OF PRESENTATION %d':",WhichInput + 1);
        MoveTo(20,20);
        MyDrawString(Legend);
        if(NumGenerators == 1)
            {
            if(NumRelators == 1)
                sprintf(Legend,"  %d GENERATOR  %d RELATOR",NumGenerators,NumRelators);
            else
                sprintf(Legend,"  %d GENERATOR  %d RELATORS",NumGenerators,NumRelators);
            }
        else
            {
            if(NumRelators == 1)
                sprintf(Legend,"  %d GENERATORS  %d RELATOR",NumGenerators,NumRelators);
            else
                sprintf(Legend,"  %d GENERATORS  %d RELATORS",NumGenerators,NumRelators);
            }
        }
    else
        {
        if(k == 1)
            sprintf(Legend,"DIAGRAM OF PRESENTATION %d:",m);
        else        
            sprintf(Legend,"DIAGRAM OF PRESENTATION %d:",WhichInput + 1);
        MoveTo(20,20);
        MyDrawString(Legend);
        if(NumGenerators == 1)
            {
            if(NumRelators == 1)
                sprintf(Legend,"  %d GENERATOR  %d RELATOR",NumGenerators,NumRelators);
            else
                sprintf(Legend,"  %d GENERATOR  %d RELATORS",NumGenerators,NumRelators);
            }
        else
            {
            if(NumRelators == 1)
                sprintf(Legend,"  %d GENERATORS  %d RELATOR",NumGenerators,NumRelators);
            else
                sprintf(Legend,"  %d GENERATORS  %d RELATORS",NumGenerators,NumRelators);
            }
        }
    MyDrawString(Legend);
    if(!Connected)
        sprintf(Legend,"  LENGTH %lu: NOT CONNECTED.",Length);
    else
        {    
        if(NonPlanar)
            sprintf(Legend,"  LENGTH %lu: NONPLANAR.",Length);
        else
            {
            if(SepPairs)
                sprintf(Legend,"  LENGTH %lu: '%c' AND '%c' SEPARATE.",Length,
                    c = LSP[WhichInput],c = LSQ[WhichInput]);
            else
                sprintf(Legend,"  LENGTH %lu",Length);
            }
        }        
    MyDrawString(Legend);
    if(k == 2)
        {
        if(F1)
            sprintf(Legend,"HIT 'P' TO PRINT THIS DGRM, 'm' FOR MORE INFO, 'n' FOR NEXT PRES, 'p' FOR PREVIOUS PRES.");
        else
            sprintf(Legend,"HIT 'P' TO PRINT THIS DGRM, 'm' FOR MORE INFO, 'n' FOR NEXT, 'p' FOR PREVIOUS, 'q' TO QUIT.");                                                  
        MoveTo(20,390);
        MyDrawString(Legend);
        }
    if(k == 1)    for(i = 0; i < Vertices; i++)
        {
        X[i] *= 1.25;
        Y[i] *= 1.46;
        }    
    for(i = 0; i < Vertices -1; i ++)
        {
        for(j = i + 1; j < Vertices; j ++)
            {
            if(A[i][j])
                {
                MoveTo(X[i],Y[i]);
                LineTo(X[j],Y[j]);
                }
            }
        }
        
    /****************************************************************************************** 
                    Draw a circle around each vertex of the graph and put the letter 
                        corresponding to that vertex in the center. 
    ******************************************************************************************/
    
    for(i = 0; i < Vertices; i ++)
        {
        if(VWG[i])
            {
            x = X[i];
            y = Y[i];
            SetRect(r,x - 8,y - 8,x + 8,y + 8);
            EraseOval(r);
            FrameOval(r);                                                        
            if(i & 1)
                {
                c = i/2 + 97;
                MoveTo(x-3,y+3);
                DrawChar(c);
                }
            else
                {
                c = i/2 + 65;
                MoveTo(x-3,y+3);
                DrawChar(c);
                }
            }
        }
    if(k == 1)
        for(i = 0; i < Vertices; i++)
            {
            X[i] *= 0.8;
            Y[i] *= 0.6849;
            }    
    if(k == 2)
        {        
        /************************************************************************************** 
                    Check whether any vertices have been superimposed and set a flag 
                    in Flags[] so the top vertex in a stack of superimposed vertices
                    will "blink" in the graph. 
        **************************************************************************************/    
        
        for(i = 0; i < Vertices; i ++) Flags[i] = 0;        
        for(i = 0; i < Vertices -1; i ++) if(VWG[i])
            {
            x = X[i];
            y = Y[i];
            for(j = Vertices - 1; j > i; j-- )
                if(VWG[j] && (abs(x - X[j]) + abs(y - Y[j]) <= 4) && !Flags[i])
                    {
                    Flags[j] ++;
                    break;
                    }
            }
        }                                        
#endif
}
     
Print_Graph(F1)
int        F1;
{
#ifndef MAC
  printf("\n************ Function Print_Graph is disabled. ************\n");
  return(0);
#else
    register int     x,
                    xx,
                    y,
                    yy;
                    
    int             Error,
                    Flag,
                    i,
                    j,
                    m,
                    MinDist,
                    PrintGraph;
                    
    long             ticks;
    
    Point             mouseLoc;
    
    Rect             *r,
                    rect;
                    
    r = &rect;
    Gauss_Seidel();    
_PLOT_GRAPH:
    printf("\f");
    Plot_Graph(2,0,F1);
    
    /******************************************************************************************
        The following code allows the user to use the mouse to redraw the graph by moving the
        vertices around the screen. And, if one vertex is superimposed on top of another
        vertex, then the code causes the top vertex to "blink". This serves to alert the user
        about superimposed vertices.
    ******************************************************************************************/
    
    ObscureCursor();    
    while(1)
        {            
        for(i = Vertices - 1; i >= 0; i --)
            {
            if(Flags[i])
                {        
                 x = X[i];
                y = Y[i];
                SetRect(r,x - 8,y - 8,x + 8,y + 8);
                FrameOval(r);
                InvertOval(r);
                }
            }
        while(Button())
            {
            GetMouse(&mouseLoc);
            MinDist    = 18;
            i        = VERTICES;
            Flag     = FALSE;
            
            /**********************************************************************************
                Find the vertex of the graph which is closest to the location of the mouse.
            **********************************************************************************/
            
            for(j = 0; j < Vertices; j++) if(VWG[j])
                {
                x = X[j];
                y = Y[j];
                if(abs(x - mouseLoc.h) <= 9 && abs(y - mouseLoc.v) <= 9 
                    && abs(x - mouseLoc.h) + abs(y - mouseLoc.v) <= MinDist)
                    {
                    MinDist = abs(x - mouseLoc.h) + abs(y - mouseLoc.v);
                    i = j;
                    }
                }
                
            if(i < VERTICES)
                {
                /******************************************************************************
                    If the closest vertex is within 9 pixels horizontally and 9 pixels
                    vertically of the location of the mouse, let the user move the closest
                    vertex.
                ******************************************************************************/
                
                x = X[i];
                y = Y[i];
                SetRect(r,x - 8,y - 8,x + 8,y + 8);
                EraseOval(r);
                while(Button())
                    {
                    GetMouse(&mouseLoc);
                    if((x != mouseLoc.h) || (y != mouseLoc.v))
                        {
                        PenPat(&qd.white);
                        for(j = 0; j < Vertices; j ++) if(A[i][j])
                            {
                            MoveTo(x,y);
                            LineTo(X[j],Y[j]);
                            }
                        PenPat(&qd.black);                
                        xx = mouseLoc.h;
                        yy = mouseLoc.v;
                        
                        /**********************************************************************
                                Don't allow users to drag a vertex off screen. Otherwise
                                they might not be able to find it to drag it back.
                        **********************************************************************/
                            
                        if((30 <= xx ) && (xx <= 550) && (30 <= yy) && (yy <= 370))
                            {
                            x = xx;
                            y = yy;
                            }
                        X[i] = x;
                        Y[i] = y;
                        for(j = 0; j < Vertices; j ++) if(A[i][j])    
                            {
                            MoveTo(x,y);
                            LineTo(X[j],Y[j]);
                            }                                
                        }
                    }        
                X[i] = x;
                Y[i] = y;
                Plot_Graph(2,0,F1);
                for(j = Vertices -1; j >= 0; j --) if(Flags[j])    
                    {        
                     xx = X[j];
                    yy = Y[j];
                    SetRect(r,xx - 8,yy - 8,xx + 8,yy + 8);
                    FrameOval(r);
                    InvertOval(r);
                    }    
                Flag = TRUE;
                }        
            if(Flag) break;
            }                        
            {
            GET_RESPONSE1:    
            switch(WaitkbHit())
                {
                case 'm':
                    PrintGraph = FALSE;
                    printf("\f");
                    Error = Diagram_Data(1,F1);
                    GET_RESPONSE2:
                    switch(WaitkbHit())
                        {
                        case 'b':
                            if(!Error && Connected && !SepPairs && NumGenerators > 1)
                                {
                                Print_Bdry_Comp_Info();
                                printf("\n");
                                printf("\nHIT 'P' TO PRINT THE DIAGRAM.");
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
                                goto GET_RESPONSE2;    
                                }
                            else
                                {
                                SysBeep(5);
                                goto GET_RESPONSE2;
                                }
                        case 'D':
                            if(!Error && Connected && !SepPairs && NumGenerators > 1)
                                {
                                printf("\n");                                
                                Print_DualRelators(F1);
                                goto GET_RESPONSE2;
                                }                                
                            else
                                {
                                SysBeep(5);
                                goto GET_RESPONSE2;
                                }                                                
                        case 'm':
                            goto _PLOT_GRAPH;
                        case 'O':
                            if(!Error && Connected && !SepPairs && NumGenerators > 1)
                                {
                                printf("\n");                                
                                Print_OutRelators(F1);
                                goto GET_RESPONSE2;
                                }                                
                            else
                                {
                                SysBeep(5);
                                goto GET_RESPONSE2;
                                }                                
                        case 'P':
                            PrintGraph = TRUE;
                            break;
                        case 'p':
                            printf("\f");
                            return(1);                                    
                        case 'v':
                            printf("\n\nPresentation %d is:\n",WhichInput + 1);
                            Print_Relators(Relators,NumRelators,stdout);
                            printf("\n");
                            printf("\nHIT 'P' TO PRINT THE DIAGRAM.");
                            printf("\n   HIT 'm' TO RETURN TO DIAGRAM %d.",WhichInput + 1);
                            printf("\n      HIT 'v' TO REVIEW PRESENTATION %d.",WhichInput + 1);
                            if(F1)
                                {
                                printf("\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
                                printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
                                if(!Error && Connected && !SepPairs && NumGenerators > 1)
                                    {
                                    printf("\n               HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");                                
                                    printf("\n                  HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
                                    printf("\n                     HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");
                                    }
                                }
                            else
                                {    
                                printf("\n         HIT 'n' TO SEE THE NEXT DIAGRAM.");
                                printf("\n            HIT 'p' TO SEE THE PREVIOUS DIAGRAM.");
                                printf("\n               HIT 'q' TO QUIT VIEWING DIAGRAMS.");
                                if(!Error && Connected && !SepPairs && NumGenerators > 1)
                                    {
                                    printf("\n                  HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");
                                    printf("\n                     HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
                                    printf("\n                        HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");
                                    }
                                }
                            goto GET_RESPONSE2;    
                        case 'n':
                            printf("\f");
                            return(0);
                        case 'q':
                            if(F1)
                                {
                                SysBeep(5);
                                goto GET_RESPONSE2;                                
                                }
                            WhichInput = NumFilled + 1;
                            printf("\f");
                            return(0);    
                        default:
                            SysBeep(5);
                            goto GET_RESPONSE2;
                            break;
                        }
                    break;                    
                case 'P':
                    PrintGraph = TRUE;    
                    break;
                case 'p':
                    printf("\f");
                    return(1);    
                case 'n':
                    printf("\f");
                    return(0);
                case 'q':
                    if(F1)
                        {
                        SysBeep(5);
                        goto GET_RESPONSE1;                            
                        }
                    WhichInput = NumFilled + 1;
                    printf("\f");
                    return(0);        
                default:
                    SysBeep(5);
                    goto GET_RESPONSE1;    
                }
            break;        
            }        
        ticks = TickCount();    
        while((TickCount() - ticks)    < 20L) ;
        }                                                    
    if(PrintGraph)
        {
        #ifdef MANUALLY_RENUMBER_PRESENTATIONS
            printf("\n\nPRINT AS PRESENTATION NUMBER ?");
            scanf("%d",&m);
        #else
            m = WhichInput + 1;
        #endif
        i = Print_Picture(m);
        printf("\f");
        if(i)
            {
            SysBeep(5);
            printf("\n                        Unable to print the diagram. Sorry!");        
            }        
        Diagram_Data(2,F1);
        printf("\n\nPRINT THIS DATA ? HIT 'y' OR 'n'.");
        if(WaitkbHit() == 'y' && Print_Diagram_Data(m))
            {
            SysBeep(5);
            printf("\n                        Unable to print this data. Sorry!");            
            }
        printf("\n\nPRINT THE PRESENTATION ? HIT 'y' OR 'n'.");
        if(WaitkbHit() == 'y' && Print_Presentation(m))
            {
            SysBeep(5);
            printf("\n                        Unable to print this presentation. Sorry!");        
            }
        }    
    printf("\f");
    return(0);
#endif
}

Diagram_Data(PrintOut,F1)
int     PrintOut,
        F1;
{
    register int     i,
                    j;
                    
    int                Error,
                    k,
                    SWhichInput;                
                    
    unsigned char     *p,
                    x,
                    y;

    unsigned int Diagram_Main();
    
    Error = FALSE;
    
if(PrintOut == 2)    
fprintf(myout,"\n\n|******************************************************************************|");    
    
    for(k = 0; k < PrintOut; k++)
        {
        if(k == 0)
            fptr = stdout;
        else
            fptr = myout;    
        fprintf(fptr,"\n\n                            Data For Diagram %d",WhichInput + 1);
        fprintf(fptr,"\n\nGenerators %d, Relators %d, Length %lu.",NumGenerators,NumRelators,Length);
        fprintf(fptr,"\n\nThe following table gives the number of edges joining each pair of vertices.\n");
        
        for(i = 0; i < Vertices - 1; i++)
            {
            for(j = i + 1; j < Vertices && !A[i][j]; j++) ;
            if(j >= Vertices) continue;
            if(i & 1)
                x = i/2 + 97;
            else
                x = i/2 + 65;
            fprintf(fptr,"\n%c: ",x);    
            for(j = i + 1; j < Vertices; j++)
                {
                if(A[i][j])
                    {
                    if(j & 1)
                        y = j/2 + 97;
                    else
                        y = j/2 + 65;
                    fprintf(fptr,"%c%u ",y,A[i][j]);    
                    }    
                }    
            }

        if(!Connected)
            {
            Error = TRUE;
            fprintf(fptr,"\n\nThe diagram is not connected. So complete data is unavailable.");
            goto END;
            }        
        if(NonPlanar)
            {
            Error = 2;
            fprintf(fptr,"\n\nThe diagram is nonplanar. So complete data is unavailable.");
            goto END;
            }
        else if(SepPairs)
            {
            Error = TRUE;
            if(UDV[WhichInput] == ANNULUS_EXISTS)
                {
                p = *SUR[WhichInput][0];
                x = *p++;
                y = *p++;                
                fprintf(fptr,"\n\nVertices '%c' and '%c' separate the diagram.",x,y);
                fprintf(fptr,"\nThe component consisting of vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    fprintf(fptr,",%c",x);
                    p++;
                    }
                fprintf(fptr,"}");    
                p++;        
                fprintf(fptr,"\nlies in an annulus which swallows the component and otherwise follows the curve:");
                fprintf(fptr,"\n%s.",p);
                }
            else
                {
                if(UDV[WhichInput] == V2_ANNULUS_EXISTS)
                    {
                    p = *SUR[WhichInput][0];
                    fprintf(fptr,"\n\nThere exists an annulus which swallows vertice(s) {%c",*p);
                    p++;
                    while((x = *p) != '@')
                        {
                        fprintf(fptr,",%c",x);
                        p++;
                        }
                    fprintf(fptr,"}");    
                    p++;        
                    fprintf(fptr,"\nand otherwise follows the curve:");
                    fprintf(fptr,"\n%s.",p);
                    }        
                else
                if(UDV[WhichInput] == SEP_PAIRS)
                        fprintf(fptr,"\n\nVertices '%c' and '%c' separate the Whitehead graph. So complete data is unavailable.",
                        x = LSP[WhichInput],y = LSQ[WhichInput]);    
                }    
            }
        else if( NumGenerators > 1 && Connected && !SepPairs && !NonPlanar)    
            {
            if(k == 0)
                {
                TestRealizability1 = TRUE;
                switch(Diagram_Main())
                    {
                    case TOO_LONG:
                        fprintf(fptr,"\n\nThe presentation is too long. So complete data is unavailable.");
                        TestRealizability1 = FALSE;
                        Error = TRUE;
                        goto END;
                    case FATAL_ERROR:
                        fprintf(fptr,"\n\nThe presentation is not realizable. So complete data is unavailable.");
                        TestRealizability1 = FALSE;
                        Error = TRUE;
                        goto END;
                    case NON_UNIQUE_1:
                        if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_1;
                        break;
                    case NON_UNIQUE_2:
                        if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_2;
                        break;
                    case NON_UNIQUE_3:
                        if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_3;
                        break;
                    case NON_UNIQUE_4:
                        if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_4;
                        break;
                    case V2_ANNULUS_EXISTS:
                        Error = TRUE;
                        break;    
                    }
                TestRealizability1 = FALSE;
                }        
        
            fprintf(fptr,"\n\nFor each (X,x) pair of vertices:");
            fprintf(fptr,"\n1) Number the edges at vertex X counter-clockwise about vertex X giving the");
            fprintf(fptr,"\n   'first-edge' at vertex X number 0.");
            fprintf(fptr,"\n2) Number the edges at vertex x clockwise about vertex x, giving the");
            fprintf(fptr,"\n   'first-edge' at vertex x the number shown in the following list.\n\n");    
    
            fprintf(fptr,"(%u",OSA[0] % VA[0]);    
            for(i = 2; i < Vertices; i += 2) fprintf(fptr,",%u",OSA[i] % VA[i/2]);
            fprintf(fptr,")");    

            SWhichInput = WhichInput;
            if(UDV[WhichInput] == DUPLICATE)
                {
                SWhichInput = Daughters[WhichInput];
                fprintf(fptr,"\n\nPresentation %d is a duplicate of presentation %d of summand %u.",
                WhichInput + 1, SWhichInput + 1,ComponentNum[SWhichInput]);
                }
                    
            switch(UDV[SWhichInput])
                {
                 case NON_UNIQUE_4:
                    fprintf(fptr,"\n\n     Note,the diagram is not unique because there is a generator which appears");
                    fprintf(fptr,"\nwith only one exponent and that exponent is greater than 6.");
                    break;
                case NON_UNIQUE_3:
                    fprintf(fptr,"\n\n     Note,the diagram is not unique because there is a generator which appears");
                    fprintf(fptr,"\nonly with exponent 5.");
                    break;
                case NON_UNIQUE_2:
                    fprintf(fptr,"\n\n     Note,the diagram is not unique because there is a generator which appears");
                    fprintf(fptr,"\nonly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
                    break;
                case NON_UNIQUE_1:
                    fprintf(fptr,"\n\n     Note,the diagram is not unique because there is a generator which appears");
                    fprintf(fptr,"\nwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
                    break;
                case ANNULUS_EXISTS:
                    p = *SUR[SWhichInput][0];
                    x = *p++;
                    y = *p++;                
                    fprintf(fptr,"\n\nVertices '%c' and '%c' separate the diagram.",x,y);
                    fprintf(fptr,"\nThe component consisting of vertice(s) { %c",*p);
                    p++;
                    while((x = *p) != '@')
                        {
                        fprintf(fptr,",%c",x);
                        p++;
                        }
                    fprintf(fptr," }");    
                    p++;        
                    fprintf(fptr,"\nlies in an annulus which swallows the component and otherwise follows the curve:");
                    fprintf(fptr,"\n%s.",p);
                    Error = TRUE;
                    break;
                case V2_ANNULUS_EXISTS:
                    p = *SUR[SWhichInput][0];
                    fprintf(fptr,"\n\nThere exists an annulus which swallows vertice(s) { %c",*p);
                    p++;
                    while((x = *p) != '@')
                        {
                        fprintf(fptr,",%c",x);
                        p++;
                        }
                    fprintf(fptr," }");    
                    p++;        
                    fprintf(fptr,"\nand otherwise follows the curve:");
                    fprintf(fptr,"\n%s.",p);
                    Error = TRUE;
                    break;                
                }            
            }
        }

END:                    
    if(PrintOut == 2)
        printf("\n\nTHESE TABLES ARE PRINTED IN THE FILE 'Heegaard_Results'.");
    if(PrintOut == 1)
        {
        printf("\n");
        printf("\nHIT 'P' TO PRINT THE DIAGRAM.");
        printf("\n   HIT 'm' TO RETURN TO DIAGRAM %d.",WhichInput + 1);
        printf("\n      HIT 'v' TO REVIEW PRESENTATION %d.",WhichInput + 1);
        if(F1)
            {
            printf("\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
            printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
            if(!Error && Connected && !SepPairs && NumGenerators > 1)
                {
                printf("\n               HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");            
                printf("\n                  HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
                printf("\n                     HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");
                }
            }
        else
            {    
            printf("\n         HIT 'n' TO SEE THE NEXT DIAGRAM.");
            printf("\n            HIT 'p' TO SEE THE PREVIOUS DIAGRAM.");
            printf("\n               HIT 'q' TO QUIT VIEWING DIAGRAMS.");
            if(!Error && Connected && !SepPairs && NumGenerators > 1)
                {
                printf("\n                  HIT 'b' FOR INFO ABOUT THE BDRY OF THIS MANIFOLD.");
                printf("\n                     HIT 'D' TO SEE THE 'DUAL' RELATORS FOR THIS DIAGRAM.");
                printf("\n                        HIT 'O' TO SEE THE 'OUT' RELATORS FOR THIS DIAGRAM.");
                }
            }
        }
    return(Error);            
} 

void Floating_Gauss_Seidel(void)
{    
    /****************************************************************************************** 
        This routine determines where the vertices of the Heegaard diagram are to be located
        in the plane. First the locations of the vertices that belong to the boundary of the
        "infinite" face of the diagram are determined. Then each remaining vertex of the
        diagram is located at the barycenter of the locations of its neighbors. Determining 
        the locations of these vertices requires solving a set of linear equations. We do this
        using the Gauss-Seidel method.
    ******************************************************************************************/
    
    register int     h,
                    i,
                    j,
                    k;
                    
    register double x,
                    y,
                    z;
                                    
    int             converge,
                    NumVert,
                    PVert,
                    Vert;
    
    double             bottom,
                    left,
                    RIC[VERTICES],
                    right,
                    top,
                    X1,
                    X2,
                    XF[VERTICES],
                    Y1,
                    Y2,
                    YF[VERTICES];
                                
    for(i = 0; i < Vertices; i++)
        {
        XF[i] = 290.0;    
        YF[i] = 200.0;
        if(VWG[i]) Flags[i] = 1;
        else Flags[i] = 0;
        }    
        
    /******************************************************************************************
        Plot the vertices of the "infinite" face of the graph. These vertices will be spaced 
        around the perimeter of a rectangle on the screen -- unless there are only 3 vertices
        in the boundary of the "infinite" face. So that the cyclic orders of the vertices are
        represented correctly, the order in which vertices appear in the boundary of the
        infinite face is such that the interior of the infinite face lies to the right-hand
        side of its oriented boundary.
    ******************************************************************************************/
    
    if((MaxLength < 2) || (MaxLength > VERTICES)) goto _DONE;        
    
    if(MaxLength == 3)
        {
        x = 290.0;
        y = 30.0;
        for(i = 0; i < MaxLength; i++)
            {
            j = SaveBdry[i];
            XF[j] = x;
            YF[j] = y;
            if(x == 290.0)
                {
                x = 30.0;
                y = 370.0;
                }
            else if(x == 30.0)
                x = 550.0;
            Flags[j] = 0;
            }    
        }
    else    
        {
        top = 31.0;
        left = 31.0;
        right = 549.0;
        bottom = 369.0;
        i = MaxLength/4;
        switch(MaxLength % 4)
            {
            case 0:
                X2 = X1 = 520.0/(double)i;
                Y2 = Y1 = 340.0/(double)i;
                break;
            case 1:
                X1 = 520.0/(double)(i+1);
                X2 = 520.0/(double)i;
                Y2 = Y1 = 340.0/(double)i;
                break;
            case 2:
                X2 = X1 = 520.0/(double)(i+1);
                if(i)
                Y2 = Y1 = 340.0/(double)i;
                break;
            case 3:
                X2 = X1 = 520.0/(double)(i+1);
                Y1 = 340.0/(double)(i+1);
                if(i)
                Y2 = 340.0/(double)i;
                break;
            }
        x = 30.0;
        y = 30.0;        
        for(i = 0; i < MaxLength; i++)
            {
            j = SaveBdry[i];
            XF[j] = x;
            YF[j] = y;
            if(y < top && x > left) x -= X1;
            else
            if(x < left && y < bottom) y += Y1;
            else
            if(y > bottom && x < right) x += X2;
            else
            if(x > right && y > top) y -= Y2;
            Flags[j] = 0;
            }    
        }            
        
    for(i = 0; i < Vertices; i++)
        {
        XF[i] -= 290.0;
        YF[i] -= 200.0;
        }
                    
    /******************************************************************************************
        Save the reciprocals of the valences of those vertices with nonzero valence in the
                                    array RIC[].
    ******************************************************************************************/
    
    if(NonPlanar || !Connected)
        {
        for(i = 0;i < Vertices; i++) if(VWG[i])
        RIC[i] = 1.0/((double)VWG[i]);
        }
    else for(i = 0;i < Vertices; i++) if(VWG[i] && Flags[i])
        {
        
        /**************************************************************************************
            Add some "virtual" edges to the graph so that the graph will have a more visually
            appealling embedding on screen.
        **************************************************************************************/
            
        NumVert = 0;
        for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
            {
            Vert = j;
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            Vert = CO[j][i];
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            Vert = i;
            do
                {
                PVert = Vert;
                Vert = CO[j][PVert];
                }
            while(Vert != i);
            Vert = PVert;    
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            }
        AJ3[i][NumVert] = VERTICES;
        z = NumVert + 2.0*VWG[i];
        RIC[i] = 1.0/z;
        }

    /******************************************************************************************
        Iterate until two successive iterations don't move any vertex by more than ERROR.
    ******************************************************************************************/
    
    if(NonPlanar || !Connected)
        {
        converge = 1;
        k = 0;
        while (1)
            {
            k ++;
            for(i = 0; i < Vertices; i++) if(Flags[i])
                {
                x = y = 0.0;
                z = RIC[i];
                for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)    
                    {
                    x += XF[j];
                    y += YF[j];
                    }
                x *= z;
                y *= z;                                
                if(converge)
                    if((fabs(x - XF[i]) + fabs(y - YF[i])) > ERROR) converge = 0;                        
                XF[i] = x;
                YF[i] = y;
                }
            if (++converge == 2 || k > 3000) break;                        
             }
         }    
    else
        {                    
        converge = 1;
        k = 0;
        while (1)
            {
            k ++;
            for(i = 0; i < Vertices; i++) if(Flags[i])
                {
                x = y = 0.0;
                z = RIC[i];
                for(h = 0; (j = AJ3[i][h]) < VERTICES; h++)    
                    {
                    if(A[i][j])
                        {
                        x += 3.0*XF[j];
                        y += 3.0*YF[j];
                        }
                    else
                        {
                        x += XF[j];
                        y += YF[j];
                        }
                    }
                x *= z;
                y *= z;                                
                if(converge)
                    if((fabs(x - XF[i]) + fabs(y - YF[i])) > ERROR) converge = 0;                        
                XF[i] = x;
                YF[i] = y;
                }
            if (++converge == 2 || k > 3000) break;                        
             }
         }    
     
    /******************************************************************************************
        Shrink the interior of the complement of the 'infinite' face of the graph slightly so
        that vertices which should not appear on the boundary of the infinite face do not
        appear there. Then translate the center of the graph to the center of the screen.
    ******************************************************************************************/
    
     for(i = 0; i < Vertices; i++)
        {
        if(Flags[i])
            {
            XF[i] *= 0.95;
            YF[i] *= 0.95;
            }
        XF[i] += 290.0;
        YF[i] += 200.0;
        }
        
_DONE:        
     /*****************************************************************************************
             Convert the floating point values in XF[] and YF[] to integers in X[] and Y[].
     *****************************************************************************************/    
     
     for(i = 0; i < Vertices; i ++) 
        {
        X[i] = XF[i] + 0.5;
        Y[i] = YF[i] + 0.5;
        }                    
}

void Gauss_Seidel(void)
{    
    /****************************************************************************************** 
        This routine determines where the vertices of the Heegaard diagram are to be located
        in the plane. First the locations of the vertices that belong to the boundary of the
        "infinite" face of the diagram are determined. Then each remaining vertex of the
        diagram is located at the barycenter of the locations of its neighbors. Determining 
        the locations of these vertices requires solving a set of linear equations. We do this
        using the Gauss-Seidel method.
    ******************************************************************************************/
    
    register int     h,
                    i,
                    j,
                    k;
                    
    register long    x,
                    y,
                    z;
                                    
    int             converge,
                    NumVert,
                    PVert,
                    Vert;
    
    long             bottom,
                    left,
                    RIC[VERTICES],
                    right,
                    top,
                    X1,
                    X2,
                    XF[VERTICES],
                    Y1,
                    Y2,
                    YF[VERTICES];
                                
    for(i = 0; i < Vertices; i++)
        {
        XF[i] = 290000;    
        YF[i] = 200000;
        if(VWG[i]) Flags[i] = 1;
        else Flags[i] = 0;
        }    
        
    /******************************************************************************************
        Plot the vertices of the "infinite" face of the graph. These vertices will be spaced 
        around the perimeter of a rectangle on the screen -- unless there are only 3 vertices
        in the boundary of the "infinite" face. So that the cyclic orders of the vertices are
        represented correctly, the order in which vertices appear in the boundary of the
        infinite face is such that the interior of the infinite face lies to the right-hand
        side of its oriented boundary.
    ******************************************************************************************/
    
    if((MaxLength < 2) || (MaxLength > VERTICES)) goto _DONE;        
    
    if(MaxLength == 3)
        {
        x = 290000;
        y = 30000;
        for(i = 0; i < MaxLength; i++)
            {
            j = SaveBdry[i];
            XF[j] = x;
            YF[j] = y;
            if(x == 290000)
                {
                x = 30000;
                y = 370000;
                }
            else if(x == 30000)
                x = 550000;
            Flags[j] = 0;
            }    
        }
    else    
        {
        top = 31000;
        left = 31000;
        right = 549000;
        bottom = 369000;
        i = MaxLength/4;
        switch(MaxLength % 4)
            {
            case 0:
                X2 = X1 = 520000/i;
                Y2 = Y1 = 340000/i;
                break;
            case 1:
                X1 = 520000/(i+1);
                X2 = 520000/i;
                Y2 = Y1 = 340000/i;
                break;
            case 2:
                X2 = X1 = 520000/(i+1);
                if(i)
                Y2 = Y1 = 340000/i;
                break;
            case 3:
                X2 = X1 = 520000/(i+1);
                Y1 = 340000/(i+1);
                if(i)
                Y2 = 340000/i;
                break;
            }
        x = 30000;
        y = 30000;        
        for(i = 0; i < MaxLength; i++)
            {
            j = SaveBdry[i];
            XF[j] = x;
            YF[j] = y;
            if(y < top && x > left) x -= X1;
            else
            if(x < left && y < bottom) y += Y1;
            else
            if(y > bottom && x < right) x += X2;
            else
            if(x > right && y > top) y -= Y2;
            Flags[j] = 0;
            }    
        }            
        
    for(i = 0; i < Vertices; i++)
        {
        XF[i] -= 290000;
        YF[i] -= 200000;
        }
                    
    /******************************************************************************************
        Save the valences of those vertices with nonzero valence in the array RIC[].
    ******************************************************************************************/
    
    if(NonPlanar || !Connected)
        {
        for(i = 0;i < Vertices; i++) if(VWG[i])
        RIC[i] = VWG[i];
        }
    else for(i = 0;i < Vertices; i++) if(VWG[i] && Flags[i])
        {
        
        /**************************************************************************************
            Add some "virtual" edges to the graph so that the graph will have a more visually
            appealling embedding on screen.
        **************************************************************************************/
            
        NumVert = 0;
        for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
            {
            Vert = j;
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            Vert = CO[j][i];
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            Vert = i;
            do
                {
                PVert = Vert;
                Vert = CO[j][PVert];
                }
            while(Vert != i);
            Vert = PVert;    
            for(k = 0; (k < NumVert) && AJ3[i][k] != Vert; k++) ;
            if(k == NumVert)
                {
                AJ3[i][NumVert] = Vert;
                NumVert ++;
                }
            }
        AJ3[i][NumVert] = VERTICES;
        z = NumVert + 2*VWG[i];
        RIC[i] = z;
        }

    /******************************************************************************************
        Iterate until two successive iterations don't move any vertex by more than ERROR.
    ******************************************************************************************/
    
    if(NonPlanar || !Connected)
        {
        converge = 1;
        k = 0;
        while (1)
            {
            k ++;
            for(i = 0; i < Vertices; i++) if(Flags[i])
                {
                x = y = 0;
                z = RIC[i];
                for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)    
                    {
                    x += XF[j];
                    y += YF[j];
                    }
                x = x/z;
                y = y/z;                                
                if(converge)
                    if((labs(x - XF[i]) + labs(y - YF[i])) > 100) converge = 0;                        
                XF[i] = x;
                YF[i] = y;
                }
            if (++converge == 2 || k > 3000) break;                        
             }
         }    
    else
        {                    
        converge = 1;
        k = 0;
        while (1)
            {
            k ++;
            for(i = 0; i < Vertices; i++) if(Flags[i])
                {
                x = y = 0;
                z = RIC[i];
                for(h = 0; (j = AJ3[i][h]) < VERTICES; h++)    
                    {
                    if(A[i][j])
                        {
                        x += 3*XF[j];
                        y += 3*YF[j];
                        }
                    else
                        {
                        x += XF[j];
                        y += YF[j];
                        }
                    }
                x = x/z;
                y = y/z;                                
                if(converge)
                    if((labs(x - XF[i]) + labs(y - YF[i])) > 100) converge = 0;                        
                XF[i] = x;
                YF[i] = y;
                }
            if (++converge == 2 || k > 3000) break;                        
             }
         }    
     
    /******************************************************************************************
        Shrink the interior of the complement of the 'infinite' face of the graph slightly so
        that vertices which should not appear on the boundary of the infinite face do not
        appear there. Then translate the center of the graph to the center of the screen.
    ******************************************************************************************/
    
     for(i = 0; i < Vertices; i++)
        {
        if(Flags[i])
            {
            XF[i] *= 95;
            XF[i] = XF[i]/100;
            YF[i] *= 95;
            YF[i] = YF[i]/100;
            }
        XF[i] += 290000;
        YF[i] += 200000;
        }
        
_DONE:        
     /*****************************************************************************************
             Convert the long integer values in XF[] and YF[] to integers in X[] and Y[].
     *****************************************************************************************/    
     
     for(i = 0; i < Vertices; i ++) 
        {
        X[i] = XF[i]/1000;
        Y[i] = YF[i]/1000;
        }                    
}

Print_Picture(m)
int m;
{
#ifndef MAC
  printf("Function Print_Picture is disabled.\n\r");
#else 
    TPrStatus        prStatus;                    
    THPrint            hPrint;
    TPPrPort        pPrPort;
    GrafPtr            savePort;
    int                Error;
    
    Error = FALSE;
    GetPort(&savePort);
    PrOpen();
    hPrint = (THPrint) NewHandle(sizeof(TPrint));
    PrintDefault(hPrint);
    PrValidate(hPrint);
    if(PrStlDialog(hPrint) == FALSE)
        {
        Error = TRUE;
        goto _END;
        }
    if(PrJobDialog(hPrint) == FALSE)
        {
        Error = TRUE;
        goto _END;
        }
    prStatus.iTotPages = 1;    
    pPrPort = PrOpenDoc(hPrint,NULL,NULL);
    if(PrError() == noErr)
        {
        PrOpenPage(pPrPort,NULL);
        TextFont(4);
        TextSize(9);
        TextFace(0);
        if(PrError() == noErr)
            Plot_Graph(1,m,0);
        else
            {
            PrClosePage(pPrPort);
            PrCloseDoc(pPrPort);
            Error = TRUE;
            goto _END;
            }    
        PrClosePage(pPrPort);
        PrCloseDoc(pPrPort);
        if((**hPrint).prJob.bJDocLoop == bSpoolLoop && PrError() == noErr)
        PrPicFile(hPrint,NULL,NULL,NULL,&prStatus);
        }
    else
        {
        PrCloseDoc(pPrPort);
        Error = TRUE;
        }    
_END:        
    PrClose();
    PrDrvrClose();
    SetPort(savePort);
    return(Error);
#endif
}    

Print_Diagram_Data(m)
int m;
{
#ifndef MAC
  printf("Function Print_Diagram_Data is disabled.\n\r");
#else
    register int     i,
                    j;
                    
    unsigned char     *p,
                    x,
                    y;
                    
    char            **HText,
                    *textptr;
                    
    long            length;                    
                    
    GrafPtr            SavePort;                                

    unsigned int Diagram_Main();
    
    HText = (char **) NewHandle(35000L);
    if((textptr = *HText) == NULL) return(TOO_LONG);
    HLock(HText);
            
    textptr += sprintf(textptr,"                                DATA FOR DIAGRAM %d",m);
    textptr --;
    textptr += sprintf(textptr,"\r\rGENERATORS %d, RELATORS %d, LENGTH %lu.",NumGenerators,NumRelators,Length);
    textptr --;
    textptr += sprintf(textptr,"\r\rTHE FOLLOWING TABLE GIVES THE NUMBER OF EDGES JOINING EACH PAIR OF VERTICES.\r");
    textptr --;
    
    for(i = 0; i < Vertices - 1; i++)
        {
        for(j = i + 1; j < Vertices && !A[i][j]; j++) ;
        if(j >= Vertices) continue;
        if(i & 1)
            x = i/2 + 97;
        else
            x = i/2 + 65;
        textptr += sprintf(textptr,"\r%c: ",x);
        textptr --;    
        for(j = i + 1; j < Vertices; j++)
            {
            if(A[i][j])
                {
                if(j & 1)
                    y = j/2 + 97;
                else
                    y = j/2 + 65;
                textptr += sprintf(textptr,"%c%u ",y,A[i][j]);
                textptr --;    
                }    
            }    
        }

    if(!Connected)
        {
        textptr += sprintf(textptr,"\r\rThe diagram is not connected, so complete data is unavailable.");    
        textptr --;
        }
    else if(SepPairs)
        {
        textptr += sprintf(textptr,"\r\rVertices '%c' and '%c' separate the Whitehead graph, so complete data is unavailable.",
        x = LSP[WhichInput],y = LSQ[WhichInput]);
        textptr --;
        }
    else if(NumGenerators > 1 && Connected && !SepPairs && !NonPlanar)    
        {
        length = textptr - *HText;
        HUnlock(HText);
        switch(Diagram_Main())
            {
            case NON_UNIQUE_1:
                if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_1;
                break;
            case NON_UNIQUE_2:
                if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_2;
                break;
            case NON_UNIQUE_3:
                if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_3;
                break;
            case NON_UNIQUE_4:
                if(UDV[WhichInput] < DONE) UDV[WhichInput] = NON_UNIQUE_4;
                break;
            }    
        
        HLock(HText);
        textptr = *HText + length;
        textptr += sprintf(textptr,"\r\rFOR EACH (X,x) PAIR OF VERTICES, IDENTIFY EDGE 0 OF THE FIRST VERTEX");
        textptr --;
        textptr += sprintf(textptr,"\rWITH THE EDGE OF THE SECOND VERTEX SHOWN IN THE FOLLOWING LIST.\r\r");
        textptr --;
        textptr += sprintf(textptr,"(%u",OSA[0] % VA[0]);
        textptr --;
        for(i = 2; i < Vertices; i += 2)
            {
            textptr += sprintf(textptr,",%u",OSA[i] % VA[i/2]);
            textptr --;
            }
        textptr += sprintf(textptr,")");
        textptr --;
    
        switch(UDV[WhichInput])
            {
            case NON_UNIQUE_4:
                textptr += sprintf(textptr,"\r\r     Note,the diagram is not unique because there is a generator which appears");
                textptr --;
                textptr += sprintf(textptr,"\rwith only one exponent and that exponent is greater than 6.");
                textptr --;
                break;
            case NON_UNIQUE_3:
                textptr += sprintf(textptr,"\r\r     Note,the diagram is not unique because there is a generator which appears");
                textptr --;
                textptr += sprintf(textptr,"\ronly with exponent 5.");
                textptr --;
                break;
            case NON_UNIQUE_2:
                textptr += sprintf(textptr,"\r\r     Note,the diagram is not unique because there is a generator which appears");
                textptr --;
                textptr += sprintf(textptr,"\ronly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
                textptr --;
                break;
            case NON_UNIQUE_1:
                textptr += sprintf(textptr,"\r\r     Note,the diagram is not unique because there is a generator which appears");
                textptr --;
                textptr += sprintf(textptr,"\rwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
                textptr --;
                break;    
            case ANNULUS_EXISTS:
                p = *SUR[WhichInput][0];
                x = *p++;
                y = *p++;                
                textptr += sprintf(textptr,"\r\rNote that, vertices '%c' and '%c' separate the diagram.",x,y);
                textptr --;
                textptr += sprintf(textptr,"\rThe component consisting of vertice(s) { %c",*p);
                textptr --;
                p++;
                while((x = *p) != '@')
                    {
                    textptr += sprintf(textptr,",%c",x);
                    textptr --;
                    p++;
                    }
                textptr += sprintf(textptr," }");
                textptr --;
                p++;        
                textptr += sprintf(textptr,"\rlies in an annulus which swallows the component and otherwise follows the curve:");
                textptr --;
                textptr += sprintf(textptr,"\r%s.",p);
                textptr --;
                break;
            case V2_ANNULUS_EXISTS:
                p = *SUR[WhichInput][0];
                textptr += sprintf(textptr,"\r\rNote that there exists an annulus which swallows vertice(s) { %c",*p);
                textptr --;
                p++;
                while((x = *p) != '@')
                    {
                    textptr += sprintf(textptr,",%c",x);
                    textptr --;
                    p++;
                    }
                textptr += sprintf(textptr," }");    
                textptr --;
                p++;        
                textptr += sprintf(textptr,"\rand otherwise follows the curve:");
                textptr --;
                textptr += sprintf(textptr,"\r%s.",p);
                textptr --;
                break;                
            }    
        }
    
    *textptr = EOS;
    textptr++;
    HUnlock(HText);
    SetHandleSize(HText,textptr - *HText);
    GetPort(&SavePort);
    PrintText(HText,GetHandleSize(HText) - 1,SavePort,StringWidth("\pmmmm"));
    DisposeHandle(HText);
    return(NO_ERROR);                    
#endif
}

Print_Presentation(m)
int m;
{
#ifndef MAC
  printf("Function Print_Presentation has been disabled.\n\r");
#else
    int             i,
                    j;
                    
    unsigned char     *p;
                    
    char            **HText,
                    *textptr;
                    
    long            length = 0L;                
                    
    GrafPtr            SavePort;                                    
    
    for(i = 1; i <= NumRelators; i++)
        length += GetHandleSize((char **) Relators[i]) + GetHandleSize((char **) Relators[i]) / 25;
    length += 20*NumRelators + 200;
        
    HText = (char **) NewHandle(length);
    if((textptr = *HText) == NULL) return(TOO_LONG);
    HLock(HText);
    
    textptr += sprintf(textptr,"\r\rPRESENTATION  %d",m);
    textptr --;
    
    textptr += sprintf(textptr,"\rGENERATORS %d, RELATORS %d, LENGTH  %ld\r",
    NumGenerators,NumRelators,Length);
    textptr --;
            
    for(i = 1; i <= NumRelators; i++)    
        {
        HLock((char **) Relators[i]);         
        textptr += sprintf(textptr,"\rR%3u)   ",i);
        textptr --;
        j = 8;
        p = *Relators[i];
        while(*textptr++ = *p++)
            {
            if(j++ >= 80 && *p)
                {
                *textptr ++ = '\r';
                *textptr ++ = '\t';
                *textptr ++ = '\t';
                j = 8;
                }
            }    
        HUnlock((char **) Relators[i]);                                                        
        }
            
    *textptr ++ = '\r';
    HUnlock(HText);
    SetHandleSize(HText,textptr - *HText);
    GetPort(&SavePort);
    PrintText(HText,GetHandleSize(HText) - 1,SavePort,StringWidth("\pmmmm"));
    DisposeHandle(HText);
    return(NO_ERROR);                    
#endif
} 
    
void Print_Bdry_Comp_Info(void)
{
    register unsigned char    x;
    
    int                        h,
                            i,
                            j,
                            k,
                            n,
                            ParallelRel;

    register unsigned char    *p;

    FILE                    *fptr;
    
    fptr = stdout;
    
    RERUN_AND_SAVE:
    
    fprintf(fptr,"\n\n                    DATA ABOUT THE BOUNDARY COMPONENTS OF DIAGRAM %d:",
        WhichInput + 1);
    
    if(fptr == myout)
        Get_Bdry_Comps(TRUE,TRUE,WhichInput);
    else
        Get_Bdry_Comps(TRUE,FALSE,WhichInput);
    
    ParallelRel = 0;    
    for(i = 0; i <= NumGenerators && BCWG[i] < BDRY_UNKNOWN; i++) if(BCWG[i])
        {
        for(j = 1; j <= NumBdryComps; j++) if(GBC[j] == i)
            {
            if(GBC[j] == 0)
                {
                for(k = 1,n = 0,h = 0; k <= NumFaces; k++) if(BCF[k] == j) n++;
                if(n == 0)
                    {
                    ParallelRel ++;
                    continue;
                    }
                }
            fprintf(fptr,"\n\nThe following faces 'form' a boundary component of genus %d.\n",i);
            for(k = 1,n = 1,h = 0; k <= NumFaces; k++) if(BCF[k] == j)
                {
                if(h >= 80)
                    {
                    fprintf(fptr,"\n");
                    h = 0;
                    }
                h += fprintf(fptr," %d) ",n);
                n++;
                p = Face[k];
                while((x = *p++) < VERTICES)
                    {
                    if(x & 1)
                        x = x/2 + 97;
                    else
                        x = x/2 + 65;    
                    h += fprintf(fptr,"%c",x);
                    }
                }
            }
        if(i == 0) switch(ParallelRel)
            {
            case 0:
                break;
            case 1:
                fprintf(fptr,"\n\nThe diagram also has a pair of 'parallel' relators which form a boundary component of genus 0.");            
                break;
            default:    
                fprintf(fptr,"\n\nThe diagram also has %d pairs of 'parallel' relators which form %d boundary components of genus 0.",
                    ParallelRel,ParallelRel);
                break;    
            }
        }
    
    if(fptr == stdout)
        {
        printf("\n\n    SAVE A COPY OF THIS DATA IN 'Heegaard_Results' ?  HIT 'y' OR 'n'.");
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
    fptr = stdout;
}    

void MyDrawString(char *p)
{
#ifndef MAC
  printf("Function MyDrawString has been disabled.\n\r");
#else
    register char    *q;
                    
    short int        len;

    q = p;
    while (*q++) ;
    len = q - p;
    DrawText(p,0,len);
#endif
}
