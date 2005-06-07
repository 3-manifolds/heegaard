#include "Heegaard.h"
#include "Heegaard_Dec.h"

#define        SQUEEZED    1

static unsigned char    *Cycle,
                        *M[MAXNUMRELATORS],
                        *N[MAXNUMRELATORS],
                        *NewSinkList,
                        *NewVertexList[VERTICES],
                        *NewNumber,
                        *Path,
                        *P_VertexList[VERTICES],
                        *RowSum,
                        *SinkList,
                        *Temp20,
                        *VertexList[VERTICES];

extern unsigned long    Num_Level_Transformations;    

Init_Find_Level_Transformations(Print)
int Print;
{
    /******************************************************************************************
        This routine should be called before Find_Level_Transformations() is called. It
        checks whether the presentation of interest has minimal length and has a connected
        Whitehead graph, and if both of these are true it gets some memory for Find_Level_
        Transformations() to use.
    ******************************************************************************************/
        
    int                        i,
                            j;
            
    register unsigned char     *p,
                            *q;            
    
    /******************************************************************************************
        Copy the presentation which we want level transformations for into Relators[].
    ******************************************************************************************/
        
    for(i = 1; i <= NumRelators; i++)
        {
        ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
        if((q = *Relators[i]) == NULL) return(TOO_LONG);            
        p = *SUR[ReadPres][i];
        while(*q++ = *p++) ;
        }
    Length = SURL[ReadPres];
        
    /******************************************************************************************
                    Check whether this presentation has minimal length.
        If this presentation does not have minimal length, we could still look 
        for level transformations but it does not seem very interesting to do so.            
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
            printf("\n\nThe program will only find level transformations for presentations ");
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
            printf("\n\nThe program will only find level transformations for presentations ");
            printf("with connected graphs!");
            }
        return(2);
        }                                                
    
    /******************************************************************************************
                    Get some memory for Find_Level_Transformations() to use.    
    ******************************************************************************************/
    
    for(i = 0; i < Vertices; i++)
        {
        if((M[i] = (unsigned char *) NewPtr(Vertices)) == NULL)
            {
            j = 4*i;
            goto END;
            }
        if((N[i] = (unsigned char *) NewPtr(Vertices)) == NULL)
            {
            j = 4*i + 1;
            goto END;
            }
        if((VertexList[i] = (unsigned char *) NewPtr(Vertices)) == NULL)
            {
            j = 4*i + 2;
            goto END;
            }
        if((NewVertexList[i] = (unsigned char *) NewPtr(Vertices)) == NULL)
             {
            j = 4*i + 3;
            goto END;
            }
        }
    if((Cycle = (unsigned char *) NewPtr(Vertices + 1)) == NULL)
        {
        j = 4*Vertices;
        goto END;
        }
    if((NewSinkList = (unsigned char *) NewPtr(Vertices)) == NULL)
        {
        j = 4*Vertices + 1;
        goto END;
        }
    if((NewNumber = (unsigned char *) NewPtr(Vertices)) == NULL)
        {
        j = 4*Vertices + 2;
        goto END;
        }
    if((Path = (unsigned char *) NewPtr(Vertices + 1)) == NULL)
        {
        j = 4*Vertices + 3;
        goto END;
        }
    if((RowSum = (unsigned char *) NewPtr(Vertices)) == NULL)
        {
        j = 4*Vertices + 4;
        goto END;
        }
    if((SinkList = (unsigned char *) NewPtr(Vertices)) == NULL)
        {
        j = 4*Vertices + 5;
        goto END;
        }
    j = 1000;

END:
    if(j < (Vertices << 2) + 6)
        {
        Free_Memory_For_Find_Level_Transformations(FALSE,j);
        printf("\n\nOut of memory. Sorry!");
        return(TOO_LONG);
        }    
    Num_Level_Transformations = 0L;    
    return(0);
}

Find_Level_Transformations(F1,F2)
int     F1,
        F2;
{
    unsigned char            InPath[VERTICES],
                            *P_VertexList[VERTICES],
                            *Temp20,
                            *UnAvailable;
    
    register unsigned char    *p,
                            *q;
    
    register int            i,
                            j,
                            k;
                    
    unsigned int             *r;
                
    int                     CycleLength,
                            FirstVertex,
                            h,
                            LastVertex,
                            Minimal,
                            NewSinks,
                            NumSinks,
                            NumVertices,
                            PathLength,
                            Source,
                            Sink,
                            Temp2,
                            TotalSinks;
    
    Minimal = TRUE;
    
    /******************************************************************************************
        Copy the presentation which we want level transformations for into Relators[].
    ******************************************************************************************/
        
    for(i = 1; i <= NumRelators; i++)
        {
        ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
        if((q = *Relators[i]) == NULL) return(5);            
        p = *SUR[ReadPres][i];
        while(*q++ = *p++) ;
        }
        
    /******************************************************************************************
        Call Fill_A(NumRelators) and ComputeValences_A() to initialize the arrays A[][],
        AJ3[][] and VA[].
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
            Find the level T-transformations corresponding to each generator in turn.
    ******************************************************************************************/
    
    for(Source = 0; Source < Vertices; Source += 2)
_BEGIN:    
        {        
        Sink = Source + 1;
        FirstVertex = Source;
        switch(Find_Level_Flow(Source))
            {
            case 0:
                break;
            case 1:
                if(Minimal)
                    {
                    printf("\n\nPresentation %d does not have minimal length.",
                    ReadPres + 1);
                    fprintf(myout,"\nPresentation %d does not have minimal length.",ReadPres + 1);
                    Minimal = FALSE;
                    }
                printf("\nThere is an automorphism which reduces the number of appearances of generator %c.",
                65 + Source/2);
                fprintf(myout,"\nThere is an automorphism which reduces the number of appearances of generator %c.",65 + Source/2);
                break;
            case 2:
                continue;
                break;
            }                
            
        /**************************************************************************************
                                     Initialize RowSum[].
        **************************************************************************************/
        
        for(i = 0; i < Vertices; i++) RowSum[i] = 0;
        
        /**************************************************************************************
                                    Initialize the VertexLists.
        **************************************************************************************/
        
        for(i = 0; i < Vertices; i++)
        for(j = 0; j < Vertices; j++) VertexList[i][j] = FALSE;
        for(i = 0; i < Vertices; i++) VertexList[i][i] = TRUE;
        
        /**************************************************************************************
                            Initialize M[][] and RowSum[]. Then restore A[][].
        **************************************************************************************/    
            
        for(i = 0; i < Vertices; i++)
        for(j = 0; j < Vertices; j++)
            {
            if(A[j][i])
                {
                M[i][j] = TRUE;
                RowSum[i] ++;
                }
            else
                M[i][j] = FALSE;
            }
        for(i = 1; i < Vertices; i++)
        for(r = AJ3[i]; (j = *r) < i; r++)
            A[i][j] = A[j][i] = (A[i][j] + A[j][i]) >> 1;
                
        /**************************************************************************************    
            The matrix M[][] specifies a directed graph determined by the flow from source
               to sink. We next want to recursively delete sinks from M[][]. We also want to
               renumber the vertices in the order in which they become sinks. We do this as
               follows. First, give the Sink Newnumber 0, and delete the sink from the graph.
               Then, if any new vertices become sinks, give them consecutive Newnumbers. Then
               delete these sinks one-by-one from the graph and if any new sinks are formed give
               them consecutive Newnumbers. Continue until this process stops creating new sinks.
        **************************************************************************************/
              
        NewSinks         = 1;
        TotalSinks         = 1;
        NewSinkList[0]     = Sink;
        NewNumber[Sink] = 0;
        NumVertices     = Vertices;
        
        while(1)
            {
            while(NewSinks)
                {
                NumSinks     = NewSinks;
                Temp20         = SinkList;
                SinkList     = NewSinkList;
                NewSinkList = Temp20;    
                for(k = NewSinks = 0; k < NumSinks; k++)
                    {
                    for(j = 0,i = SinkList[k]; j < NumVertices; j++) if(M[j][i])
                        {
                        RowSum[j] --;
                        if(RowSum[j] == 0)
                            {
                            NewSinkList[NewSinks ++] = j;
                            NewNumber[j] = TotalSinks;
                            TotalSinks ++;
                            if(TotalSinks >= NumVertices) goto _RELABEL_VERTEX_LISTS;
                            }
                        }
                    }        
                }
                
            /**********************************************************************************
                 If TotalSinks = NumVertices, then the directed graph determined by the flow is
                acyclic. This is the desired condition. If TotalSinks is less than NumVertices,
                then we go looking for cycles in the graph. Our aim is to delete all of the
                cycles in the graph by successively identifying all of the vertices in a cycle
                to a single vertex.     
                    At this point, the graph has no sinks. We now look for a cycle. Any vertex,
                which is not a sink, can serve as the initial vertex of a path which will
                yield a cycle.                                                                
            **********************************************************************************/
        
            Path[0]      = LastVertex = FirstVertex;
            PathLength     = 1;
               for(j = 0; j < NumVertices; j++) InPath[j] = EOS;
               InPath[FirstVertex] = 1;
               do
                   {
                   for(j = 0; j < NumVertices; j++) if(M[LastVertex][j] && RowSum[j])
                       {
                       Path[PathLength ++] = j;
                       LastVertex = j;
                       break;
                       }
                InPath[LastVertex] ++;
                   }
               while(InPath[LastVertex] < 2);

               for(j = 0; j < PathLength - 1; j++) if(Path[j] == LastVertex)
                   {
                   CycleLength = PathLength - j - 1;
                   for(k = j, i = 0; k < PathLength - 1; k++, i++) Cycle[i] = Path[k];
                   break;
                   }
                                  
               /**********************************************************************************
                   If NumVertices < CycleLength + 3, then there are no non-trivial level
                   transformations corresponding to the current choice of generator. If there
                   are generators which have not been checked go to the beginning, otherwise
                   we are done.
               **********************************************************************************/
                   
               if(NumVertices < CycleLength + 3)
                   {
                   Source += 2;
                   if(Source < Vertices) goto _BEGIN;
                   goto _END;
                   } 
               
               /**********************************************************************************
                    Next, we want to delete this cycle from the graph and we want to simplify
                   the graph. We do this by identifying all of the vertices in this cycle to
                   a single vertex. As we do this, we swap the row and column of one of the
                   vertices being identified with the last row and column of M[][], being careful
                   to update all of the relevant data, and then we decrement the number of rows
                   and columns of M[][]. 
               **********************************************************************************/
               
               CycleLength --;           
               while(CycleLength)
                   {
                   j = Cycle[CycleLength];
                   k = NumVertices - 1;
                   if(j != k)
                       {
                       /**************************************************************************
                                     If FirstVertex = k, we need to update FirstVertex.        
                       **************************************************************************/
                       
                       if(FirstVertex == k) FirstVertex = j;
                       
                       /**************************************************************************
                                   Swap the jth entry and the kth entry of NewNumber[].            
                       **************************************************************************/
                       
                       Temp2             = NewNumber[j];
                       NewNumber[j]     = NewNumber[k];
                       NewNumber[k]     = Temp2;
                       
                       /**************************************************************************
                                Swap the jth column of M[][] with the last column of M[][].
                       **************************************************************************/
                       
                       for(i = 0; i < NumVertices; i++)
                           {
                           Temp2     = M[i][j];
                           M[i][j] = M[i][k];
                           M[i][k] = Temp2;
                           }
                           
                       /**************************************************************************
                                    Swap the jth row of M[][] with the last row of M[][]. 
                       **************************************************************************/
                           
                       Temp20     = M[j];
                       M[j]    = M[k];
                       M[k]     = Temp20;
                       
                       /**************************************************************************
                                   Swap the RowSums of row j and the last row of M[][].
                       **************************************************************************/
                       
                       Temp2         = RowSum[j];
                       RowSum[j]     = RowSum[k];
                       RowSum[k]     = Temp2;
                       
                       /**************************************************************************
                                   Swap the VertexLists corresponding to j and k.
                       **************************************************************************/
                       
                       Temp20             = VertexList[j];
                       VertexList[j]     = VertexList[k];
                       VertexList[k]     = Temp20;
                       
                       /**************************************************************************
                                   If necessary swap entries in Cycle[].
                       **************************************************************************/
                       
                       for(i = 0; i <= CycleLength; i++)
                           {
                           if(Cycle[i] == j)
                               {
                               Cycle[i] = k;
                               continue;
                               }
                           if(Cycle[i] == k) Cycle[i] = j;
                           }        
                       }
                       
                   NumVertices --;
                   j = Cycle[CycleLength - 1];
                   
                   /******************************************************************************
                           Update the entries of column j and row j of M[][] together with
                           RowSum[j].                                                            
                   ******************************************************************************/
                       
                   for(i = 0; i < NumVertices; i++)
                       {
                       if(M[i][j] && M[i][k]) RowSum[i] --;
                       if(!M[i][j] && M[i][k]) M[i][j] = TRUE;
                       }
                   if(M[j][j])
                       {
                       M[j][j] = FALSE;
                       RowSum[j] --;
                       }
                   for(i = 0; i < NumVertices; i++) if(i != j && !M[j][i] && M[k][i])
                       {
                       M[j][i] = TRUE;
                       if(RowSum[i]) RowSum[j] ++;
                       }
                   
                   /******************************************************************************
                                       Add VertexList[k] to VertexList[j].
                       (Note that this step handles the possibility that Vertex k may itself be
                       the union of vertices in a number of cycles.)                
                   ******************************************************************************/
                   
                   for(i = 0; i < Vertices; i++) VertexList[j][i] += VertexList[k][i];
                   
                   CycleLength --;                    
                   }
               
               /**********************************************************************************
                   After all of the vertices in a cycle have been identified to a single vertex,
                   namely the vertex Cycle[0], we want to check whether the vertex Cycle[0] has
                   become a new sink. This will be true iff RowSum[Cycle[0]] is now 0.
               **********************************************************************************/
                        
               if(RowSum[Cycle[0]] == 0)
                   {
                   /******************************************************************************
                       The vertex Cycle[0] has become a sink. Set new sinks to 1, and set things
                       up for the loop which looks for sinks and renumbers vertices at the top of
                       this while loop we are in.
                   ******************************************************************************/
                       
                   NewSinks             = 1;
                   NewSinkList[0]         = Cycle[0];
                   NewNumber[Cycle[0]] = TotalSinks;
                   TotalSinks ++;
                   if(TotalSinks >= NumVertices) goto _RELABEL_VERTEX_LISTS;
                   }                            
               }

_RELABEL_VERTEX_LISTS:
            
        /**************************************************************************************
                Relabel the VertexLists[] to reflect the relabeling of the vertices.             
        **************************************************************************************/
                        
        for(i = 0; i < NumVertices; i++) P_VertexList[NewNumber[i]] = VertexList[i];    
                
        /**************************************************************************************
                Use M[][] to create N[][] according to the permutation of the vertices
                induced by relabeling. 
                (The array N[][] corresponds to an acyclic graph obtained from M[][] by
                contracting cycles in the graph correspond to M[][] to a single vertex 
                and also renumbering vertices in the order in which they become sinks.)                                                            
        **************************************************************************************/    
        
        for(i = 0; i < NumVertices; i++)
        for(j = 0; j < NumVertices; j++)
            {
            if(M[i][j])
                N[NewNumber[i]][NewNumber[j]] = TRUE;
            else
                N[NewNumber[i]][NewNumber[j]] = FALSE;
            }
                                        
        /**************************************************************************************
            Set up the NewVertexLists. NewVertexList[i] lists all of the relabeled vertices
            of N[][] that are ancestors of vertex i. 
                (A vertex x is an ancestor of vertex y if x = y, or there is a directed path 
            from vertex x to vertex y.)                            
        **************************************************************************************/
         
        for(i = 0; i < NumVertices; i++)
        for(j = 0; j < NumVertices; j++) NewVertexList[i][j] = FALSE;
        for(i = 0; i < NumVertices; i++)
             {
             NewVertexList[i][i] = TRUE;
             NewSinks             = 1;
            NewSinkList[0]         = i;
            while(NewSinks)
                {
                NumSinks     = NewSinks;
                Temp20         = SinkList;
                SinkList     = NewSinkList;
                NewSinkList = Temp20;    
                for(k = NewSinks = 0; k < NumSinks; k++)
                    {
                    for(j = 0,h = SinkList[k]; j < NumVertices; j++)
                        if(N[j][h] && !NewVertexList[i][j])
                            {
                            NewSinkList[NewSinks ++] = j;
                            NewVertexList[i][j] = TRUE;
                            }
                    }        
                }
             }
         
         /*************************************************************************************
             Now we are ready to read off and output the level-transformations. First, we 
             initialize the UnAvailable list, PathLength and Path[0].
                 Any partition of the vertices of our graph which separates Source and Sink, 
             can be used to determine a T-transformation. In our case, the vertices are 
             represented by the indices of the UnAvailable list, and the partition of the 
             vertices is determined by whether the entries of the UnAvailable list are zero 
             or nonzero. 
                 It is not hard to show that a partition of the vertices of the graph
             will induce a level-transformation iff the set of entries in the UnAvailable list
             which are nonzero is the union of a set of vertices together with all of the
             ancestors of these vertices. NewVertexlist[] holds just such information since
             NewVertexList[i] lists all of the relabeled vertices of N[][] that are ancestors 
             of vertex i.
                 The array Path[] is used to keep track of the current NewVertexlist[]s whose
             union partitions the entries of the UnAvailable list. Path[] lists the current 
             NewVertesList[]s in monotonically increasing order, so that we can keep track of
             the sets of NewVertexLists[]s used lexicographically, and thus ensure that each
             distinct partition of the vertices producing a level-transformation occurs once
             and only once.
         *************************************************************************************/
         
         UnAvailable = Cycle;
         for(i = 0; i < NumVertices; i++) UnAvailable[i] = FALSE;
         PathLength = 0;
         Path[0]     = 0;
         while(1)
             {
             /**********************************************************************************
                 AA) Look for the first vertex with number greater than Path[PathLength] which 
                     is not on the UnAvailable list. 
                     A)     If such a vertex exists:
                         1) Set Path[PathLength] equal to this vertex then increment PathLength 
                         and set the next entry of Path[] equal to the number of this vertex. 
                         (This ensures that the vertices that appear in the list Path[] will 
                         always appear in strict monotonically increasing order.)
                         2) Output a level-transformation corresponding to the vertices on the 
                         Unavailable list.
                     B)     If no such vertex exists:
                         1)     Decrement PathLength by one, update the Unavailable list and goto 
                             step AA) above, or
                         2)     If PathLength becomes negative, stop, we are done listing 
                             level-transformations for the current generator, goto the next
                             generator.
             **********************************************************************************/
                 
             for(i = Path[PathLength] + 1; i < NumVertices && UnAvailable[i]; i++) ;
             if(i < NumVertices)
                 {
                 Path[PathLength] = i;
                 PathLength ++;
                 Path[PathLength] = i;
                 
                 /******************************************************************************
                         Increment the number of appearances of vertex i and all of its
                         ancestors on the UnAvailable list. We can do this by adding a copy
                         of the entries of NewVertexList[i][j] to UnAvailable[j] since 
                         NewVertexList[i][j] is 1 iff vertex j is an ancestor of vertex i and
                         otherwise NewVertexList[i][j] is 0.
                 ******************************************************************************/
                 
                 for(j = i,p = NewVertexList[i]; j < NumVertices; j++)
                     UnAvailable[j] += p[j];
                     
                 /******************************************************************************
                                     Setup the array ZZ[] for Do_Aut().
                 ******************************************************************************/
                 
                 for(i = 0; i < Vertices; i++) ZZ[i] = FALSE;
                 for(i = 0; i < NumVertices; i++) if(!UnAvailable[i])
                     {
                     for(j = 0,p = P_VertexList[i]; j < Vertices; j++) if(p[j])
                         ZZ[j] = TRUE;
                     }    
                 
                 /******************************************************************************
                     Count the number of nonzero entries of ZZ[]. If this number equals 1 or
                     Vertices - 1, then the induced level-transformation is the identity.
                 ******************************************************************************/
                     
                 for(i = j = 0; i < Vertices; i++) if(ZZ[i]) j++;
                 
                 if((j > 1) && (j < Vertices - 1))
                     {
                         
                     /**************************************************************************
                         Call Do_Aut() to perform the level-transformation specified in ZZ[].
                         Increment Num_Level_Transformations.
                     **************************************************************************/
                     
                     if(Micro_Print)
                         {
                         printf("\n\nPerformed the following level-transformation:");
                         if(Micro_Print_F)
                              fprintf(myout,"\n\nPerformed the following level-transformation:");
                         }
                 
                     if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);
                     Num_Level_Transformations ++;
                     
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
                     
                     /**************************************************************************
                         Call In_File() or On_File() to see whether this new presentation is
                         already on file. If not, save the new presentation.
                     **************************************************************************/
                     
                     Length = SURL[ReadPres];
                     if(F2)
                         {
                         Length = SURL[ReadPres];
                         j = In_File();
                         if(j == TOO_LONG) return(5);
                         if(j == NumFilled - 1)
                             {
                             if(NumFilled >= MAX_SAVED_PRES - 3) return(5);
                             if(BytesUsed > BytesAvailable)
                                 {
                                 if(UserSaidQuit) return(5);
                                 if(UserSaidQuit = User_Says_Quit()) return(5);
                                 }
                             }
                         }
                     else
                         {            
                         j = On_File();
                         if(j == NumFilled)
                             {
                             if(BytesUsed > BytesAvailable)
                                 {
                                 if(UserSaidQuit) return(5);
                                 if(UserSaidQuit = User_Says_Quit()) return(5);
                                 }
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
                             if(NumFilled >= MAX_SAVED_PRES - 3) return(5);
                             if(F1)
                                 {
                                 for(i = 1; i <= NumRelators; i++)
                                     {
                                     if(GetHandleSize((char **) Relators[i]) <= 3L) return(3);
                                     }
                                 }
                             }
                         }    
                 
                     /**************************************************************************
                                     Restore the original presentation.
                         This is necessary because Do_Aut changes the presentation in
                         Relators[] into its image under the level-transformation and we do
                         not want to compose level-transformations at this point. In particular,
                         the data in the arrays used above is invalid unless Relators[] has been
                         restored.            
                     **************************************************************************/
                     
                     for(i = 1; i <= NumRelators; i++)
                        {
                        ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
                        if((q = *Relators[i]) == NULL) return(5);                
                        p = *SUR[ReadPres][i];
                        while(*q++ = *p++) ;
                        }
                     }
                 }
             else
                 {
                 /******************************************************************************
                     Decrement PathLength and then decrement the number of times vertex
                     Path[PathLength] and its ancestors appear on the UnAvailable list.
                     When PathLength becomes negative we have found all of the level
                     transformations corresponding to the current choice of source and sink.                            
                 ******************************************************************************/
                     
                 if(--PathLength < 0) break;
                 i = Path[PathLength];
                 for(j =    i,p = NewVertexList[i]; j < NumVertices; j++)
                     UnAvailable[j] -= p[j];    
                 }    
             }                            
        }    

_END:
    return(0);            
}

Find_Level_Flow(Source)
unsigned int     Source;
{
    /******************************************************************************************
        This routine is a variant of Find_Flow_A(). Unlike Find_Flow_A(), this
        routine only finds a maximal flow in the Whitehead graph between a single pair of
        vertices, namely the source and sink passed to it.
    ******************************************************************************************/
            
    register unsigned int     i,
                            j,
                            max,
                            min,
                            *p,
                            *q,
                            *r;
                            
    unsigned int             Flow,
                            k,
                            MaxFlow,
                            S[VERTICES],
                            Sink;
                
    Sink = Source + 1;
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

Free_Memory_For_Find_Level_Transformations(Print,NumPtrs)
int     Print,
        NumPtrs;
{
    /******************************************************************************************
        This routine should be called when Find_Level_Transformations has terminated or been
        interrupted by the user. Its major job is to free the memory allocated for
        Find_Level_Transformations().
    ******************************************************************************************/
        
    int i;
    
    if(Print)
        {
        printf("\n\nThe program found %lu level transformations(s).",
            Num_Level_Transformations);
        fprintf(myout,"\n\nThe program found %lu level transformations(s).",
            Num_Level_Transformations);
        }

    /******************************************************************************************
             Dispose of the memory that we allocated for Find_Level_Transformations().
    ******************************************************************************************/
            
    for(i = 0; i < Vertices && NumPtrs > 0; i++)
        {
        if(NumPtrs-- > 0)
            DisposePtr((char *) M[i]);
        if(NumPtrs-- > 0)
            DisposePtr((char *) N[i]);
        if(NumPtrs-- > 0)
            DisposePtr((char *) VertexList[i]);
        if(NumPtrs-- > 0)
            DisposePtr((char *) NewVertexList[i]);
        }
    if(NumPtrs-- > 0)    
        DisposePtr((char *) Cycle);
    if(NumPtrs-- > 0)
        DisposePtr((char *) NewSinkList);
    if(NumPtrs-- > 0)
        DisposePtr((char *) NewNumber);
    if(NumPtrs-- > 0)
        DisposePtr((char *) Path);
    if(NumPtrs-- > 0)
        DisposePtr((char *) RowSum);
    if(NumPtrs-- > 0)
        DisposePtr((char *) SinkList);
}

/* This seems to be a duplicate -- MC */
#ifdef MAC
Connected_(i,k)
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

#endif
       
Test_LT_For_Pseudo_Min()
{
    /******************************************************************************************
        This is a rather crude routine which takes the set of presentations in the orbit of a
        presentation under level transformations, and checks them for pairs of separating
        vertices, reducibility of their dual diagrams and pseudo-minimality. It does
        essentially no error checking!!!!
    ******************************************************************************************/
        
    unsigned char     *p,
                    *q,
                    **Temp,
                    AA[MAX_SAVED_PRES];
    
    int                Realizable,
                    SNumGenerators,
                    SNumRelators;
                    
    unsigned int    Flag,
                    h,
                    i,
                    j,
                    k,
                    Num_Pseudo_Min,
                    Num_Reducible,
                    Num_Sep_Vert,
                    Num_Standard;
                    
    unsigned long    SLength;                
                    
    unsigned int Whitehead_Graph();
    
    DrawingDiagrams = TRUE;
    TestRealizability1 = TRUE;
    Realizable = -1;
    Num_Pseudo_Min = Num_Reducible = Num_Sep_Vert = Num_Standard = 0;
    printf("\n");
        
    for(ReadPres = 0; ReadPres < NumFilled; ReadPres++)
        {
        putc('*',stdout);
        
        /**************************************************************************************
                    Copy the presentation that we want to test into Relators[].
        **************************************************************************************/
        
        NumGenerators = NG[ReadPres];
        NumRelators = NR[ReadPres];
        Vertices = 2*NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
            if((q = *Relators[i]) == NULL) goto NEXT_PRES;            
            p = *SUR[ReadPres][i];
            while(*q++ = *p++) ;
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
            if(SQUEEZED)
                {
                SNumGenerators = NumGenerators;
                SNumRelators = NumRelators;
                NumRelators = NumGenerators;
                NumGenerators = SNumRelators;
                Vertices = 2*NumGenerators;
                }
            if(NumGenerators == NumRelators || SQUEEZED)
                {
                SLength = SURL[ReadPres];
                for(i = 1; i <= NumRelators; i++)
                    { 
                    Temp                 = Relators[i];
                    Relators[i]         = DualRelators[i];
                    DualRelators[i]     = Temp;    
                    }
                if(Freely_Reduce() == TOO_LONG) continue;
                Length = OrigLength;
                
            /* 3/27/95 Changed following code to check whether dual presenations are standard. */

            /*    if(Length < SLength)
                    {
                    Num_Reducible ++;
                    AA[ReadPres] = 3;
                    }
                else
                    {        
                    switch(Find_Flow_A(NORMAL,TRUE))
                        {
                        case 0:    
                            Num_Pseudo_Min ++;
                            AA[ReadPres] = 2;
                            break;
                        case 2:
                            Num_Reducible ++;
                            AA[ReadPres] = 3;
                            break;
                        }
                    }    */
                    
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
                if(SQUEEZED)
                    {
                    NumRelators = SNumRelators;
                    NumGenerators = SNumGenerators;
                    Vertices = 2*NumGenerators;
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
                    AA[ReadPres] = EOS;
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
        NEXT_PRES:
        continue;            
        }

    DrawingDiagrams = FALSE;
    TestRealizability1 = FALSE;
    
    printf("\n\nTotal presentations with no pairs of separating vertices:  %u.",
        NumFilled - Num_Sep_Vert);
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
            printf("\nThe program found both realizable and unrealizable presentations in the orbit!!!??");
            break;            
        }
    
    fptr = stdout;
    PRINT_NO_SEP_VERT:
    if(Num_Sep_Vert < NumFilled)
        {
        j = k = 0;
        j += fprintf(fptr,"\n\nPresentations with no pairs of separating vertices:    ");
        for(h = 0; h < NumFilled; h++) if(AA[h])
            {
            if(++k < NumFilled - Num_Sep_Vert)
                j += fprintf(fptr,"{%3d,",h+1);
            else
                j += fprintf(fptr,"{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i])
            {
            if(++k < NumFilled - Num_Sep_Vert)
                j += fprintf(fptr,"%3d,",i+1);
            else
                j += fprintf(fptr,"%3d}.",i+1);    
            if(j > 83)
                {
                j = 0;
                fprintf(fptr,"\n");
                }
            }
        if(fptr == stdout)
            {
            GET_RESPONSE1:
            printf("\n\nSAVE THIS LIST IN THE FILE 'Heegaard_Results' ?     HIT 'y' OR 'n'.");
            switch(WaitkbHit())
                {
                case 'y':
                    fptr = myout;
                    goto PRINT_NO_SEP_VERT;
                case 'n':
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE1;    
                }
            }                            
        }    
    
    fptr = stdout;
    PRINT_PSEUDO_MINIMAL:
    if(Num_Pseudo_Min)
        {
        j = k = 0;
        j += fprintf(fptr,"\n\nPresentations which are pseudo-minimal:    ");
        for(h = 0; h < NumFilled; h++) if(AA[h] == 2)
            {
            if(++k < Num_Pseudo_Min)
                j += fprintf(fptr,"{%3d,",h+1);
            else
                j += fprintf(fptr,"{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i] == 2)
            {
            if(++k < Num_Pseudo_Min)
                j += fprintf(fptr,"%3d,",i+1);
            else
                j += fprintf(fptr,"%3d}.",i+1);    
            if(j > 83)
                {
                j = 0;
                fprintf(fptr,"\n");
                }
            }
        if(fptr == stdout)
            {
            GET_RESPONSE2:
            printf("\n\nSAVE THIS LIST IN THE FILE 'Heegaard_Results' ?     HIT 'y' OR 'n'.");
            switch(WaitkbHit())
                {
                case 'y':
                    fptr = myout;
                    goto PRINT_PSEUDO_MINIMAL;
                case 'n':
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE2;    
                }
            }                            
        }

    fptr = stdout;
    PRINT_REDUCIBLE:
    if(Num_Reducible)
        {
        j = k = 0;
        j += fprintf(fptr,"\n\nPresentations which are reducible:     ");
        for(h = 0; h < NumFilled; h++) if(AA[h] >= 3)
            {
            if(++k < Num_Reducible)
                j += fprintf(fptr,"{%3d,",h+1);
            else
                j += fprintf(fptr,"{%3d}.",h+1);    
            break;
            }
        for(i = h + 1; i < NumFilled; i++) if(AA[i] >= 3)
            {
            if(++k < Num_Reducible)
                j += fprintf(fptr,"%3d,",i+1);
            else
                j += fprintf(fptr,"%3d}.",i+1);
            if(j > 83)
                {
                j = 0;
                fprintf(fptr,"\n");
                }
            }
        if(Num_Standard)
            {
            j = k = 0;
            j += fprintf(fptr,"\n\nPresentations whose duals are Standard:     ");
            for(h = 0; h < NumFilled; h++) if(AA[h] == 4)
                {
                if(++k < Num_Standard)
                    j += fprintf(fptr,"{%3d,",h+1);
                else
                    j += fprintf(fptr,"{%3d}.",h+1);    
                break;
                }
            for(i = h + 1; i < NumFilled; i++) if(AA[i] == 4)
                {
                if(++k < Num_Standard)
                    j += fprintf(fptr,"%3d,",i+1);
                else
                    j += fprintf(fptr,"%3d}.",i+1);
                if(j > 83)
                    {
                    j = 0;
                    fprintf(fptr,"\n");
                    }
                }
            }
        if(fptr == stdout)
            {
            GET_RESPONSE3:
            printf("\n\nSAVE THIS LIST IN THE FILE 'Heegaard_Results' ?     HIT 'y' OR 'n'.");
            switch(WaitkbHit())
                {
                case 'y':
                    fptr = myout;
                    goto PRINT_REDUCIBLE;
                case 'n':
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE3;    
                }
            }                    
        }
    
    printf("\n\nPLEASE SELECT ONE OF THE FOLLOWING ALTERNATIVES."); 
    printf("\n    0) Save no   'level' presentations.");
    printf("\n    1) Save only 'level' presentations which have no pairs of separating vertices.");
    printf("\n    2) Save only 'level' presentations which are pseudo-minimal."); 
    printf("\n    3) Save only 'level' presentations which are reducible.");
    printf("\n    4) Save all  'level' presentations.");
    printf("\nAny presentations saved will appear in 'Heegaard_Results'. ");
    printf("HIT 0,1,2,3,OR 4.        ");
    GET_RESPONSE4:
    switch(WaitkbHit())
        {
        case '0':
            break;
        case '1':
            if(Num_Sep_Vert < NumFilled)
                Report(0,0,0,0,1,1,0,1,1,AA);
            break;
        case '2':
            if(Num_Pseudo_Min)
                Report(0,0,0,0,1,1,0,2,1,AA);
            break;
        case '3':
            if(Num_Reducible)
                Report(0,0,0,0,1,1,0,3,1,AA);
            break;
        case '4':
            Report(0,0,0,0,1,1,0,0,1,0);
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE4;                    
        }            
}

unsigned int In_File(void)
{    
    register unsigned char    *p,
                            *q;
    
    unsigned char            *r;
                            
    int                        i,
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
                         {
                         printf("  %d -> %d.",ReadPres + 1,NumFilled + 1);
                         fprintf(myout,"  %d -> %d.",ReadPres + 1,NumFilled + 1);
                         }                                        
                    if(Save_Pres(ReadPres,0,Length,1,75,0,0,0)) return(TOO_LONG);
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
                     {
                     printf("  %d -> %d.",ReadPres + 1,Node + 1);
                     fprintf(myout,"  %d -> %d.",ReadPres + 1,Node + 1);
                     }
                return(Node);
            case -1:
                if(Right[Node] == INFINITE)
                    {
                    if(Compute_Stabilizers)
                         {
                         printf("  %d -> %d.",ReadPres + 1,NumFilled + 1);
                         fprintf(myout,"  %d -> %d.",ReadPres + 1,NumFilled + 1);
                         }                    
                    if(Save_Pres(ReadPres,0,Length,1,75,0,0,0)) return(TOO_LONG);
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
