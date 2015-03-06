#include "Heegaard.h"
#include "Heegaard_Dec.h"                                                                                                                                                                             

/****************************** function prototypes *****************************************
L   28 Diagram_Main(void)
L   85 Diagram_1(void)
L  189 Diagram_2(void)
L  280 Diagram_3(void)
L  485 Diagram_4(void)
L  733 Diagram_5(void)
L  758 Offset(void)
L 1208 OffsetSub(unsigned int i,unsigned int p,unsigned int q,unsigned int a,unsigned int b,unsigned int c)
L 1248 Diagram_6(void)
L 1287 Diagram_7(void)
L 1462 Compare(unsigned char *p)
L 1487 Compare_Str(unsigned char *S1,unsigned char *S2,unsigned long length)
L 1537 Get_Bdry_Comps(int Print,int Where,unsigned int WhichPres)
L 1664 Delete_Redundant_Relators(void)
L 1798 MG_Bdry_Comp_Data(unsigned int WhichPres)
L 1986 Sep_Surface(void)
L 2040 Fill_DRA(void)
L 2092 GCD(unsigned int p,unsigned int q)
L 2141 Inverse(register unsigned char *p)
L 2169 Non_Unique(void)
L 2221 Debug(void)		
********************************************************************************************/

unsigned int Diagram_Main()
{
    /******************************************************************************************
        Given a presentation which has been reduced to minimal length and whose reduced
        Whitehead graph does not have any separating pairs of vertices, the routines in this
        file are used to determine how the edges of the Heegaard diagram meeting each pair of
        inverse vertices must by identified.
    ******************************************************************************************/
         
    unsigned int Diagram_5();
    unsigned int Offset();
    unsigned int Valence_Two();
    unsigned int Diagram_7();
    
    if(Diagram_5()) return(TOO_LONG);
    if(Diagram_1()) return(FATAL_ERROR);    
    if(Diagram_2()) return(FATAL_ERROR);                                        
    Diagram_3();                                        
    Diagram_4();            
    switch(Offset())
        {
        case NO_ERROR: 		break;
        case NON_UNIQUE_1: 	return(NON_UNIQUE_1);
        case NON_UNIQUE_2: 	return(NON_UNIQUE_2);
        case NON_UNIQUE_3: 	return(NON_UNIQUE_3);
        case NON_UNIQUE_4: 	return(NON_UNIQUE_4);
        case FATAL_ERROR:  	return(FATAL_ERROR);
        case TOO_LONG:      return(TOO_LONG);    
        }
    Diagram_6();
    
    /******************************************************************************************
        If NGV2 is > 0, there exist generators for which both vertices in the reduced Whitehead
        graph have valence two and we have not been able to determine the proper
        identifications of the edges at these vertices. In this case, Valence_Two() is
        called to determine how the edges at such pairs of vertices should be identified.
    ******************************************************************************************/

    if(NGV2) switch(Valence_Two(FALSE))
        {
        case NO_ERROR:
            break;
        case FATAL_ERROR:
            return(FATAL_ERROR);
        case TOO_LONG:
            return(TOO_LONG);
        case V2_ANNULUS_EXISTS:
            if(TestRealizability3) return(V2_ANNULUS_EXISTS);
            if(DrawingDiagrams || TestRealizability2)
                break;
            else
                return(V2_ANNULUS_EXISTS);                
        }  
                 
    return(Diagram_7());                                    
}

int Diagram_1()
{
    /******************************************************************************************
        This routine does two things. First, it determines the absolute value of the
    exponents with which each generator appears in the relations. It leaves these exponents
    in the three columns of the ith row of the array EXP[][].
        (Note that if a presentation corresponds to a Heegaard diagram, the presentation has
    minimal length and the reduced Whitehead graph of the presentation has no pairs of
    separating vertices, then any generator that appears in the presentation can appear with
    at most three distinct exponents up to sign. Furthermore, these exponents must be pairwise
    relatively prime and of the form p,q and p+q.)
        Second, the routine counts the number of appearances of each exponent and leaves
    these counts in the array NEX[][].
    ******************************************************************************************/    
    
    register unsigned char  i,
                            *p,
                            x,
                            y;
                            
    register unsigned int   ex,
                            sex,
                            *q;
                            
    unsigned char           j;                        
                            
    unsigned long           length1,
                            length2;                            
                            
    for(i = NumGenerators; i > 0; )
    for(x = 0,i--; x < 3; x++) NEX[i][x] = EXP[i][x] = 0;
    for(i = NumRelators; i >= 1; i--) LR[i] = GetHandleSize((char **) Relators[i]) - 1;
    for(j = 1; j <= NumRelators; j++)
        {
        p = *Relators[j];
        x = *p++;
        if(!x) continue;
        ex = 1;
        while(*p == x)
            {
            ex++;
            p++;
            }
        sex = ex;
        ex = 0;
        if(*p) x = *p;
        while( (y = *p++) )
            {
            if(x == y)
                ex++;
            else
                {
                if(x < 95) i = x - 65;
                else i = x - 97;
                for(x = 0,q = EXP[i]; x < 3 && ex != q[x] && q[x]; x++) ;
                q[x] = ex;
                NEX[i][x]++;
                ex = 1;    
                }
            x = y;
            }
        y = **Relators[j];        
        if(x == y)
            {
            ex += sex;
            if(x < 95) i = x - 65;
            else i = x - 97;
            for(x = 0,q = EXP[i]; x < 3 && ex != q[x] && q[x]; x++) ;
            q[x] = ex;
            NEX[i][x]++;
            }
        else
            {
            if(x < 95) i = x - 65;
            else i = x - 97;
            for(x = 0,q = EXP[i]; x < 3 && ex != q[x] && q[x]; x++) ;
            q[x] = ex;
            NEX[i][x]++;
            ex = sex;
            x = y;
            if(x < 95) i = x - 65;
            else i = x - 97;
            for(x = 0,q = EXP[i]; x < 3 && ex != q[x] && q[x]; x++) ;
            q[x] = ex;
            NEX[i][x]++;
            }
        }
    
    if(TestRealizability1)
        {
        for(i = 0,length1 = 0L; i < NumGenerators; i++)
        for(x = 0; x < 3; x++)
        length1 += EXP[i][x]*NEX[i][x];
        for(i = 1,length2 = 0L; i <= NumRelators; i++) length2 += LR[i];
        if(length1 != length2)
            {
            if(DrawingDiagrams == FALSE)
                printf("\n\nThere is a generator which appears with more than three distinct exponents!");
            return(1);
            }
        }
    return(0);            
}                                

int Diagram_2()
{        
    /******************************************************************************************
        This routine takes the array EXP[][] as determined by Diagram_1() and permutes the 
        entries of the three columns of each row so that the exponents with which each
        generator appears are in ascending numerical order. It also applies the same 
        permutation to the corresponding columns of NEX[][] so that if an exponent appears in 
        column j then the number of times that exponent appears in the relators is found in 
        column j of NEX[][]. Permuting these columns makes things somewhat simpler for other
        routines that use EXP[][] and NEX[][]. 
    ******************************************************************************************/                 
    
    register int            i;
    
    register unsigned int   *p,
                            *q,
                            temp;
    
    unsigned int            GCD();
        
    for(i = NumGenerators - 1; i >= 0; i--) 
        {
        p = EXP[i];
        q = NEX[i];
        if(p[0] > p[1]) 
            {
            temp = p[1];
            p[1] = p[0];
            p[0] = temp;
            temp = q[1];
            q[1] = q[0];
            q[0] = temp;
            }
        if(p[1] > p[2]) 
            {
            temp = p[2];
            p[2] = p[1];
            p[1] = temp;
            temp = q[2];
            q[2] = q[1];
            q[1] = temp;
            }
        if(p[0] > p[1]) 
            {
            temp = p[1];
            p[1] = p[0];
            p[0] = temp;
            temp = q[1];
            q[1] = q[0];
            q[0] = temp;
            }
        }
        
    if(TestRealizability1) for(i = NumGenerators - 1; i >= 0; i--)
        {
        p = EXP[i];
        if(p[0])
            {
            if(GCD(p[0],p[1]) != 1)
                {
                if(DrawingDiagrams == FALSE)
                    {
                    i += 'A';
                    printf("\n\nGenerator '%c' appears with exponents which are not relatively prime!",i);
                    }
                return(1);
                }
            if(p[2] != p[0] + p[1])
                {
                if(DrawingDiagrams == FALSE)
                    {
                    printf("\n\nGenerator '%c' appears with exponents %u, %u, and %u, but %u != %u + %u.",
                        i + 'A',p[0],p[1],p[2],p[2],p[0],p[1]);
                    }
                return(1);
                }
            }
        else
            if(p[1] && GCD(p[1],p[2]) != 1)
                {
                if(DrawingDiagrams == FALSE)
                    {
                    i += 'A';
                    printf("\n\nGenerator '%c' appears with exponents which are not relatively prime!",i);
                    }
                return(1);
                }
        }
    return(0);            
}

void Diagram_3(void)
{
    /******************************************************************************************
        This routine sets some flags that appear in the array DF[] as detailed below, and also
        sets up a table T[][] of subwords of the relators.
    ******************************************************************************************/
        
    register int        i,
                        j,
                        k,
                        sj,
                        sk;
                        
    unsigned int        max,
                        smax;
    
    /******************************************************************************************
            Set the flags that appear in the array DF[]. The default value is 0. If the
        maximal exponent with which generator i appears is 2 or less and both vertex 2i
        and vertex 2i+1 have valence 2, then we set DF[i] equal to 1. Otherwise, the maximal
        exponent with which generator i appears is greater than 2 and we set DF[i] equal to 2
        unless both vertex 2i and vertex 2i+1 have valence 2, in which case, we set DF[i]
        equal to 3.            
    ******************************************************************************************/    
    
    for(i = 0; i < NumGenerators; i++)
        {
        if(EXP[i][2] > 2)
            {
            if(VWG[(i << 1)] > 2 || VWG[(i << 1) + 1] > 2)
                DF[i] = 2;
            else
                DF[i] = 3;
            }
        else
            {
            if(VWG[(i << 1)] > 2 || VWG[(i << 1) + 1] > 2)
                DF[i] = 0;
            else
                DF[i] = 1;
            }            
        }
            
    /******************************************************************************************
                        Next, set up the table of chars T[][].
    ******************************************************************************************/        
    
    for(i = 0; i < NumGenerators; i++) switch(DF[i])
        {
        case 0:
        case 1:
            /**********************************************************************************
                If DF[i] is 0 or 1, then we want to find a band of parallel edges meeting
                vertex 2i (vertex 2i+1), which contains as many edges as possible. If we then
                read through the relators and determine what vertices the paths containing
                these edges proceed to after leaving vertex 2i+1 (vertex 2i), we will have
                enough information to determine how to identify the edges at vertex 2i and
                vertex 2i+1 in case 0. When DF[i] = 1, there will be no more than two 
                possible ways to identify these edges; which possible identification is
                correct will be resolved later.
            **********************************************************************************/    
            
            max = 0;
            k = 2*i;
            sk = k;    
            j = FV[k];
            do
                {
                if(A[k][j] > max)
                    {
                    max = A[k][j];
                    sj = j;
                    }
                j = CO[k][j];
                }
            while(j != FV[k]);        
            k++;
            if(VWG[k] == 2) max = 0;
            
            /**********************************************************************************
                If one of the two inverse vertices has valence two, we want to take one of the
                bands of edges meeting this vertex as our maximal band.
            **********************************************************************************/
                
            smax = max;
            j = FV[k];
            do
                {
                if(A[k][j] > max)
                    {
                    max = A[k][j];
                    sj = j;
                    }
                j = CO[k][j];
                }
            while(j != FV[k]);    
            if(max > smax) sk = k;    
            TV[i] = sk;
            SV[i] = sj;                
            if(sj & 1)
                { 
                T[i][0] = (sj >> 1) + 65;
                T[i][2] = (sj >> 1) + 97;    
                }
            else
                {
                T[i][0] = (sj >> 1) + 97;
                T[i][2] = (sj >> 1) + 65;
                }
            if(sk & 1)
                {
                T[i][1] = (sk >> 1) + 97;
                T[i][3] = (sk >> 1) + 65;
                }
            else
                {
                T[i][3] = (sk >> 1) + 97;
                T[i][1] = (sk >> 1) + 65;
                }
            break;        
        
        case 2:
            /**********************************************************************************
                In this case, generator i appears with exponent greater than two, and either
                vertex 2i or vertex 2i+1 has valence greater than two in the reduced Whitehead
                graph. (We know at this point that the exponents with which generator i
                appears are of the form {p,q,p+q}. What we need to determine is whether they
                appear in the order (p,p+q,q) or in the order (q,p+q,p) as one proceeds around
                vertex 2i.) 
                    Note: There is an edge joining vertex 2i and vertex 2i+1 in the reduced 
                graph. Suppose that vertex 2i has valence greater than two. Then we want to
                determine which exponents appear in words corresponding to paths that come
                from the first vertex following vertex 2i+1 in counter-clockwise cyclic order
                about vertex 2i and to also determine which exponents appear in words
                corresponding to paths that come from the vertex preceeding vertex 2i+1
                in counter-clockwise cyclic order about vertex 2i. This gives us the
                information we need.
            **********************************************************************************/     
            
            k = 2*i;                                                                    
            if(VWG[k] > 2)
                {
                j = CO[k][k+1];
                if(j & 1)
                    { 
                    T[i][0] = (j >> 1) + 65;
                    T[i][3] = (j >> 1) + 97;
                    }
                else
                    {
                    T[i][0] = (j >> 1) + 97;
                    T[i][3] = (j >> 1) + 65;
                    }
                T[i][1] = (k >> 1) + 65;
                T[i][2] = (k >> 1) + 97;
                while(CO[k][j] != (k+1)) j = CO[k][j];
                if(j & 1)
                    { 
                    T[i][4] = (j >> 1) + 65;
                    T[i][7] = (j >> 1) + 97;
                    }
                else
                    {
                    T[i][4] = (j >> 1) + 97;
                    T[i][7] = (j >> 1) + 65;
                    }
                T[i][5] = (k >> 1) + 65;
                T[i][6] = (k >> 1) + 97;
                }
            else        
                {
                j = CO[k+1][k];
                if(j & 1)
                    { 
                    T[i][0] = (j >> 1) + 65;
                    T[i][3] = (j >> 1) + 97;
                    }
                else
                    {
                    T[i][0] = (j >> 1) + 97;
                    T[i][3] = (j >> 1) + 65;
                    }
                T[i][1] = (k >> 1) + 97;
                T[i][2] = (k >> 1) + 65;
                while(CO[k+1][j] != k) j = CO[k+1][j];
                if(j & 1)
                    { 
                    T[i][4] = (j >> 1) + 65;
                    T[i][7] = (j >> 1) + 97;
                    }
                else
                    {
                    T[i][4] = (j >> 1) + 97;
                    T[i][7] = (j >> 1) + 65;
                    }
                T[i][5] = (k >> 1) + 97;
                T[i][6] = (k >> 1) + 65;
                }
            break;    
        
        case 3:
            break;
        }
}

void Diagram_4(void)
{     
    /******************************************************************************************
        This routine counts the number of appearances of certain subwords of the relator
        Relators[k] of the form xy......yz where x != y and y != z. The subwords of interest
        are determined by the entries in the character table T[][].
    ******************************************************************************************/
    
    register unsigned char  i,
                            *p,
                            *q,
                            x,
                            y,
                            z;
                            
    register unsigned int   e,
                            *r;
    
    unsigned char           *Rk;
                            
    unsigned int            k;
    
    for(i = 0; i < NumGenerators; i++)
        {
        if(DF[i] < 2)
            for(z = 0,r = ED[i]; z < Vertices; z++) r[z] = 0;
        else
            for(z = 0; z < 3; z++) EXL[i][z] = EXR[i][z] = 0;
        }        
    for(i = 0; i < NumGenerators; i++) PG[i] = 0;                    
    for(k = 1; k <= NumRelators; k++) if(LR[k])
        {
        e = 0;
        Rk = p = *Relators[k];
        x = *p++;
        while(x == *p) p++;
        
        if(*p == EOS)
            {
            /**********************************************************************************
                In this case, Relators[k] is a power. Increment PG[i] if generator i is the
                generator that appears in Relators[k].    
            **********************************************************************************/    
            
            if(x < 95) i = x - 65;
            else  i = x - 97;
            switch(DF[i])
                {
                case 0:
                case 1:
                    PG[i]++;
                    q = T[i];
                    if(x == q[0] &&  x == q[1])
                        {
                        if(x < 95) ED[i][(x << 1) -130] += LR[k];
                        else ED[i][(x << 1) -193] += LR[k];
                        }
                    else        
                    if(x == q[2] &&  x == q[3])
                        {
                        if(x < 95) ED[i][(x << 1) -129] += LR[k];
                        else ED[i][(x << 1) -194] += LR[k];
                        }
                    break;    
                case 2:
                    PG[i]++;
                    break;    
                case 3:
                    PG[i]++;
                    break;        
                }    
            continue;
            }
                
        q = Rk + LR[k] - 1; 
        if(x != *q) 
            {
            x = *q;
            p = Rk;
            }    
        y = *p;
        
        while(1)
            {
            while(y == *p) 
                {
                p++;
                e++;
                }
                                
            if(*p == EOS)
                {
                p = Rk;
                while(y == *p)    
                    {
                    p++; 
                    e++;
                    }
                z = *p;                    
                if(y < 95) i = y - 65;
                else  i = y - 97;
                q = T[i];
                switch(DF[i])
                    {
                    case 0:
                    case 1:
                        if(x == q[0] && y == q[1])
                            {
                            if( e > 1)
                                {
                                if(y < 95) ED[i][(y << 1) -130]++;
                                else ED[i][(y << 1) -193]++;
                                }
                            else
                                {
                                if(z < 95) ED[i][(z << 1) -130]++;
                                else ED[i][(z << 1) -193]++;
                                }
                            }
                        else        
                        if(z == q[2] && y == q[3])
                            {
                            if( e > 1)
                                {
                                if(y < 95) ED[i][(y << 1) -129]++;
                                else ED[i][(y << 1) -194]++;
                                }
                            else
                                {
                                if(x < 95) ED[i][(x << 1) -129]++;
                                else ED[i][(x << 1) -194]++;
                                }
                            }
                        else    
                        if(e > 1)
                            {    
                            if(y == q[0] && y == q[1])
                                {
                                if(z < 95) ED[i][(z << 1) -130]++;
                                else ED[i][(z << 1) -193]++;
                                }
                            else    
                            if(y == q[2] && y == q[3])
                                {
                                if(x < 95) ED[i][(x << 1) -129]++;
                                else ED[i][(x << 1) -194]++;
                                }
                            }                    
                        break;
                            
                    case 2:
                        if((x == q[0] && y == q[1]) || (z == q[3] && y == q[2]))
                            {
                            for(z = 0; z < 3 && (EXP[i][z] != e); z++) ;
                            EXR[i][z] = TRUE;
                            }
                        else        
                        if((x == q[4] && y == q[5]) || (z == q[7] && y == q[6]))
                            {
                            for(z = 0; z < 3 && (EXP[i][z] != e); z++) ;
                            EXL[i][z] = TRUE;
                            }
                        break;    
                    
                    case 3:
                        break;            
                    }                
                break;
                }
            else
                {
                z = *p;
                if(y < 95) i = y - 65;
                else  i = y - 97;
                q = T[i];
                switch(DF[i])
                    {
                    case 0:
                    case 1:
                        if(x == q[0] && y == q[1])
                            {
                            if( e > 1)
                                {
                                if(y < 95) ED[i][(y << 1) -130]++;
                                else ED[i][(y << 1) -193]++;
                                }
                            else
                                {
                                if(z < 95) ED[i][(z << 1) -130]++;
                                else ED[i][(z << 1) -193]++;
                                }
                            }
                        else        
                        if(z == q[2] && y == q[3])
                            {
                            if( e > 1)
                                {
                                if(y < 95) ED[i][(y << 1) -129]++;
                                else ED[i][(y << 1) -194]++;
                                }
                            else
                                {
                                if(x < 95) ED[i][(x << 1) -129]++;
                                else ED[i][(x << 1) -194]++;
                                }
                            }
                        else
                        if(e > 1)
                            {    
                            if(y == q[0] && y == q[1])
                                {
                                if(z < 95) ED[i][(z << 1) -130]++;
                                else ED[i][(z << 1) -193]++;
                                }
                            else    
                            if(y == q[2] && y == q[3])
                                {
                                if(x < 95) ED[i][(x << 1) -129]++;
                                else ED[i][(x << 1) -194]++;
                                }
                            }                    
                        break;
                    
                    case 2:
                        if((x == q[0] && y == q[1]) || (z == q[3] && y == q[2]))
                            {
                            for(z = 0; z < 3 && (EXP[i][z] != e); z++) ;
                            EXR[i][z] = TRUE;
                            }
                        else        
                        if((x == q[4] && y == q[5]) || (z == q[7] && y == q[6]))
                            {
                            for(z = 0; z < 3 && (EXP[i][z] != e); z++) ;
                            EXL[i][z] = TRUE;
                            }
                        break;
                    
                    case 3:
                        break;    
                    }
                x = y;
                y = *p;
                e = 0;
                }
            }    
        }    
}            

unsigned int Diagram_5()
{
    /******************************************************************************************
        This routine computes the valences of the vertices of the Heegaard diagram and enters
    	the valence of vertex i in the array V at V[i]. It also checks for a possible overflow
    	condition. An entry in the array B[] can be as large as the sum of the lengths of two
    	dual relators. Since B[] holds unsigned ints which must be no larger than 2^16, an
    	entry in B[] could overflow if there are dual relators longer than 2^15 characters.
    ******************************************************************************************/
    
    register unsigned char	h,
                            i,
                            j;
                            
    register unsigned long  t;

    for(i = Vertices; i > 0; )
        {                                    
        for(h = 0,t = 0L,i -= 2; (j = AJ1[i][h]) < VERTICES; h++) t += A[i][j];
        if(t > 32765) return(TOO_LONG);
        V[i] = V[i+1] = t;
        }                    
    return(NO_ERROR);        
}

unsigned int Offset()
{
    /******************************************************************************************
        The routine Offset() determines how the edges at vertex 2i and vertex 2i+1 are
        identified in the Heegaard diagram, for i in the range: 0 <= i < Numgenerators. 
            If v is a vertex of the diagram, the first vertex in the link of vertex v -- found
        in FV[v] -- is that vertex in the link of vertex v which appears first in the ordering
        A,a,B,b,C,c etc. Then the edges that meet vertex v in the diagram are numbered 0,1,2 ...
        in counter-clockwise order about v, with edge number 0 the first edge in the band of
        edges joining v and FV[v].
            Offset() computes an integer, which we call the offset, with the property that if
        we enter vertex i at edge e, then we leave vertex 2i+1 at edge (offset - e). It is a
        pleasant fact that the offset is a symmetric function of the two vertices 2i and 2i+1.
            If both vertex 2i and vertex 2i+1 have valence two in the reduced Whitehead graph,
        Offset() may not be able to completely determine the offset. In this case however, if
        the presentation is realizable, then there are only two possible values that offset
        could assume. Offset computes both of these possibilities, and leaves one possible
        value of offset in the array OSA[] at OSA[i] and the other in the array OSB[] at
        OSB[i]. NGV2 is then incremented to indicate the existence of this ambiguity and GV2[i]
        set to TRUE. The routine Valence_Two() will be called later to attempt to resolve the
        ambiguous offset.     
    ******************************************************************************************/
        
    register unsigned int   a,
                            b,
                            c,
                            p,
                            q;
                            
    unsigned int            d,
                            i,
                            Nexp346;
    
    unsigned int            OffsetSub();
    
    NGV2 = 0;
    Nexp346 = 0;            
    for(i = 0; i < NumGenerators; i++) switch(DF[i])
        {
        case 0:
            /**********************************************************************************
                This is the case if generator i appears with maximal exponent 2 or less, and
                either vertex 2i or vertex 2i+1 has valence greater than two in the reduced
                Whitehead graph.
            **********************************************************************************/
            
            /**********************************************************************************
                The following check determines whether a special case holds in genus two which
                allows Heegaard to find a Heegaard diagram even though a maximal band
                contains more than half of the edges meeting a vertex.
            **********************************************************************************/
            
            if(NumGenerators == 2 && NumRelators == 2)
                {
                a = SV[i];
                b = TV[i];
                if(a & 1)
                    c = a - 1;
                else
                    c = a + 1;
                d = b + 1;    
                if((A[a][b] > A[b][d] + A[c][b]) &&
                    (ED[i][a] == A[d][a]) &&
                    (ED[i][b] == A[d][b]))
                return(TOO_LONG);
                }
                    
            /**********************************************************************************
                Next take care of the special case where generator i appears in a relator
                which is a power of generator i.
            **********************************************************************************/    
            
            if(PG[i])
                {
                c = i << 1;
                a = c + 1;
                b = FV[c];
                d = 0;
                while(b != a)
                    {
                    d += A[c][b];
                    b = CO[c][b];                    
                    }
                a = c;
                c ++;
                b = FV[c];
                while(b != a)
                    {
                    d += A[c][b];
                    b = CO[c][b];
                    }
                d += PG[i] - 1;
                d = d % V[c] + V[c];    
                OSA[i << 1] = OSA[(i << 1) + 1] = d;
                GV2[i] = FALSE;        
                break;
                }        
            
            a = SV[i];
            c = TV[i];
            b = FV[c];
            d = 0;
            while(b != a)
                {
                d += A[c][b];
                b = CO[c][b];
                }        
            if(c & 1) c--;
            else c++;
            a = b = FV[c];
                            
            do
                {
                d += A[c][b];
                if(ED[i][b] == A[c][b])
                    {
                    while(ED[i][b] == A[c][b])
                        {
                        b = CO[c][b];
                        d += ED[i][b];
                        }
                    break;
                    }
                if(ED[i][b] == 0)
                    {
                    while(ED[i][b] == 0)
                        {
                        b = CO[c][b];        
                        d += A[c][b];
                        }            
                    d += A[SV[i]][TV[i]];
                    d -= ED[i][b];
                    break;
                    }    
                b = CO[c][b];
                }
            while(a != b); 
            
            d = d % V[c] + V[c] - 1;
            OSA[i << 1] = OSA[(i << 1) + 1] = d;
            GV2[i] = FALSE;        
            break;        
        
        case 1:
            /**********************************************************************************
                This is the case where generator i appears with maximal exponent 2 or less,
                and both vertex 2i and vertex 2i+1 have valence two in the reduced Whitehead
                graph.
            **********************************************************************************/    
            
            /**********************************************************************************
                First take care of the special case where generator i appears in a relator
                which is a power of generator i.
            **********************************************************************************/    
            
            if(PG[i])
                {
                c = i << 1;
                a = c + 1;
                b = FV[c];
                d = 0;
                while(b != a)
                    {
                    d += A[c][b];
                    b = CO[c][b];                    
                    }
                a = c;
                c ++;
                b = FV[c];
                while(b != a)
                    {
                    d += A[c][b];
                    b = CO[c][b];
                    }
                d += PG[i] - 1;
                d = d % V[c] + V[c];    
                OSA[i << 1] = OSA[(i << 1) + 1] = d;
                GV2[i] = FALSE;        
                break;
                }            
            
            a = SV[i];
            c = TV[i];
            b = FV[c];
            d = 0;
            if(b != a) d += A[c][b];
            if(c & 1) c --;
            else c ++;
            b = FV[c];
                            
            if(ED[i][b] == A[c][b])
                {
                d += A[c][b];
                d = d % V[c] + V[c] - 1;
                OSA[i << 1] = OSA[(i << 1) + 1] = d;
                GV2[i] = FALSE;
                }
            else
                {
                if(ED[i][b] == 0)
                    {
                    d += A[c][b];
                    b = CO[c][b];
                    d += A[c][b];
                    d = d % V[c] + V[c] - 1;
                    OSA[i << 1] = OSA[(i << 1) + 1] = d;
                    GV2[i] = FALSE;
                    }
                else
                    {
                    d += ED[i][b];
                    OSA[i << 1] = OSA[(i << 1) + 1] = d % V[c] + V[c] - 1;
                    d -= ED[i][b];
                    d += A[c][b];
                    b = CO[c][b];
                    d += ED[i][b];
                    OSB[i << 1] = OSB[(i << 1) + 1] = d % V[c] + V[c] - 1;
                    NGV2 ++;
                    GV2[i] = TRUE;
                    }
                }
            break;                                            
        
        case 2:
            GV2[i] = FALSE;
            if(!EXP[i][1])
                {
                /******************************************************************************
                    This is the case where generator i appears only with one exponent and 
                    that exponent is greater than two.
                ******************************************************************************/
                
                p = EXP[i][2];
                q = 1;
                a = NEX[i][2];
                b = c = 0;
                OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                if(!DrawingDiagrams && !TestRealizability2) switch(p)
                    {
                    case 3:
                    case 4:
                    case 6:
                        if(++Nexp346 > 1) return(NON_UNIQUE_2);
                        if(Non_Unique()) return(NON_UNIQUE_1);
                        break;
                    case 5:
                        return(NON_UNIQUE_3);
                    default:
                        return(NON_UNIQUE_4);
                    }
                if(DrawingDiagrams && UDV[WhichInput] == 0) switch(p)
                    {
                    case 3:
                    case 4:
                    case 6:
                        if(++Nexp346 > 1)
                            {
                            UDV[WhichInput] = NON_UNIQUE_2;
                            break;
                            }
                        if(Non_Unique())
                            {
                            UDV[WhichInput] = NON_UNIQUE_1;
                            break;
                            }
                        break;
                    case 5:
                        UDV[WhichInput] = NON_UNIQUE_3;
                        break;
                    default:
                        UDV[WhichInput] = NON_UNIQUE_4;
                        break;
                    }                    
                }
                
            else
            if(!EXP[i][0])
                {
                /******************************************************************************
                    This is the case where generator i appears with only two distinct 
                        exponents and the largest exponent is greater then two.
                ******************************************************************************/
                
                a = b = c = d = 0;
                if(EXR[i][1]) a = 1;
                if(EXR[i][2]) b = 1;
                if(EXL[i][1]) c = 1;
                if(EXL[i][2]) d = 1;
                switch(a*d - b*c)
                    {
                    case 1:
                        p = EXP[i][1];
                        q = EXP[i][2];
                        a = NEX[i][1];
                        b = NEX[i][2];
                        c = 0;
                        OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                        break;
                    case -1:
                        p = EXP[i][2];
                        q = EXP[i][1];
                        a = NEX[i][2];
                        b = NEX[i][1];
                        c = 0;
                        OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                        break;
                    default:
                        return(FATAL_ERROR);
                    }
                }
            else
                {
                /******************************************************************************
                    This is the remaining case. Generator i appears with three distinct 
                            exponents and the largest of these is greater than two.
                ******************************************************************************/
                
                a = b = c = d = 0;
                if(EXR[i][0]) a = 1;
                if(EXR[i][1]) b = 1;
                if(EXL[i][0]) c = 1;
                if(EXL[i][1]) d = 1;
                switch(a*d - b*c)
                    {
                    case 1:
                        p = EXP[i][0];
                        q = EXP[i][1];
                        a = NEX[i][0];
                        b = NEX[i][1];
                        c = NEX[i][2];
                        OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                        break;
                    case -1:
                        p = EXP[i][1];
                        q = EXP[i][0];
                        a = NEX[i][1];
                        b = NEX[i][0];
                        c = NEX[i][2];
                        OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                        break;
                    default:
                        return(FATAL_ERROR);
                    }
                }
            break;        
        
        case 3:
            if(!EXP[i][1])
                {
                /******************************************************************************
                    This is the case where generator i appears only with one exponent, that
                    exponent is greater than two, and both vertex 2i and vertex 2i+1 have
                    valence two.
                ******************************************************************************/
                
                p = EXP[i][2];
                q = 1;
                a = NEX[i][2];
                b = c = 0;
                OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                GV2[i] = FALSE;
                if(!DrawingDiagrams && !TestRealizability2) switch(p)
                    {
                    case 3:
                    case 4:
                    case 6:
                        if(++Nexp346 > 1) return(NON_UNIQUE_2);
                        if(Non_Unique()) return(NON_UNIQUE_1);
                        break;
                    case 5:
                        return(NON_UNIQUE_3);
                    default:
                        return(NON_UNIQUE_4);
                    }
                if(DrawingDiagrams && UDV[WhichInput] == 0) switch(p)
                    {
                    case 3:
                    case 4:
                    case 6:
                        if(++Nexp346 > 1)
                            {
                            UDV[WhichInput] = NON_UNIQUE_2;
                            break;
                            }
                        if(Non_Unique())
                            {
                            UDV[WhichInput] = NON_UNIQUE_1;
                            break;
                            }
                        break;
                    case 5:
                        UDV[WhichInput] = NON_UNIQUE_3;
                        break;
                    default:
                        UDV[WhichInput] = NON_UNIQUE_4;
                        break;
                    }                
                }
            else
            if(!EXP[i][0])
                {
                /******************************************************************************
                    This is the case where generator i appears with only two distinct 
                    exponents, the largest exponent is greater then two, and both vertex 2i
                    and vertex 2i+1 have valence two.
                ******************************************************************************/
                
                p = EXP[i][1];
                q = EXP[i][2];
                a = NEX[i][1];
                b = NEX[i][2];
                c = 0;
                OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c); 
                p = EXP[i][2];
                q = EXP[i][1];
                a = NEX[i][2];
                b = NEX[i][1];
                c = 0;
                OSB[i << 1] = OSB[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c); 
                NGV2 ++;
                GV2[i] = TRUE;
                }
            else
                {
                /******************************************************************************
                    This is the remaining case. Generator i appears with three distinct 
                    exponents, the largest of these is greater than two, and both vertex 2i
                    and vertex 2i+1 have valence two.
                ******************************************************************************/
                
                p = EXP[i][0];
                q = EXP[i][1];
                a = NEX[i][0];
                b = NEX[i][1];
                c = NEX[i][2];
                OSA[i << 1] = OSA[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                p = EXP[i][1];
                q = EXP[i][0];
                a = NEX[i][1];
                b = NEX[i][0];
                c = NEX[i][2];
                OSB[i << 1] = OSB[(i << 1) + 1] = OffsetSub(i,p,q,a,b,c);
                NGV2 ++;
                GV2[i] = TRUE;        
                }    
            break;    
        }            
    return(NO_ERROR);                
}

unsigned int OffsetSub(unsigned int i,unsigned int p,unsigned int q,unsigned int a,unsigned int b,
	unsigned int c)
{
    /******************************************************************************************
        OffsetSub() is a routine called by Offset() in those cases where a generator appears
        with maximal exponent greater than two. The generator appears with exponents p,q and
        p+q and these exponents appear respectively a,b and c times in the relators.
    ******************************************************************************************/
        
    register long   d;
    
    unsigned int 	GCD();
    
    GCD(p,q);
    d = (Recip_Q - 1)*(a + c) - (Recip_P + 1)*(b + c) + c;        
    c = i << 1;
    if(d < 0L)
        {
        d = - d;
        d = d % V[c];
        d = V[c] - d;
        }    
    b = FV[c];
    d += PG[i];
    while(b != c + 1)
        {
        d += A[c][b];
        b = CO[c][b];
        }            
    c++;
    b = FV[c];
    while(b != c - 1)
        {
        d += A[c][b];
        b = CO[c][b];
        }
    d = d % V[c] + V[c] - 1;
    return(d);
}

void Diagram_6(void)
{    
    /******************************************************************************************
            This routine sets up the array B[][]. This information is used when the curves of
        the Heegaard diagram are traced. More specifically, suppose we are tracing curves in
        the Heegaard diagram. And suppose we leave vertex i at edge e, with edge e belonging to
        a band of edges joining vertex i and vertex j. Then we will enter vertex j at the edge
        B[i][j] - e. It happens that B[][] is symmetric, so B[i][j] = B[j][i].    
    ******************************************************************************************/            
    
    register int            h;
    
    register unsigned int   a,
                            b,
                            i,
                            j;
                            
    for(i = 0; i < Vertices - 1; i++)
    for(h = VWG[i] - 1; h >= 0; h --)
        {
        if((j = AJ1[i][h]) <= i) break;
        b = 0;
        a = FV[i];
        while(a != j)
            {
            b += A[i][a];
            a = CO[i][a];
            }
        a = FV[j];
        while(a != i)
            {
            b += A[j][a];
            a = CO[j][a];
            }
        b += A[j][a] -1;
        B[i][j] = B[j][i] = b;            
        }
}    

unsigned int Diagram_7()
{
    /******************************************************************************************
        Diagram_7() uses the information developed, by the other routines in this file, to 
        determine whether the presentation in Relators[] is realizable by a Heegaard diagram.
        The routine "runs" around the diagram and reads off what should be the realizations of
        the relators. It leaves these in the array OutRelators[]. While it is reading off the
        relators, the routine also determines what the dual relators are. It leaves the dual
        relators in the array DualRelators[].
            Finally, when TestRealizability1 is TRUE, Diagram_7() compares the relators in
        OutRelators[] with the relators in Relators[]. If these lists agree, up to permutation,
        taking cyclic conjugates, and forming inverses, then the presentation is realizable.
        Otherwise the presentation is not realizable.    
    ******************************************************************************************/
            
    register unsigned char  i,
                            *p,
                            *q,
                            *r,
                            v,
                            w;
                            
    register unsigned int   d,
                            e;
    
    unsigned char           *DR[VERTICES];
                            
    int                     j;
                            
    unsigned int            length,
                            max;
    
    unsigned long           TotLength1,
                            TotLength2;
    
    for(i = 1,TotLength1 = 0L; i <= NumGenerators; i++)
        {
        max = V[(i-1) << 1];
        TotLength1 += max;
        max ++;
        if(DualRelators[i] != NULL) DisposeHandle((char **) DualRelators[i]);
        DualRelators[i] = (unsigned char **) NewHandle(max);
        if(DualRelators[i] == NULL) Mem_Error();
        p = *DualRelators[i];
        for(e = V[(i-1) << 1]; e > 0; e--) *p++ = '@';
        *p = EOS;
        }    
    for(i = 1,max = 0; i <= NumRelators; i++) if(LR[i] > max) max = LR[i];
    if(Temp16 != NULL) DisposeHandle((char **) Temp16);
    Temp16 = (unsigned char **) NewHandle(max + 2);
    if(Temp16 == NULL) Mem_Error();    
        
    for(i = 0; i < Vertices; i++) DR[i] = *DualRelators[(i >> 1) + 1];
    
    SRError = 3;        /* This tells other routines that OutRelators[] is dirty. */
        
    for(j = 1,i = 0,TotLength2 = 0L; j <= NumGenerators; j++) 
        {
        p = *DualRelators[j];
        while(*p)
            {
            if(*p == '@')
                {
                if(i >= NumRelators) return(FATAL_ERROR);
                r = *Temp16;
                length = 0;
                v = (j - 1) << 1;
                e = p - *DualRelators[j];
                while(1)
                    {
                    if(v & 1)
                        {
                        *r = (v >> 1) + 97;
                        e = OSA[v] - e;
                        if(e >= V[v]) e -= V[v];
                        q = DR[v] + e;
                        *q = i + 97;
                        w = v - 1;
                        }
                    else 
                        {
                        *r = (v >> 1) + 65;
                        q = DR[v] + e;
                        if(*q != '@') break;
                        *q = i + 65;
                        e = OSA[v] - e;
                        if(e >= V[v]) e -= V[v];
                        w = v + 1;
                        }                        
                    r++;    
                    if(++length > max)
                        {
                        *r = EOS;
                        return(FATAL_ERROR);
                        }                    
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
                
                if(OutRelators[i+1] != NULL) DisposeHandle((char **) OutRelators[i+1]);
        		OutRelators[i+1] = (unsigned char **) NewHandle(length + 1);
        		if(OutRelators[i+1] == NULL) Mem_Error();
        		q = *OutRelators[i+1];
        		r = *Temp16;
        		while( (*q++ = *r++) ) {}
                TotLength2 += length;
                i++;
                if(TotLength2 >= TotLength1) goto END;
                }
            p++;
            }    
        }

    END:
    
    if(Batch == 51)
    	{
    	printf("\n The Dual Relators are:");
    	Print_Relators(DualRelators,NumGenerators);
    	if(H_Results != NULL)
    		{
    		fprintf(H_Results,"\n\n%s",PresName);
    		Print_Relators2(DualRelators,NumGenerators);
    		}
    	return(NO_ERROR);
    	}
    
    /******************************************************************************************
        Make sure that, in all cases, the number of OutRelators allocated equals NumRelators.
    ******************************************************************************************/
            
    for(j = i + 1; j <= NumRelators; j++)
        {
        if(OutRelators[j] != NULL) DisposeHandle((char **) OutRelators[j]);
        OutRelators[j] = (unsigned char **) NewHandle(sizeof(char));
        if(OutRelators[j] == NULL) Mem_Error();
        r = *OutRelators[j];
        *r = EOS;        
        }
            
    if(TestRealizability1 || TestRealizability4)
        {
        for(i = 1; i <= NumRelators; i++)
            {
            LR[0] = GetHandleSize((char **) OutRelators[i]) - 1;
            j = Compare(*OutRelators[i]);
            if(j == 0)
                {
                for(i = 1; i <= NumRelators; i++)
                    {
                    p = *Relators[i];
                    if(*p > 124) *p -= 125;
                    }
                return(FATAL_ERROR);
                }
            if(j < 0) j = -j;
            p = *Relators[j];
            *p += 125;
            }
        for(i = 1; i <= NumRelators; i++)
            {
            p = *Relators[i];
            if(*p > 124) *p -= 125;
            }            
        }        
    return(NO_ERROR);    
}

int Compare(unsigned char *p)
{    
    /******************************************************************************************
            This routine returns i if the string p is a cyclic conjugate of the 
            string *Relators[i]. It returns -i if the string p is a cyclic conjugate
            of the inverse of *Relators[i] and otherwise returns 0.
    ******************************************************************************************/
    
    register int i;
    
    for(i = 1; i <= NumRelators; i++) if(LR[i] == LR[0] && Compare_Str(*Relators[i],p,LR[0]))
        return(i);
    Inverse(p);    
    for(i = 1; i <= NumRelators; i++) if(LR[i] == LR[0] && Compare_Str(*Relators[i],p,LR[0]))
        {
        Inverse(p);
        return(-i);
        }
    Inverse(p);                    
    return(0);
}

#define PRIME_P         65521
#define PRIME_P_X_128     8386688

int Compare_Str(unsigned char *S1,unsigned char *S2,unsigned long length)      
{    
    /******************************************************************************************
        Given strings S1 and S2 of equal length, this routine returns TRUE if string S2 is a
        cyclic conjugate of string S1 and otherwise returns FALSE. The routine uses a
        modified version of a method developed by Rabin and Karp to efficiently compare S1
        and S2. This routine determines a hash value for a string by forming a polynomial in x,
        whose ith coefficient is the ASCII decimal value of the ith character in the string,
        and then the routine evaluates the polynomial at x = 128. The hash value of the
        string is the value of this polynomial taken mod 65521. Note that this method lends
        itself to testing cyclic conjugates because it is easy to compute the new hash value
        when a character is moved from the beginning to the end of a string. The running time
        of the routine should be on the order of 3 times the length of the strings.
    ******************************************************************************************/    
    
    register unsigned char  C,
                            *p,
                            *q,
                            *r;
    
    register unsigned long  D,
                            HP,
                            HT;                        

    if(length == 0) return(TRUE);
    length --;
    for(HT = D = 1L; D < length; HT++) D = (D << 1);
    for(D = 1L; HT > 0; )
        {
        D = (D*D) % PRIME_P;
        HT--;
        if(length & (1L << HT)) D = (D << 7) % PRIME_P;
        }
    
    for(p = S1,HT = 0L; (C = *p); p++) HT = ((HT << 7) + C) % PRIME_P;
    for(q = S2,HP = 0L; (C = *q); q++) HP = ((HP << 7) + C) % PRIME_P;

    for(p = S1; (C = *p); p++)
        {
        if(HT == HP) for(q = p,r = S2;  ;q++,r++)
            {
            if(*q == EOS) q = S1;
            if(*r == EOS) return(TRUE);
            if(*q != *r) break;
            }
        HT = ((((HT + PRIME_P_X_128) - C*D) << 7) + C) % PRIME_P;
        }
    return(FALSE);    
}

void Get_Bdry_Comps(int Print,int Where,unsigned int WhichPres)
{
    unsigned char   *p,
                    x;
    
    int             CompNum,
                    Genus;
    
    unsigned int    Edge,
                    i,
                    j,
                    TheBdryComp,
                    v,
                    V2,
                    V3;

    Planar(FALSE,TRUE);
    Sep_Surface();
    
    for(i = 1; i <= NumBdryComps; i++)
        {
        NEBC[i] = EOS;
        NFBC[i] = EOS;
        }
    
    for(i = 1; i <= NumFaces; i++)
        {
        V2 = Face[i][1];
        V3 = Face[i][2];
        for(Edge = 0, v = FV[V2]; v != V3; v = CO[V2][v]) Edge += A[V2][v];
        if(V2 & 1)
            {
            Edge = OSA[V2] - Edge;
            Edge ++;
            while(Edge >= V[V2]) Edge -= V[V2];    
            }
        p = *DualRelators[(V2 >> 1) + 1] + Edge;
        x = *p;
        x = x << 1;
        if(x < 194) x -= 130;
        else x -= 193;
        TheBdryComp = zz[x];
        BCF[i] = TheBdryComp;
        NFBC[TheBdryComp] ++;
        for(p = Face[i],j = 0; *p < VERTICES; p++,j++) ;
        NEBC[TheBdryComp] += j;
        }
    
    for(i = 0; i <= NumGenerators; i++) BCWG[i] = EOS;

    for(i = 1; i <= NumBdryComps; i++)
        {
        Genus = NFBC[i] - NEBC[i]/2 + NRBC[i];
        if(Genus < 0)
            {
            Genus = - Genus;
            Genus /= 2;
            }
        else
            {
            Genus /= 2;
            Genus = - Genus;
            }
        Genus ++;
        BCWG[Genus] ++;
        GBC[i] = Genus;
        }
    
    for(i = NumGenerators + 1; i > 0; i--) if(BCWG[i-1])
        {
        BCWG[i] = BDRY_UNKNOWN;
        break;
        }
        
    if(Print && Batch == FALSE)
        {
        CompNum = ComponentNum[WhichPres];        
        j = NumBdryComps - BCWG[0];
        switch(j)
            {
            case 0:
                if(TotalComp == 1)
                    printf("\n\n                    This manifold is closed.\n");                
                else
                    printf("\n\n                    'Summand' %d of M is closed.\n",CompNum);
                break;
            case 1:
                for(i = 1; i <= NumGenerators; i++) if(BCWG[i])
                    {
                    if(TotalComp == 1)
                        printf("\n\n                    This manifold has one boundary component of genus %u.\n",
                            i);                    
                    else
                        printf("\n\n                    'Summand' %d of M has one boundary component of genus %u.\n",
                            CompNum,i);
                    break;
                    }
                break;        
            default:
                for(i = 1; BCWG[i] < BDRY_UNKNOWN; i++) if(BCWG[i] == j)
                    {
                    if(TotalComp == 1)
                        printf("\n\n                    This manifold has %u boundary components of genus %u.\n",
                            j,i);                    
                    else
                        printf("\n\n                    'Summand' %d of M has %u boundary components of genus %u.\n",
                            CompNum,j,i);
                    return;
                    }
                if(TotalComp == 1)
                    printf("\n\n                    This manifold has the following boundary components:\n");                
                else    
                    printf("\n\n                    'Summand' %d of M has the following boundary components:\n",
                        CompNum);
                for(i = 1; (j = BCWG[i]) < BDRY_UNKNOWN; i++) if(j)
                    {
                    if(j == 1)
                        printf("\n                         %2u component  of genus %2u.",j,i);
                    else
                        printf("\n                         %2u components of genus %2u.",j,i);
                    }
                printf("\n");    
                break;    
            }
        }        
}

int Delete_Redundant_Relators()
{
    /******************************************************************************************
        This routine is called by Heegaard to delete redundant relators from the Heegaard
        diagram. It looks for relators which separate a 2-sphere boundary component of the
        manifold from another boundary component of the manifold. Such relators can be
        deleted. There may be several possible ways in which this can be done. Topologically,
        this is irrelevant, but which relator is deleted can have an influence on which
        presentations Heegaard finds.
    ******************************************************************************************/    
    
    unsigned char           DRL[MAXNUMRELATORS + 1],
                            **Temp;
    
    int                     NumS2BC,
                            OrigNumRelators;
    
    unsigned int            i,
                            j,
                            k,
                            N1,
                            N2,
                            RL[MAXNUMRELATORS];
    
    unsigned long           TempLR;
                                    
    for(i = 1; i <= NumRelators; i++) DRL[i] = EOS;
    for(i = 0; i < NumRelators; i++) RL[i] = i;
    
    NumS2BC = BCWG[0];
        
    TRY_AGAIN:
    
    /******************************************************************************************
                        Randomize the order in which relators will be tested.
    ******************************************************************************************/
    
    for(i = NumRelators; i > 1; )
        {
        j = abs(rand()) % i;
        i--;
        k = RL[i];
        RL[i] = RL[j];
        RL[j] = k;
        }    
    
    for(k = 0; k < NumRelators; k++)
        {
        i = RL[k];
        i += i;
        N1 = zz[i];
        N2 = zz[i+1];
        if(N1 == N2) continue;
        if(GBC[N1] && GBC[N2]) continue;
        if(GBC[N1])
            {
            N1 = N2;
            N2 = zz[i];
            }
        
        /**************************************************************************************
            Amalgamate boundary component N1, with boundary component N2, decrement the number
            of boundary components which are 2-spheres, decrement the total number of boundary
            components, and flag Relator i/2 + 1 so we can delete it later.
        **************************************************************************************/
                 
        for(j = 0; j < 2*NumRelators; j++) if(zz[j] == N1) zz[j] = N2;
        NumBdryComps --;
        DRL[i/2 + 1] = TRUE;
        NumS2BC --;
        if(NumS2BC) goto TRY_AGAIN;
        break;
        }
    
    OrigNumRelators = NumRelators;
    
    for(i = 1; i <= OrigNumRelators; i++) if(DRL[i])
        {    
        LR[0] = GetHandleSize((char **) OutRelators[i]) - 1;
        j = abs(Compare(*OutRelators[i]));
        if(j == 0) return(TOO_LONG);
        Temp = Relators[NumRelators];
        Relators[NumRelators] = Relators[j];
        Relators[j] = Temp;
        TempLR = LR[NumRelators];
        LR[NumRelators] = LR[j];
        LR[j] = TempLR;
        NumRelators --;
        Length -= LR[0];
        }

    if(NumRelators == OrigNumRelators) return(FALSE);

    j = OrigNumRelators - NumRelators;
        
    if(Micro_Print)
        {
        printf("\n\nDeleted %d redundant relator(s) from presentation %d to get:\n",
            j,ReadPres + 1);
        Print_Relators(Relators,NumRelators);
        }
                    
    if(Find_Flow_A(NORMAL,FALSE) == TOO_LONG) return(TOO_LONG);
                    
    if(On_File() == NumFilled)
        {
        if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);    
        if(Dup_On_File < INFINITE)
        	{
            if(Save_Pres(ReadPres,Dup_On_File,Length,1,1,1,0,0)) return(TOO_LONG);
            ER[NumFilled - 1] = 0;                                
            Mark_As_Duplicate(Dup_On_File);             
            }
         else
            {
            if(Save_Pres(ReadPres,0,Length,1,1,1,0,0)) return(TOO_LONG);             
            if(BCWG[1] == BDRY_UNKNOWN)    
                BDY[NumFilled - 1] = FALSE;
            else
                BDY[NumFilled - 1] = TRUE;    
            UDV[NumFilled - 1] = 0;
            ER[NumFilled - 1] = 0;        
            }    
        }
    
    if(ER[ReadPres] == 0) ER[ReadPres] = -1;
    
    /* Try setting UDV[ReadPres] = DONE so we won't process it more than once. */
    
    UDV[ReadPres] = DONE;
        
    return(TRUE);            
}

int MG_Bdry_Comp_Data(unsigned int WhichPres)
{
    unsigned char   PS1[MAXNUMGENERATORS + 1],
                    PS2[MAXNUMGENERATORS + 2];
                    
    int             Ambiguity = FALSE,
                    BdryComps1,
                    BdryComps2,
                    CNumPres1,
                    CNumPres2,
                    CNumPres3,
                    i,
                    j,
                    k,
                    NumEmtyHand,
                    Num1H,
                    NumS1XD2,
                    Pres1,
                    Pres2,
                    Pres3,
                    TGenus1,
                    TGenus2;
    
    if(ComponentNum[WhichPres] == 1) return(FALSE);
    
    /******************************************************************************************
        Follow the entries of FR[] backwards until we find the presentation which "SPLIT" and
        created this summand. Then check whether the "SPLIT" was created by the routine
        Missing_Gen(). If the "SPLIT" was not created by Missing_Gen(), return.
    ******************************************************************************************/    
        
    for(Pres1 = FR[WhichPres]; UDV[Pres1] != SPLIT; Pres1 = FR[Pres1]) ;
    Pres2 = Daughters[Pres1];
    Pres3 = Pres2 + 1;
    if(PRIM[Pres3] != 30 && PRIM[Pres3] != 130) return(FALSE);
    CNumPres1 = ComponentNum[Pres1];
    CNumPres2 = ComponentNum[Pres2];
    CNumPres3 = ComponentNum[Pres3];
    if(CS[CNumPres3] != 3) return(FALSE);
    if(CBC[CNumPres1][0] < BDRY_UNKNOWN)
        {
        NumEmtyHand = NG[Pres3];
        if(CBC[CNumPres1][1] == BDRY_UNKNOWN)
            CBC[CNumPres1][0] = 1;
        else
            CBC[CNumPres1][0] = 0;    
        for(i = BdryComps1 = TGenus1 = 0; (j = CBC[CNumPres1][i]) < BDRY_UNKNOWN; i++)
            {
            BdryComps1 += j;
            TGenus1 += i*j;
            }
        for(i = BdryComps2 = TGenus2 = 0; (j = CBC[CNumPres2][i]) < BDRY_UNKNOWN; i++)
            {
            BdryComps2 += j;
            TGenus2 += i*j;
            }
        if(BdryComps2 >= BdryComps1 && BdryComps1 >= 1)
            Num1H = BdryComps2 - BdryComps1;
        else
            {
            CS[CNumPres3] = 4;
            return(FALSE);
            }
        if(Num1H > NumEmtyHand)
            {
            CS[CNumPres3] = 4;
            return(FALSE);            
            }
        NumS1XD2 = NumEmtyHand - Num1H;
        if(TGenus1 - TGenus2 != NumS1XD2)
            {
            CS[CNumPres3] = 4;
            return(FALSE);            
            }        
                    
        N1H[CNumPres3]         = Num1H;
        NS1XS2[CNumPres3]     = 0;
        NS1XD2[CNumPres3]     = NumS1XD2;
        
        if(Num1H == 0)
            {                
            for(i = j = 0; CBC[CNumPres1][i] < BDRY_UNKNOWN; i++)
                {
                j += CBC[CNumPres1][i];
                PS1[i] = j;
                }
            for(k = j = 0; CBC[CNumPres2][k] < BDRY_UNKNOWN; k++)
                {
                j += CBC[CNumPres2][k];
                PS2[k] = j;
                }
            for( ; k <= i; k++) PS2[k] = j;    
            for(i = 0, Ambiguity = FALSE; (j = CBC[CNumPres1][i]) < BDRY_UNKNOWN; i++)
                {
                if(PS2[i] < PS1[i])
                    {
                    CS[CNumPres3] = 4;
                    return(FALSE);    
                    }
                if(j && PS2[i] > PS1[i])
                    {
                    Ambiguity = TRUE;
                    break;
                    }
                }    
            }
        if(Ambiguity)
            {    
            CBC[CNumPres3][0]      = Num1H;
            CBC[CNumPres3][1]      = NumS1XD2;
            CBC[CNumPres3][2]      = BDRY_UNKNOWN;
            UDV[Pres3]             = MISSING_GEN_DONE2;
            CS[CNumPres3]          = 2;
            return(TRUE);            
            }
        
        for(i = 0; i <= MAXNUMGENERATORS; i++) CBC[CNumPres3][i] = EOS;
        for(i = k = 0; (j = CBC[CNumPres1][i]) < BDRY_UNKNOWN; i++) if(j)
        for( ; k <= i; k++) if(CBC[CNumPres2][k] < BDRY_UNKNOWN)
            CBC[CNumPres3][i-k] += CBC[CNumPres2][k] - CBC[CNumPres1][k];        
    
        for(i = MAXNUMGENERATORS + 1; i > 0; i--) if(CBC[CNumPres3][i-1])
            {
            CBC[CNumPres3][i] = BDRY_UNKNOWN;
            break;
            }
        if(i == 0)
            {
            CBC[CNumPres3][0] = 1;
            CBC[CNumPres3][1] = BDRY_UNKNOWN;
            }            
        
        CS[CNumPres3] = 2;        
        UDV[Pres3] = MISSING_GEN_DONE1;
        return(TRUE);
        }
        
    /******************************************************************************************
        Check whether Heegaard is done processing all of the presentations associated with
        the summand that "SPLIT". If it is, then we will make an attempt to determine the
        number and types of "handles" associated with this manifold.
    ******************************************************************************************/
    
    for(i = 0; i < NumFilled; i++) if(ComponentNum[i] == CNumPres1 && UDV[i] < DONE) break;
    if(i >= NumFilled)
        {
        NumEmtyHand = NG[Pres3];
        BdryComps1 = 1;
        for(i = BdryComps2 = 0; (j = CBC[CNumPres2][i]) < BDRY_UNKNOWN; i++)
            BdryComps2 += j;
        if(BdryComps2 >= BdryComps1)
            Num1H = BdryComps2 - BdryComps1;
        else
            {
            CS[CNumPres3] = 4;
            return(FALSE);            
            }    
        if(Num1H > NumEmtyHand)
            {
            CS[CNumPres3] = 4;
            return(FALSE);            
            }            
        NumS1XD2 = NumEmtyHand - Num1H;
        
        NS1XS2[CNumPres3] = 0;
        
        if(Num1H)
            {
            N1H[CNumPres3]         = 0;
            NS1XD2[CNumPres3]     = NumEmtyHand;
            CBC[CNumPres3][0]      = BDRY_UNKNOWN;        
            UDV[Pres3]            = MISSING_GEN_DONE2;
            }
        else
            {
            N1H[CNumPres3]         = 0;
            NS1XD2[CNumPres3]     = NumS1XD2;
            for(i = 0; i < NumS1XD2; i++) CBC[CNumPres3][i] = EOS;
            CBC[CNumPres3][i] = 1;
            CBC[CNumPres3][i+1] = BDRY_UNKNOWN;
            UDV[Pres3] = MISSING_GEN_DONE1;
            }
        CS[CNumPres3] = 2;                
        return(TRUE);
        }
    return(FALSE);            
}

int Sep_Surface()
{
    /******************************************************************************************
        This routine determines the components which arise when the curves representing the
        relators are deleted from the Heegaard surface. It calls Fill_DRA(), which gets the
        Whitehead graph of the dual presentation. The number of components of this graph gives
        the number of components which arise when the curves corresponding to the relators are
        deleted from the Heegaard surface, and it also allows us to determine which curves
        form the boundaries of each of these components. This information is returned in the
        global array zz[]. Relators which bound the same component of the Heegaard surface
        have equal values in zz[].
    ******************************************************************************************/    
    
    register unsigned int   h,
                            i,
                            j,
                            k,
                            *p,
                            *r;
                            
    int                     Vertices;                        
                                
    Fill_DRA();        
    Vertices = 2*NumRelators;    
    for(i = 0, p = zz; i < Vertices; i++,p++) *p = 0;
    for(k = 0, NumBdryComps = 0; k < Vertices; )
        {
        NumBdryComps ++;
        for(h = 0; h < Vertices && zz[h]; h++) ;
        if(h >= Vertices) goto FOUND_COMPS;
        zz[h] = NumBdryComps;
        k++;
        NRBC[NumBdryComps] = 1;
        for(r = UpDate,*r = h,p = r + 1; r < p; r++)
            {
            i = *r;
            for(j = 0; j < Vertices; j++) if(DRA[i][j])
                {
                if(i == j) continue;
                if(zz[j] == 0)
                    {
                    zz[j] = NumBdryComps;
                    *p++ = j;
                    NRBC[NumBdryComps] ++;
                    if(++k >= Vertices) goto FOUND_COMPS;
                    }
                }
            }
        }
        
    FOUND_COMPS:
    return(NumBdryComps);
}                                

void Fill_DRA(void)
{                        
    /******************************************************************************************
        This routine takes the cyclic relator pointed to by "p" and examines each pair of 
        consecutive letters that appear in that relator, including the pair comprised of the
        last letter and the first letter. For each such pair of letters, it adds the 
        appropriate edge to the array DRA[][]. Note that for this to work correctly in the
        case where NumRelators != NumGenerators, we must redefine Vertices to be 2*NumRelators.
    ******************************************************************************************/
    
    register unsigned char  i,
                            j,
                            *p,
                            Vertices,
                            x;
                            
    Vertices = 2*NumRelators;    
    for(i = 0; i < Vertices; i++)
    for(j = 0; j < Vertices; j++) DRA[i][j] = 0;
    for(j = 1; j <= NumGenerators; j++)
        {
        p = *DualRelators[j];
        if(*p == EOS) continue;
        x = *p << 1;
        if(x < 194) x -= 129;
        else x -= 194;
        i = x;
        p++;
        while( (x = *p) )
            {
            x = x << 1;
            if(x < 194) x -= 130;
            else x -= 193;
            DRA[i][x]++;
            if(x & 1) i = x - 1;
            else i = x + 1;
            p++;
            }
        x = **DualRelators[j];
        x = x << 1;            
        if(x < 194) x -= 130;
        else x -= 193;
        DRA[i][x]++;
        }
    for(i = 0; i < Vertices - 1; i++)
    for(j = i + 1; j < Vertices; j++)
        {
        DRA[i][j] += DRA[j][i];
        DRA[j][i] = DRA[i][j];
        }        
}

unsigned int GCD(unsigned int p,unsigned int q)
{
    /******************************************************************************************
        Given a pair of nonnegative integers p and q, this routine computes the GCD
        of p and q and returns GCD(p,q) as an unsigned integer. The routine simultaneously
        computes a pair of integers Recip_P and Recip_Q such that the relation
                    GCD(p,q) = p*Recip_P + q*Recip_Q  holds.
    ******************************************************************************************/    
    
    register long   d,
                    t,
                    u2,
                    v1,
                    v2;
                    
    long            u1;
        
    u1 = 0L;
    u2 = p;
    v1 = 1L;
    v2 = q;
    while(v2)
        {
        d = u2/v2;
        t = u1 - v1*d;
        u1 = v1;
        v1 = t;
        t = u2 - v2*d;
        u2 = v2;
        v2 = t;
        }
    if(p)
        {
        t = u2 - q*u1;
        if(t < 0)
            {
            t = - t;
            Recip_P = t/p;
            Recip_P = -Recip_P;
            }
        else    
            Recip_P = t/p;
        }
    else
        Recip_P = 0L;        
    Recip_Q = u1;
    return(u2);
}

void Inverse(register unsigned char *p)
{
    /******************************************************************************************
        This routine transforms the string pointed to by p into its inverse by first replacing
        each char in the string by its inverse and then writing the string backwards.
    ******************************************************************************************/
    
    register unsigned char  *q,
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

int Non_Unique()
{
    /******************************************************************************************
        Generally, if a generator appears in the relators with only one exponent, the Heegaard
        diagram is not unique. There are exceptions to this however. These exceptions occur
        when the "reduced" Whitehead graph is homeomorphic to a circle and there are not too
        many edges in the Heegaard diagram. Non_Unique() is called to check whether we have
        uniqueness in this special set of circumstances. For example, without this special
        check, Heegaard would declare that the diagram corresponding to the relator AAABB
        was not unique because generator A appears only with exponent 3. However this diagram
        is unique because there is an involution of the diagram which exchanges the two major
        faces of the diagram.
    ******************************************************************************************/
        
    unsigned int    i,
                    j,
                    k;
    
    /******************************************************************************************
        If more than one generator appears with maximal exponent greater than two, then the
        diagram is not unique.
    ******************************************************************************************/
    
    for(i = j = 0; i < NumGenerators; i++) if(EXP[i][2] > 2 && (++j > 1)) return(TRUE);
    
    /******************************************************************************************
        If there is a vertex which is joined to more than two other vertices, then the diagram
        is not unique.
    ******************************************************************************************/
    
    for(i = 0; i < Vertices; i++)
    for(j = 0,k = 0; j < Vertices; j++) if(A[i][j] && (++k > 2)) return(TRUE);        
            
    /******************************************************************************************
        If i and j are vertices, which are not inverses, and i and j are joined by more than
        one edge in the Whitehead graph, then the diagram is not unique.
    ******************************************************************************************/
    
    for(i = 0; i < Vertices; i++)
        {
        if(i & 1)
            {
            for(j = i + 1; j < Vertices; j++) if(A[i][j] > 1) return(TRUE);
            }
        else
            {
            for(j = i + 2; j < Vertices; j++) if(A[i][j] > 1) return(TRUE);
            }
        }
    return(FALSE);
}

void Debug(void)
{

    register int    i,
                    j;
    
    /******************************************************************************************
                                Print the array CO[i][j]
           CO[i][j] is the vertex following vertex j in counterclockwise order about vertex i.                            
    ******************************************************************************************/                            
                    printf("\n\nThis is the array CO[i][j].\n");         
                    for(i = 0; i < Vertices; i++)
                        {
                        printf("\n");
                        for(j = 0; j < Vertices; j++)
                            {
                            if(i == j)
                                printf("   *");
                            else
                                printf("%4u",CO[i][j]);
                            }     
                        }
                    printf("\n");        
                                                        
    /******************************************************************************************
                                    Print the array LR[i].
                            LR[i] gives the length of the ith relator.
    ******************************************************************************************/
    
                    printf("\n\nThis is the array LR[i].\n\n");    
                    for(i = 1; i <= NumRelators; i++) printf("%4lu",LR[i]);
                    printf("\n");
                            
    /******************************************************************************************
                                    Print the array EXP[i][j].
        The three columns of EXP[i][j] give the exponents in ascending order with which
        generator i appears in the set of relators.
    ******************************************************************************************/
    
                    printf("\n\nThis is the array EXP[i][j].\n");    
                    for(i = 0; i < NumGenerators; i++)
                        {
                        printf("\n");
                        for(j = 0; j < 3; j++) printf("%4u",EXP[i][j]);
                        }
                    printf("\n");
                    
    /******************************************************************************************
                                    Print the array NEX[i][j].
        NEX[i][j] gives the number of times that the exponent EXP[i][j] occurs in the relators.
    ******************************************************************************************/
    
                    printf("\n\nThis is the array NEX[i][j].\n");    
                    for(i = 0; i < NumGenerators; i++)
                        {
                        printf("\n");
                        for(j = 0; j < 3; j++) printf("%4u",NEX[i][j]);
                        }
                    printf("\n");
                            
    /******************************************************************************************
                                    Print the array DF[i].
            DF[] is an array of flags. The default value is 0. If the maximal exponent
        with which generator i appears is 2 or less and both vertex 2i and vertex 2i+1 have
        valence 2, then DF[i] is 1. Otherwise, the maximal exponent with which generator i
        appears is greater than 2 and DF[i] is 2 unless both vertex 2i and vertex 2i+1 have
        valence 2, in which case, DF[i] is 3.                                    
    ******************************************************************************************/
    
                    printf("\n\nThis is the array DF[i].\n\n");    
                    for(i = 0; i < NumGenerators; i++) printf("%4u",DF[i]);
                    printf("\n");
    
    /******************************************************************************************
                                    Print the array PG[i].                        
    ******************************************************************************************/
    
                    fprintf(stdout,"\n\nThis is the array PG[i].\n\n");    
                    for(i = 0; i < NumGenerators; i++) fprintf(stdout,"%4u",PG[i]);
                    fprintf(stdout,"\n");
                                    
    /******************************************************************************************
                                    Print the array T[i][j].
    ******************************************************************************************/    
    
                    printf("\n\nThis is the array T[i][j].\n");
                    for(i = 0; i < NumGenerators; i++)
                        {
                        printf("\n");
                        for(j = 0; j < 8; j++) printf("%c ",T[i][j]);
                        }
                    printf("\n");
                    
    /******************************************************************************************                                        
                                    Print the array ED[i][j].
    ******************************************************************************************/    
                                                        
                    printf("\n\nThis is the array ED[i][j].\n");
                    for(i = 0; i < NumGenerators; i++)
                        {
                        printf("\n");
                        for(j = 0; j < Vertices; j++)
                            printf("%4u",ED[i][j]);
                        }
                    printf("\n");
                                                        
    /******************************************************************************************
                                    Print the array A[i][j].
    ******************************************************************************************/
        
                    printf("\n\nThis is the array A[i][j].\n");
                    for(i = 0; i < Vertices; i++)
                        {
                        printf("\n");
                        for(j = 0; j < Vertices; j++) printf("%4u",A[i][j]);
                        }
                    printf("\n");
                    
    /******************************************************************************************
                                    Print the array AJ1[i][j].
    ******************************************************************************************/    
    
                    printf("\n\nThis is the array AJ1[i][j].\n");
                    for(i = 0; i < Vertices; i++)
                        {
                        printf("\n");
                        for(j = 0; AJ1[i][j] < VERTICES; j++) printf("%4u",AJ1[i][j]);
                        }
                    printf("\n");                    
    
    /******************************************************************************************
                                    Print the array AJ2[i][j].
    ******************************************************************************************/    
    
                    printf("\n\nThis is the array AJ2[i][j].\n");
                    for(i = 0; i < Vertices; i++)
                        {
                        printf("\n");
                        for(j = 0; AJ2[i][j] < VERTICES; j++) printf("%4u",AJ2[i][j]);
                        }
                    printf("\n");                    
                                    
    /******************************************************************************************
                                    Print the array FV[i].
        FV[i] gives the first vertex in the cyclic ordering of the vertices about vertex i
        in the graph. We have taken this to be the first vertex, in the usual lexicographic
        ordering of the generators, that is joined to vertex i by an edge. 
    ******************************************************************************************/
    
                    printf("\n\nThis is the array FV[i].\n\n");
                    for(i = 0; i < Vertices; i++)
                        printf("%4u",FV[i]);
                    printf("\n");                
    
    /******************************************************************************************                        
                                    Print the array OSA[i].
    ******************************************************************************************/
        
                    printf("\n\nThis is the array OSA[i].\n\n");    
                    for(i = 0; i < Vertices; i++) printf("%4u",OSA[i]);
                    printf("\n");
                    
    /******************************************************************************************                        
                                    Print the array OSB[i].
    ******************************************************************************************/
        
                    printf("\n\nThis is the array OSB[i].\n\n");    
                    for(i = 0; i < Vertices; i++) printf("%4u",OSB[i]);
                    printf("\n");
                            
    /******************************************************************************************
                                    Print the array B[i][j].
    ******************************************************************************************/
                    printf("\n\nThis is the array B[i][j].\n");
                    for(i = 0; i < Vertices; i++)
                        {
                        printf("\n");
                        for(j = 0; j < Vertices; j++)
                            {
                            if(i == j)
                                printf("   *");
                            else
                                printf("%4u",B[i][j]);
                            }    
                        }
                    printf("\n");            
}
