#include "Heegaard.h"
#include "Heegaard_Dec.h"            

/****************************** function prototypes *****************************************
L   24 Reduce_Genus(Input,InitPres,F1)
L  222 Defining_Relator(h,Input,InitPres,F1)
L  603 Proper_Power(void)
L  978 Check_For_Primitives(int Input,int MyNumRelators)
L 1073 Lens_Space_D(int Prim)
L 1310 Lens_Space(void)
L 1510 Transverse(unsigned char *ptr)
L 1817 Recip_Mod_P(p,q)
L 1851 Find_Flow_B(unsigned int Source)
L 1972 Get_Connected_Components(void)
L 2011 LGCD(p,q)
L 2061 Split_At_Empty_Relators(F1)
L 2472 Split_At_Empty_Relators_Sub1(void)
L 2537 Split_At_Empty_Relators_Sub2(int NumGen,int NumNewPres)
L 2648 Splitting_Pres_On_File(int WhoCalled,int NumNewPres)
********************************************************************************************/

int InitialNumGenerators;
    
unsigned int Reduce_Genus(int Input,int InitPres,int F1)
{    
    /******************************************************************************************
        Reduce_Genus() first calls Defining_Relator() which looks for defining relators among
        the relators. Defining relator removes defining relators, making appropriate
        substitutions and updates, and returns when there are no more defining relators among
        the relators. Reduce_Genus() then calls Check_For_Primitives() which looks for relators
        which are primitives or proper powers. If any of these are found they are removed and
        appropriate reductions made to the number of generators and to the number of relators.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q,
                            **Temp;
                            
    int                     k,
                            InitialNumGenerators,
                            SMicro_Print;
                            
    unsigned int            j;                        
    
    unsigned long           STotalAuts;
    
    unsigned int Lens_Space();
    unsigned int Proper_Power();
    
    /******************************************************************************************
        A) If Do_Not_Reduce_Genus is true, Heegaard is in a mode where it will not
        look for primitives, defining relators, proper powers of generators or lens spaces.
        B) If Delete_Only_Short_Primitives is true, Heegaard is in a mode where
        it will delete only relators of length one or two from the presentation.
    ******************************************************************************************/    

    if(Do_Not_Reduce_Genus) return(0);
    
    InitialNumGenerators = NumGenerators;
    FoundPrimitive = FoundPower = LensSpace = EmtyRel = FALSE;
    
    if(Input == BANDSUM)
        FoundPrimitive = Defining_Relator(1,Input,InitPres,F1);
    else
        FoundPrimitive = Defining_Relator(0,Input,InitPres,F1);
        
    switch(FoundPrimitive)
        {
        case TOO_LONG:
            return(TOO_LONG);
        case CAN_NOT_DELETE:
            return(CAN_NOT_DELETE);
        }
        
    if(FoundPrimitive || LensSpace || EmtyRel || NumRelators == 1) return(0);
            
    if(Delete_Only_Short_Primitives) return(0);
    
    SMicro_Print = Micro_Print;
    Micro_Print = FALSE;
    STotalAuts = TotalAuts;                
    k = Check_For_Primitives(Input,NumRelators);
    TotalAuts = STotalAuts;
    Micro_Print = SMicro_Print;
        
    if(k == TOO_LONG) return(TOO_LONG);
    
    if(k > 0) 
        {
        if(Micro_Print)
            {
            printf("\n\nRelator %d is a primitive.",k);
            if(k > 1)
                printf(" Swapped Relator %d and Relator 1.",k);
            printf("\n");   
            }
            
        FoundPrimitive = TRUE;
        
        if(k > 1)
            {
            Temp = Relators[1];
            Relators[1] = Relators[k];
            Relators[k] = Temp;
            }
            
        if(NumGenerators == 2 && NumRelators == 2)
            {
            /******************************************************************************
                If there are two generators and two relators, this may be a lens space.
                If InitialNumGenerators = NumGenerators and Input = DUALIZE then this is a
                lens space. Call Lens_Space_D() to determine which lens space it is.
                Otherwise, call Lens_Space() to see if we have a lens space and, if so,
                determine which lens space it is. But first,save copies of the relators.
            ******************************************************************************/    
            
            if(InitialNumGenerators == NumGenerators && Input == DUALIZE)
                {
                if(k > 1)
                    {
                    Temp = Relators[1];
                    Relators[1] = Relators[k];
                    Relators[k] = Temp;
                    }
                if(Lens_Space_D(k) == TOO_LONG) return(TOO_LONG);
                LensSpace = TRUE;
                return(0);    
                }
            
            if(Relators[3] != NULL) DisposeHandle((char **) Relators[3]);
            Relators[3] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[1]));  
            if(Relators[3] == NULL) Mem_Error();
            p = *Relators[3];
            q = *Relators[1];
            while( (*p++ = *q++) ) ;
            if(Relators[4] != NULL) DisposeHandle((char **) Relators[4]);
            Relators[4] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[2]));
            if(Relators[4] == NULL) Mem_Error();
            p = *Relators[4];
            q = *Relators[2];
            while( (*p++ = *q++) ) ;
            
            switch(Lens_Space())
                {
                case NO_ERROR:
                    LensSpace = TRUE;                     /*** We found the Lens Space. ***/
                    return(0);
                case NOT_CONNECTED:
                    LensSpace = NOT_CONNECTED;        /*** The diagram is not connected. ***/
                    return(0);
                default:
                	if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
                	Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[3]));
                    if(Relators[1] == NULL) Mem_Error();
                    p = *Relators[1];
                    q = *Relators[3];
                    while( (*p++ = *q++) ) ;
                    if(Relators[2] != NULL) DisposeHandle((char **) Relators[2]);
                    Relators[2] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[4]));
                    if(Relators[2] == NULL) Mem_Error();
                    p = *Relators[2];
                    q = *Relators[4];
                    while( (*p++ = *q++) ) ;
                    
                    j = Defining_Relator(0,Input,InitPres,F1);    
                    if(j == TOO_LONG) return(TOO_LONG);
                    if(j == CAN_NOT_DELETE) return(CAN_NOT_DELETE);
                    if(j) return(0);
                }    
            }
                
        if(NumGenerators > 1)
            {
            if(Find_Primitives(1) == TOO_LONG) return(TOO_LONG);
            j = Defining_Relator(1,Input,InitPres,F1);
            if(j == TOO_LONG) return(TOO_LONG);
            if(j == CAN_NOT_DELETE) return(CAN_NOT_DELETE);
            if(EmtyRel) return(0);
            }                
        }
        
    if(k < 0)
        {
        k = -k;
        
        if(Micro_Print)
            {
            printf("\n\nRelator %d is a proper power of a free generator.",k);
            if(k > 1)
                printf(" Swapped Relator %d and Relator 1.",k);
            }
            
        FoundPower = TRUE;
            
        if(k > 1)
            {
            Temp = Relators[1];
            Relators[1] = Relators[k];
            Relators[k] = Temp;
            }
            
        if(NumGenerators > 1)
            {
            if(Find_Primitives(3) == TOO_LONG) return(TOO_LONG);
            switch(Proper_Power())
                {
                case TOO_LONG:
                    return(TOO_LONG);
                case CAN_NOT_DELETE:
                    return(CAN_NOT_DELETE);    
                case FATAL_ERROR:
                    return(FATAL_ERROR);
                case NO_ERROR:
                    break;        
                }
            }            
        }
                
    return(0);    
}
        
int Defining_Relator(int h,int Input,int InitPres,int F1)
{        
    /******************************************************************************************
        This routine checks the relators for defining relators. If it finds a defining
        relator, it makes the appropriate substitutions in the remaining relators, reduces
        the number of generators and the number of relators and if h = 0, continues until
        there are no defining relators present or EmtyRel is TRUE. If h != 0, the
        routine only checks whether Relators[1] is a defining relator. (This saves some time
        when the presentation differs from a preceeding presentation by a single bandsum,
        which by convention, leaves the new relator in Relators[1].) The routine returns TRUE
        if it finds a defining relator and otherwise returns FALSE.
    ******************************************************************************************/
    
    register unsigned char  s,
                            t,
                            x,
                            y,
                            z,
                            *p,
                            *q,
                            *r;
                            
    unsigned char           **Temp;
    
    int                     PassNum,
                            SNumGenerators,
                            SSNumGenerators;
                            
    unsigned int            C[125],
                            f,
                            i,
                            j,
                            k,
                            RL[MAXNUMRELATORS + 1];
                            
    unsigned long           length,
                            M;
    
    unsigned int            Lens_Space(),
                            Split_At_Empty_Relators();
    
    PassNum = 1;
    if(F1) PassNum = 2;
    SNumGenerators = NumGenerators;
    
RE_ENTER:
    
    /******************************************************************************************
                        If h = 0, test the relators in random order.
    ******************************************************************************************/
        
    for(i = 1; i <= NumRelators; i++) RL[i] = i;
    if(h == 0) for(i = NumRelators; i > 1; i--)
        {
        j = abs(rand()) % i;
        j++;
        k = RL[i];
        RL[i] = RL[j];
        RL[j] = k;
        }    
                
    for(f = 1;f <= NumRelators ; f++)
        {
        if(NumGenerators <= 1 || NumRelators <= 1) break;
                                                          
        /**************************************************************************************
                If there is only one generator or only one relator we return because there is 
            nothing to do. (If there is a single relator and it is a primitive, that fact
            will be detected elsewhere.)
        **************************************************************************************/
            
        if(h && PassNum == 2 && f > 1) break;
        i = RL[f];
        if((Delete_Only_Short_Primitives || PassNum != 2) && GetHandleSize((char **) Relators[i]) > 3)
            continue;    
        
        /**************************************************************************************            
                Look for a generator which appears only once in Relators[i].
        **************************************************************************************/
        
        for(x = 'A',y = 'a'; x < 'A' + NumGenerators; x++,y++)    C[x] = C[y] = 0;
        p = *Relators[i];
        while( (z = *p++) ) C[z]++;
        for(x = 'A',y = 'a',(j = 0); x < 'A' + NumGenerators; x++,y++) if(C[x] + C[y] == 1) j++;
        if(j == 0) continue;
        k = abs(rand()) % j;
        k++;
        for(x = 'A',y = 'a',j = 0; x < 'A' + NumGenerators; x++,y++)
            if(C[x] + C[y] == 1 && ++j == k) break;
            
        if(j)
            {
            /**********************************************************************************             
                                Relators[i] is a defining relator.
                    Next, find the location of a defining generator in Relators[i].
                If there are two generators and two relators, this is probably a lens
                space. If InitialNumGenerators = NumGenerators and Input = DUALIZE, then
                this is indeed a lens space. Call Lens_Space_D() to determine which lens
                space it is. Otherwise, call Lens_Space() to see if we have a lens space,
                and, if so, determine which lens space it is. But first, save copies
                of the relators.        
            **********************************************************************************/
            
            if(NumGenerators == 2 && NumRelators == 2)
                {
                if(InitialNumGenerators == NumGenerators && Input == DUALIZE)
                    {
                    if(Lens_Space_D(i) == TOO_LONG) return(TOO_LONG);
                    LensSpace = TRUE;
                    return(FALSE);
                    }
                            
                if(i == 2)
                    {
                    Temp = Relators[2];
                    Relators[2] = Relators[1];
                    Relators[1] = Temp;
                    }
                
                if(Relators[3] != NULL) DisposeHandle((char **) Relators[3]); 
                Relators[3] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[1]));   
                if(Relators[3] == NULL) Mem_Error();
                p = *Relators[3];
                q = *Relators[1];
                while( (*p++ = *q++) ) ;
                if(Relators[4] != NULL) DisposeHandle((char **) Relators[4]);
                Relators[4] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[2]));
                if(Relators[4] == NULL) Mem_Error();
                p = *Relators[4];
                q = *Relators[2];
                while( (*p++ = *q++) ) ;
                    
                switch(Lens_Space())
                    {
                    case NO_ERROR:
                        LensSpace = TRUE;                 /*** We found the Lens Space. ***/
                        return(FALSE);
                    case NOT_CONNECTED:
                        LensSpace = NOT_CONNECTED;    /*** The diagram is not connected. ***/
                        return(FALSE);
                    }
                
                if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
                Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[3]));  
                if(Relators[1] == NULL) Mem_Error();
                p = *Relators[1];
                q = *Relators[3];
                while( (*p++ = *q++) ) ;
                if(Relators[2] != NULL) DisposeHandle((char **) Relators[2]);
                Relators[2] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[4]));
                if(Relators[2] == NULL) Mem_Error();
                p = *Relators[2];
                q = *Relators[4];
                while( (*p++ = *q++) ) ;
                        
                if(i == 2)
                    {
                    Temp = Relators[2];
                    Relators[2] = Relators[1];
                    Relators[1] = Temp;
                    }
                }
            
            /**********************************************************************************
                Save a copy of the current relators in Copy_Of_Rel_1[]. We will need this
                copy of the current set of relators if an empty relator appears.
            **********************************************************************************/
                
                for(j = 1; j <= NumRelators; j++)
                    {
                    if(Copy_Of_Rel_1[j] != NULL) DisposeHandle((char **) Copy_Of_Rel_1[j]);
                    Copy_Of_Rel_1[j] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[j]));
                    if(Copy_Of_Rel_1[j] == NULL) Mem_Error();
                    p = *Copy_Of_Rel_1[j];
                    q = *Relators[j];
                    while( (*p++ = *q++) ) ;                    
                    }
                                        
            /**********************************************************************************     
                    Create two strings which will be substituted for the defined
                                 generator and its inverse.                                            
            **********************************************************************************/

            if(Micro_Print)
                {
                printf("\n\nRelator %d is a defining relator which was used to eliminate generator %c.",i,x);
                if(i != NumRelators)
                    printf("\nSwapped Relator %d with Relator %d and reduced the number of relators.",i,NumRelators);
                }
                
            p = *Relators[i];
            while( (z = *p++) ) if(z == x || z == y) break;
            if(Temp6 != NULL) DisposeHandle((char **) Temp6);
            Temp6 = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]) -1);
            if(Temp6 == NULL) Mem_Error();
            q = *Temp6;
            while(*p) *q++ = *p++;
            p = *Relators[i];
            while( (z != *p) ) *q++ = *p++;
            *q = EOS;
            p = *Temp6;
            if(Temp7 != NULL) DisposeHandle((char **) Temp7);
            Temp7 = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]) -1);
            if(Temp7 == NULL) Mem_Error();
            q = *Temp7;
            while( (*q++ = *p++) ) ;
            if(z == x) Inverse(*Temp6);
            else Inverse(*Temp7);
                
            /**********************************************************************************
                 Next, if i < NumRelators, exchange Relator[NumRelators] and Relator[i]. Also
                 exchange Copy_Of_Rel_1[NumRelators] and Copy_Of_Rel_1[i].
             *********************************************************************************/
             
            if(i < NumRelators)
                {
                Temp = Relators[i];
                Relators[i] = Relators[NumRelators];
                Relators[NumRelators] = Temp;
                Temp = Copy_Of_Rel_1[i];
                Copy_Of_Rel_1[i] = Copy_Of_Rel_1[NumRelators];
                Copy_Of_Rel_1[NumRelators] = Temp;                
                }
                
            /**********************************************************************************
                            Make the appropriate substitutions in the relators.         
            **********************************************************************************/
            
            for(j = 1; j < NumRelators;    j++)
                {
                k = 0;    
                p = *Relators[j];
                while( (z = *p++) ) if(z == x || z == y) k++;
                if(k)
                    {
                    M = GetHandleSize((char **) Relators[j]) + k*GetHandleSize((char **) Relators[NumRelators]);
                    M -= 3*k;
                    if(M > MAXLENGTH)
                        {
                        if(Batch == FALSE)
                        printf("\nA relator derived from Pres %d is too long after deleting a primitive.",ReadPres + 1);
                        NumRelTooLong ++;
                        return(TOO_LONG);
                        }
                    if(Temp4 != NULL) DisposeHandle((char **) Temp4);
                    Temp4 = (unsigned char **) NewHandle(M);
                    if(Temp4 == NULL) Mem_Error();
                    q = *Temp4;    
                    p = *Relators[j];
                    while( (z = *p++) )
                        {
                        *q++ = z;
                        if(z == x)
                            {
                            q--;
                            r = *Temp6;
                            while( (*q++ = *r++) ) ;
                            q--;
                            }
                        if(z == y)
                            {
                            q--;
                            r = *Temp7;
                            while( (*q++ = *r++) ) ;
                            q--;
                            }
                        }
                    *q = EOS;    
                    Temp = Relators[j];
                    Relators[j] = Temp4;
                    Temp4 = Temp;                    
                    }        
                }
                
            /**********************************************************************************
                        If x is not the largest generator, replace the largest generator
                        throughout with x and its inverse with y.                                
            **********************************************************************************/
                     
            if(x < NumGenerators + 64)
                {    
                s = NumGenerators + 64;
                t = NumGenerators + 96;
                for(j = 1; j < NumRelators; j++)
                    {
                    p = *Relators[j];
                    while( (z = *p) )
                        {
                        if(z == s) *p = x;
                        if(z == t) *p = y;
                        p++;
                        }
                    }
                if(Micro_Print)
                    printf("\n\nRewrote the relators by replacing %c with %c.",s,x);
                }    
                
            /**********************************************************************************
                        Update the variables NumGenerators, NumRelators and Vertices.
            **********************************************************************************/
            
            NumGenerators --;
            NumRelators --;
            Vertices -= 2;
            
            /**********************************************************************************
                                    Freely reduce the new relators.
            **********************************************************************************/
            
            if(Micro_Print)
                {
                for(k = 1,length = 0L; k <= NumRelators; k++)
                    length += GetHandleSize((char **) Relators[k]);
                length -= NumRelators;
                }    
                
            k = Freely_Reduce();
            if(k == TOO_LONG) return(TOO_LONG);
            Length = OrigLength;
            
            /**********************************************************************************
                We may have created some empty relators. This can only occur when there are
                relators which are consequences of the primitive relator we just deleted. If
                such relators exist, we call Split_At_Empty_Relators(). This routine attempts
                to find a Heegaard diagram which will allow us to determine where these empty
                relators lie on the Heegaard surface.
            **********************************************************************************/
            
            if(k & !InitPres)
                {
                if(BDY[ReadPres] == FALSE && NumGenerators == NumRelators) return(TRUE);
                SSNumGenerators = NumGenerators;
                if(Split_At_Empty_Relators(FALSE))
                    {
                    ReadPres = SReadPres;
                    DrawingDiagrams = FALSE;
                    TestRealizability3 = FALSE;
                    EmtyRel = FALSE;
                    if(Micro_Print)
                        {
                        printf("\n\nDeleting the previous primitive created some empty relators. This means the manifold");
                        printf("\nis a connected sum. However, Heegaard could not determine what the summands are.");
                        }                    
                    return(CAN_NOT_DELETE);
                    }
                ReadPres = SReadPres;    
                DrawingDiagrams = FALSE;
                TestRealizability3 = FALSE;                    
                if(TestRealizability1 && !TestRealizability2)
                    {
                    if(SSNumGenerators < SNumGenerators) return(TRUE);
                    return(FALSE);
                    }
                EmtyRel = TRUE;                    
                }
            if(k & InitPres) EmtyRel = TRUE;   
            
            if(Micro_Print && length == Length && !EmtyRel)
                {
                printf("\n\nThe presentation is currently:\n");
                Print_Relators(Relators,NumRelators);
                }
            
            h = 0;    
            if(EmtyRel) return(TRUE);
            if(InitPres && PassNum == 1) return(TRUE);
            if(PassNum != 2) goto RE_ENTER;
            if(PassNum == 2)
                {
                if(InitPres) return(TRUE);
                PassNum = 3;
                goto RE_ENTER;
                }                                                     
            }    
        }
    if(InitPres && PassNum == 2) return(FALSE);    
    if(++PassNum < 4) goto RE_ENTER;
    if(NumGenerators < SNumGenerators) return(TRUE);
    return(FALSE);        
}        

unsigned int Proper_Power()
{        
    /******************************************************************************************
        This routine is called when Relators[1] is a proper power. In this case, the routine
        deletes appearances of Relators[1] and its inverse from the remaining relators. If the
        presentation is one which Heegaard knows is realizable, then this process is fairly
        simple. If the presentation is not known to be realizable, then Heegaard must run
        some further checks in order to insure that it does not transform an unrealizable
        presentation into a realizable presentation.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q,
                            x,
                            y,
                            z;
                            
    register unsigned int   j,
                            length;                                            
                            
    unsigned char           SaveMinV,
                            SaveMinZ,
                            **Temp,
                            v;
        
    unsigned int            h,
                            i,
                            MaxExp,
                            MinExp,
                            VL,
                            VLI;
    
    int                     NumComps,
                            OrigNumComps;
                            
    unsigned long           HS,
    						jj,
   							longlength;                                
    
    unsigned int            Split_At_Empty_Relators();
    
    longlength = GetHandleSize((char **) Relators[1]) - 1;
    if(longlength > MAXLENGTH) return(TOO_LONG);
    length = longlength;
                
    /******************************************************************************************         
                    Look for the generator which appears in Relators[1].
                    Suppose this generator is G and g is its inverse.
    ******************************************************************************************/
            
    p = *Relators[1];
    z = *p;
    if(z < 'a')
        {
        x = z;
        y = x + 32;
        }
    else
        {
        y = z;
        x = y - 32;
        }
                
    /******************************************************************************************
                         Exchange Relators[NumRelators] and Relators[1].
    ******************************************************************************************/
        
    Temp = Relators[1];
    Relators[1] = Relators[NumRelators];
    Relators[NumRelators] = Temp;
    
    if(Micro_Print)
        printf("\n\nSwapped Relator 1 and Relator %d.",NumRelators);

    /******************************************************************************************
        Save a copy of the current relators in Copy_Of_Rel_1[]. We will need this copy of
        the current set of relators if an empty relator appears.
    ******************************************************************************************/
        
    for(j = 1; j <= NumRelators; j++)
        {
        if(Copy_Of_Rel_1[j] != NULL) DisposeHandle((char **) Copy_Of_Rel_1[j]);
        Copy_Of_Rel_1[j] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[j]));
        if(Copy_Of_Rel_1[j] == NULL) Mem_Error();
        p = *Copy_Of_Rel_1[j];
        q = *Relators[j];
        while( (*p++ = *q++) ) ;                    
        }
                
    /******************************************************************************************
            Note that if we already know that the given presentation is realizable, then we
            can immediately delete all appearances of the generator G and its inverse g from
            the remaining relators and we are guaranteed that the resulting presentation is
            again realizable. Otherwise, it is necessary to run the following section of code.
    ******************************************************************************************/    
    
    if(TestRealizability1)
        {
        /**************************************************************************************
            Call Find_Flow_B() to see if there are T-transformations which reduce the total
            number of appearances of G and g in the relators.
        **************************************************************************************/    

        i = x - 65;
        if(Find_Flow_B(2*i) == TOO_LONG) return(TOO_LONG);
        if(VA[i] == length) return(NO_ERROR);
        if((VA[i] % length) != 0) return(FATAL_ERROR);
        Get_Matrix();
        for(j = 0; j < Vertices; j++) ZZ[j] = 0;
        OrigNumComps = Get_Connected_Components();
        for(j = 0; j < Vertices; j++) ZZ[j] = 0;
        ZZ[2*i] = ZZ[2*i + 1] = INFINITE;
        NumComps = Get_Connected_Components();
        }
        
    while(TestRealizability1)
        {
        /**************************************************************************************
            Find the smallest exponent and the largest exponent with which G appears in the
            relators and also find a corresponding component of the separation if G and g are
            a separating pair of vertices.
        **************************************************************************************/    
        
        for(i = 1,MinExp = INFINITE,MaxExp = 0; i < NumRelators; i++)
            {
            p = *Relators[i];
            h = 0;
            while(*p == x || *p == y)
                {
                p++;
                h++;
                }
            if(*p == EOS)
                {
                if(h != length)
                    return(FATAL_ERROR);
                else
                    {
                    if(length < MinExp)
                        MinExp = length;
                    if(length > MaxExp)
                        MaxExp = length;    
                    continue;
                    }
                }        
            while( (z = *p) )
                {
                p++;
                if(z != x && z != y)
                    {
                    if(*p == EOS && h > 0)
                        {
                        if(h < MinExp)
                            {
                            MinExp = h;
                            SaveMinV = z;
                            SaveMinZ = **Relators[i];
                            }
                        if(h > MaxExp)
                            MaxExp = h;
                        }
                    if(*p == x || *p == y)
                        {
                        j = 0;
                        v = z;
                        z = *p;
                        while(z == *p)
                            {
                            j++;
                            p++;
                            }
                        if(*p)
                            {
                            if(j < MinExp)
                                {
                                MinExp = j;
                                SaveMinV = v;
                                SaveMinZ = z;
                                }
                            if(j > MaxExp)
                                MaxExp = j;
                            }
                        else
                            {        
                            if(j + h < MinExp)
                                {
                                MinExp = j + h;
                                SaveMinV = v;
                                SaveMinZ = z;
                                }
                            if(j + h > MaxExp)
                                MaxExp = j + h;
                            }    
                        }
                    }    
                }
            }
            
        /**************************************************************************************
            Note that if the presentation is realizable, and the pair of vertices G and g do
            not separate the Whitehead graph of the presentation, then it is necessary that
            MinExp = MaxExp = length.
        **************************************************************************************/        
        
        if(MinExp > length || MaxExp > length) return(FATAL_ERROR);
        if(MinExp == length) break;
        if(NumComps == OrigNumComps)
            {
            if(MinExp != length || MaxExp != length) return(FATAL_ERROR);
            break;
            }

        /**************************************************************************************
            Since NumComps > OrigNumComps, the pair of vertices G and g separate the Whitehead
            graph of the presentation. Since we also have MinExp < length, if the presentation
            is realizable, there must exist a component of this separation which we can
            "slide around" so that it amalgamates with another component -- thus reducing
            NumComps. It follows that: if the presentation is realizable, then after a number
            of iterations of this while loop, MinExp must increase and eventually come to equal
            length. 
            So, find the component of the separation corresponding to SaveMinV and set up the
            array ZZ[] for Do_Aut(). Then call Do_Aut() MinExp times.
        **************************************************************************************/
        
        if(SaveMinZ < 'a')
            VL = 2*SaveMinZ - 130;
        else    
            VL = 2*SaveMinZ - 193;     
            
        if(VL & 1)
            VLI = VL - 1;
        else
            VLI = VL + 1;
            
        if(SaveMinV < 'a')
            j = ZZ[2*SaveMinV - 129];
        else
            j = ZZ[2*SaveMinV - 194];
        
        for(i = 0; i < Vertices; i++)
            {
            if(ZZ[i] == j)
                ZZ[i] = 1;
            else
                ZZ[i] = 0;
            }
        ZZ[VL] = 1;
                    
        if(!(VL & 1))
            for(i = 0; i < Vertices; i++) 
                {
                if(ZZ[i])
                    ZZ[i] = 0;
                else
                    ZZ[i] = 1;
                }
        
        if(Micro_Print)
            printf("\n\nPerforming level-transformations in order to eliminate a proper power.");
                                    
        if(VL & 1)
            {
            if(MinExp > 2)
                Do_Auts(VLI,MinExp,NumRelators);
            else
                Do_Aut(VLI,MinExp,NumRelators);
            }
        else
            {
            if(MinExp > 2)
                Do_Auts(VL,MinExp,NumRelators);
            else
                Do_Aut(VL,MinExp,NumRelators);
            }    
    
        /**************************************************************************************
            If the presentation is realizable, then the number of components of the
            separation must decrease by one after we have dragged this component around. So
            recompute NumComps and check whether it has decreased. If NumComps has not
            decreased, then the presentation is not realizable, and we return FATAL_ERROR.
        **************************************************************************************/
        
        i = x - 65;    
        for(j = 0; j < Vertices; j++) ZZ[j] = 0;
        ZZ[2*i] = ZZ[2*i + 1] = INFINITE;
        Fill_A(NumRelators);
        Get_Matrix();
        j = Get_Connected_Components();
        if(j >= NumComps) return(FATAL_ERROR);
        NumComps = j;            
        }
                    
    /******************************************************************************************
            Delete all appearances of G and g from Relators[i] for 1 <= i < NumRelators.                         
    ******************************************************************************************/
    
    for(i = 1; i < NumRelators;    i++)
        {
        HS = GetHandleSize((char **) Relators[i]);
        if(HS <= 1) continue;
        if(Temp4 != NULL) DisposeHandle((char **) Temp4);
        Temp4 = (unsigned char **) NewHandle(HS);
        if(Temp4 == NULL) Mem_Error();
        q = *Temp4;      
        p = *Relators[i];
        z = *p;
        jj = 1;
        while( (z = *p++) )
            {
            if(z != x && z != y)
                {    
                *q++ = z;
                jj++;
                }
            }
        *q = EOS;
        if(jj < HS)
        	{
        	if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
        	Relators[i] = (unsigned char **) NewHandle(jj);
        	if(Relators[i] == NULL) Mem_Error();
        	p = *Relators[i];
        	q = *Temp4;
        	while( (*p++ = *q++) ) ;
        	}                           
        }
            
    /******************************************************************************************
                                Freely reduce the new relators.
    ******************************************************************************************/
    
    if(Micro_Print)
        {
        printf("\n\nThe presentation is currently:\n");
        Print_Relators(Relators,NumRelators);
        }
        
    i = Freely_Reduce();
    Length = OrigLength;        
    if(i == TOO_LONG) return(TOO_LONG);

    
    /******************************************************************************************
            We may have created some empty relators; i.e. relators which were consequences of
            the proper power relator. If so, call Split_At_Empty_Relators(TRUE) to remove any
            empty relators we have created.
    ******************************************************************************************/
                
    if(i)
        {
        if(BDY[ReadPres] == FALSE && NumGenerators == NumRelators)
            return(NO_ERROR);
        if(Split_At_Empty_Relators(TRUE))
            {
            EmtyRel = FALSE;
            DrawingDiagrams = FALSE;
            TestRealizability3 = FALSE;                
            ReadPres = SReadPres;
            if(Micro_Print)
                {
                printf("\n\nDeleting the previous proper power created some empty relators. This means the manifold");
                printf("\nis a connected sum. However, Heegaard could not determine what the summands are.");
                }                            
            return(CAN_NOT_DELETE);
            }
        ReadPres = SReadPres;    
        DrawingDiagrams = FALSE;
        TestRealizability3 = FALSE;            
        if(TestRealizability1 && !TestRealizability2) return(CAN_NOT_DELETE);            
        EmtyRel = TRUE;                    
        }
            
    return(NO_ERROR);    
}            
    
int Check_For_Primitives(int Input,int MyNumRelators)
{        
    /******************************************************************************************
        If Input = BANDSUM, this routine checks whether the relator Relators[1] is a primitive
        or a proper power of a free generator. Otherwise the routine takes each relator in turn
        and calls CheckPrimitivity() to determine if that relator is a primitive or proper
        power of a free generator. It returns i if Relators[i] is primitive, returns -i if
        Relators[i] is a proper power and otherwise returns 0 or TOO_LONG.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q;
                            
    register int            i,
                            j;         
    
    int                     k,
                            RL[MAXNUMRELATORS + 1];
    
    if(NumGenerators == 1) return(0);

    /******************************************************************************************
              Since Relators[1] gets trashed by Find_Primitives() etc., we need to make a copy.
    ******************************************************************************************/    
    
    if(Temp11 != NULL) DisposeHandle((char **) Temp11);
    Temp11 = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[1]));
    if(Temp11 == NULL) Mem_Error();
    q = *Temp11;
    p = *Relators[1];
    while( (*q++ = *p++) ) ;    
    
    if(Input == BANDSUM)
        {
        j = CheckPrimitivity();
        if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
        Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Temp11));
        if(Relators[1] == NULL) Mem_Error();
        p = *Relators[1];
        q = *Temp11;
        while( (*p++ = *q++) ) ;
        if(j == TOO_LONG) return(TOO_LONG);
        if(j > 0) return(1);
        if(j < 0) return(-1);
        return(0);                        
        }
    
    /******************************************************************************************
                        If Input != 1 test the relators in random order. 
    ******************************************************************************************/    
    
    for(i = 1; i <= MyNumRelators; i++) RL[i] = i;
            
    if(Input != 1) for(i = MyNumRelators; i > 1; i--)
        {
        j = abs(rand()) % i;
        j++;
        k = RL[i];
        RL[i] = RL[j];
        RL[j] = k;
        }
    
    for(i = 1; i <= MyNumRelators ; i++)
        {
        if(RL[i] == 1)
            {
            if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
            Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Temp11));
            if(Relators[1] == NULL) Mem_Error();
            p = *Temp11;            
            }
        else
            { 
            if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
            Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[RL[i]]));
            if(Relators[1] == NULL) Mem_Error();
            p = *Relators[RL[i]];
            }
        q = *Relators[1];
        while( (*q++ = *p++) ) ;
        if( (j = CheckPrimitivity()) ) break;    
        }
    
    if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
    Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Temp11));  
    if(Relators[1] == NULL) Mem_Error();
    p = *Relators[1];
    q = *Temp11;
    while( (*p++ = *q++) ) ;
    if(j == TOO_LONG) return(TOO_LONG);
    if(j > 0) return(RL[i]);
    if(j < 0) return(-RL[i]);        
    return(0);                            
}

int Lens_Space_D(int Prim)
{
    /******************************************************************************************
        Lens_Space_D() is called when we have a two generator, two relator presentation of a
        closed manifold and one of the dualrelators is a primitive. Thus the manifold is a lens
        space. This routine can then determine which lens space we have by a procedure
        similar to that used by Lens_Space(). (Note that Relators[] and DualRelators[] have
        been swapped here since this presentation comes from dualizing!)
    ******************************************************************************************/
        
    register unsigned char  x,
                            *p,
                            *q;
    
    unsigned char           **Temp;
    
    int                     SMicro_Print;
                            
    unsigned int            edge,
                            ee,
                            s,
                            u,
                            v,
                            vertex,
                            w;
    
    long                    a,
                            b,
                            c,
                            d,
                            e,
                            f,
                            *ptr = NULL,
                            t;

    unsigned long           Recip_Mod_P();
    unsigned int            Whitehead_Graph();
    
    /******************************************************************************************
        Look for a simple closed curve on the Heegaard surface which is transverse to
        the "non-primitive" dual relator at a single point and which is disjoint from the
        "primitive" dual relator. But first, because some free reductions may have been
        performed, we need to go back to the original diagram and recover a set of dual
        relators which have not been freely reduced.
    ******************************************************************************************/    

    for(v = 1; v <= 2; v++)
        {
        Temp = Relators[v];
        Relators[v] = DualRelators[v];
        DualRelators[v] = Temp;
        }    
    
    for(v = 1,Length = 0L; v <= 2; v++) Length += GetHandleSize((char **) Relators[v]);
    Length -= 2;
        
    SMicro_Print = Micro_Print;
    Micro_Print = FALSE;            
    if(Find_Flow_A(NORMAL,FALSE))
        {
        Micro_Print = SMicro_Print;
        return(TOO_LONG);
        }
    if(Whitehead_Graph())
        {
        Micro_Print = SMicro_Print;        
        return(TOO_LONG);
        }
    Micro_Print = SMicro_Print;    
    
    for(v = 1; v <= 2; v++)
        {
        Temp = Relators[v];
        Relators[v] = DualRelators[v];
        DualRelators[v] = Temp;
        }    
        
    /*****************************************************************************************
        The DualRelators now have the same length as the Relators i.e. no free reductions
        have been performed.
    ******************************************************************************************/            
    
    if(Prim == 1)
        vertex = 2;
    else
        vertex = 0;
    
    edge = 0;
    s = FV[vertex];
    while(1)
        {
        ee = edge;
        w = s;    
        v = vertex;
        while(1)
            {
            ee = B[w][v] - ee;
            if(w == vertex) break;
            if(w == (vertex + 1)) goto END;
            u = CO[w][v];
            v = w;
            w = u;
            ee++;
            if(ee >= V[v]) ee -= V[v];        
            }
        edge += A[vertex][s];
        s = CO[vertex][s];    
        }
END:        
    ee = OSA[w] - ee;
    if(ee >= V[w]) ee -= V[w];
    
    if(edge > ee)
        {
        s = ee;
        ee = edge;
        edge = s;
        }
    
    s = ee - edge;
    if(Temp13 != NULL) DisposeHandle((char **) Temp13);
    Temp13 = (unsigned char **) NewHandle(s + 1);
    if(Temp13 == NULL) Mem_Error();
    q = *Temp13;
    if(Prim == 2)
        p = *Relators[1];
    else
        p = *Relators[2];
    p += edge;
    for(v = 0; v < s; v++) *q++ = *p++;
    *q = EOS;        
    
    if(Prim == 2)
        {
        Temp = Relators[1];
        Relators[1] = Relators[2];
        Relators[2] = Temp;
        }
                                                
    /******************************************************************************************
        Once we have found the transverse curve, we abelianize it and the two relators R1 and
        R2. Then we can reduce R1 to the canonical form (1,0) or (0,1) by performing column
        operations on the presentation matrix given by R1 and R2 while carrying the abelianized
        transverse curve along for the ride. The image of the transverse curve then gives
        enough information to allow us to determine which particular lens space we have.
    ******************************************************************************************/
    
    ptr = (long *)NewPtr(sizeof(long)*100);
    if(ptr == NULL) Mem_Error();
    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = *Relators[1];
    while( (x = *p++) ) ptr[x] ++;
    a = ptr['A'] - ptr['a'];
    b = ptr['B'] - ptr['b'];
    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = *Relators[2];
    while( (x = *p++) ) ptr[x] ++;    
    c = ptr['A'] - ptr['a'];
    d = ptr['B'] - ptr['b'];
    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = *Temp13;
    while( (x = *p++) ) ptr[x] ++;
    e = ptr['A'] - ptr['a'];
    f = ptr['B'] - ptr['b'];
    DisposePtr((char *) ptr);
    
    if(a < 0)
        {
        a = - a;
        c = - c;
        e = - e;
        }
    if(b < 0)
        {
        b = - b;
        d = - d;
        f = - f;
        }        
    while(1)
        {
        if(b == 0)
            {
            if(d < 0) d = - d;
            if(f < 0) f = - f;
            P = d;    
            Q = f;
            Q = Recip_Mod_P(P,Q);
            break;
            }
        if(a >= b)
            {
            t = a/b;
            a -= t*b;
            c -= t*d;
            e -= t*f;
            }    
        if(a == 0)
            {
            if(c < 0) c = - c;
            if(e < 0) e = - e;
            P = c;    
            Q = e;
            Q = Recip_Mod_P(P,Q);
            break;
            }        
        if(b > a)
            {
            t = b/a;
            b -= t*a;
            d -= t*c;
            f -= t*e;
            }
        }
    
    if(Prim == 2)
        {
        Temp = Relators[1];
        Relators[1] = Relators[2];
        Relators[2] = Temp;
        }
    for(v = 1; v <= 2; v++)
        {
        Temp = Relators[v];
        Relators[v] = DualRelators[v];
        DualRelators[v] = Temp;
        }
    
    if(Micro_Print)
        {
        printf("\n\nThe current presentation is a two generator, two relator presentation of a closed");
        printf("\nmanifold M, for which one of the relators is a primitive. Thus M is a lens space.");
        printf("\nIn particular, M is the lens space L(%lu,%lu).",P,Q);
        printf("\nSwapped the current relators and dual relators.");
        }    
    return(0);                
}

unsigned int Lens_Space()
{
    /******************************************************************************************
        Lens_Space() is called when we have a two generator two relator presentation and one of
        the relators is a primitive. Usually, such a presentation yields a lens space.
        This routine tries to determine whether this is the case, and if we do indeed have a
        lens space, the routine tries to determine which lens space it is.
    ******************************************************************************************/
        
    register int		i,
						j;
                            
    register unsigned char 	x,
	                        *ptr1 = NULL,
							*p;
    
    unsigned char         	**Temp;
    
    int                   	LTRV;
                            
    unsigned int          	r;
    
    long                  	a,
                            b,
                            c,
                            d,
                            e,
                            f,
                            *ptr = NULL,
                            t;

    unsigned long         	Recip_Mod_P();
    unsigned int            Whitehead_Graph();
    
    /******************************************************************************************
         Call Find_Flow to reduce the relators to minimal length before we try to find the
         diagram. Then call Whitehead_Graph() to actually find the diagram. If Whitehead_Graph()
         doesn't return an error, then call Sep_Surface() and check whether this is a diagram
         of a closed manifold.
     *****************************************************************************************/
                             
    if(Find_Flow_A(NORMAL,FALSE)) return(1);
    
_GET_DIAGRAM:
    
    Saved_Vertices = 0;
    switch(r = Whitehead_Graph())
        {
        case NO_ERROR:
            break;
        case SEP_PAIRS:
            Num_Saved_LPres = 0;
            NotNewPres = 0;
            LTRV = Level_Transformations(FALSE,TRUE,FALSE);
            for(i = 0; i < Num_Saved_LPres; i++)
            for(j = 0; j <= NumRelators; j++) if(SLR[i][j] != NULL)
            	{
            	DisposeHandle((char **) SLR[i][j]);
            	SLR[i][j] = NULL;
            	}
            if(Micro_Print)
                {
                printf("\n\nThe current Presentation is:\n");
                Print_Relators(Relators,NumRelators);
                }
                        
            switch(LTRV)
                {
                case 0:
                case 1:
                    return(SEP_PAIRS);
                case 2:
                    /**************************************************************************
                        Lens_Space() expects that the non-primitive relator is Relators[2].
                        The relators may have been swapped by Canonical_Rewrite(). So we need
                        to call Check_For_Primitives() to find out which relator is the 
                        primitive.
                    **************************************************************************/    
                        
                    i = Check_For_Primitives(1,2);
                    if(i == 2)
                        {
                        Temp = Relators[1];
                        Relators[1] = Relators[2];
                        Relators[2] = Temp;
                        }
                    if(i == TOO_LONG) return(TOO_LONG);    
                    goto _GET_DIAGRAM;    
                default:
                    return(SEP_PAIRS);
                }    
        default:
            return(r);
        }                                            
    if(Sep_Surface() > 1) return(1);
                                                                    
    /******************************************************************************************
            We have found the diagram, and the manifold is closed. So this really is a lens
            space. Now look for a curve which is transverse to the "non-primitive" relator at
            a single point and which is disjoint from the "primitive" relator. This curve is a
            core for one of the solid tori composing the lens space.
    ******************************************************************************************/
    
    a = LR[1];
    if(LR[2] > a) a = LR[2];
    a *= 2;
    ptr1 = (unsigned char *) NewPtr(a + 1);
    if(ptr1 == NULL) Mem_Error();
    if(Transverse(ptr1))
        {                            
        DisposePtr((char *) ptr1);
        return(1);
        }
                                            
    /******************************************************************************************
        Once we have found the transverse curve, we abelianize it and the two relators R1 and
        R2. Then we can reduce R1 to the canonical form (1,0) or (0,1) by performing column
        operations on the presentation matrix given by R1 and R2 while carrying the abelianized
        transverse curve along for the ride. The image of the transverse curve then gives
        enough information to allow us to determine which particular lens space we have.
    ******************************************************************************************/
    
    ptr = (long *)NewPtr(sizeof(long)*100);
    if(ptr == NULL) Mem_Error();

    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = *Relators[1];
    while( (x = *p++) ) ptr[x] ++;
    a = ptr['A'] - ptr['a'];
    b = ptr['B'] - ptr['b'];
    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = *Relators[2];
    while( (x = *p++) ) ptr[x] ++;    
    c = ptr['A'] - ptr['a'];
    d = ptr['B'] - ptr['b'];
    ptr['A'] = ptr['B'] = ptr['a'] = ptr['b'] = 0;
    p = ptr1;
    while( (x = *p++) ) ptr[x] ++;
    e = ptr['A'] - ptr['a'];
    f = ptr['B'] - ptr['b'];
    DisposePtr((char *) ptr);
    
    if(a < 0)
        {
        a = - a;
        c = - c;
        e = - e;
        }
    if(b < 0)
        {
        b = - b;
        d = - d;
        f = - f;
        }        
    while(1)
        {
        if(b == 0)
            {
            if(d < 0) d = - d;
            if(f < 0) f = - f;
            P = d;    
            Q = f;
            Q = Recip_Mod_P(P,Q);
            break;
            }
        if(a >= b)
            {
            t = a/b;
            a -= t*b;
            c -= t*d;
            e -= t*f;
            }    
        if(a == 0)
            {
            if(c < 0) c = - c;
            if(e < 0) e = - e;
            P = c;    
            Q = e;
            Q = Recip_Mod_P(P,Q);
            break;
            }        
        if(b > a)
            {
            t = b/a;
            b -= t*a;
            d -= t*c;
            f -= t*e;
            }
        }
                                                                            
    DisposePtr((char *) ptr1);
    if(Micro_Print)
        {
        printf("\n\nThe current presentation is a two generator, two relator presentation of a closed");
        printf("\nmanifold M, in which one of the relators is a primitive. Thus M is a lens space.");
        printf("\nIn particular, M is the lens space L(%lu,%lu).",P,Q);
        }
    return(NO_ERROR);                            
}

int Transverse(ptr)
unsigned char *ptr;
{
    /******************************************************************************************
        This routine is called by Lens_Space(). The routine tries to find a simple closed
        curve on the genus-two Heegaard surface that is transverse to the curve representing
        the "non-primitive" relator and disjoint from the "primitive" relator. If successful,
        it returns the word representing the transverse curve in the string pointed to by ptr.
    ******************************************************************************************/
    
    register unsigned char 	*p,
                            *q,
                            *r;
                            
    register unsigned int  	d,
                            e,
                            v,
                            vv;
    
    unsigned char          	x,
                            y,
                            z,
                            w;
                                                                            
    int                     i,
                            j;
    
    unsigned int          	edge,
                            edge2,
                            vertex;
    
    unsigned long          	max;
    
    max = GetPtrSize((char *) ptr) - 1;
    
    LR[0] = GetHandleSize((char **) OutRelators[1]) - 1;    
    if(abs(Compare(*OutRelators[1])) == 1)
        {
        x = 'B';
        y = 'b';
        z = 'A';
        w = 'a';
        }
    else
        {
        x = 'A';
        y = 'a';
        z = 'B';
        w = 'b';
        }
            
    for(i = 1; i <= 2; i++)
        {
        if(GetHandleSize((char **) DualRelators[i]) < 3) continue;
        for(p = *DualRelators[i], q = p + 1; *p; p++,q++)
            {
            if(!*q) q = *DualRelators[i];
            if((*p == x || *p == y) && *p == *q)
                {
                edge = p - *DualRelators[i];
                e = q - *DualRelators[i];
                if(i == 1)
                    vertex = 0;
                else
                    vertex = 2;    
                v = vertex;
                r = ptr;
                do
                    {
                    if(v & 1)
                        {
                        *r = (v >> 1) + 'a';
                        vv = v - 1;
                        }
                    else
                        {
                        *r = (v >> 1) + 'A';
                        vv = v + 1;
                        }
                    e = OSA[v] - e;
                    if(e >= V[v]) e -= V[v];
                    r++;
                    if(r - ptr > max )
                        {
                        r--;
                        *r = EOS;
                        return(5001);
                        }
                    v = FV[vv];
                    d = A[vv][v];
                    while(d <= e)
                        {
                        v = CO[vv][v];
                        d += A[vv][v];
                        }
                    e = B[vv][v] - e;
                    }
                while(v != vertex || e != edge);
                *r = EOS;
                return(NO_ERROR);
                }
            }
        }
    
    /******************************************************************************************
        If execution gets to this point, then the "non-primitive" relator appears in the dual
        relators only with exponent 1. We next look for an occurance of the "primitive" relator
        in the dual relators with exponent greater than 1. If we find such an occurance, then
        we can find the desired transverse curve by next locating an occurance of the "non-
        primitive relator flanked on either side by appearances of the "primitive" relator.
    ******************************************************************************************/
    
    for(i = 1; i <= 2; i++)
        {
        for(p = *DualRelators[i],q = p + 1; *p; p++,q++)
            {
            if(!*q) q = *DualRelators[i];
            if((*p == z || *p == w) && *p == *q)
                {
                edge = p - *DualRelators[i];
                if(i == 1)
                    vertex = 0;
                else
                    vertex = 2;
                for(j = 1; j <= 2; j++)
                    {
                    if(GetHandleSize((char **) DualRelators[j]) < 4) continue;
                    for(p = *DualRelators[j],r = p + 1, q = p + 2; *p; p++,r++,q++)
                        {
                        if(!*q) q = *DualRelators[j];
                        if(!*r) r = *DualRelators[j];
                        if((*p == z || *p == w) && *p == *q && (*r == x || *r == y))
                            {
                            e = p - *DualRelators[j];
                            edge2 = q - *DualRelators[j];
                            if(j == 1)
                                v = 0;
                            else
                                v = 2;    
                            r = ptr;
                            do
                                {
                                if(v & 1)
                                    {
                                    *r = (v >> 1) + 'a';
                                    vv = v - 1;
                                    }
                                else
                                    {
                                    *r = (v >> 1) + 'A';
                                    vv = v + 1;
                                    }
                                e = OSA[v] - e;
                                if(e >= V[v]) e -= V[v];
                                r++;
                                if(r - ptr > max )
                                    {
                                    r--;
                                    *r = EOS;
                                    return(5001);
                                    }
                                v = FV[vv];
                                d = A[vv][v];
                                while(d <= e)
                                    {
                                    v = CO[vv][v];
                                    d += A[vv][v];
                                    }
                                e = B[vv][v] - e;
                                }
                            while(v != vertex || e != edge);
                            
                            e++;
                            if(e >= V[v]) e -= V[v];
                            edge = edge2;
                            if(j == 1)
                                vertex = 0;
                            else
                                vertex = 2;
                            max *= 2;    
                            do
                                {
                                if(v & 1)
                                    {
                                    *r = (v >> 1) + 'a';
                                    vv = v - 1;
                                    }
                                else
                                    {
                                    *r = (v >> 1) + 'A';
                                    vv = v + 1;
                                    }
                                e = OSA[v] - e;
                                if(e >= V[v]) e -= V[v];
                                r++;
                                if(r - ptr > max )
                                    {
                                    r--;
                                    *r = EOS;
                                    return(5001);
                                    }
                                v = FV[vv];
                                d = A[vv][v];
                                while(d <= e)
                                    {
                                    v = CO[vv][v];
                                    d += A[vv][v];
                                    }
                                e = B[vv][v] - e;
                                }
                            while(v != vertex || e != edge);
                            *r = EOS;
                            return(NO_ERROR);
                            }
                        }
                    }
                }
            }
        } 

    /******************************************************************************************
        If execution gets to this point, then neither relator appears in the dual relators
        with exponent greater than 1. We next look for an occurance of the "non-primitive" 
        relator lying between two oppositely oriented edges of the "primitive" relator.
        If this situation occurs, a "wave" to the "primitive" relator will intersect the
        "non-primitive" relator once. 
        	Note that this situation will not arise when the so-called "primitive" relator 
        is a "positive" relator. Including the following code allows this subroutine to be
        used to locate transverse curves in more general situations.
    ******************************************************************************************/
       
    for(i = 1; i <= 2; i++)
		{
		if(GetHandleSize((char **) DualRelators[i]) < 4) continue;
		for(p = *DualRelators[i],r = p + 1, q = p + 2; *p; p++,r++,q++)
			{
			if(!*q) q = *DualRelators[i];
            if(!*r) r = *DualRelators[i];
			if((*r == x || *r == y) && ((*p == z && *q == w) || (*p == w && *q == z)))
				{
				e = p - *DualRelators[i];
                edge = q - *DualRelators[i];
                if(i == 1)
                    {
					v = 0;
					vertex = 1;
					}
                else
                    {
					v = 2;
					vertex = 3;
					}				
				edge = OSA[v] - edge;
				if(edge >= V[v]) edge -= V[v];	 
                r = ptr;
				
				do
					{
					if(v & 1)
						{
						*r = (v >> 1) + 'a';
						vv = v - 1;
						}
					else
						{
						*r = (v >> 1) + 'A';
						vv = v + 1;
						}
					e = OSA[v] - e;
					if(e >= V[v]) e -= V[v];
					r++;
					if(r - ptr > max )
						{
						r--;
						*r = EOS;
						return(5001);
						}
					v = FV[vv];
					d = A[vv][v];
					while(d <= e)
						{
						v = CO[vv][v];
						d += A[vv][v];
						}
					e = B[vv][v] - e;
					}
				while(v != vertex || e != edge);
                        
				/****** Freely reduce the string in ptr. ******/
				
				r--;
				for(p = ptr + 1; abs(*p - *r) == 32 && p < r; p++, r--) ;
				r++;
				*r = EOS;
				
				/****** Move the freely reduced string in ptr so it starts at ptr. *****/
				
				q = ptr;
				while( (*q++ = *p++) ) ;
				return(NO_ERROR);
				}
			}
		}
                           
    return(1);
}        
    
unsigned long Recip_Mod_P(p,q)
register unsigned long  p,
                        q;
{
    /******************************************************************************************
        If L(p,q) is a lens space, with 0 <= q < p, 0 <= q' < p and qq' congruent to 1 mod p,
        then L(p,q), L(p,p - q), L(p, q') and L(p, p - q') are homeomorphic. This routine is
        used to find that element of the set {q, p - q, q', p - q'} which is minimal. This 
        means that a user can directly compare the parameters of two lens spaces returned by
        Heegaard and determine whether the lens spaces are homeomorphic.
    ******************************************************************************************/
        
    unsigned long     min;
    unsigned long LGCD();
    
    if(p <= 1) return(0);
    if(q == 0) return(0);             	/** Error!!! **/
    if(LGCD(p,q) != 1) return(0);        	/** Error!!! **/
    q = q % p;
    min = q;
    if(p - q < min) min = p - q;
    if(Recip_Q < 0)
        {
        Recip_Q = -Recip_Q;
        Recip_Q = Recip_Q % p;
        Recip_Q = p - Recip_Q;
        }
    else
        Recip_Q = Recip_Q % p;        
    if(Recip_Q < min) min = Recip_Q;
    if(p - Recip_Q < min) min = p - Recip_Q;    
    return(min);            
}

int Find_Flow_B(unsigned int Source)
{
    /******************************************************************************************
        Let G be the generator which corresponds to the "Source". This variant of Find_Flow_A()
        looks for T-transformations which reduce the number of appearances of G and its
        inverse g, in the relators. Find_Flow_B() is called by the routine which looks for
        relators which are proper powers of free generators.
    ******************************************************************************************/
            
    register unsigned int  	i,
                            j,
                            max,
                            min,
                            *p,
                            *q,
                            *r;
                            
    unsigned int           	Flow,
                            k,
                            MaxFlow,
                            S[VERTICES],
                            Sink;
    
    Fill_A(NumRelators);        
    if(ComputeValences_A()) return(TOO_LONG);                    
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
    Automorphisms = 0L;            
    while(1)
        {
        Sink = Source + 1;
        MaxFlow = VA[Source/2];
        Flow = A[Source][Sink];
        A[Source][Sink] = 0;
        A[Sink][Source] *= 2;                
        while(MaxFlow)
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
            if(Flow >= MaxFlow) break;
            if(max == 0) break;
            i = Source;
            while(i != Sink)
                {
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
        **************************************************************************************/
        
        if(MaxFlow > Flow)    
            {    
            Automorphisms ++;
            if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(TOO_LONG);                                                
            Fill_A(NumRelators);        
            Length -= MaxFlow - Flow;
            VA[Source/2] = Flow;
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
            return(NO_ERROR);
            }                                                        
        }
}

int Get_Connected_Components(void)
{    
    /******************************************************************************************
        This routine finds the components of the graph specified by the adjacency lists AJ1[].
        The array ZZ[] is initialized by the calling routine which sets the entries of vertices
        which should be deleted from the adjacency lists to the value INFINITE.
            The routine returns an integer equal to the number of components of the graph
        determined by: the original adjacency lists less the vertices deleted from these lists.
    ******************************************************************************************/    
    
    int                    	CompNum;
     
    register unsigned int  	h,
                            i,
                            j,
                            *p,
                            *r;
    CompNum = 0;
    while(1)
        {
        for(i = 0; ZZ[i] && (i < Vertices); i++) ;
        if(i == Vertices) return(CompNum);
        CompNum ++;                        
        ZZ[i] = CompNum;
        for(r = UpDate,*r = i,p = r + 1; r < p; r++)
            {
            i = *r;
            for(h = 0; (j = AJ1[i][h]) < VERTICES; h++)
                {
                if(ZZ[j] == 0)
                    {
                    ZZ[j] = CompNum;
                    *p++ = j;
                    }
                }
            }            
        }        
}

unsigned long 	LGCD(p,q)
unsigned long  	p,
                q;
{
    /******************************************************************************************
        Given a pair of nonnegative long integers p and q, this routine computes the GCD
        of p and q and returns LGCD(p,q) as an unsigned long integer. The routine
        simultaneously computes a pair of long integers Recip_P and Recip_Q such that the
        following equation holds.
                        LGCD(p,q) = p*Recip_P + q*Recip_Q.
    ******************************************************************************************/    
    
    register long  	d,
                    t,
                    u2,
                    v1,
                    v2;
                    
    long           	u1;
        
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
    Recip_Q = u1;
    return(u2);
}

unsigned int Split_At_Empty_Relators(F1)
int        F1;
{
    register unsigned char  *p,
                            *q,
                            *r;

    unsigned char           BCL1[MAXNUMRELATORS + 2],
                            BCL2[MAXNUMRELATORS + 2],
                            CGBC[MAXNUMRELATORS + 2],
                            ERL1[MAXNUMRELATORS + 1],
                            ERL2[MAXNUMRELATORS + 1],
                            **Temp;
                    
    int                    	i,
                            j,
                            k,
                            N1,
                            N2,
                            NumEmtyGens,
                            NumEmtyRel,
                            NumNewPres,
                            OrigNumBdryComps,
                            SaveCS,
                            SaveNumGens,
                            SaveNumRelators;
        
    unsigned int            Error,
                            SaveUDV;

    unsigned long          	HS;
    
    Error = FALSE;
    
    /***************************************************************************************
        Swap Relators[] and Copy_Of_Rel_1[]. We are going to try and find a Heegaard
        diagram corresponding to the set of relators we saved in Copy_Of_Rel_1[].
    ***************************************************************************************/        

    if(!F1)
        {
        NumRelators ++;
        NumGenerators ++;
        Vertices += 2;
        }
    
    for(i = 1; i <= NumRelators; i++)
        { 
        Temp               	= Relators[i];
        Relators[i]         = Copy_Of_Rel_1[i];
        Copy_Of_Rel_1[i]    = Temp;    
        }
    
    for(i = 1,Length = 0L; i <= NumRelators; i++) Length += GetHandleSize((char **) Relators[i]);
    Length -= NumRelators;

    if(Micro_Print)
        {
        printf("\n\nSwapped the presentation with a copy of the previous presentation to get:\n");
        Print_Relators(Relators,NumRelators);
        }

    if( (Error = Find_Flow_A(NORMAL,FALSE)) ) return(Error);
    
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
        
    DrawingDiagrams = TRUE;
    TestRealizability3 = TRUE;
    WhichInput = ReadPres;
    Error = Whitehead_Graph();
    if(!Connected) return(NOT_CONNECTED);    
    switch(Error)
        {
        case NO_ERROR:
            break;
        case NON_UNIQUE_1:
        case NON_UNIQUE_2:
        case NON_UNIQUE_3:
        case NON_UNIQUE_4:
            Error = FALSE;
            break;                /* Ignore these  non-unique "errors" ? */
        case V2_ANNULUS_EXISTS:    
        case REDUCE_GENUS:
        case NON_PLANAR:
        case TOO_MANY_COMPONENTS:
        case FATAL_ERROR:        
        case TOO_LONG:        
            return(Error);            
        case SEP_PAIRS:
            Num_Saved_LPres = 0;
            NotNewPres = 0; 
            Error = Level_Transformations(FALSE,TRUE,TRUE);
            for(i = 0; i < Num_Saved_LPres; i++)
            for(j = 0; j <= NumRelators; j++) if(SLR[i][j] != NULL)
            	{
            	DisposeHandle((char **) SLR[i][j]);
            	SLR[i][j] = NULL;
            	}
            if(Micro_Print)
                {
                printf("\n\nThe current Presentation is:\n");
                Print_Relators(Relators,NumRelators);
                }                        
            if(Error != 2)
                return(TRUE);
            Fill_A(NumRelators);
            DrawingDiagrams = FALSE;
            if(Whitehead_Graph()) return(TRUE);
            DrawingDiagrams = TRUE;        
            break;        
        }

    /***************************************************************************************
        Heegaard has found a diagram corresponding to the presentation in Relators[].
        If TestRealizability1 is TRUE and TestRealizability2 is FALSE, we can return
        because Heegaard is only looking for an initial diagram which is realizable,
        and we have found one.
    ***************************************************************************************/
    
    if(TestRealizability1 && !TestRealizability2) return(NO_ERROR);    

    /***************************************************************************************
        Determine which relators are consequences of the primitive relator, or proper power,
        Relators[NumRelators]. And, if this routine is not being called when a proper power
        is present, then flag the relator Relators[NumRelators] as TRUE.
    ***************************************************************************************/
                            
    for(i = 1; i <= NumRelators; i++) ERL1[i] = EOS;
    NumEmtyRel = 0;
    for(i = 1; i < NumRelators; i++) if(GetHandleSize((char **) Copy_Of_Rel_1[i]) == 1)
        {
        ERL1[i] = TRUE;
        NumEmtyRel ++;
        }
    if(F1 == FALSE)    ERL1[NumRelators] = TRUE;    

    /***************************************************************************************
            Determine which relators in OutRelators[] correspond to which relators in
                                        Relators[].
    ***************************************************************************************/
        
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
            return(TRUE);
            }
        if(j < 0) j = -j;
        p = *Relators[j];
        *p += 125;
        ERL2[i] = j;
        }
    for(i = 1; i <= NumRelators; i++)
        {
        p = *Relators[i];
        if(*p > 124) *p -= 125;
        }        

    Split_At_Empty_Relators_Sub1();
    
    for(i = 1; i <= NumBdryComps; i++) CGBC[i] = GBC[i];
    
    OrigNumBdryComps = NumBdryComps;
    
    for(i = 1; i <= NumBdryComps; i++) BCL1[i] = TRUE;
    
    /*****************************************************************************************
        Compute the genus of each handlebody that results when we cut open along the empty
        relators.
    *****************************************************************************************/
                
    for(i = 0; i < NumRelators; i++)
        {
        if(ERL1[ERL2[i+1]]) continue;
        j = 2*i;
        N1 = zz[j];
        N2 = zz[j+1];
        if(N1 == N2)
            {
            GBC[N1] ++;
            NRBC[N1] -= 2;
            continue;
            }
        if(N2 < N1)
            {
            j = N1;
            N1 = N2;
            N2 = j;
            }
        for(j = 0; j < 2*NumRelators; j++) if(zz[j] == N2) zz[j] = N1;
        GBC[N1]  += GBC[N2];
        NRBC[N1] += NRBC[N2] - 2;
        NFBC[N1] += NFBC[N2];
        NEBC[N1] += NEBC[N2];
        BCL1[N2]  = FALSE;
        BCL2[N2]  = N1;
        NumBdryComps --;
        }        
    
    /*****************************************************************************************
        Check whether the presentation that is "splitting" is already on file, if it is, we
        don't want to save another copy of it.
    *****************************************************************************************/    
    
    switch(Splitting_Pres_On_File(12,NumBdryComps + 1))
        {
        case TOO_LONG:
            return(TOO_LONG);
        case TOO_MANY_COMPONENTS:
            return(TOO_MANY_COMPONENTS);
        case TRUE:
            return(TRUE);    
        case NO_ERROR:
            break;
        }
    
    SaveUDV = UDV[ReadPres];
    switch(SaveUDV)
        {
        case SPLIT:
        case ANNULUS_EXISTS:
        case V2_ANNULUS_EXISTS:
            return(TRUE);
        }    
    UDV[ReadPres]          	= SPLIT;                
    Daughters[ReadPres]     = NumFilled;
    NCS[ReadPres]          	= 0;    
    SaveCS                 	= CS[CurrentComp];
    if(!CS[CurrentComp])    CS[CurrentComp] = TRUE;
    
    SaveNumRelators 	= NumRelators;
    SaveNumGens     	= NumGenerators;
    NumNewPres         	= 0;
    NumEmtyGens     	= 0;

	/***************************************************************************************** 
	The line below sets UDV[] == DONE for each presentation in the component that just split
	so it won't be run again. Comment out the line below here and in Heegaard3.c and 
	Heegaard8.c to let Heegaard rerun presentations that have split.
	******************************************************************************************/ 
	
	for(i = 0; i < NumFilled; i++) 
		if(ComponentNum[i] == ComponentNum[ReadPres] && UDV[i] < DONE) UDV[i] = DONE;
    
    for(i = 0; i <= SaveNumGens; i++) BCWG[i] = EOS;
    
    /*****************************************************************************************
        Partition the nonempty relators of the original presentation into disjoint sets
        corresponding to each summand of the splitting induced by the empty relators.
    *****************************************************************************************/    
        
    for(i = 1; i <= OrigNumBdryComps; i++)
        {
        if(!BCL1[i]) continue;
        NumRelators = 0;
        
        for(j = 0; j < SaveNumRelators; j++) if(zz[2*j] == i)
            {
            k = ERL2[j+1];
            if(ERL1[k]) continue;
            NumRelators ++;
            HS = GetHandleSize((char **) Copy_Of_Rel_1[k]);
            if(Relators[NumRelators] != NULL) DisposeHandle((char **) Relators[NumRelators]);
            Relators[NumRelators] = (unsigned char **) NewHandle(HS);
            if(Relators[NumRelators] == NULL) Mem_Error();
            p = *Relators[NumRelators];
            q = *Copy_Of_Rel_1[k];
            while( (*p++ = *q++) ) ;
            }
            
        if(NumRelators == 0)
            {
            if(GBC[i])
                {
                NumEmtyGens += GBC[i];
                BCWG[GBC[i]] ++;
                }
            continue;
            }
            
        if(Split_At_Empty_Relators_Sub2(GBC[i],NumNewPres) == TOO_LONG)
            {
            if(NumNewPres)
                {
                for(i = NumFilled - NumNewPres; i < NumFilled; i++)
                    {
                    BytesUsed -= SURL[i];
                    UDV[i] = PRIM[i] = 0;
                    CBC[TotalComp][0] = BDRY_UNKNOWN;
                    MLC[TotalComp][NG[i]] = BIG_NUMBER;
                    OnStack -= 2*NG[i];
                    TotalComp --;
                    }
                NumFilled -= NumNewPres;
                }
            UDV[ReadPres] = SaveUDV;
            CS[ComponentNum[ReadPres]] = SaveCS;                    
            return(TOO_LONG);
            }
        
        /*************************************************************************************
            Set up the information about the boundary components of each new summand.
        *************************************************************************************/
                
        for(j = 0; j <= SaveNumGens; j++) CBC[TotalComp][j] = EOS;
        for(j = 1; j <= OrigNumBdryComps; j++)
            {
            k = j;
            while(BCL1[k] == FALSE) k = BCL2[k];
            if(k == i) CBC[TotalComp][CGBC[j]] ++;
            }
            
        for(j = SaveNumGens + 1; j > 0; j--) if(CBC[TotalComp][j-1])
            {
            CBC[TotalComp][j] = BDRY_UNKNOWN;
            break;
            }
        if(j == 0) CBC[TotalComp][0] = BDRY_UNKNOWN;    
                    
        NumNewPres ++;
        NCS[ReadPres] ++;
        }
    
    /*****************************************************************************************
        If there are any empty summands or any S1 X S2s that exist, save a presentation
                corresponding to the union of these empty summands and or S1 X S2s.
    *****************************************************************************************/
    
    NumEmtyRel -= NumBdryComps - 1;
    BCWG[0] += NumEmtyRel;
    NumEmtyGens += NumEmtyRel;
    for(i = SaveNumGens + 1; i > 0; i--) if(BCWG[i-1])
        {
        BCWG[i] = BDRY_UNKNOWN;
        break;
        }
    if(i == 0) BCWG[0] = BDRY_UNKNOWN;
    
    if(BCWG[0] < BDRY_UNKNOWN)
        {
        if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);        
        TotalComp                  	++;
        NCS[ReadPres]              	++;
        BDY[NumFilled]            	= 3;
        ComponentNum[NumFilled]     = TotalComp;            
        ER[NumFilled]              	= 0;
        FR[NumFilled]             	= ReadPres;
        MLC[TotalComp][NumEmtyGens] = 0L;    
        NG[NumFilled]              	= NumEmtyGens;
        NR[NumFilled]              	= NumEmtyRel;
        PRIM[NumFilled]             = 12;
        SURL[NumFilled]             = 0L;
        TP[NumFilled]              	= 0;
        for(i = 1; i <= NumEmtyRel; i++)
            {
            HS = sizeof(char);
            if(SUR[NumFilled][i] != NULL) DisposeHandle((char **) SUR[NumFilled][i]);
            SUR[NumFilled][i] = (unsigned char **) NewHandle(HS); 
            if(SUR[NumFilled][i] == NULL) Mem_Error();         
            q = *SUR[NumFilled][i];
            r = q;
            *q++ = EOS; 
        	if((q-r) != HS) 
        		{
        		NumErrors ++;
        		printf("\n\n4) Error in Presentation %u! |Relator[%d]| = 0, HS = %lu.",NumFilled + 1,i,HS);
        		}        	
            }
        NumNewPres ++;
        UDV[NumFilled]     	= MISSING_GEN_DONE1;
        CS[TotalComp]      	= 2;
        N1H[TotalComp]     	= 0;
        NS1XS2[TotalComp]  	= NumEmtyRel;
        NS1XD2[TotalComp]  	= NumEmtyGens - NumEmtyRel;
        for(j = 0; (CBC[TotalComp][j] = BCWG[j]) < BDRY_UNKNOWN; j++) ;
        
        if(Micro_Print)
            {
            printf("\nThe presentation of summand %d is:\n",NumNewPres);
            Print_Relators(Relators,NumRelators);
            printf("\n\nSaved this presentation as: Presentation %u\n",NumFilled + 1);
            }
                        
        NumFilled ++;
        SaveMinima = TRUE;
        printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
        printf("Gen%3d  Rel%3d  Length%6lu  From%6d  ER",
            NumEmtyGens,NumEmtyRel,0L,ReadPres + 1);        
        }
        
    return(Error);        
}

void Split_At_Empty_Relators_Sub1(void)
{
    unsigned char  	*p,
                    x;
    
    int           	Genus;
    
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
        GBC[i] = Genus;
        }    
}

int Split_At_Empty_Relators_Sub2(int NumGen,int NumNewPres)
{
    register unsigned char 	*p,
                            *q,
                            *r;
                            
    int                    	h,
                            i,
                            j;
                            
    unsigned long			HS;
    
    /*****************************************************************************************
        Call Rewrite_Input() to rewrite the relators corresponding to the current summand.
        Save the rewritten relators in SUR[]. Increase TotalComp by 1 and update Daughters[].
    *****************************************************************************************/
    
    if(NumFilled >= MAX_SAVED_PRES - 3) return(TOO_LONG);
    
    if(Micro_Print)
        {
        printf("\n\nThe presentation of summand %d is currently:\n",NumNewPres + 1);
        Print_Relators(Relators,NumRelators);
        }
                
    Rewrite_Input();
    for(i = 1, Length = 0L; i <= NumRelators; i++)
        Length += GetHandleSize((char **) Relators[i]);
    Length -= NumRelators;
    
    if(Length && Find_Flow_A(NORMAL,FALSE) == TOO_LONG) return(TOO_LONG);

    Canonical_Rewrite(Relators,FALSE,FALSE);
    
    /*****************************************************************************************
        Reset the number of generators of this presentation so that it equals the genus of
        the handlebody summand to which the relators are attached. This number may be
        greater than the number of generators which appear explicitly in the relators.
    *****************************************************************************************/
            
    NumGenerators                  	= NumGen;
            
    TotalComp                      	++;
    ComponentNum[NumFilled]         = TotalComp;    
    ER[NumFilled]                  	= -4;
    FR[NumFilled]                 	= ReadPres;        
    MLC[TotalComp][NumGenerators]  	= Length;
    NG[NumFilled]                  	= NumGenerators;
    NR[NumFilled]                  	= NumRelators;
    PRIM[NumFilled]                	= 12;
    SURL[NumFilled]                 = Length;
    UDV[NumFilled]                 	= 0;
    TP[NumFilled]                 	= NumRelators;
    BDY[NumFilled]                 	= 3;
    OnStack                        	+= 2*NumGenerators;
    
    for(i = 1; i <= NumRelators; i++)
        {
        HS = GetHandleSize((char **) Relators[i]);
        if(SUR[NumFilled][i] != NULL) DisposeHandle((char **) SUR[NumFilled][i]);
        SUR[NumFilled][i] = (unsigned char **) NewHandle(HS);            
        if(SUR[NumFilled][i] == NULL) Mem_Error();
        p = *Relators[i];
        q = *SUR[NumFilled][i];    
        r = q;
        while( (*q++ = *p++) ) ; 
        if((q-r) != HS) 
        	{
        	NumErrors ++;
        	printf("\n\n5) Error in Presentation %u! |Relator[%d]| = %lu, HS = %lu.",NumFilled + 1,i,q-r-1,HS);
        	}
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
    
    if(Micro_Print)
        {
        printf("\nThe presentation of summand %d is:\n",NumNewPres + 1);
        Print_Relators(Relators,NumRelators);
        printf("\nSaved this presentation as: Presentation %u\n",NumFilled + 1);
        }    
                                                             
    NumFilled ++;
    
    printf("\nPres%6u  ToDo%6u  Summand%3d  ",NumFilled,OnStack,TotalComp);
    printf("Gen%3d  Rel%3d  Length%6lu  From%6u  ER",
        NG[NumFilled-1],NR[NumFilled-1],SURL[NumFilled-1],ReadPres + 1);
    
    return(NO_ERROR);
}

int Splitting_Pres_On_File(int WhoCalled,int NumNewPres)
{
    register int    i,
                    j;

    /**************************************************************************************
        Check whether the presentation that is "splitting" is already on file, if it is,
        we don't want to save another copy of it.
    **************************************************************************************/    
    
    Canonical_Rewrite(Relators,FALSE,FALSE); 
    
/*          
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
                 return(TRUE);    
             break;        
             }    
         }
*/
 /* Tried changing this so it only searches the current component. 2/18/15 */
 
     for(i = 0,Dup_On_File = INFINITE; i < NumFilled; i++)
        if(ComponentNum[i] == CurrentComp && SURL[i] == Length && NG[i] == NumGenerators && NR[i] == NumRelators)
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
                 return(TRUE);    
             break;        
             }    
         }
         
    if(i == NumFilled && Dup_On_File < INFINITE)
        {
         WhoCalled = 12;
         if((NumFilled >= MAX_SAVED_PRES - 3) || 
             Save_Pres(ReadPres,Dup_On_File,Length,1,WhoCalled,1,0,0)) return(TOO_LONG);                
         Mark_As_Duplicate(Dup_On_File);
         return(TRUE);
         }
    if(i < NumFilled && UDV[ReadPres] == SPLIT) return(TRUE);
    
    if(Micro_Print)
        {
        if(WhoCalled == 12)
            {
            printf("\n\nRelator %d of the preceding presentation is a primitive or proper power whose",
                NumRelators);
            printf("\ndeletion creates empty relators in the presentation.");    
            printf("\nThis means the manifold 'splits' as a connected sum.");    
            printf("\nHeegaard will look for presentations corresponding to each summand.\n");
            }
        if(WhoCalled == 40)
            {
            printf("\n\nA relator in the preceding dual presentation freely reduces to an empty relator.");
            printf("\nThis means the manifold 'splits' as a connected sum.");    
            printf("\nHeegaard will look for presentations corresponding to each summand.\n");
            }
        }
        
    if(i == NumFilled)
        {    
        /**********************************************************************************
            The presentation that is "splitting" is new, so we want to save a copy.
        **********************************************************************************/

        WhoCalled = 12;
        if((NumFilled >= MAX_SAVED_PRES - 3) || 
            Save_Pres(ReadPres,0,Length,1,WhoCalled,1,0,0)) return(TOO_LONG);        
        BDY[NumFilled - 1] = BDY[ReadPres];
        UDV[NumFilled - 1] = 0;
        ReadPres = NumFilled - 1;
        }    
        
    /**************************************************************************************
         If Heegaard already has as many summands as it can handle, flag any other
         presentations corresponding to this summand so that we will quit processing them.
     **************************************************************************************/    
         
     if(TotalComp + NumNewPres > MAXNUMCOMPONENTS - 3)
         {
         ReadPres = SReadPres;
         Mark_As_Found_Elsewhere(CurrentComp);
         if(Batch == FALSE) SysBeep(5);
         printf("\n\nStopping because Heegaard cannot deal with any more summands. Sorry!"); 
         if(NumErrors == 1)
			printf("\nOne error was detected. Scroll back for details.");
		 if(NumErrors > 1)
			printf("\n%lu errors were detected. Scroll back for details.",NumErrors); 
         printf("\n\nRerunning using Depth-First Search may help.");
         Too_Many_Components_ALert();   
         return(TOO_MANY_COMPONENTS);
        }
    return(NO_ERROR);    
}
