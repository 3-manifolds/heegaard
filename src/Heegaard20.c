#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L  10 Is_Knot_Relator(void)
L  87 Genus_Two_Meridian_Reps(int OrbitNum, int Flag)
l 310 Genus_Two_Meridian_Reps_Sub(unsigned char* ptr1, unsigned char* ptr2)
********************************************************************************************/

int Is_Knot_Relator(void)
{
	char			x;
	
	unsigned char 	*p,
					*q;
					
	int				NumA,
					NumB,
					Numa, 
					Numb;
	
	unsigned int	AP, 
					AQ;													
	
	if(Batch == 2 && H_Results != NULL) fprintf(H_Results,"\n\n%s",PresName);
	
	NumA = NumB = Numa = Numb = 0;
	p = *Relators[1];
	while( (x = *p++) ) switch(x)
		{
		case 'A':
			NumA ++;
			break;
		case 'a':
			Numa ++;
			break;
		case 'B':
			NumB ++;
			break;
		case 'b':
			Numb ++;
			break;
		default:
			printf("\nThe relator R must contain only the characters 'A', 'B', 'a' and 'b'!");
			return(1);
		}
	
	if(NumA >= Numa)
		AP = NumA - Numa;
	else
		AP = Numa - NumA;
	if(NumB >= Numb)
		AQ = NumB - Numb;
	else
		AQ = Numb - NumB;
	if(AP == 0 && AQ == 0)
		printf("\n\nThe input relator has trivial abelianization.");
	else
		{
		if(1 != AlexanderPolynomial_GCD(AP,AQ))
		printf("\n\nThe input relator is not a homology generator, hence not a knot relator.");
		}	

	
	if(Genus_Two_Meridian_Reps(1,1))
		{
		if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);
		Relators[1] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[1]));
		if(Relators[1] == NULL) Mem_Error();
		p = *Copy_Of_Input[1];
        q = *Relators[1];
        while( (*q++ = *p++) ) ;
        NumGenerators = 2;
		NumRelators = 1;
		Vertices = 4;
		printf("\n\nRestored the original initial presentation.");
		printf("\nSo R1 is: %s",*Relators[1]);
		return(1);
		}
	
	NumGenerators = 2;
	NumRelators = 2;
	Vertices = 4;	
	return(0);		
}

int Genus_Two_Meridian_Reps(int OrbitNum, int Flag)
{
    /******************************************************************************************
        Genus_Two_Meridian_Reps() can be used to find the canonical representatives of the
        'meridian' of a manifold M = H[R] obtained by adding a 2-handle to a genus two 
        handlebody H along a nonseparating simple closed curve R. 
  			The routine won't work when H[R] is not uniquely determined by R. For instance, 
  		when one of the generators of H appears in R only with exponents having only one 
  		absolute value E with E greater than two. Thus, for example, M usually can't be the 
  		exterior of a torus or cable knot in S^3.	
    ******************************************************************************************/
                                                
    unsigned char 	*ptr1 = NULL,
                    *ptr2 = NULL,
                    *p,
                    *q,
                    x;
    
    int           	i,
                    j,
                    LTRV;
                                           
    unsigned int  	r;
    
    unsigned long   HS,
    				SLength;

    unsigned int   	Whitehead_Graph();
    
    if(NumGenerators != 2) return(11);
    if(NumRelators != 1)  return(12);
    
    /*******************************************************************************************
        Call Find_Flow_A() to reduce Relators[1] to minimal length before we try to find the 
        Heegaard diagram. Then call Whitehead_Graph() to actually find the diagram. If 
        Whitehead_Graph() doesn't return an error, then call Sep_Surface() to see if 
        Relators[1] is realized by a non-separating simple closed curve.
     ******************************************************************************************/
    
    SLength = Length;
    if(Flag) 
    	{
    	printf("\n");
    	Micro_Print = TRUE;
    	}                         
    if(Find_Flow_A(NORMAL,FALSE))
    	{
    	if(Flag) Micro_Print = FALSE;
    	return(1);
    	}
    if(Flag) Micro_Print = FALSE;
    if(Length < SLength)
    	{
    	printf("\n R 1 did not have minimal length. The above automorphism(s) reduced R 1 to:");
    	printf("\n %s",(char *) *Relators[1]);
    	}
    	    
_GET_DIAGRAM:
    
    Saved_Vertices = 0;
    TestRealizability1 = TRUE;
    r = Whitehead_Graph();
    switch(r)
        {
        case NO_ERROR:
            break;
        case NON_PLANAR:
            printf("\n\n                    The Whitehead graph is nonplanar.");
            TestRealizability1 = FALSE;
            return(r);
        case FATAL_ERROR:
            Fatal_Error(); 
            TestRealizability1 = FALSE;       
            return(r);
        case TOO_LONG:
            printf("\n\n                    This presentation may be too long!");
            TestRealizability1 = FALSE;
            return(r);
        case TOO_MANY_COMPONENTS:
        	TestRealizability1 = FALSE;
            return(r);        
        case NON_UNIQUE_1:
        case NON_UNIQUE_2:
        case NON_UNIQUE_3:
        case NON_UNIQUE_4:
         	printf("\n\nCan't find meridian reps unless H[R] is uniquely determined by R. Here R is:");
        	printf("\n%s",(char *) *Relators[1]);
        	TestRealizability1 = FALSE;
            return(r);
        case V2_ANNULUS_EXISTS:    
			printf("\n\nA V2_Annulus exists.");
			TestRealizability1 = FALSE;
			return(r);    
        case NOT_CONNECTED:
            printf("\n\n                    The Whitehead graph is not connected.");
            printf("\n\n                    Please check each summand separately.");
            TestRealizability1 = FALSE;
            return(r);
        case SEP_PAIRS:
            Num_Saved_LPres = 0;
            NotNewPres = 0;
            LTRV = Level_Transformations(FALSE,TRUE,TRUE);
 	        for(i = 0; i < Num_Saved_LPres; i++)
            for(j = 1; j <= NumRelators; j++) if(SLR[i][j] != NULL)
            	{
            	DisposeHandle((char **) SLR[i][j]);
            	SLR[i][j] = NULL;
            	}
            switch(LTRV)
                {
                case 0:
                case 1:
                	printf("\n\nCan't find meridian reps unless H[R] is uniquely determined by R. Here R is:");	
                	printf("\n%s",(char *) *Relators[1]);
                	TestRealizability1 = FALSE;	
                    return(SEP_PAIRS);
                case 2: 
                	{
                	if(Micro_Print)
                		printf("\n\nPerformed some level-transformations on the input presentation.");	   
                	goto _GET_DIAGRAM;
                	}
                default:
                	printf("\n\nCan't find meridian reps unless H[R] is uniquely determined by R. Here R is:");       	
                	printf("\n%s",(char *) *Relators[1]);
                	TestRealizability1 = FALSE;				
                    return(SEP_PAIRS);
                }    
        default:
        	printf("\n\nCan't find meridian reps unless H[R] is uniquely determined by R. Here R is:");
        	printf("\n%s",(char *) *Relators[1]);
        	TestRealizability1 = FALSE;
            return(r);
        }
    
    TestRealizability1 = FALSE;
    
    if(LTRV == 2)
    	{
    	printf("\n\n Note: The Whitehead graph of R 1 has a separating pair of vertices.");
    	printf("\n Performed some level-transformations on R 1 to get R 1':");
    	printf("\n R 1' = %s",(char *) *Relators[1]);
    	}
                                                    
    if(Sep_Surface() > 1) 		/***** H[R] has two boundary components. *****/
    	{
    	if(LTRV == 2)
    		printf("\n\nR 1' is a separating curve. Can't find M1 and M2 when R 1' is separating.");
    	else	
			printf("\n\nR 1 is a separating curve. Can't find M1 and M2 when R 1 is separating.");
    	return(67);    
    	}
                                                                   
    /******************************************************************************************
            	We have found the Heegaard diagram, and H[R] has torus boundary. 
    ******************************************************************************************/
    
	if(Micro_Print)
		{
		printf("\n\nThe current Presentation is:\n");
		Print_Relators(Relators,NumRelators);
		}       

    HS = GetHandleSize((char **) Relators[1]);
    ptr1 = (unsigned char *) NewPtr(HS);
    if(ptr1 == NULL) Mem_Error();
    ptr2 = (unsigned char *) NewPtr(HS);
    if(ptr2 == NULL) Mem_Error();
    
    if(Genus_Two_Meridian_Reps_Sub(ptr1,ptr2))
        {                            
        DisposePtr((char *) ptr1);
   		DisposePtr((char *) ptr2);
        return(1);
        }
    
    if(Flag)
    	{
    	if(LTRV == 2) 
    		{
    		if(Batch == FALSE) printf("\n\n");
    		printf("\n\nCanonical 'meridian' reps M1 and M2 of R 1'.");
    		}
    	else
    		{
    		if(Batch == FALSE) printf("\n\n");
    		printf("Canonical 'meridian' reps M1 and M2 of R 1.");
    		}
    	}
    else	
    	{
    	printf("\n\nCanonical 'meridian' reps M1 and M2 of the Rep Pres R of Orbit %d:", OrbitNum);
    	printf("\nR  = %s L = %lu",(char *) *Relators[1], Length);
    	}
    printf("\nM1 = %s",ptr1);
    if(Batch == 2 && H_Results != NULL) fprintf(H_Results,"\n%s",ptr1);
    printf("\nM2 = %s",ptr2);
    if(Batch == 2 && H_Results != NULL) fprintf(H_Results,"\n%s",ptr2);
    if(Batch) printf("\n\n");
    if(Flag)
    	{
	for(p = ptr1, r = 1; (x = *p); p++) r++;
    	if(Relators[1] != NULL) DisposeHandle((char **) Relators[1]);   	
    	Relators[1] = (unsigned char **) NewHandle(r);
    	if(Relators[1] == NULL) Mem_Error();
    	
    	for(p = ptr2, r = 1; (x = *p); p++) r++;
    	if(Relators[2] != NULL) DisposeHandle((char **) Relators[2]);
    	Relators[2] = (unsigned char **) NewHandle(r); 
    	if(Relators[2] == NULL) Mem_Error();
    	
    	p = ptr1;
    	q = *Relators[1]; 
    	while( (*q++ = *p++) ) ;
    	p = ptr2;
    	q = *Relators[2];
    	while( (*q++ = *p++) ) ; 	
    	}
    if(ptr1) DisposePtr((char *) ptr1);
    if(ptr2) DisposePtr((char *) ptr2);
    return(NO_ERROR);                            
}

int Genus_Two_Meridian_Reps_Sub(unsigned char* ptr1, unsigned char* ptr2)
{
    /******************************************************************************************
        This routine is called by Genus_Two_Meridian_Reps_Sub(). The routine tries to find a 
        simple closed curve on the genus-two Heegaard surface that is disjoint from a wave 
        based at Relators[1]. If successful, it returns the word representing the disjoint 
        curve in the string pointed to by ptr.
    ******************************************************************************************/
    
    register unsigned char  *p,
                            *q,
                            *r;
                            
    register unsigned int   d,
                            e,
                            v,
                            vv;
    
    int                     i;
    
    unsigned int            edge,
                            initial_e,
                            initial_v,
                            vertex;
    
    unsigned long           max;
    
    max = GetPtrSize((char *) ptr1) - 1;
    

    /******************************************************************************************
 				We first check if Relators[1] is non-positive, and handle that case.
    ******************************************************************************************/
       
    for(i = 1; i <= 2; i++)
      	{
	  	if(GetHandleSize((char **) DualRelators[i]) < 3) continue;
	  	for(p = *DualRelators[i], q = p + 1; *p; p++,q++)
			{
			if(!*q) q = *DualRelators[i];
			if(*p != *q)
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
				initial_e = e;
				initial_v = v;	 
				r = ptr1;
		
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
				  if(r - ptr1 > max )
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
				
				/****** Freely reduce the string in ptr1. ******/
		
				r--;
				for(p = ptr1 + 1; abs(*p - *r) == 32 && p < r; p++, r--) ;
				r++;
				*r = EOS;
		
		/****** Move the freely reduced string in ptr1 so it starts at ptr1. *****/
				
		q = ptr1;
		while( (*q++ = *p++) ) ;
		
		/***** Set things up to find the second representative of the meridian. *****/
		
		e = edge;
		v = vertex;
		edge = initial_e;
		vertex = initial_v;				
		r = ptr2;

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
			if(r - ptr2 > max)
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

		/****** Freely reduce the string in ptr2. ******/
		
		r--;
		for(p = ptr2 + 1; abs(*p - *r) == 32 && p < r; p++, r--) ;
		r++;
		*r = EOS;
		
		/****** Move the freely reduced string in ptr2 so it starts at ptr2. *****/
		
		q = ptr2;
		while( (*q++ = *p++) ) ;
		
		return(NO_ERROR);
	    }
	  }
   }
		
    /******************************************************************************************  
		If execution gets to this point, Relators[1] is positive, and we deal with that case.
    ******************************************************************************************/
		
    e = 0;
    v = 0;
    if(FV[2] == 0)
      	{
		vertex = 3;
		edge = A[3][1];
      	}
    if(FV[2] == 1)
      	{
		vertex = 2;
		edge = A[2][1];
      	}		
		r = ptr1;
		
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
		    if(r - ptr1 > max )
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
		
		/***** Set things up to find the second representative of the meridian. *****/
		
		edge = 0;
		vertex = 0;
		if(FV[2] == 0)
		  	{
		    v = 3;
		    e = A[3][1];
		  	}
		if(FV[2] == 1)
		  	{
		    v = 2;
		    e = A[2][1];
		  	}
		r = ptr2;	
		
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
		    if(r - ptr2 > max )
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
		while(v != vertex || e != edge) ;	
		*r = EOS;
		
	return(NO_ERROR);	
}
