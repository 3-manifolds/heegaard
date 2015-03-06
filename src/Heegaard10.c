#include "Heegaard.h"
#include "Heegaard_Dec.h"
#include <ctype.h>

/****************************** function prototypes *****************************************
L   14 Wirtinger(void)
L  669 Try_Exponent_Surgery(void)
L 1051 Try_Cutting_Disk_Surgery(void)
L 1311 Cutting_Disk_Surgery_Diagram_7(void)
L 1432 Stabilize(void)
L 1667 Heegaard_Splash_Screen(void)
********************************************************************************************/

int Wirtinger()
{
	/*****************************************************************************************
		Wirtinger() checks the presentation in Relators[] to see if it has the form of a 
		Wirtinger presentation of a knot or link. If Wirtinger() decides the presentation
		in Relators[] is a Wirtinger presentation of a knot or link, Wirtinger() offers the
		user the opportunity to perform Dehn surgery on each component of the knot or link.
			(Note that Wirtinger() expects a Wirtinger presentation obtained from an 
		n-crossing projection of a knot or link to have n relators, even though any one of 
		the n relators is superfluous.)
	*****************************************************************************************/
	
	unsigned char	GL[MAXNUMGENERATORS][4],
					*p,
					*q,
					*r,
					*ptr = NULL,
					**Temp,
					*T166 = NULL,
					w,
					x,
					y,
					z;
	
	int				i,
					j,
					L,
					M,
					Mark1,
					Mark2,
					NumComponents,
					Rnum,
					Rnum1,
					Rnum2,
					SavedPresentation,
					TheCase;
							
	unsigned int	h;
	
	long			NRL;						
	
	unsigned int	GCD();
	
	/******************************************************************************************
								See if each relator has length 4.
	******************************************************************************************/	
	
	for(i = 1; i <= NumRelators; i++) if(GetHandleSize((char **) Relators[i]) != 5)
		return(1);

	/******************************************************************************************
		Save a copy of the initial relators so that the user can perform additional surgeries
		if the presentation is a Wirtinger presentation of a knot of link.
	******************************************************************************************/
	
	if(Knot_Or_Link == FALSE) for(i = 1; i <= NumRelators; i++)
		{
		if(KorLRelators[i] != NULL) DisposeHandle((char **) KorLRelators[i]);
		KorLRelators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(KorLRelators[i] == NULL) Mem_Error();
		p = *KorLRelators[i];	
		q = *Relators[i];
		while( (*p++ = *q++) ) ;
		}
	
	NumKnot_Or_Link_Rel = NumRelators;	
	Knot_Or_Link = FALSE;	
	
	/******************************************************************************************
		Each relator has length 4. Next check if each relator has the form XYxW or XYZy.
	******************************************************************************************/	
	
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p++;
		w = *p;
		if(abs(x-z) != 32 && abs(y-w) != 32) return(2);
		if(abs(y-w) != 32)
			{
			p = *Relators[i];
			*p++ = y;
			*p++ = z;
			*p++ = w;
			*p = x;
			}
		}
	
	/******************************************************************************************
		If necessary, replace a relator with a cyclic conjugate of itself so that the first
		letter in each relator is uppercase and the third letter of each relator is lowercase.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p++;
		w = *p;
		if((x < 'a' && z < 'a') || (x > 'Z' && z > 'Z')) return(3);
		if(x > 'Z')
			{
			p = *Relators[i];
			*p++ = z;
			*p++ = w;
			*p++ = x;
			*p = y;
			}
		}
	
	/******************************************************************************************
		Sort the relators. Put all the relators of the form AAaa or AaaA last. Then put all
		of the relators of the form ABab or AbaB next to last.
	******************************************************************************************/
	
	for(i = j = NumRelators; i > 0; i--)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p;
		if(x == z) return(4);	/* Relator has the form ABAb. Not a knot or link! */
		if(abs(x-z) == 32)		/* Relator has the form ABab or AbaB or AAaa or AaaA. */
			{
			if(x == y || y == z) /* Relator has the form AAaa or AaaA. */
				{
				if(i < j)
					{
					Temp = Relators[j];
					Relators[j] = Relators[i];
					Relators[i] = Temp;
					}
				j--;	
				}
			}
		}
	Mark2 = j;
	
	for(i = j; i > 0; i--)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p;
		if(abs(x-z) == 32)		/* Relator has the form ABab or AbaB */
			{
			if(i < j)
				{
				Temp = Relators[j];
				Relators[j] = Relators[i];
				Relators[i] = Temp;
				}
			j--;	
			}
		}
	Mark1 = j;
	
	for(i = 0; i < MAXNUMGENERATORS; i++) GL[i][0] = GL[i][1] = GL[i][2] = GL[i][3] = 0;
	
	/****************************************************************************************** 
								Set up the array GL[][].
		Check the relators which have the form ABcb with c != a. If this is a Wirtinger 
		presentation of a knot or link, each generator should appear as the first or third 
		char in exactly two of these relators. 
			Check Relators[j]. Increment GL[i][0] by 50 if generator i appears in Relators[j] 
		as first or third char. Store j in GL[i][1] if GL[i][0] = 50. Store j in GL[i][2] if 
		GL[i][0] = 100.	Return if GL[i][0] > 100. 	
	******************************************************************************************/
	
	for(i = 1; i <= Mark1; i++)
		{
		p = *Relators[i];
		x = *p;
		x -= 65;
		z = *(p+2);
		z -= 97;
		if(GL[x][0] == 100 || GL[z][0] == 100) return(5); /* Not a knot or link! */	
		if(GL[x][0] == 0)
			GL[x][1] = i;
		else
			GL[x][2] = i;
		GL[x][0] += 50;
		GL[x][3]++;
		if(GL[z][0] == 0)
			GL[z][1] = i;
		else
			GL[z][2] = i;
		GL[z][0] += 50;
		GL[z][3]++;
		}
		
	/*****************************************************************************************
		Check that each generator appears in two relators in the range 0 < i <= Mark1.
		If a generator appears in only one relator, this may be a Wirtinger presentation
		with n relators of a knot or link with more than n crossings. Wirtinger() expects
		that a presentation of an n crossing knot or link has n relators.
	*****************************************************************************************/	
		
		for(i = 0; i < MAXNUMGENERATORS; i++) if(GL[i][0] == 50) return(6);
			
	/******************************************************************************************
			Deal with those relators which are of the form ABab with A and B distinct.
	******************************************************************************************/
		
	for(i = Mark1 + 1; i <= Mark2; i++)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p++;
		w = *p;	
		x -= 65;
		GL[x][3]++;
		if(y < 'a')
			{
			y -= 65;
			GL[y][3]++;
			if((GL[x][0] == 100) && (GL[y][0] == 100)) return(7); /* Not a knot or link! */
			if((GL[x][0] == 0) && (GL[y][0] == 0))   TheCase = 1;
			if((GL[x][0] == 0) && (GL[y][0] == 100)) TheCase = 2;
			if((GL[x][0] == 100) && (GL[y][0] == 0)) TheCase = 3;
			if(TheCase == 1) /* There is a sublink which is a Hopf link. */				
				{
				GL[x][1] = GL[x][2] = i;
				GL[x][0] = 100;
				}	
			if(TheCase == 2)
				{
				GL[x][1] = GL[x][2] = i;
				GL[x][0] = 100;
				}			
			if(TheCase == 3)
				{
				GL[y][1] = GL[y][2] = i;
				GL[y][0] = 100;
				
				/* Rotate Relator[i] so that 'B' becomes the leading char i.e. ABab --> BabA. */
				
				p = *Relators[i];
				x = *p++;
				y = *p++;
				z = *p++;
				w = *p;	
				p = *Relators[i];
				*p++ = y;
				*p++ = z;
				*p++ = w;
				*p = x;
				}			
			}
		else
			{
			y -= 97;
			GL[y][3]++;
			if((GL[x][0] == 100) && (GL[y][0] == 100)) return(8); /* Not a knot or link! */
			if((GL[x][0] == 0) && (GL[y][0] == 0))   TheCase = 1;
			if((GL[x][0] == 0) && (GL[y][0] == 100)) TheCase = 2;
			if((GL[x][0] == 100) && (GL[y][0] == 0)) TheCase = 3;
			if(TheCase == 1) /* There is a sublink which is a Hopf link. */
				{
				GL[x][1] = GL[x][2] = i;
				GL[x][0] = 100;
				}	
			if(TheCase == 2)
				{
				GL[x][1] = GL[x][2] = i;
				GL[x][0] = 100;
				}			
			if(TheCase == 3)
				{
				GL[y][1] = GL[y][2] = i;
				GL[y][0] = 100;
				
				/* Rotate Relator[i] so that 'B' becomes the leading char i.e. AbaB --> BAba. */
				
				p = *Relators[i];
				x = *p++;
				y = *p++;
				z = *p++;
				w = *p;	
				p = *Relators[i];
				*p++ = w;
				*p++ = x;
				*p++ = y;
				*p = z;
				}			
			}	
		}	
	
	/******************************************************************************************
			Check that each generator appears twice in GL[][3]. This ensures Wirtinger() does 
		not state that ABab is a knot presentation.
	******************************************************************************************/	
		for(i = 0; i < MAXNUMGENERATORS; i++) 
			if(GL[i][3] == 1) return(10); /* Not a Wirtinger presentation. */
			
	/******************************************************************************************
					Deal with those relators which are of the form AAaa or AaaA.
	******************************************************************************************/
	
	for(i = Mark2 + 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *p;	
		x -= 65;
		if(GL[x][0]) return(9); /* Not a knot or link! */
		GL[x][1] = GL[x][2] = i;
		GL[x][0] = 100;	
		}	

	/******************************************************************************************
		Find each component of the link and get a meridian and longitude for each component.			
	******************************************************************************************/
	
	T166 = (unsigned char *) NewPtr(sizeof(char)*(2*MAXNUMGENERATORS + 1));
	if(T166 == NULL) Mem_Error();
	
	NumComponents = 0;
	while(1)
		{
		for(i = 0; i < MAXNUMGENERATORS; i++) if(GL[i][0] == 100) break;
		if(i == MAXNUMGENERATORS) break;
		NumComponents ++;
		GL[i][0] = NumComponents;		
		if(WirtingerM[NumComponents] != NULL) DisposeHandle((char **) WirtingerM[NumComponents]);
		WirtingerM[NumComponents] = (unsigned char **) NewHandle(2*sizeof(char));
		if(WirtingerM[NumComponents] == NULL) Mem_Error();
		p = *WirtingerM[NumComponents];	
		*p++ = i + 'A';
		*p = EOS;

		r = T166;
		j = i;
		do
			{
			if(j == i)
				{
				x = j + 'A';
				Rnum1 = GL[j][1];
				Rnum2 = GL[j][2];
				p = *Relators[Rnum1];
				q = *Relators[Rnum2];
				if((*p == x) && (*q == x)) TheCase = 1;
				if((*p == x) && (*q != x)) TheCase = 2;
				if((*p != x) && (*q == x)) TheCase = 3;
				if((*p != x) && (*q != x)) TheCase = 4;
				switch(TheCase)
					{
					case 1:
						{
						if(Rnum1 <= Rnum2) Rnum = Rnum1;
						else Rnum = Rnum2;
						break;
						}
					case 2:
						{
						Rnum = Rnum1;
						break;
						}
					case 3:
						{
						Rnum = Rnum2;
						break;
						}
					case 4:
						{
						if(Rnum1 <= Rnum2) 
							{
							Rnum = Rnum1;
							q = p;
							x = *p++;
							y = *p++;
							z = *p++;
							w = *p;
							*q++ = z - 32;
							*q++ = w;
							*q++ = x + 32;
							*q++ = y;
							}
						else
							{
							Rnum = Rnum2;
							p = q;
							x = *p++;
							y = *p++;
							z = *p++;
							w = *p;
							*q++ = z - 32;
							*q++ = w;
							*q++ = x + 32;
							*q++ = y;
							}	
						break;
						}			
					}
				}
			else
				{
				Rnum1 = GL[j][1];
				Rnum2 = GL[j][2];
				if(Rnum1 == Rnum) Rnum = Rnum2;
				else Rnum = Rnum1;
				p = q = *Relators[Rnum];
				x = j + 'A';
				if(*p != x)
					{
					x = *p++;
					y = *p++;
					z = *p++;
					w = *p;
					*q++ = z - 32;
					*q++ = w;
					*q++ = x + 32;
					*q++ = y;					
					}				
				}	
			p = *Relators[Rnum];
			x = *p++ - 'A';
			GL[x][0] = NumComponents;
			*r++ = *p++;
			j = *p - 97;
			}
		while(j != i);
		j = 0;
		q = T166;
		while(q < r)
			{
			x = *q++;
			if(x < 'a')
				{
				x -= 'A';
				if(GL[x][0] == NumComponents) j++;
				}
			else
				{
				x -= 'a';
				if(GL[x][0] == NumComponents) j--;
				}	
			}
		if(j < 0)
			{
			x = i + 'A';
			j = -j;
			}
		else
			x = i + 'a';
		while(j)
			{
			*r++ = x;
			j--;
			}
		*r = EOS;		
		q = T166;		
		if(WirtingerL[NumComponents] != NULL) DisposeHandle((char **) WirtingerL[NumComponents]);
		WirtingerL[NumComponents] = (unsigned char **) NewHandle(sizeof(char)*(r - q + 1));
		if(WirtingerL[NumComponents] == NULL) Mem_Error();
		r = *WirtingerL[NumComponents];		
		while( (*r++ = *q++) ) ;
		}
	
	DisposePtr((unsigned char *) T166);
	
	if(NumComponents == 0) return(14);
	
	Knot_Or_Link = TRUE;
	
	if(NumComponents == 1)
		{
		printf("\n\n		This may be a Wirtinger presentation of a knot.");
		printf("\n\n		TRY DEHN-FILLING ON THIS HYPOTHETICAL KNOT ?  HIT 'y' OR 'n'.    ");
		}
	else
		{	
		printf("\n\n		This may be a Wirtinger presentation of a %d-component link.",
			NumComponents);
		printf("\n\n		TRY DEHN-FILLING(S) ON THIS HYPOTHETICAL LINK ? HIT 'y' OR 'n'.    ");
		}
	GET_RESPONSE1:	
	switch(WaitkbHit())
		{
		case 'y':
			break;
		case 'n':
			printf("\n");
			
			/* Restore the original set of relators. */
			
			NumRelators = NumKnot_Or_Link_Rel;
			for(i = 1; i <= NumRelators; i++)
				{
				if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
				Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) KorLRelators[i]));
				if(Relators[i] == NULL) Mem_Error();
				p = *Relators[i];	
				q = *KorLRelators[i];
				while( (*p++ = *q++) ) ;
				}	
			Knot_Or_Link = FALSE;	
			return(17);
		default:
			if(Batch == FALSE) SysBeep(5);
			goto GET_RESPONSE1;
		}
	
	if(NumComponents > 1)
	printf("\n\n		NOTE: TO SKIP FILLING A PARTICULAR BOUNDARY COMPONENT HIT 'return'.");
	else
	printf("\n\n		NOTE: TO SKIP DOING THIS FILLING HIT 'return'.");
	
	SavedPresentation = FALSE;
	ptr = (unsigned char *) NewPtr(200);
	if(ptr == NULL) Mem_Error();
	for(i = 1; i <= NumComponents; i++)
		{
		if(NumRelators >= MAXNUMRELATORS)
			{
			printf("\n\nThis presentation now has %d relators.",NumRelators);
			printf("\nHeegaard cannot accept more than %d relators.",NumRelators);
			if(i == 1)
				printf("\nSo no fillings are possible. Sorry!");
			else
				printf("\nSo no more fillings are possible. Sorry!");
			break;
			}
		printf("\n\n				Component %d:",i);		
		printf("\nMeridian:  %s",*WirtingerM[i]);
		printf("\nLongitude: %s",*WirtingerL[i]);
		
GET_COEFFICIENTS:
		printf("\n\nPlease specify the slope M/L of the surgery on this component or hit 'return'.");
REGET_COEFFICIENT_M:
		printf("\n1) ENTER THE VALUE OF M AND HIT 'return'.    ");	
		ReadString((char *)ptr, GetPtrSize(ptr));
		p = ptr;
		x = *p++;
		if(x == EOS || x == '\n') continue;
		if(!isdigit(x) && x != '+' && x != '-')
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\nThe first character must be a digit, '+' or '-' sign!");
			goto REGET_COEFFICIENT_M;
			}
		while( (x = *p++) ) if(!isdigit(x))
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\nEach character must be a digit or a leading '+' or '-' sign!");
			goto REGET_COEFFICIENT_M;
			}
		sscanf((char *)ptr,"%d",&M);
REGET_COEFFICIENT_L:
		printf("2) ENTER THE VALUE OF L AND HIT 'return'.    ");
		ReadString((char *)ptr, GetPtrSize(ptr));
		p = ptr;
		x = *p++;
		if(!isdigit(x) && x != '+' && x != '-')
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\nThe first character must be a digit, '+' or '-' sign!\n");
			goto REGET_COEFFICIENT_L;
			}
		while( (x = *p++) ) if(!isdigit(x)) 
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
			goto REGET_COEFFICIENT_L;
			}
		sscanf((char *)ptr,"%d",&L);
		if(M == 0 && L == 0)
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\nBoth M and L are zero!");
			printf("\nPlease specify different surgery coefficients for this component!");
			goto GET_COEFFICIENTS;
			}
		h = GCD(abs(M),abs(L));
		if(L < 0)
			{
			M = -M;
			L = -L;
			}	
		L = L/h;
		if(M < 0)
			{
			M = -M;
			M = M/h;
			M = -M;
			}
		else	
			M = M/h;
			
		/* Since meridians always have length one, surgery relators have length given by: */
					
		NRL = abs(M) + abs(L)*(GetHandleSize((char **) WirtingerL[i]) - 1);			
		if(NRL > MAXLENGTH)
			{
			if(Batch == FALSE) SysBeep(5);
			printf("\n\nThese coefficients give a new relator which will probably be too long!");
			printf("\n\nPlease specify new surgery coefficients for this component.");
			goto GET_COEFFICIENTS;
			}
		if(SavedPresentation == FALSE)
			{
			SavedPresentation = TRUE;
			printf("\n\nTried Dehn-filling on: %s\n",PresName);
			Print_Relators(Relators,NumRelators);
			}			
		printf("\nTried %d/%d surgery with meridian: %s and longitude: %s.",M,L,
			*WirtingerM[i],*WirtingerL[i]);							
		NumRelators ++;	
		if(Relators[NumRelators] != NULL) DisposeHandle((char **) Relators[NumRelators]);
		Relators[NumRelators] = (unsigned char **) NewHandle(sizeof(char)*(NRL + 1));
		if(Relators[NumRelators] == NULL) Mem_Error();
		LR[NumRelators] = NRL;
		q = *Relators[NumRelators];	
		r = q;
		if(M < 0)
			{
			Inverse(*WirtingerM[i]);
			M = -M;
			}	
		h = 0;
		do
			{
			if(h < M)
				{
				p = *WirtingerM[i];
				while( (*q++ = *p++) ) ;
				q--;
				}
			else
				{
				p = *WirtingerL[i];
				while( (*q++ = *p++) ) ;
				q--;
				}
			if(h >= L)
				h -= L;
			else
				h += M;	
			}
		while(h);
		if((q-r) != NRL)
			{
			printf("\n q-r = %lu, NRL = %ld in Wirtinger.",q-r,NRL);
			printf("%s",*Relators[NumRelators]);
			NumErrors ++;
			}	
		}
	DisposePtr((char *) ptr);
	return(0);		
}

int Try_Exponent_Surgery()
{
	register unsigned char	*p,
							*q,
							*ptr = NULL,
							x,
							y;
							
	char					Change[125];
	
	unsigned char 			**Temp;
	
	int			e,
				i,
				j,
				k,
				P,
				Q,
				NewEXP[MAXNUMGENERATORS][3];
	
	unsigned int		gcd;
							
	unsigned long		length;
	
	unsigned int 		GCD(),
						Whitehead_Graph();

	if(Find_Flow_A(NORMAL,FALSE))
		{
		printf("\n\nExponent surgery is not possible for this presentation. Sorry!\n");
		return(1);
		}
	if(Automorphisms && Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG)
		{
		printf("\n\nExponent surgery is not possible for this presentation. Sorry!\n");
		return(1);
		}
	Fill_A(NumRelators);		
	Saved_Vertices = 0;
	DrawingDiagrams = TRUE;
	WhichInput = 0;
	if(Whitehead_Graph())
		{
		printf("\n\nExponent surgery is not possible for this presentation. Sorry!\n");
		DrawingDiagrams = FALSE;
		return(1);
		}
	DrawingDiagrams = FALSE;	
	
	for(i = 0; i < NumGenerators; i++) if(EXP[i][2] > 2) break;
	if(i == NumGenerators)
		{
		printf("\n\nExponent surgery is not possible for this presentation. Sorry!\n");
		return(1);
		}
	
	switch(Automorphisms)
		{
		case 0:
			break;
		case 1:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\none automorphism reduced its length to %lu.",Length);
			break;
		default:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\n%lu automorphisms reduced its length to %lu.",Automorphisms,Length);
			break;
		}
		
	printf("\n\nThe presentation is currently:\n");
	Print_Relators(Relators,NumRelators);

	/******************************************************************************************
		Save a copy of the initial relators so that the user can perform additional exponent
		surgeries.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		if(Exp_Surgery_Rel[i] != NULL) DisposeHandle((char **) Exp_Surgery_Rel[i]);
		Exp_Surgery_Rel[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(Exp_Surgery_Rel[i] == NULL) Mem_Error();
		p = *Exp_Surgery_Rel[i];	
		q = *Relators[i];
		while( (*p++ = *q++) ) ;
		}
		
	NumExp_Sur_Rel = NumRelators;
	
	GET_NEW_EXPONENTS:	
	ptr = (unsigned char *) NewPtr(100);
	if(ptr == NULL) Mem_Error();
	for(i = 0; i < NumGenerators; i++) Change[i+'A'] = Change[i+'a'] = FALSE;	
	for(i = 0; i < NumGenerators; i++) if(EXP[i][2] > 2)
		{
		Change[i+'A'] = Change[i+'a'] = TRUE;
		if(EXP[i][0])
			{
			printf("\n\nGENERATOR '%c' APPEARS WITH EXPONENTS P = %u, Q = %u, AND P+Q = %u.",
					i + 'A',EXP[i][0],EXP[i][1],EXP[i][2]);
			GET_NEW_EXPONENTS1P:		
			printf("\nREPLACE P WITH:     ");
			P = EXP[i][0];
			ReadString((char *)ptr, GetPtrSize(ptr));
			q = ptr;			
			x = *q++;
			if(!isdigit(x) && x != '+' && x != '-')
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\nThe first character must be a digit, '+' or '-' sign!");
				goto GET_NEW_EXPONENTS1P;
				}
			while( (x = *q++) ) if(!isdigit(x))
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\nEach character must be a digit or a leading '+' or '-' sign!");
				goto GET_NEW_EXPONENTS1P;
				}
			sscanf((char *) ptr,"%d",&P);
			GET_NEW_EXPONENTS1Q:
			printf("REPLACE Q WITH:     ");
			Q = EXP[i][1];
			ReadString((char *)ptr, GetPtrSize(ptr));
			q = ptr;
			x = *q++;
			if(!isdigit(x) && x != '+' && x != '-')
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\nThe first character must be a digit, '+' or '-' sign!\n");
				goto GET_NEW_EXPONENTS1Q;
				}
			while( (x = *q++) ) if(!isdigit(x))
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
				goto GET_NEW_EXPONENTS1Q;
				}		
			sscanf((char *) ptr,"%d",&Q);
			if(P == 0 && Q == 0)
				{
				if(Batch == FALSE) SysBeep(5);
				printf("\n\nBoth P and Q are zero!");
				goto GET_NEW_EXPONENTS1P;
				}
			gcd = GCD(abs(P),abs(Q));
			if(gcd != 1)
				{
				printf("\n\nP and Q are not relatively prime! Dividing P and Q by their GCD %u.",gcd);
				if(P < 0)
					{
					P = -P;
					P /= gcd;
					P = -P;
					}
				else
					P /= gcd;
				if(Q < 0)
					{
					Q = -Q;
					Q /= gcd;
					Q = -Q;
					}
				else
					Q /= gcd;
				printf("\nSetting P = %d, Q = %d.",P,Q);				
				}
			NewEXP[i][0] = P;
			NewEXP[i][1] = Q;
			NewEXP[i][2] = P + Q;	
			}
		else
			{
			if(EXP[i][1])
				{
				printf("\n\nGENERATOR '%c' APPEARS WITH EXPONENTS P = %u, AND Q = %u.",
					i + 'A',EXP[i][1],EXP[i][2]);
				GET_NEW_EXPONENTS2P:	
				printf("\nREPLACE P WITH:     ");
				P = EXP[i][1];		
				ReadString((char *)ptr, GetPtrSize(ptr));
				q = ptr;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!");
					goto GET_NEW_EXPONENTS2P;
					}
				while( (x = *q++) ) if(!isdigit(x))
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!");
					goto GET_NEW_EXPONENTS2P;
					}				
				sscanf((char *) ptr,"%d",&P);
				GET_NEW_EXPONENTS2Q:
				printf("REPLACE Q WITH:     ");
				Q = EXP[i][2];
				ReadString((char *)ptr, GetPtrSize(ptr));
				q = ptr;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!\n");
					goto GET_NEW_EXPONENTS2Q;
					}
				while( (x = *q++) ) if(!isdigit(x))
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
					goto GET_NEW_EXPONENTS2Q;
					}								
				sscanf((char *) ptr,"%d",&Q);
				if(P == 0 && Q == 0)
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\n\nBoth P and Q are zero!");
					goto GET_NEW_EXPONENTS2P;
					}
				gcd = GCD(abs(P),abs(Q));
				if(gcd != 1)
					{
					printf("\n\nP and Q are not relatively prime! Dividing P and Q by their GCD.");
					if(P < 0)
						{
						P = -P;
						P /= gcd;
						P = -P;
						}
					else
						P /= gcd;
					if(Q < 0)
						{
						Q = -Q;
						Q /= gcd;
						Q = -Q;
						}
					else
						Q /= gcd;
					printf("\nSetting P = %d, Q = %d.",P,Q);									
					}
				NewEXP[i][0] = 0;
				NewEXP[i][1] = P;
				NewEXP[i][2] = Q;		
				}
			else
				{
				printf("\n\nGENERATOR '%c' APPEARS ONLY WITH EXPONENT P = %u.",
					i + 'A',EXP[i][2]);
				GET_NEW_EXPONENTS3P:	
				printf("\nREPLACE P WITH:     ");
				P = EXP[i][2];
				ReadString((char *)ptr, GetPtrSize(ptr));
				q = ptr;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!");
					goto GET_NEW_EXPONENTS3P;
					}
				while( (x = *q++) ) if(!isdigit(x))
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!");
					goto GET_NEW_EXPONENTS3P;
					}							
				sscanf((char *) ptr,"%d",&P);
				NewEXP[i][0] = 0;
				NewEXP[i][1] = 0;
				NewEXP[i][2] = P;	
				}	
			}	
		}
	DisposePtr((char *) ptr);
	
	if(Temp16 != NULL) DisposeHandle((char **) Temp16);
	Temp16 = (unsigned char **) NewHandle(MAXLENGTH + 2);
	if(Temp16 == NULL) Mem_Error();
	q = *Temp16;
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		q = *Temp16;
		length = 0L;
		while( (x = *p) )
			{
			if(Change[x])
				{
				e = 0;
				while(*p == x)
					{
					e++;
					p++;
					}
				y = x;
				if(x < 'a')
					x -= 65;
				else
					x -= 97;	
				for(j = 0; j < 3; j++) if(EXP[x][j] == e) break;
				k = NewEXP[x][j];
				if(k < 0)
					{
					k = -k;
					if(y < 'a')
						y += 32;
					else
						y -= 32;	
					}
				length += k;					
				if(length > MAXLENGTH)
					{
					if(Batch == FALSE) SysBeep(5);
					*q = EOS;
					printf("\n\nRelator %d will be too long!",i);
					printf("\n\nTRY DIFFERENT EXPONENTS ?  HIT 'y' OR 'n'.");
					GET_RESPONSE1:					
					switch(WaitkbHit())
						{
						case 'y':
							goto GET_NEW_EXPONENTS;
						case 'n':
							return(1);
						default:
							if(Batch == FALSE) SysBeep(5);
							goto GET_RESPONSE1;
						}
					}							
				for(j = 0; j < k; j++) *q++ = y;		
				}
			else
				{
				length ++;
				if(length > MAXLENGTH)
					{
					if(Batch == FALSE) SysBeep(5);
					*q = EOS;
					printf("\n\nRelator %d will be too long!",i);
					printf("\n\nTRY DIFFERENT EXPONENTS ?  HIT 'y' OR 'n'.");
					GET_RESPONSE2:					
					switch(WaitkbHit())
						{
						case 'y':
							goto GET_NEW_EXPONENTS;
						case 'n':
							return(1);
						default:
							if(Batch == FALSE) SysBeep(5);
							goto GET_RESPONSE2;
						}
					}				
				*q++ = x;
				p++;
				}	
			}
		*q = EOS;
		
		if(OutRelators[i] != NULL) DisposeHandle((char **) OutRelators[i]);
		OutRelators[i] = (unsigned char **) NewHandle(length + 1);
		if(OutRelators[i] == NULL) Mem_Error();
		q = *OutRelators[i];	
		p = *Temp16;
		while( (*q++ = *p++) ) ;
		}
		
	for(i = 1; i <= NumRelators; i++)
		{
		Temp = Relators[i];
		Relators[i] = OutRelators[i];
		OutRelators[i] = Temp;
		LR[i] = GetHandleSize((char **)Relators[i]) - 1;
		}
		
	printf("\n\nThe Exponent-Surgered presentation is:\n");
	Print_Relators(Relators,NumRelators);
		
	return(0);			
}

int Try_Cutting_Disk_Surgery()
{
	/******************************************************************************************
		Try_Cutting_Disk_Surgery() is a routine which provides a way to obtain different 
	Heegaard diagrams from a diagram with a given Whitehead graph by allowing the user to
	change how the ends of edges meeting pairs of inverse vertices are identified. 
		Suppose (X,x) is a pair of inverse vertices of the Whitehead graph, the identification
	of X and x is specifed by an integer Vx in the range 0 <= Vx < V[X], where V[X] is the 
	valence of X and x in the Whitehead graph. By choosing different values for the offsets, 
	you can obtain a Heegaard diagram of any 3-manifold which whose Whitehead graph is 
	isomorphic to the Whitehead graph of the initial diagram. 
		(Note however, that the new presentation can't have more than MAXNUMRELATORS relators.)
	*******************************************************************************************/
	register unsigned char	*p,
							*q,
							*ptr = NULL,
							x;
	
	unsigned char 			**Temp;
	
	int						i,j;
	
	unsigned int			NOSA;
	
	long					Scratch;						
	
	unsigned int 			Whitehead_Graph();

	if(Find_Flow_A(NORMAL,FALSE))
		{
		printf("\n\nSurgery is not possible for this presentation. Sorry!\n");
		return(1);
		}
	if(Automorphisms && Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG)
		{
		printf("\n\nSurgery is not possible for this presentation. Sorry!\n");
		return(1);
		}
	Fill_A(NumRelators);		
	Saved_Vertices = 0;
	DrawingDiagrams = TRUE;
	WhichInput = 0;
	if(Whitehead_Graph())
		{
		printf("\n\nSurgery is not possible for this presentation. Sorry!\n");
		DrawingDiagrams = FALSE;
		return(1);
		}
	DrawingDiagrams = FALSE;	
	
	switch(Automorphisms)
		{
		case 0:
			break;
		case 1:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\none automorphism reduced its length to %lu.",Length);
			break;
		default:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\n%lu automorphisms reduced its length to %lu.",Automorphisms,Length);
			break;
		}
		
	/******************************************************************************************
		Save a copy of the initial relators so that the user can perform additional surgeries.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		if(CD_Surgery_Rel[i] != NULL) DisposeHandle((char **) CD_Surgery_Rel[i]);
		CD_Surgery_Rel[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(CD_Surgery_Rel[i] == NULL) Mem_Error();
		p = *CD_Surgery_Rel[i];	
		q = *Relators[i];
		while( (*p++ = *q++) ) ;
		}
		
	NumCuttingDiskSurgeryRel = NumRelators;
	
	printf("\n\n	Cutting disk surgery provides a way to obtain different Heegaard diagrams");
	printf("\n	by changing how the ends of edges meeting pairs of inverse vertices are");
	printf("\n  identified. By choosing different values for the offsets, you can obtain");
	printf("\n  a Heegaard diagram of any 3-manifold which has a Heegaard diagram whose");
	printf("\n  Whitehead graph is isomorphic to the Whitehead graph of the initial diagram.");
	printf("\n	(Note however, that the new presentation can't have more than %d relators.)",MAXNUMRELATORS);		

	printf("\n\nThe presentation is currently:\n");
	Print_Relators(Relators,NumRelators);

	RETRY:
	
	/******************************************************************************************						
						Print the array V[] and the array OSA[] % V[].
	******************************************************************************************/
	
	printf("\n\n        The following arrays give the 'valence' of each generator and the");
	printf("\n        current values of the offsets.\n");
	
	do
		{
		printf("\nGens:");	
		for(i = 1,j = 5; i <= NumGenerators && j < 85; i++) j += printf("%5d",i);
		printf("\nVals:");	
		for(i = 0,j = 5; i < Vertices && j < 85; i += 2) j += printf("%5u",V[i]);	
		printf("\nOSet:");	
		for(i = 0 ,j = 5; i < Vertices && j < 85; i += 2) j += printf("%5u",OSA[i] % V[i]);
		}
	while(i < Vertices);
		
	printf("\n\nPlease enter new non-negative values for the offsets of each generator.\n\n");	
							
	ptr = (unsigned char *) NewPtr(100);
	if(ptr == NULL) Mem_Error();
	for(i = 0; i < Vertices; i += 2)
		{
		printf("Generator %2d, Valence %5u, Current offset %5u. Change offset to ?    ",
			i/2 + 1,V[i],OSA[i] % V[i]);
		GET_NEW_OFFSET:
		ReadString((char *)ptr, GetPtrSize(ptr));
		q = ptr;
		while( (x = *q++) ) if(!isdigit(x))
			{
			if(Batch == FALSE) SysBeep(5);
			goto GET_NEW_OFFSET;
			}
		OSB[i] = OSA[i];
		OSB[i+1] = OSA[i+1];	
		sscanf((char *) ptr,"%u",&NOSA);
		NOSA = NOSA % V[i];
		NOSA += V[i];	
		OSA[i] = OSA[i+1] = NOSA;				
		}
	DisposePtr((char *) ptr);
	
	if(Cutting_Disk_Surgery_Diagram_7())
		{
		printf("\n\nTRY ANOTHER SURGERY ON THIS PRESENTATION ? HIT 'y' OR 'n'.    ");
		GET_RESPONSE1:		
		switch(WaitkbHit())
			{
			case 'y':
				goto RETRY;
			case 'n':
				Did_Cutting_Disk_Surgery = FALSE;
				return(1);
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE1;		
			}
		}
		
	for(i = 1; i <= NumRelators; i++)
		{
		Temp = Relators[i];
		Relators[i] = OutRelators[i];
		OutRelators[i] = Temp;
		LR[i] = GetHandleSize((char **)Relators[i]) - 1;
		}
		
	/**************************************************************************************
		Echo the surgered relators to the output so we will have a copy of them. Then call
		Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
		which will serve as the initial presentation for Heegaard. 
	**************************************************************************************/	
	
	fprintf(stdout,"\n\nThe surgered presentation is:\n");
	Print_Relators(Relators,NumRelators);
	for(i = 1,Scratch = 0L; i <= NumRelators; i++)
		Scratch += GetHandleSize((char **) Relators[i]);
	Scratch -= NumRelators;
	printf("\n\nThis presentation has length %ld ",Scratch);	
	if(Freely_Reduce() == TOO_LONG)
		{
		printf("\n\nThis presentation is too long!!");
		if(Batch == FALSE) SysBeep(5);
		printf("\n\nTRY ANOTHER SURGERY ON THIS PRESENTATION ? HIT 'y' OR 'n'.    ");
		GET_RESPONSE2:		
		switch(WaitkbHit())
			{
			case 'y':
				goto RETRY;
			case 'n':
				Did_Cutting_Disk_Surgery = FALSE;
				return(1);
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE2;		
			}
		}
	if(Scratch > OrigLength)
		{
		printf("and freely reduces to length %lu.",OrigLength);
		Scratch = OrigLength;
		}
	else
		printf("and is freely reduced.");
	if(Rewrite_Input())
		{
		printf("\n\nThere must be at least one non-empty relator!!");
		if(Batch == FALSE) SysBeep(5);
		printf("\n\nTRY ANOTHER SURGERY ON THIS PRESENTATION ? HIT 'y' OR 'n'.    ");
		GET_RESPONSE3:		
		switch(WaitkbHit())
			{
			case 'y':
				goto RETRY;
			case 'n':
				Did_Cutting_Disk_Surgery = FALSE;
				return(1);
			default:
				if(Batch == FALSE) SysBeep(5);
				goto GET_RESPONSE3;		
			}
		}
	Length = OrigLength;	
									
	/**************************************************************************************
			Call Canonical_Rewrite() to rewrite the presentation in canonical form.
	**************************************************************************************/
			
	Canonical_Rewrite(Relators,FALSE,FALSE);					
	printf("\n\nNumRelators = %d, NumGenerators = %d\n",NumRelators, NumGenerators);
	
	printf("\n\nThe rewritten initial presentation is:\n");
	Print_Relators(Relators,NumRelators);

	printf("\n\nUSE THIS PRESENTATION ? HIT 'y' OR 'n'.    ");
	GET_RESPONSE4:	
	switch(WaitkbHit())
		{
		case 'y':
			break;
		case 'n':
			goto RETRY;
		default:
			if(Batch == FALSE) SysBeep(5);
			goto GET_RESPONSE4;		
		}
	
	/**************************************************************************************
		Save a copy of the initial set of relators in case we want to refer to them later.
	**************************************************************************************/
		
	CopyNumRelators 	= NumRelators;
	CopyNumGenerators 	= NumGenerators;
	for(i = 1; i <= NumRelators; i++)
		{
		if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
		Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
		if(Copy_Of_Input[i] == NULL) Mem_Error();
		p = *Copy_Of_Input[i];	
		q = *Relators[i];
		while( (*p++ = *q++) ) ;
		}
	
	Did_Cutting_Disk_Surgery = TRUE;	
	return(0);			
}

int Cutting_Disk_Surgery_Diagram_7()
{
	/******************************************************************************************
		Cutting_Disk_Surgery_Diagram_7() "runs" around the diagram and reads off relators. It
		leaves these relators in the array OutRelators[]. While it is reading off relators,
		the routine also determines what the dual relators are. It leaves the dual relators in
		the array DualRelators[].
	******************************************************************************************/
			
	register unsigned char 	i,
							*p,
							*q,
							*r,
							v,
							w;
							
	register unsigned int 	d,
							e;
							
	int 					j,
							SNumRelators;
							
	unsigned int 			length,
							max;
	
	for(i = 1; i <= NumGenerators; i++)
		{
		max = V[(i-1) << 1] + 1;
		if(DualRelators[i] != NULL) DisposeHandle((char **) DualRelators[i]);
		DualRelators[i] = (unsigned char **) NewHandle(max);
		if(DualRelators[i] == NULL) Mem_Error();
		p = *DualRelators[i];	
		for(e = V[(i-1) << 1]; e > 0; e--) *p++ = '@';
		*p = EOS;
		}	

	max = MAXLENGTH;
	SNumRelators = NumRelators;
	NumRelators = 0;
	
	if(Temp16 != NULL) DisposeHandle((char **) Temp16);	
	Temp16 = (unsigned char **) NewHandle(max + 2);
	if(Temp16 == NULL) Mem_Error();
		
	for(j = 1,i = 0; j <= NumGenerators; j++) 
		{
		p = *DualRelators[j];
		while(*p)
			{
			if(*p == '@')
				{
				NumRelators ++;
				if(NumRelators > MAXNUMRELATORS)
					{
					if(Batch == FALSE) SysBeep(5);
					printf("\n\nThis presentation now has %d relators.",NumRelators);
					printf("\nHeegaard cannot accept more than %d relators. Sorry!",MAXNUMRELATORS);
					NumRelators = SNumRelators;
					return(TOO_LONG);
					}
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
						q = *DualRelators[(v >> 1) + 1] + e;
						*q = i + 97;
						w = v - 1;
						}
					else 
						{
						*r = (v >> 1) + 65;
						q = *DualRelators[(v >> 1) + 1] + e;
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
						if(Batch == FALSE) SysBeep(5);
						printf("\n\nRelator %d of this presentation is longer than %d.",
							NumRelators, MAXLENGTH);
						printf("\nThis is too long for Heegaard to handle. Sorry!");
						NumRelators = SNumRelators;	
						return(TOO_LONG);
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
				while( (*q++ = *r++) ) ;
				i++;
				}
			p++;
			}	
		}
	return(NO_ERROR);	
}

int Stabilize()
{
	unsigned char	c,
					d,
					*p,
					*q,
					*r,
					X,
					x,
					Y,
					y,
					Z,
					z;
								
	unsigned char 	*ptr,
					**Temp;
	
	int				i,
					j;
									
	long			MaxHS;

	if(NumGenerators >= MAXNUMGENERATORS)
		{
		printf("\n\nCan't stabilize. This presentation already has the maximum allowed number of generators.");
		return(1);
		}
		
	if(NumRelators >= MAXNUMRELATORS)
		{
		printf("\n\nCan't stabilize. This presentation already has the maximum allowed number of relators.");
		return(1);
		}			

	if(Freely_Reduce() == TOO_LONG)
		{
		printf("\n\nThis presentation is too long!!");
		if(Batch == FALSE) SysBeep(5);
		return(1);
		}	
									
	if(Rewrite_Input())
		{
		printf("\n\nThere must be at least one non-empty relator!!");
		if(Batch == FALSE) SysBeep(5);
		return(1);
		}	
	
	Vertices = 2*NumGenerators;
	Fill_A(NumRelators);
	Get_Matrix();								

	if(Find_Flow_A(NORMAL,FALSE))
		{
		printf("\n\nStabilization is not possible for this presentation. Sorry!\n");
		return(1);
		}
	
	i = Input;
	Input = BEGIN;	
	if(Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG)
		{
		printf("\n\nStabilization is not possible for this presentation. Sorry!\n");
		Input = i;
		return(1);
		}
	Input = i;	

	switch(Automorphisms)
		{
		case 0:
			break;
		case 1:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\none automorphism reduced its length to %lu.",Length);
			break;
		default:
			printf("\n\nThe initial presentation did not have minimal length;");
			printf("\n%lu automorphisms reduced its length to %lu.",Automorphisms,Length);
			break;
		}
	
	if(Automorphisms)
		{	
		printf("\n\nThe presentation is currently:\n");
		Print_Relators(Relators,NumRelators);
		Vertices = 2*NumGenerators;
		Fill_A(NumRelators);
		Get_Matrix();	
		}
		
	for(i = 0; i < Vertices; i++) ZZ[i] = 0;
	if(Connected_(0,0) == FALSE)
		{
		printf("\n\nCan't stabilize because the Whitehead graph of this presentation is not connected.");
		return(1);
		}
	
	if(Planar(TRUE,FALSE))
		{
		printf("\n\nCan't stabilize because the Whitehead graph of this presentation is not planar.");
		return(1);		
		}
		
	if(Sep_Pairs(0,0,1) && (Level_Transformations(FALSE,FALSE,TRUE) != 2))
		{
		printf("\n\nCan't stabilize! Could not remove all separating pairs of vertices.");
		return(1);		
		}	
		
	/******************************************************************************************
		Choose a random cyclic relator of length greater than one and a random consecutive 
		pair of letters in that relator.
	******************************************************************************************/

	j = NumRelators;
	do
		{
		i = abs(rand()) % j;
		i++;
		MaxHS = GetHandleSize((char **) Relators[i]);
		j--;
		if(j == 0) break;
		}
	while(MaxHS <= 2);		
	
	MaxHS --;
	j = abs(rand()) % MaxHS;
	
	p = *Relators[i];
	X = *(p + j);
	Y = *(p + j + 1);
	if(Y == EOS) Y = *p;
	if(X >= 97) 
		x = X - 32;
	else
		x = X + 32;
	if(Y >= 97) 
		y = Y - 32;
	else
		y = Y + 32;	
							
	/*******************************************************************************************
			Let Z be the new generator which we will be adding. Then Stabilize() will replace 
		each appearance of XY in the relators with XZY, and it will replace each appearance of 
		yx in the relators with yzx. Since we have verified that the Whitehead graph of the 
		incoming presentation P is planar, connected, and has no separating pairs of vertices, 
		all of the edges of the Heegaard diagram of P which connect vertices x and Y are
		parallel in the Heegaard surface. This means that the new presentation P' produced by 
		Stabilize() will be realizable if P is realizable.
	*******************************************************************************************/						

	Z = 'A' + NumGenerators;
	z = 'a' + NumGenerators;		
	
	for(i = 1, MaxHS = 0; i <= NumRelators; i++)
		{
		if(GetHandleSize((char **) Relators[i]) > MaxHS)
		MaxHS = GetHandleSize((char **) Relators[i]);
		}
	
	ptr = (unsigned char *) NewPtr(2*MaxHS);
	if(ptr == NULL) Mem_Error();
	
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		r = ptr;
		while( (c = *p) )
			{
			*r++ = c;
			if(c == X)
				{
				d = *(p+1);
				if(d == EOS) d = **Relators[i];
				if(d == Y) *r++ = Z;
				}
			if(c == y)
				{
				d = *(p+1);
				if(d == EOS) d = **Relators[i];
				if(d == x) *r++ = z;
				}	
			p++;
			}
		*r = EOS;		
		p = ptr;	
		if(r - p > MAXLENGTH)
			{
			printf("\n\nCan't stabilize because a stabilized relator will exceed the current maximum allowed length.");
			DisposePtr((char *) ptr);
			return(1); 
			}
				
		if(OutRelators[i] != NULL) DisposeHandle((char **) OutRelators[i]);
        OutRelators[i] = (unsigned char **) NewHandle(r - p + 1);
        if(OutRelators[i] == NULL) Mem_Error();
        p = *OutRelators[i];
        q = ptr;
        while( (*p++ = *q++) ) ;
		}
			
	DisposePtr((char *) ptr);
	
	for(i = 1; i <= NumRelators; i++)
		{
		Temp = Relators[i];
		Relators[i] = OutRelators[i];
		OutRelators[i] = Temp;
		}
		
	NumGenerators ++;	
	NumRelators ++;
	Vertices = 2*NumGenerators;
	if(Relators[NumRelators] != NULL) DisposeHandle((char **) Relators[NumRelators]);
    Relators[NumRelators] = (unsigned char **) NewHandle(2);
    if(Relators[NumRelators] == NULL) Mem_Error();
    p = *Relators[NumRelators];
    *p++ = Z;
    *p = EOS;

	for(i = 1,Length = 0; i <= NumRelators; i++) Length += GetHandleSize((char **)Relators[i]);
	Length -= NumRelators;
	
	printf("\n\nThe stabilized presentation is:\n");
	Print_Relators(Relators,NumRelators);
	if(Batch == 13 && H_Results != NULL) 
		{
		fprintf(H_Results,"\n\n%s",PresName);
		Print_Relators2(Relators,NumRelators);
		}
		
	return(0);			
}

void Heegaard_Splash_Screen()
{
	fprintf(Gvizdata, "graph G{layout = neato; model = circuit; size = \04210.0,8.0\042; ratio = fill;");
	fprintf(Gvizdata, "\nlabel = \042A genus 21 Heegaard diagram.\042"); 
	fprintf(Gvizdata, "\nnode [shape = circle, fontsize = 10, height = 0.1, style = white]");
	fprintf(Gvizdata, "\nA [pos = \04230,30!\042]; a [pos = \04299,63!\042]; B [pos = \042149,93!\042]; b [pos = \04283,350!\042];");
	fprintf(Gvizdata, "\nC [pos = \04297,381!\042]; c [pos = \042550,264!\042]; D [pos = \042550,342!\042]; d [pos = \042412,337!\042];"); 
	fprintf(Gvizdata, "\nE [pos = \042550,420!\042]; e [pos = \042236,175!\042]; F [pos = \042303,110!\042]; f [pos = \042290,30!\042];"); 
	fprintf(Gvizdata, "\nG [pos = \042376,30!\042]; g [pos = \04230,342!\042]; H [pos = \04230,264!\042]; h [pos = \04230,108!\042];");
	fprintf(Gvizdata, "\nI [pos = \04230,186!\042]; i [pos = \04230,420!\042]; J [pos = \042134,420!\042]; j [pos = \042142,213!\042];"); 
	fprintf(Gvizdata, "\nK [pos = \04292,257!\042]; k [pos = \042446,420!\042]; L [pos = \042342,420!\042]; l [pos = \042469,94!\042];");
	fprintf(Gvizdata, "\nM [pos = \042452,136!\042]; m [pos = \042463,30!\042]; N [pos = \042550,30!\042]; n [pos = \042550,108!\042];"); 
	fprintf(Gvizdata, "\nO [pos = \042506,92!\042]; o [pos = \042381,101!\042]; P [pos = \04275,315!\042]; p [pos = \042377,377!\042];"); 
	fprintf(Gvizdata, "\nQ [pos = \042378,336!\042]; q [pos = \042550,186!\042]; R [pos = \042319,395!\042]; r [pos = \042507,218!\042];");
	fprintf(Gvizdata, "\nS [pos = \042341,281!\042]; s [pos = \042238,420!\042]; T [pos = \042203,30!\042]; t [pos = \04276,170!\042];"); 	
	fprintf(Gvizdata, "\nU [pos = \042107,188!\042]; u [pos = \042116,30!\042];"); 
	fprintf(Gvizdata, "\nedge [fontsize = 10]; { A -- a ; A -- h ; A -- u ; a -- B ; a -- u ; B -- e ; B -- h ;"); 
	fprintf(Gvizdata, "\nB -- T ; B -- u ; b -- C ; b -- P ; C -- i ; C -- J ; c -- D ; c -- q ; c -- r ; D -- E ;"); 
	fprintf(Gvizdata, "\nD -- r ; d -- E ; d -- Q ; d -- S ; E -- k ; e -- F ; e -- j ; e -- S ; F -- f ; F -- o ;");
	fprintf(Gvizdata, "\nf -- G ; f -- T ; G -- m ; G -- o ; g -- H ; g -- i ; g -- P ; H -- I ; H -- K ; h -- I ;");
	fprintf(Gvizdata, "\nh -- t ; I -- t ; i -- J ; J -- s ; j -- K ; j -- U ; K -- P ; k -- L ; k -- p ; L -- R ;");
	fprintf(Gvizdata, "\nL -- s ; l -- M ; l -- m ; l -- O ; M -- o ; M -- r ; m -- N ; N -- n ; N -- O ; n -- O ;");
	fprintf(Gvizdata, "\nn -- q ; p -- Q ; p -- R ; Q -- S ; q -- r ; R -- s ; T -- u ; t -- U ; }}");
	fclose(Gvizdata);
}
