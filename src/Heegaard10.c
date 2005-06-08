#include "Heegaard.h"
#include "Heegaard_Dec.h"
#include <ctype.h>


Wirtinger()
{
	unsigned char 			GL[MAXNUMGENERATORS][4],
							*p,
							*q,
							*r,
							**Temp,
							w,
							x,
							y,
							z;
	
	int						i,
							j,
							L,
							M,
							Mark1,
							Mark2,
							NumComponents,
							SavedPresentation;
							
	unsigned int			h;
	
	long					NRL;						
	
	unsigned int			GCD();
	
	/******************************************************************************************
		See if each relator has length 4. If each relator has length 4, see if each relator
		has the form XYxW or XYZy.
	******************************************************************************************/	
	
	for(i = 1; i <= NumRelators; i++) if(GetHandleSize((char **) Relators[i]) != 5L)
		return(1);

	/******************************************************************************************
		Save a copy of the initial relators so that the user can perform additional surgeries
		if the presentation is a Wirtinger presentation of a knot of link.
	******************************************************************************************/
	
	Knot_Or_Link = TRUE;
	for(i = 1; i <= NumRelators; i++)
		{
		ReallocateHandle((char **) KorLRelators[i],GetHandleSize((char **) Relators[i]));				
		if((p = *KorLRelators[i]) == NULL)
			{
			Knot_Or_Link = FALSE;
			break;
			}
		q = *Relators[i];
		while(*p++ = *q++) ;
		}
		
	NumKnot_Or_Link_Rel = NumRelators;		
		
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p++;
		w = *p;
		if(abs(x-z) != 32 && abs(y-w) != 32) return(1);
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
		if((x < 'a' && z < 'a') || (x > 'Z' && z > 'Z')) return(1);
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
		if(x == z) return(1);	/* Relator has the form ABAb. Not a knot or link! */
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
	
	for(i = 0; i < MAXNUMGENERATORS; i++) GL[i][0] = GL[i][2] = GL[i][3] = 0;
	
	/****************************************************************************************** 
		Look for the two crossings which mark the ends of each overpass, in those relators
		which have the form ABCb where C != a.
	******************************************************************************************/
		
	for(i = 1; i <= Mark1; i++)
		{
		p = *Relators[i];
		x = *p;
		z = *(p + 2);
		x -= 65;
		if(GL[x][0]) return(1); /* Not a knot or link! */
		GL[x][0] ++;
		GL[x][3] = 100;
		GL[x][1] = i;
		z -= 97;
		if(GL[z][2]) return(1); /* Not a knot or link! */
		GL[z][2] ++;
		GL[z][3] = 100;
		}
			
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
		if(GL[x][0])
			{
			p = *Relators[i];
			x = *p;
			if(y < 'a')
				{
				*p++ = y;
				*p++ = z;
				*p++ = w;
				*p++ = x;	
				}
			else
				{
				*p++ = w;
				*p++ = x;	
				*p++ = y;
				*p++ = z;
				}
			p = *Relators[i];
			x = *p++;
			y = *p++;
			z = *p++;
			w = *p++;
			x -= 65;
			if(GL[x][0]) return(1);	/* Not a knot or link! */		
			}
		GL[x][0] ++;
		GL[x][2] ++;
		GL[x][1] = i;
		GL[x][3] = 100;			
		}	
	
	/******************************************************************************************
					Deal with those relators which are of the form AAaa or AaaA.
	******************************************************************************************/
	
	for(i = Mark2 + 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *p++;
		y = *p++;
		z = *p++;
		w = *p;	
		x -= 65;
		if(GL[x][3]) return(1); /* Not a knot or link! */
		if(GL[x][0]) return(1); /* Not a knot or link! */
		GL[x][0] ++;
		GL[x][2] ++;
		GL[x][1] = i;
		GL[x][3] = 100;			
		}

	/******************************************************************************************
						Check that each overpass has two ends listed.
	******************************************************************************************/
	
	for(i = 0; i < MAXNUMGENERATORS; i++)
		if(GL[i][3] && (GL[i][0] != 1 || GL[i][2] != 1)) return(1);
		
	/******************************************************************************************
				Check that each generator now appears on the generator list.
	******************************************************************************************/		
	
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		x = *(p+1);
		if(x < 'a')
			x -= 65;
		else
			x -= 97;
		if(GL[x][3] != 100) return(1); /* Not a knot or link! */
		}
		
	/******************************************************************************************
		Find each component of the link and get a meridian and longitude for each component.			
	******************************************************************************************/
	
	NumComponents = 0;
	while(1)
		{
		for(i = 0; i < MAXNUMGENERATORS; i++) if(GL[i][3] == 100) break;
		if(i == MAXNUMGENERATORS) break;
		NumComponents ++;
		GL[i][3] = NumComponents;
		ReallocateHandle((char **) DualRelators[NumComponents],2L);
		if((p = *DualRelators[NumComponents]) == NULL) return(1);
		*p++ = i + 'A';
		*p = EOS;
		ReallocateHandle((char **) OutRelators[NumComponents],2*MAXNUMGENERATORS + 1);
		if((q = *OutRelators[NumComponents]) == NULL) return(1);
		j = i;
		do
			{
			p = *Relators[GL[j][1]];
			x = *p++ - 'A';
			GL[x][3] = NumComponents;
			*q++ = *p++;
			j = *p - 97;
			}
		while(j != i);
		j = 0;
		r = *OutRelators[NumComponents];
		while(r < q)
			{
			x = *r++;
			if(x < 'a')
				{
				x -= 'A';
				if(GL[x][3] == NumComponents) j++;
				}
			else
				{
				x -= 'a';
				if(GL[x][3] == NumComponents) j--;
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
			*q++ = x;
			j--;
			}
		*q = EOS;
		r = *OutRelators[NumComponents];
		SetHandleSize((char **) OutRelators[NumComponents],q - r + 1);	
		}
	
	if(NumComponents == 0) return(1);
	
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
			return(0);
		default:
			SysBeep(5);
			goto GET_RESPONSE1;
		}
	
	if(NumComponents > 1)
	printf("\n\n		NOTE: TO SKIP FILLING A PARTICULAR BOUNDARY COMPONENT HIT 'return'.");
	else
	printf("\n\n		NOTE: TO SKIP DOING THIS FILLING HIT 'return'.");
	
	SavedPresentation = FALSE;
	r = (unsigned char *) NewPtr(200L);
	for(i = 1; i <= NumComponents; i++)
		{
		if(NumRelators >= MAXNUMRELATORS)
			{
			printf("\n\nThis presentation now has %d relators.",NumRelators);
			printf("\nThe program cannot accept more than %d relators.",NumRelators);
			if(i == 1)
				printf("\nSo no fillings are possible. Sorry!");
			else
				printf("\nSo no more fillings are possible. Sorry!");
			break;
			}
		printf("\n\n				Component %d:",i);
		printf("\nMeridian:  %s",*DualRelators[i]);
		printf("\nLongitude: %s",*OutRelators[i]);
		GET_COEFFICIENTS:
		printf("\n\nPlease specify the slope M/L of the surgery on this component or hit 'return'.");
		REGET_COEFFICIENT_M:
		printf("\n1) ENTER THE VALUE OF M AND HIT 'return'.    ");	
		ReadString((char *)r, GetPtrSize(r));
		p = r;
		x = *p++;
		if(x == EOS || x == '\n') continue;
		if(!isdigit(x) && x != '+' && x != '-')
			{
			SysBeep(5);
			printf("\nThe first character must be a digit, '+' or '-' sign!");
			goto REGET_COEFFICIENT_M;
			}
		while(x = *p++) if(!isdigit(x))
			{
			SysBeep(5);
			printf("\nEach character must be a digit or a leading '+' or '-' sign!");
			goto REGET_COEFFICIENT_M;
			}
		sscanf((char *)r,"%d",&M);
		REGET_COEFFICIENT_L:
		printf("2) ENTER THE VALUE OF L AND HIT 'return'.    ");
		ReadString((char *)r, GetPtrSize(r));
		p = r;
		x = *p++;
		if(!isdigit(x) && x != '+' && x != '-')
			{
			SysBeep(5);
			printf("\nThe first character must be a digit, '+' or '-' sign!\n");
			goto REGET_COEFFICIENT_L;
			}
		while(x = *p++) if(!isdigit(x)) 
			{
			SysBeep(5);
			printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
			goto REGET_COEFFICIENT_L;
			}
		sscanf((char *)r,"%d",&L);
		if(M == 0 && L == 0)
			{
			SysBeep(5);
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
		NRL = abs(M) + abs(L)*(GetHandleSize((char **) OutRelators[i]) - 1);
		if(NRL > MAXLENGTH)
			{
			SysBeep(5);
			printf("\n\nThese coefficients give a new relator which will probably be too long!");
			printf("\n\nPlease specify new surgery coefficients for this component.");
			goto GET_COEFFICIENTS;
			}
		if(SavedPresentation == FALSE)
			{
			SavedPresentation = TRUE;
			fprintf(myout,"\n\nTried Dehn-filling on: %s\n",PresName);
			Print_Relators(Relators,NumRelators,myout);
			}
		fprintf(myout,"\nTried %d/%d surgery with meridian: %s and longitude: %s.",M,L,
			*DualRelators[i],*OutRelators[i]);		
		NumRelators++;	
		ReallocateHandle((char **) Relators[NumRelators],NRL + 1);
		if((q = *Relators[NumRelators]) == NULL)
			{
			NumRelators--;
			SysBeep(5);
			printf("\n\nThese coefficients give a new relator which is too long!");
			printf("\n\nPlease specify new surgery coefficients for this component.");
			goto GET_COEFFICIENTS;
			}
		if(M < 0)
			{
			HLock((char **) DualRelators[i]);
			Inverse(*DualRelators[i]);
			HUnlock((char **) DualRelators[i]);
			M = -M;
			}	
		h = 0;
		do
			{
			if(h < M)
				{
				p = *DualRelators[i];
				while(*q++ = *p++) ;
				q--;
				}
			else
				{
				p = *OutRelators[i];
				while(*q++ = *p++) ;
				q--;
				}
			if(h >= L)
				h -= L;
			else
				h += M;	
			}
		while(h);
		q++;
		*q = EOS;	
		}
	DisposePtr((char *) r);	
	return(0);		
}

Try_Exponent_Surgery()
{
	register unsigned char	*p,
							*q,
							x,
							y;
							
	char					Change[125];
	
	unsigned char 			**Temp;
	
	int						e,
							i,
							j,
							k,
							P,
							Q,
							NewEXP[MAXNUMGENERATORS][3];
	
	unsigned int			gcd;
							
	unsigned long			length;
	
	long					Scratch;						
	
	unsigned int 			GCD(),
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
	Print_Relators(Relators,NumRelators,stdout);

	/******************************************************************************************
		Save a copy of the initial relators so that the user can perform additional exponent
		surgeries.
	******************************************************************************************/
	
	for(i = 1; i <= NumRelators; i++)
		{
		ReallocateHandle((char **) Exp_Surgery_Rel[i],GetHandleSize((char **) Relators[i]));				
		if((p = *Exp_Surgery_Rel[i]) == NULL)
			{
			Did_Exponent_Surgery = FALSE;
			break;
			}
		q = *Relators[i];
		while(*p++ = *q++) ;
		}
		
	NumExp_Sur_Rel = NumRelators;
	
	GET_NEW_EXPONENTS:	
	p = (unsigned char *) NewPtr(100L);
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
			ReadString((char *)p, GetPtrSize(p));
			q = p;			
			x = *q++;
			if(!isdigit(x) && x != '+' && x != '-')
				{
				SysBeep(5);
				printf("\nThe first character must be a digit, '+' or '-' sign!");
				goto GET_NEW_EXPONENTS1P;
				}
			while(x = *q++) if(!isdigit(x))
				{
				SysBeep(5);
				printf("\nEach character must be a digit or a leading '+' or '-' sign!");
				goto GET_NEW_EXPONENTS1P;
				}
			sscanf((char *) p,"%d",&P);
			GET_NEW_EXPONENTS1Q:
			printf("REPLACE Q WITH:     ");
			Q = EXP[i][1];
			ReadString((char *)p, GetPtrSize(p));
			q = p;
			x = *q++;
			if(!isdigit(x) && x != '+' && x != '-')
				{
				SysBeep(5);
				printf("\nThe first character must be a digit, '+' or '-' sign!\n");
				goto GET_NEW_EXPONENTS1Q;
				}
			while(x = *q++) if(!isdigit(x))
				{
				SysBeep(5);
				printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
				goto GET_NEW_EXPONENTS1Q;
				}		
			sscanf((char *) p,"%d",&Q);
			if(P == 0 && Q == 0)
				{
				SysBeep(5);
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
				ReadString((char *)p, GetPtrSize(p));
				q = p;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!");
					goto GET_NEW_EXPONENTS2P;
					}
				while(x = *q++) if(!isdigit(x))
					{
					SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!");
					goto GET_NEW_EXPONENTS2P;
					}				
				sscanf((char *) p,"%d",&P);
				GET_NEW_EXPONENTS2Q:
				printf("REPLACE Q WITH:     ");
				Q = EXP[i][2];
				ReadString((char *)p, GetPtrSize(p));
				q = p;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!\n");
					goto GET_NEW_EXPONENTS2Q;
					}
				while(x = *q++) if(!isdigit(x))
					{
					SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!\n");
					goto GET_NEW_EXPONENTS2Q;
					}								
				sscanf((char *) p,"%d",&Q);
				if(P == 0 && Q == 0)
					{
					SysBeep(5);
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
				ReadString((char *)p, GetPtrSize(p));
				q = p;
				x = *q++;
				if(!isdigit(x) && x != '+' && x != '-')
					{
					SysBeep(5);
					printf("\nThe first character must be a digit, '+' or '-' sign!");
					goto GET_NEW_EXPONENTS3P;
					}
				while(x = *q++) if(!isdigit(x))
					{
					SysBeep(5);
					printf("\nEach character must be a digit or a leading '+' or '-' sign!");
					goto GET_NEW_EXPONENTS3P;
					}							
				sscanf((char *) p,"%d",&P);
				NewEXP[i][0] = 0;
				NewEXP[i][1] = 0;
				NewEXP[i][2] = P;	
				}	
			}	
		}
	DisposePtr((char *) p);
	
	for(i = 1; i <= NumRelators; i++)
		{
		p = *Relators[i];
		ReallocateHandle((char **) OutRelators[i],MAXLENGTH + 2);
		if((q = *OutRelators[i]) == NULL)
			{
			printf("\n\nOut of memory! Can't do exponent surgery.");
			return(1);
			}
		length = 0L;
		while(x = *p)
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
					SysBeep(5);
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
							SysBeep(5);
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
					SysBeep(5);
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
							SysBeep(5);
							goto GET_RESPONSE2;
						}
					}				
				*q++ = x;
				p++;
				}	
			}
		*q = EOS;
		SetHandleSize((char **) OutRelators[i],length + 1);	
		}
		
	for(i = 1; i <= NumRelators; i++)
		{
		Temp = Relators[i];
		Relators[i] = OutRelators[i];
		OutRelators[i] = Temp;
		}
		
	/**************************************************************************************
		Echo the surgered relators to the output so we will have a copy of them. Then call
		Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
		which will serve as the initial presentation for the program. 
	**************************************************************************************/	
	
	ObscureCursor();
	fprintf(stdout,"\n\nThe surgered presentation is:\n");
	Print_Relators(Relators,NumRelators,stdout);
	fprintf(myout,"\n\nThe surgered presentation is:\n");
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
		return(1);
		}
	if(Scratch > OrigLength)
		{
		printf("and freely reduces to length %lu.",OrigLength);
		fprintf(myout,"and freely reduces to length %lu.",OrigLength);
		Scratch = OrigLength;
		}
	else
		{
		printf("and is freely reduced.");
		fprintf(myout,"and is freely reduced.");
		}	
	if(Rewrite_Input())
		{
		printf("\n\nThere must be at least one non-empty relator!!");
		SysBeep(5);
		return(1);
		}
	Length = OrigLength;	
									
	/**************************************************************************************
			Call Canonical_Rewrite() to rewrite the presentation in canonical form.
	**************************************************************************************/
			
	Canonical_Rewrite(Relators,FALSE,FALSE);					
	printf("\n\nNumRelators = %d, NumGenerators = %d\n",NumRelators, NumGenerators);
	
	printf("\n\nThe rewritten initial presentation is:\n");
	Print_Relators(Relators,NumRelators,stdout);
	fprintf(myout,"\n\nThe rewritten initial presentation is:\n");
	Print_Relators(Relators,NumRelators,myout);
	
	/**************************************************************************************
		Save a copy of the initial set of relators in case we want to refer to them later.
	**************************************************************************************/
		
	CopyNumRelators 	= NumRelators;
	CopyNumGenerators 	= NumGenerators;
	for(i = 1; i <= NumRelators; i++)
		{
		ReallocateHandle((char **) Copy_Of_Input[i],GetHandleSize((char **) Relators[i]));				
		if((p = *Copy_Of_Input[i]) == NULL)
			{
			printf("\n\n	Memory Error. Sorry!");
			return(1);
			}
		q = *Relators[i];
		while(*p++ = *q++) ;
		}
		
	return(0);			
}

Try_Cutting_Disk_Surgery()
{
	register unsigned char	*p,
							*q,
							x;
	
	unsigned char 			**Temp;
	
	int						i,
							j;
	
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
		ReallocateHandle((char **) CD_Surgery_Rel[i],GetHandleSize((char **) Relators[i]));				
		if((p = *CD_Surgery_Rel[i]) == NULL)
			{
			Did_Cutting_Disk_Surgery = FALSE;
			printf("\n\nMemory error. Sorry!");
			return(1);
			}
		q = *Relators[i];
		while(*p++ = *q++) ;
		}
		
	NumCuttingDiskSurgeryRel = NumRelators;
	
	printf("\n\n        Cutting disk surgery gives a way to obtain different Heegaard diagrams");
	printf("\n        by changing the way in which the ends of edges meeting pairs of inverse");
	printf("\n        vertices are identified. By choosing values for the offsets, you can");
	printf("\n        obtain the diagrams of all 3-manifolds which have Heegaard diagrams with");
	printf("\n        Whitehead graphs isomorphic to the Whitehead graph of the initial diagram.");

	printf("\n\nThe presentation is currently:\n");
	Print_Relators(Relators,NumRelators,stdout);

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
							
	p = (unsigned char *) NewPtr(100L);	
	for(i = 0; i < Vertices; i += 2)
		{
		printf("Generator %2d, Valence %5u, Current offset %5u. Change offset to ?    ",
			i/2 + 1,V[i],OSA[i] % V[i]);
		GET_NEW_OFFSET:
		ReadString((char *)p, GetPtrSize(p));
		q = p;
		while(x = *q++) if(!isdigit(x))
			{
			SysBeep(5);
			goto GET_NEW_OFFSET;
			}
		OSB[i] = OSA[i];
		OSB[i+1] = OSA[i+1];	
		sscanf((char *) p,"%u",&NOSA);
		NOSA = NOSA % V[i];
		NOSA += V[i];	
		OSA[i] = OSA[i+1] = NOSA;				
		}
	DisposePtr((char *) p);
	
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
				SysBeep(5);
				goto GET_RESPONSE1;		
			}
		}
		
	for(i = 1; i <= NumRelators; i++)
		{
		Temp = Relators[i];
		Relators[i] = OutRelators[i];
		OutRelators[i] = Temp;
		}
		
	/**************************************************************************************
		Echo the surgered relators to the output so we will have a copy of them. Then call
		Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
		which will serve as the initial presentation for the program. 
	**************************************************************************************/	
	
	ObscureCursor();
	fprintf(stdout,"\n\nThe surgered presentation is:\n");
	Print_Relators(Relators,NumRelators,stdout);			
	fprintf(myout,"\n\nThe surgered presentation is:\n");
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
				SysBeep(5);
				goto GET_RESPONSE2;		
			}
		}
	if(Scratch > OrigLength)
		{
		printf("and freely reduces to length %lu.",OrigLength);
		fprintf(myout,"and freely reduces to length %lu.",OrigLength);
		Scratch = OrigLength;
		}
	else
		{
		printf("and is freely reduced.");
		fprintf(myout,"and is freely reduced.");
		}	
	if(Rewrite_Input())
		{
		printf("\n\nThere must be at least one non-empty relator!!");
		SysBeep(5);
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
				SysBeep(5);
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
	Print_Relators(Relators,NumRelators,stdout);

	printf("\n\nUSE THIS PRESENTATION ? HIT 'y' OR 'n'.    ");
	GET_RESPONSE4:
	switch(WaitkbHit())
		{
		case 'y':
			break;
		case 'n':
			goto RETRY;
		default:
			SysBeep(5);
			goto GET_RESPONSE4;		
		}
		
	fprintf(myout,"\n\nThe rewritten initial presentation is:\n");
	Print_Relators(Relators,NumRelators,myout);
	
	/**************************************************************************************
		Save a copy of the initial set of relators in case we want to refer to them later.
	**************************************************************************************/
		
	CopyNumRelators 	= NumRelators;
	CopyNumGenerators 	= NumGenerators;
	for(i = 1; i <= NumRelators; i++)
		{
		ReallocateHandle((char **) Copy_Of_Input[i],GetHandleSize((char **) Relators[i]));				
		if((p = *Copy_Of_Input[i]) == NULL)
			{
			printf("\n\n	Memory Error. Sorry!");
			Did_Cutting_Disk_Surgery = FALSE;
			return(1);
			}
		q = *Relators[i];
		while(*p++ = *q++) ;
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
		ReallocateHandle((char **) DualRelators[i],max);
		if((p = *DualRelators[i]) == NULL)
			{
			printf("\n\nMemory error. Sorry!");
			return(TOO_LONG);
			}
		for(e = V[(i-1) << 1]; e > 0; e--) *p++ = '@';
		*p = EOS;
		}	

	max = MAXLENGTH;
	SNumRelators = NumRelators;
	NumRelators = 0;	
		
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
					SysBeep(5);
					printf("\n\nThis presentation now has %d relators.",NumRelators);
					printf("\nThe program cannot accept more than %d relators. Sorry!",MAXNUMRELATORS);
					NumRelators = SNumRelators;
					return(TOO_LONG);
					}
				ReallocateHandle((char **) OutRelators[i+1],max + 2);
				if(*OutRelators[i+1] == NULL)
					{
					SysBeep(5);
					printf("\n\nMemory error. Sorry!");
					NumRelators = SNumRelators;
					return(TOO_LONG);
					}
				r = *OutRelators[i+1];
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
						SysBeep(5);
						printf("\n\nRelator %d of this presentation is longer than %d.",
							NumRelators, MAXLENGTH);
						printf("\nThis is too long for the program to handle. Sorry!");
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
				SetHandleSize((char **) OutRelators[i+1],length + 1);
				i++;
				}
			p++;
			}	
		}
	return(NO_ERROR);	
}

