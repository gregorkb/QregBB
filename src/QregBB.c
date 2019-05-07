#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void ETBBgetweights( 	int *n,
			int *l,
			int *B,
			int *blkpts,
                   	double *f_tilde, 
                  	int *seed,
			double *weights,
			double *c)
{


	int i,t,j,*intp,k,blks,max;
	double *ip,tmp,accum,scale,weightnorm1;
    
    	GetRNGstate();
    
	max = *n - *l + 1;
	blks = floor( *n / *l );

	ip = &weights[0];
	weightnorm1 = 0;

	for(i=0;i<*l;i++)
	{
		tmp = ((i+1) - .5)/ *l;

		if(tmp <= *c) *ip = tmp / *c;
		else if(tmp >= *c && tmp < 1 - *c) *ip = 1;
		else *ip = (1-tmp)/ *c;
		
		weightnorm1 += *ip;
		ip++;
	}

	scale = 1 / ( blks * weightnorm1 );


    for(i=0;i<*B;i++)
    { 

	ip = &f_tilde[i * (*n)]; // initialize pointer appropriate position in f_tilde vector

	intp = &blkpts[0];
        
	for(j=0;j<blks;j++) //  choose blocks
        {
            	*intp = floor(unif_rand()*max) ;            
		intp++;
        }

	for(t=0;t<*n;t++) // go through each observation to define f_tilde[t] values
	{
		accum = 0;
		for(j=0;j<blks;j++) // over all the blocks...
		{		

			for(k=0;k<fmin(*l,t+1);k++) // ...get # times each obs is chosen; weight it according to its position in each block.
			{
				if(blkpts[j]==t-k) accum += weights[k];  
		
			}
		}

		*ip = accum * scale  ; // save f_tilde value
		ip++;	// advance pointer to next position in f_tilde vector
	}
    }

    PutRNGstate();
}


void MBBgetweights( 	int *n,
			int *l,
			int *B,
			int *blkpts,
                   	double *f_tilde, 
                  	int *seed)
{


	int i,t,j,*intp,k,blks,max;
	double *ip,accum, doublen;
    
    	GetRNGstate();
    
	max = *n - *l + 1;
	blks = floor( *n / *l);



    for(i=0;i<*B;i++)
    { 

	ip = &f_tilde[i * (*n)]; // initialize pointer appropriate position in f_tilde vector

	intp = &blkpts[0];
        
	for(j=0;j<blks;j++) //  choose blocks
        {
            	*intp = floor(unif_rand()*max) ;            
		intp++;
        }

	doublen = *n;

	for(t=0;t<*n;t++) // go through each observation to define f_tilde[t] values
	{
		accum = 0;
		for(j=0;j<blks;j++) // over all the blocks...
		{		

			for(k=0;k<fmin(*l,t+1);k++) // ...get # times each obs is chosen; weight it according to its position in each block.
			{
				if(blkpts[j]==t-k) accum += 1 / doublen;  
		
			}
		}

		*ip = accum ; // save f_tilde value
		ip++;	// advance pointer to next position in f_tilde vector
	}
    }

    PutRNGstate();
}



void AllBBgetweights(	int *n,
                   	int *l,
                  	int *B,
			int *blkpts,
                   	double *f_MBB,
			double *f_ETBB, 
                  	int *seed,
			double *weights,
			double *c,
			double *m_l)
{
   
	int 	i,
		t,
		j,
		*ip,
		k,
		blks,
		max;
	double 	*dp1,
		*dp2,
		tmp,
		accumMBB,
		accumETBB,
		scale,
		weightnorm1,
		weightnorm2;
    
	GetRNGstate();
    
	max = *n - *l + 1;
	blks = floor( *n / *l );

	dp1 = &weights[0];
	weightnorm1 = 0;
	weightnorm2 = 0;

	for(i=0;i<*l;i++) // compute weights for taper window as well as their sum and the sum of their squares
	{
		tmp = ((i+1) - .5)/ *l;

		if(tmp <= *c) *dp1 = tmp / *c;
		else if(tmp >= *c && tmp < 1 - *c) *dp1 = 1;
		else *dp1 = (1-tmp)/ *c;
		
		weightnorm1 += *dp1;
		weightnorm2 += (*dp1)*(*dp1);
		dp1++;
	}

	scale = 1 / ( blks * weightnorm1 );

	*m_l = weightnorm1 * weightnorm1 / ( (*l) * weightnorm2);
	
// BEGIN BOOTSTRAP LOOP //

    for(i=0;i<*B;i++)
    { 

	dp1 = &f_MBB[i * (*n)]; // initialize pointers to appropriate positions in f_MBB and f_ETBB vectors
	dp2 = &f_ETBB[i * (*n)];

	ip = &blkpts[i * blks];
        
	for(j=0;j<blks;j++) //  choose blocks
        {
            	*ip = floor(unif_rand()*max) ;            
		ip++;
        }

	for(t=0;t<*n;t++) // go through each observation to define f_tilde[t] values
	{
		accumMBB = 0 ;
		accumETBB = 0 ;
		for(j=0;j<blks;j++) // over all the blocks...
		{		

			for(k=0;k<fmin(*l,t+1);k++) // ...get # times each obs is chosen; weight it according to its position in each block.
			{
				if(blkpts[(i * blks) + j]==t-k)
				{
				 	accumMBB += 1;  
					accumETBB += weights[k];  
				}
			}
		}

		*dp1 = accumMBB / ( blks * (*l) ) ; // save f_MBB value
		*dp2 = accumETBB * scale  ; // save f_ETBB value
		
		dp1++; // advance pointers to next positions in f_MBB and f_ETBB vectors
		dp2++;
	}

    }

    PutRNGstate();
}



//////////////////////////////////////////////////////////////////

