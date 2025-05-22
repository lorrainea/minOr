#include <vector>
#include <limits>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <complex>
#include <divsufsort64.h>
#include <cstring>
#include <sdsl/bit_vectors.hpp>                                
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp> 
#include <set>


using namespace std;
using namespace sdsl;
using Complex = complex<double>;

#ifdef _USE_64
#include <divsufsort64.h>
typedef int64_t INT;                                     
#endif

#ifdef _USE_32
#include <divsufsort.h>
typedef int32_t INT;                                      	
#endif


INT sum_Xi(vector<INT> Xi, INT i, INT n, INT w)
{
    INT sum = 0;
    for (INT j = 0; j < n; j++)
    {
        sum += Xi[ (i+j) % w];
    }
    return sum;
}



/* Kasai et al algorithm for O(n)-time LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}


INT ordering_t( unsigned char * seq, INT w, INT t, INT n, INT * t_order )
{
	unordered_set<INT> minimizers;
	 	
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT * rank;
	
	rank = 	( INT * ) malloc( ( n ) * sizeof( INT ) );
	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	
	/* Compute suffix array */
	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	

	for ( INT i = 0; i < n; i ++ )
	{
	        invSA [SA[i]] = i;        
	}

	/* Compute the LCP array */
	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
	{
	        fprintf(stderr, " Error: LCP computation failed.\n" );
	        exit( EXIT_FAILURE );
	}


	INT rank_count =  0;	
	rank[SA[0]] = rank_count;
	
	
	/* Compute the ranks */
	for(INT j =1; j<n; j++)
	{
		if( LCP[j] >= t )
		{
			
			rank[SA[j]] = rank_count;
		}
		else 
		{
			rank_count++;
			rank[SA[j]] = rank_count;
		}
		
	}	
	
	INT window = w + t - 1;
	
	for(INT a = 0; a<=n-window; a++)
	{
		INT min = n;
		INT min_i = 0;
		for(INT b = a; b<a+window; b++)
		{
			if( rank[b] < min )
			{
				min = rank[b];
				min_i = b;
			}
		}
		
		t_order[a] = min_i;
	}
	
return 0;
}

/* Computes the minimizers of a string of length n in O(n) time */
INT standard(  INT * k_order, INT n, INT window, INT k, unordered_set<INT> &minimizers  )	
{
	
	deque<pair<INT,INT>> min_rank = {};
	INT w = window + k - 1;
	
	// minimum rank in first window
   	for (INT j = 0; j < w - 1 ; j++) 
   	{
 		while ( !min_rank.empty() && k_order[j] < min_rank.back().first )
 			min_rank.pop_back();
				
		min_rank.push_back(std::make_pair(k_order[j], j));
    	}
    	
	minimizers.insert( min_rank.at(0).second );

	
	// minium rank in remaining windows
	for( INT i = 0; i<=n-w; i++ )
	{
		while (!min_rank.empty() && min_rank.back().first > k_order[i+w-1])
			min_rank.pop_back();
	
		min_rank.push_back(std::make_pair(k_order[i+w-1], i+w-1));
		
	
		while( !min_rank.empty() && min_rank.front().second < i )
		{
			min_rank.pop_front();
		}	
		
			
		minimizers.insert( min_rank.at(0).second );
	}
	
	return 0;
}

INT miniception( unsigned char * seq, INT * k_order, INT n, INT w, INT k, unordered_set<INT> &minimizers )
{

	
	INT * t_order = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	INT t = 8;
	
	ordering_t( seq, w, t, n, t_order );
	
	INT w_prime = k-t;
	INT p = 0;
	for(INT j = 0; j<=n-w; j++)
	{
		std::pair<INT, INT> o_min = {1, std::numeric_limits<int>::max()};
		
		for(INT i = 0; i<w+k-t; i++)
		{
			INT p_prime = t_order[j+i];
			
			std::pair<INT, INT> o;
			
			if( p_prime == 0 or p_prime == w_prime )
				o = {0, k_order[j+i]};
			else o = {1, k_order[j+i]};
			
			if( o < o_min )
			{
				o_min = o;
				p = i;
			}

		}
		minimizers.insert(j+p);
	}

return 0;
}

INT rotational( unsigned char * seq, INT * k_order, INT len, INT w, INT k, unordered_set<INT> &minimizers )
{

    	INT n = k/w;
    	
    	vector<INT> sum_0(n);
    	vector<INT> sum_j(n);
    	vector<INT> Xi(k);
	memcpy(Xi.data(), seq, k * sizeof(INT));
    	
    	sum_0[0] = sum_Xi(Xi, 0, n, k);
    	sum_j[0] = sum_Xi(Xi, 0, n, k);
    	
	for(INT i = 1; i<len; i++)
	{
	
		Xi.erase(Xi.begin()); 
		Xi.push_back( seq[i+k-1] );
		
		sum_0.push_back(sum_Xi(Xi, 0, n, k));
    		sum_j.push_back(sum_Xi(Xi, i, n, k));
	}
	
    	
    			
	
	for (INT j = 0; j <= len - w ; j++) 
	{
		INT o_min = numeric_limits<int>::max();
		INT p = 0;

		for (INT i = 0; i<w; i++)
		{
			
			INT v0 = sum_0[j+i];
			bool valid = true;

			for (INT z = 1; z < w; z++)
			{
			    INT vj = sum_j[j+z];
			    if (vj > v0+4) 
			    {
			    	valid == false;
			    	break;
			    }
			}
			
			if( valid == false )
				continue;
				 
			INT o = k_order[j+i];
			
			if (o < o_min)
			{
				o_min = o;
				p = i;
			}
		}
		
		minimizers.insert(j+p);
	}



return 0;
}

/* Computes minimizers using mod-sampling scheme */
INT mod_sampling( unsigned char * seq, INT * k_order, INT n, INT w, INT k, unordered_set<INT> &minimizers )	
{
	
	INT * t_order = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	INT t = k/2;
	
	ordering_t( seq, w, t, n, t_order );
	
	
	INT window_size = w + k - t;
	deque<INT> dq;

	for(INT j = 0; j<window_size-1; j++)
		dq.push_back( j );
		
	for (INT j = 0; j <= n - window_size; j++)
	{
	
		while (!dq.empty() && dq.front() < j)
			dq.pop_front();
  
		INT new_idx = j + window_size - 1;
		while (!dq.empty() && t_order[dq.back()] > t_order[new_idx])
			dq.pop_back();
        
		dq.push_back(new_idx);

		INT x = dq.front() - j;	
		INT p = x % w;
		
		minimizers.insert(j + p);
		
	}
	
	return 0;
}


/* Computes orderings based on k */
INT compute_minimizers(  unsigned char * seq, INT * k_order, INT k, INT w, unordered_set<INT> &minimizers, char m, INT a )
{
	INT n = strlen ( (char*) seq );
	

		
	if( m == 's')
		standard( k_order, n, w, k, minimizers );		
	else if( m == 'm' )
		mod_sampling( seq, k_order, n, w, k, minimizers );
	else if( m == 'r' )
		rotational( seq, k_order, n, w, k, minimizers );	
	else if( m == 'c')
		miniception( seq, k_order, n, w, k, minimizers );
	
	
return 0;
	
}


int main(int argc, char **argv)
{

	if( argc < 6 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./minOr <text> <ordering> <w> <k> <mode>\n";
 		exit(-1);
 	}
 	
 	ifstream is;
 	is.open(argv[1], ios::in | ios::binary);

	std::string str3(argv[3]);
	std::string str4(argv[4]);
	std::string str5(argv[5]);

 	INT w;
 	std::stringstream(str3)>>w;
 	
 	INT k;
 	std::stringstream(str4)>>k;
 	
 	char m;
 	std::stringstream(str5)>>m;
 	
	
 	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
	INT text_file_size = in_file.tellg();
  	
  	
   	unsigned char * text = ( unsigned char * ) malloc (  ( text_file_size + 1 ) * sizeof ( unsigned char ) );	
   	
  	char c = 0;
  	INT text_size = 0;
	for (INT i = 0; i < text_file_size; i++)
	{	
		is.read(reinterpret_cast<char*>(&c), 1);
		if( (unsigned char) c == '\n')
			continue;
		else
		{
			text[text_size] = (unsigned char) c ;
			text_size++;
		}
		
	}
	in_file.close();
	text[text_size] = '\0';
	
	if( text_size < w )
	{
		fprintf( stderr, " Error: Number of windows (w) cannot be larger than sequence length!\n");
		return ( 1 );
	}

	if( m != 'm' and m != 's' and m != 'r' and m != 'c' )
	{
		fprintf( stderr, " Error: Mode <mode> must be s for standard minimizer computation, m for mod-sampling, r for rotational minimizers and c for miniception.\n");
		return ( 1 );
	}
	
	INT n = strlen( (char*) text );
	
	INT * k_order = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	
	INT num;
	
	ifstream fin(argv[2]);
	INT a = 0;
	while( fin >> num )
	{
		k_order[a] = num;
		a++;
	}	
	
		
	unordered_set<INT> minimizers;

	compute_minimizers( text, k_order, k, w, minimizers, m, a);
	
	cout<<"Minimizers: "<<minimizers.size()<<endl;
	
	free( text );
	free( k_order );
			
	return 0;
}
