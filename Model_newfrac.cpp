#include <stdio.h>
#include <vector>
#include <math.h>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <Eigen/Sparse>
//---
using namespace Eigen;
typedef SparseMatrix<double> SparseMatrixXd;
using namespace std;
//--- utility funcs
void initialize( vector<int> *svect, const unsigned int nsize )
{
	for( unsigned int i = 0; i < nsize; i++ )
		( *svect )[ i ] = 0;	
}
void initialize_double( vector<double> *svect, const unsigned int nsize )
{
    for( unsigned int i = 0; i < nsize; i++ )
        ( *svect )[ i ] = 0.0;
}
//---
bool isFailed( vector<double> *tau, vector<double> *S, vector<double> *rho, unsigned int nsize )
{
	for( unsigned int i = 0; i < nsize; i++ )
	{
		if( (*tau)[ i ] + (*rho)[ i ]  > (*S)[ i ] )
			return true;
	}
	return false;
}
//---
inline void getIndex( const unsigned int n, const unsigned int i, const unsigned int j, unsigned int *index )
{
	*index = i * n + j;
}
//---
vector <double> ReadInDoubles( const string& filename )
  {
  vector <double> result;
  double          d;
  std::ifstream        inf( filename.c_str() );
  while (inf >> d)
    {
    result.push_back( d );
    }
  return result;
  }
//---
//--- main program
//---
int main()
{
	const unsigned int ngrid = 501;
	const double L = 1.0; //--- Distance between the center of a site to it's neighbor
    const double A = 1.0; //--- Area of interface between two neighbors
    const double V = 1.0; //--- Volume of a site
    const double K = 1e-16; //--- Characteristic permeability (unit : m^2)
    const double alpha = 1e-8; //--- effective compressibility (unit :  1/Pa )
	const double Viscosity = 1e-4; //---  unit( Pa.sec ) water viscosity in 100 Celsius degrees
	const double dt =  16.0 ; //--- time discretization (dimensionless)
	FILE *fptr;
	//--- solid part
	const unsigned int index_injection = ( ngrid / 2 ) * ngrid + ngrid / 2;
	const unsigned int nsize = ngrid * ngrid;
	//--- declare vectors
    VectorXd force_shear( nsize ); //--- with non zero elements associated with unstable sites
	VectorXd d_tau( nsize ); //--- change in shear_stress
	vector<int> failure( nsize ); //--- state variable 0: stable 1: unstable
	vector<double> tau( nsize ); //--- shear stress
    vector<double> S( nsize ); //--- friction and cohesion are all absorbed in this quantitiy
	vector<double> rho( nsize ); //--- pressure
	vector<double> rho_updated( nsize ); //--- updated pressure
	vector<double> permeability_x; //--- Contains permeability of bonds on the X axis
    vector<double> permeability_y; //--- Contains permeability of bonds on the Y axis
    vector<double> porosity_x; //--- Contains porosity of bonds on the X axis
    vector<double> porosity_y; //--- Contains porosity of bonds on the Y axis
    vector<int> K_x_change( nsize ); //--- state variable of 0: it did not get increased   1: it's already increased
    vector<int> K_y_change( nsize ); //--- state variable of 0: it did not get increased   1: it's already increased
    //---
    permeability_x = ReadInDoubles( "K_x.txt" );
    permeability_y = ReadInDoubles( "K_y.txt" );
    porosity_x = ReadInDoubles( "Phi_x.txt" );
    porosity_y = ReadInDoubles( "Phi_y.txt" );
    //--- parse Shear stress and Yield function
    tau = ReadInDoubles( "shear.txt" ); // uniform distribution  (0.0,1.0)
    S = ReadInDoubles( "yield.txt" ); // uniform distribution (S_min , S_max)
    //---
    FILE *f_time_series = fopen( "avalanche.txt", "w" );
    fprintf( f_time_series, "#ITIME\tDURATION\tArea\tSIZE\tEnergy\n" );
    
    FILE *f_pressure = fopen("pressure.txt","w");
    fprintf( f_pressure, "#ITIME\tMSD\n" );
    
    FILE *f_fracture = fopen("fracture.txt","w");
    fprintf( f_fracture, "#ITIME\tRad\tMSD\tRMS\n" );
    //--
    //--
    initialize( &K_x_change, nsize ); //-- this will set all the values of this equall to zero
    initialize( &K_y_change, nsize ); //-- this will set all the values of this equall to zero
    //--
	const double dimensionless_number = (A * dt * K )/( V * L * alpha * Viscosity );
	printf( "Dimensionless_number = %e \n", dimensionless_number);
	//--- sparse matrix for redistributing tau
	SparseMatrixXd smat( nsize, nsize );
	for( int  i = 0; i < ngrid; i++ )
	{
		for(  int j = 0; j < ngrid; j++ )
		{
			int index_center = i * ngrid + j;
			assert( index_center < nsize );
			smat.coeffRef( index_center, index_center ) = -1.0;
			int index_right = i * ngrid + ( j + 1 ) % ngrid;
			assert( index_right < nsize );
			smat.coeffRef( index_right, index_center ) = 0.25;
			int index_left = i * ngrid + ( j - 1  < 0 ? j - 1 + ngrid : j - 1 );
			assert( index_left < nsize );
			smat.coeffRef( index_left, index_center ) = 0.25;
			int index_top = ( ( i + 1 ) % ngrid ) * ngrid + j;
			assert( index_top < nsize );
			smat.coeffRef( index_top, index_center ) = 0.25;
			int index_bot = ( i - 1  < 0 ? i - 1 + ngrid : i - 1 ) * ngrid + j;
			assert( index_bot < nsize );
			smat.coeffRef( index_bot, index_center ) = 0.25;
		}
	}
    unsigned int Number_ofevents = 0;
    unsigned int itime = 0;
    while (Number_ofevents != 1000) {
        if (itime % 1000 == 0) {
            printf( "itime = %d\n", itime);
        }
        //--- BC
		rho[ index_injection ] = 1.0;
		rho_updated[ index_injection ] = 1.0;
		//---
		unsigned int index = 0, index_W, index_E, index_N, index_S;
		for( int node_i = 0; node_i < ngrid; node_i++ )
		{
			for( int node_j = 0; node_j < ngrid; node_j++ )
			{
                
                getIndex( ngrid, node_i, node_j - 1  < 0 ? node_j - 1 + ngrid : node_j - 1 , &index_W );
                getIndex( ngrid, node_i, ( node_j + 1 ) % ngrid, &index_E );
                getIndex( ngrid, ( node_i + 1 ) % ngrid, node_j, &index_S );
                getIndex( ngrid, ( node_i - 1  < 0 ? node_i - 1  + ngrid : node_i - 1 ), node_j, &index_N );
                
				if( index == index_injection )
				{
                    permeability_y[ index ]   = 20.0 ;
                    permeability_x[ index ]   = 20.0 ;
                    permeability_y[ index_N ] = 20.0 ;
                    permeability_x[ index_W ] = 20.0 ;
					index++;
					continue;
				}
				
                double p_E = ( rho[ index_E ] - rho[ index ] );
                double p_N = ( rho[ index_N ] - rho[ index ] );
                double p_S = ( rho[ index_S ] - rho[ index ] );
                double p_W = ( rho[ index_W ] - rho[ index ] );
                //-- Ratio of permeability to porosity of the bonds attached to the site ij.
                //-- Keep in mind that for the eastern and southern neighbours, the index of bond is the same as the index of the site ij.
                double K_E = (permeability_x[ index ] / porosity_x[ index ]) ;
                double K_N = (permeability_y[ index_N ] / porosity_y[ index_N ]) ;
                double K_S = (permeability_y[ index ] / porosity_y[ index ]) ;
                double K_W = (permeability_x[ index_W ] / porosity_x[ index_W ]) ;
                //---
                double q_E =  K_E * p_E;
                double q_N =  K_N * p_N;
                double q_S =  K_S * p_S;
                double q_W =  K_W * p_W;
				//---
				rho_updated[ index ] = rho[ index ] + (dimensionless_number) * (q_E + q_N + q_S + q_W);
				//---
				index++;
        
			}
		}
        //--- initialize
        for( unsigned int index = 0; index < nsize; index++ )
            rho[ index ] = rho_updated[ index ] ;
        //--- yield function
        unsigned int time_internal = 0, avalanche_area = 0, avalanche_size = 0;
        double  energy_rel = 0.0;
        initialize( &failure, nsize );
        while( isFailed( &tau, &S, &rho, nsize) )
        {
            //printf( "time_internal = %d\n", time_internal);
            unsigned int  i = 0 , index_w, index_e, index_n, index_s;
            for( int Node_i = 0; Node_i < ngrid; Node_i++ )
            {
                for( int Node_j = 0; Node_j < ngrid; Node_j++ )
                {
                    if( tau[ i ] + rho[ i ] > S[ i ] ){
                        failure[ i ] = 1;
                        avalanche_size += 1;
                        energy_rel += tau[i];
                        force_shear[ i ] = tau[ i ]; //--- effective force;
                        //-- if the site fail, just for one time, the permeability corresponding to that site (four bonds) would increase by a factor 200
                        getIndex( ngrid, Node_i, Node_j - 1  < 0 ? Node_j - 1 + ngrid : Node_j - 1 , &index_w );
                        getIndex( ngrid, Node_i, ( Node_j + 1 ) % ngrid, &index_e );
                        getIndex( ngrid, ( Node_i + 1 ) % ngrid, Node_j, &index_s );
                        getIndex( ngrid, ( Node_i - 1  < 0 ? Node_i - 1  + ngrid : Node_i - 1 ), Node_j, &index_n );
                        //--- Basically the thing that is happening is that each site, that is going to fail will look at
                        //--- at all of it's neighbours which are failing as well, if the permeability of the bond between them
                        //-- did not increase, the permeability of the bond between those two sites will be multipled by 100
                        if ( tau[ index_n ] + rho[ index_n ] > S[ index_n ]) { // if the north neighbor fails as well, permeability of the bond
                            if (K_y_change[ index_n ] == 0){                    // between them would increase
                                permeability_y[ index_n ] = 1 * permeability_y[ index_n ] ;
                                K_y_change[ index_n ] = 1;
                            }
                        }
                        if (tau[ index_e ] + rho[ index_e ] > S[ index_e ]) { //--for the eastern neighbour, the index of bond is the same as the index of the i site
                            if (K_x_change[ i ] == 0 ){
                                permeability_x[ i ] = 1 * permeability_x[ i ] ;
                                K_x_change[ i ] = 1;
                            }
                        }
                        if ( tau[ index_s ] + rho[ index_s ] > S[ index_s ]) { //-- Also for the southern neighbour, the index of bond is the same as the index of the i site
                            if ( K_y_change[ i ] == 0 ){
                                permeability_y[ i ] = 1 * permeability_y[ i ] ;
                                K_y_change[ i ] = 1;
                            }
                        }
                        if ( tau[ index_w ] + rho[ index_w ] > S[ index_w ]) {
                            if (K_x_change[ index_w ] == 0){
                                permeability_x[ index_w ] = 1 * permeability_x[ index_w ] ;
                                K_x_change[ index_w ] = 1;
                            }
                        }
                    }
                    else{
                        force_shear[ i ] = 0.0;
                    }
                i++;
              }
            }
            d_tau = smat * force_shear;
            for( unsigned int i = 0; i < nsize; i++ ){
                tau[ i ] = tau[ i ] + d_tau[ i ] ;
                }
            time_internal++;
        } //--- end of while loop
        
        // if an event happened, add an event number to the number of events
        if (avalanche_size != 0) {
            Number_ofevents += 1 ;
            printf( "No.event = %d\n", Number_ofevents );
        }
        //---
        //--- output
        for( unsigned int i = 0; i < nsize; i++ )
            avalanche_area += ( failure[ i ] == 1 );
        
        if ( avalanche_area > 0) {
            fprintf( f_time_series, "%d\t%d\t%d\t%d\t%f\n", itime, time_internal, avalanche_area, avalanche_size, energy_rel );
            fflush( f_time_series );
        }
        //--- /Distance to event centroid, MSD(mean squared distance), RMS(root mean squared) of events
        if( avalanche_area > 0 ){
            double X_com = 0.0 ;
            double Y_com = 0.0 ;
            double Number = 0.0 ;
            double r3 = 0.0 ;
            double denominatorR = 0.0 ;
            unsigned int IND = 0;
            for( int xx = 0; xx < ngrid; xx++ )
            {
                for( int yy = 0; yy < ngrid; yy++ )
                {
                    if (failure[ IND ] == 1 ){
                        int XXX = ((ngrid - 1 ) / 2 ) - (xx);
                        int YYY = ((ngrid - 1 ) / 2 ) - (yy);
                        X_com += XXX;
                        Y_com += YYY;
                        float R = sqrt(XXX * XXX + YYY * YYY);
                        r3 += R * R * R ;
                        denominatorR += R ;
                        Number += 1.0;
                    }
                    IND++;
                }
            }
            fprintf( f_fracture, "%d\t%f\t%e\t%e\n", itime,sqrt((X_com * X_com)/(Number * Number)+(Y_com * Y_com)/(Number * Number)),r3/denominatorR,sqrt(r3/denominatorR));
            fflush( f_fracture );
        }
        //--- Mean-squared-Displacement for pressure profile
        if( avalanche_area > 0  ){
            double r3p = 0.0 ;
            double denominator = 0.0 ;
            unsigned int index =0;
            for( int XX = 0; XX < ngrid; XX++ )
            {
                for( int YY = 0; YY < ngrid; YY++ )
                {
                    if( rho[ index ]  > 0.0 )
                    {
                        int X = ( (ngrid - 1 ) / 2) - (XX);
                        int Y = ( (ngrid - 1 ) / 2) - (YY);
                        double r = sqrt(X * X + Y * Y);
                        r3p += r * r * r * rho[ index];
                        denominator += r * rho[ index];
                    }
                    index++;
                }
            }
            fprintf( f_pressure, "%d\t%e\n",itime, r3p / denominator );
            fflush( f_pressure );
        }
        //--- plot fluid pressure
        if( avalanche_area > 0  ){
            fptr = fopen( "vmd.xyz", "a" );
            fprintf( fptr, "%d\n", nsize );
            fprintf( fptr, "ITIME=%d\n", itime );
            for( unsigned int index = 0; index < nsize; index++ )
                
                fprintf( fptr, "%d\t%d\t%d\t%e\n", index, index / ngrid , index % ngrid, fabs( rho[ index ] ) > 0.0 ? rho[ index ] : 0.0 );
            fclose( fptr );
        }
        //--- plot fractured sites
        if( avalanche_area > 0 ){
            fptr = fopen( "fracture.xyz", "a" );
            if ( avalanche_area == 0 ) {
                fprintf( fptr, "%d\n", avalanche_area + 1 );
                fprintf( fptr, "ITIME=%d\n", itime);
                fprintf( fptr, "%d\t%d\n", (ngrid - 1 )/2 , (ngrid-1)/2 );
            } else {
                fprintf( fptr, "%d\n", avalanche_area );
                fprintf( fptr, "ITIME=%d\n", itime);
                for( unsigned int index = 0; index < nsize; index++ )
                    if(failure[ index ] == 1){
                        fprintf( fptr, "%d\t%d\t%d\n", index ,index / ngrid, index % ngrid );
                    }
            }
            fclose( fptr );
        }
        itime += 1 ;
    }
    fclose( f_time_series );
    fclose( f_pressure );
    return 0;
}
