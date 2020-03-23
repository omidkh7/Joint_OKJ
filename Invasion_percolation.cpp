 #include <stdio.h>
#include <vector>
#include <math.h>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <Eigen/Sparse>
#include <algorithm>
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
//--
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
    const unsigned int ngrid = 2001;
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
    
    vector<double> permeability_x; //--- Contains permeability of bonds on the X axis
    vector<double> permeability_y; //--- Contains permeability of bonds on the Y axis
    vector<double> porosity_x; //--- Contains porosity of bonds on the X axis
    vector<double> porosity_y; //--- Contains porosity of bonds on the Y axis
    vector<double> rho( nsize ); //--- pressure
    vector<int> K_x_change( nsize ); //--- state variable of 0: it did not get increased   1: it's already increased
    vector<int> K_y_change( nsize ); //--- state variable of 0: it did not get increased   1: it's already increased
    
    vector<int> Fluid( nsize ); //--- state variable 0: No fluid, 1: neighbor to filled with fluid, 2: filled with fluid
    vector<double> Delta_up (nsize); // This contains (mu_ij + p_ij) to maximize it and find the next invading site.
    vector<double> Delta_down (nsize); // This contains (mu_ij + p_ij) to maximize it and find the next invading site.
    vector<double> Delta_left (nsize); // This contains (mu_ij + p_ij) to maximize it and find the next invading site.
    vector<double> Delta_right (nsize); // This contains (mu_ij + p_ij) to maximize it and find the next invading site.
    //--- reading porosity and permeability fields
    permeability_x = ReadInDoubles( "K_x.txt" );
    permeability_y = ReadInDoubles( "K_y.txt" );
    porosity_x = ReadInDoubles( "Phi_x.txt" );
    porosity_y = ReadInDoubles( "Phi_y.txt" );
    //--
    //--- parse Shear stress and Yield function
    tau = ReadInDoubles( "shear.txt" ); // uniform distribution  (0.0,1.0)
    S = ReadInDoubles( "yield.txt" ); // uniform distribution (S_min , S_max)
    
    //---
    FILE *f_fracture = fopen("Data.txt","w");
    fprintf( f_fracture, "#ITIME\tRad\tRMS\tMSD\tMSD_rad\tArea\tSIZE\tEnergy\tX_hypo\tY_hypo\n" );
    FILE *f_pressure = fopen("pressure.txt","w");
    FILE *f_Tau = fopen("shear_new.txt","w");
    //--
    initialize( &K_x_change, nsize ); //-- this will set all the values of this equall to zero which means all bonds on the x axis are not broken
    initialize( &K_y_change, nsize ); //-- this will set all the values of this equall to zero which means all bonds on the y axis are not broken
    
    initialize( &Fluid, nsize ); //-- this will set all the values equall to zero, which means all sites don't contain fluids
    
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
    //--
    unsigned int itime = 0;
    unsigned int Number_ofevents = 0;
    unsigned Total_events= 30000 ;
    while (Number_ofevents != Total_events) {
        
        if (itime % 10000 == 0 ) {
            printf( "itime = %d\n", itime);
        }
        
        //-- injecting fluid to the center
        rho[ index_injection ] = 1.0;
        Fluid[ index_injection ] = 2; // This will set the injection point to contain fluids
        
        //--the for loop below would set all the neighbors of the sites filled with fluid to 1.
        unsigned int IndeX=0,index_west, index_north,index_east,index_south;
        for( int sitex = 0; sitex < ngrid; sitex++ )
        {
            for( int sitey = 0; sitey < ngrid; sitey++ )
            {
            if ( Fluid[ IndeX ] == 2 ){

                getIndex( ngrid, sitex, sitey - 1  < 0 ? sitey - 1 + ngrid : sitey - 1 , &index_west );
                getIndex( ngrid, ( sitex - 1  < 0 ? sitex - 1  + ngrid : sitex - 1 ), sitey, &index_north );
                getIndex( ngrid, sitex, ( sitey + 1 ) % ngrid, &index_east );
                getIndex( ngrid, ( sitex + 1 ) % ngrid, sitey, &index_south );
                if( Fluid[ index_west ] == 0){
                    Fluid[ index_west ] = 1;
                    }
                if( Fluid[ index_north ] == 0){
                    Fluid[ index_north ] = 1;
                    }
                if( Fluid[ index_south ] == 0){
                    Fluid[ index_south ] = 1;
                    }
                if( Fluid[ index_east ] == 0){
                    Fluid[ index_east ] = 1;
                    }
                }
            IndeX++;
            }
        }
        unsigned int indeX = 0 ,index_LEFT, index_UP,index_RIGHT,index_DOWN;
        //-- this will set all the values equall to zero
        for( unsigned int iI = 0; iI < nsize; iI++ ){

                Delta_up[ iI ] = 0.0;

                Delta_down[ iI ] = 0.0;

                Delta_left[ iI ] = 0.0;

                Delta_right[ iI ] = 0.0;
        }
        
        for( int Node_i = 0; Node_i < ngrid; Node_i++ )
        {
            for( int Node_j = 0; Node_j < ngrid; Node_j++ )
            {
                if ( Fluid[ indeX ] == 1) {
                    getIndex( ngrid, Node_i, Node_j - 1  < 0 ? Node_j - 1 + ngrid : Node_j - 1 , &index_LEFT );
                    getIndex( ngrid, ( Node_i - 1  < 0 ? Node_i - 1  + ngrid : Node_i - 1 ), Node_j, &index_UP );
                    getIndex( ngrid, Node_i, ( Node_j + 1 ) % ngrid, &index_RIGHT );
                    getIndex( ngrid, ( Node_i + 1 ) % ngrid, Node_j, &index_DOWN );
                    
                    if ( Fluid[ index_UP ] == 2 ){
                        Delta_up[indeX] =  permeability_y[ index_UP ] * rho[index_UP]; //up
                    }
                    
                    if ( Fluid[ index_DOWN ] == 2 ){
                        Delta_down[indeX] = permeability_y[ indeX ] * rho[index_DOWN]; //down
                    }
                    
                    if ( Fluid[ index_LEFT ] == 2 ){
                        Delta_left[indeX] = permeability_x[ index_LEFT ] * rho[index_LEFT];  //left
                    }
                    
                    if ( Fluid[ index_RIGHT ] == 2 ){
                        Delta_right[indeX] = permeability_x[ indeX ] * rho[index_RIGHT];  //right
                    }
                }
                indeX++;
            }
        }
        // Now we want to find the maximum (most permeable bond)
        int maxIndex_up,maxIndex_down,maxIndex_left,maxIndex_right;
        double Max_up,Max_down,Max_left,Max_right;
        int i_up=0;
        int i_down=0;
        int i_right=0;
        int i_left=0;
        maxIndex_up = i_up ;
        maxIndex_down = i_down;
        maxIndex_left = i_left ;
        maxIndex_right = i_right;
        Max_up=Delta_up[i_up];
        Max_down=Delta_down[i_down];
        Max_right=Delta_right[i_right];
        Max_left=Delta_left[i_left];
        int iii=0;
        while(iii < nsize ){
            if (Max_up < Delta_up[i_up]){
                Max_up = Delta_up[i_up];
                maxIndex_up = i_up;
                }
            i_up++;

            if (Max_down < Delta_down[i_down]){
                Max_down = Delta_down[i_down];
                maxIndex_down = i_down;
            }
            i_down++;
            
            if (Max_left < Delta_left[i_left]){
                Max_left = Delta_left[i_left];
                maxIndex_left = i_left;
            }
            i_left++;
            
            if (Max_right < Delta_right[i_right]){
                Max_right = Delta_right[i_right];
                maxIndex_right = i_right;
            }
            i_right++;
            iii++;
        }
        // invading the maximum permeable site
        // let's say when a fluid invade to a site, it doeas have a minimum pressure = p_min
        if ( (Max_up > Max_left) && (Max_up > Max_right) && (Max_up > Max_down) ) { //the fluid will go to the upper site
            Fluid[maxIndex_up] = 2;
            }
            
        if ( (Max_down > Max_left) && (Max_down > Max_right) && (Max_down > Max_up) ) { //the fluid will go to the southern site
            Fluid[maxIndex_down] = 2;
            }
            
        if ( (Max_left > Max_up) && (Max_left > Max_right) && (Max_left > Max_down) ) {//the fluid will go to the left site
            Fluid[maxIndex_left] = 2;
            }
        if ( (Max_right > Max_up) && (Max_right > Max_left) && (Max_right > Max_down) ) {//the fluid will go to the right site
            Fluid[maxIndex_right] = 2;
            }
        
        unsigned int index = 0, index_W, index_E, index_N, index_S;
        for( int node_i = 0; node_i < ngrid; node_i++ )
        {
            for( int node_j = 0; node_j < ngrid; node_j++ )
            {
                
                if( Fluid[index] == 2 )
                {
                    if( index == index_injection )
                    {
                        index++;
                        continue;
                    }
                    
                    int X_i = ((ngrid - 1 ) / 2 ) - (node_i);
                    int Y_j = ((ngrid - 1 ) / 2 ) - (node_j);
                    float R2 = X_i * X_i + Y_j * Y_j ;
                    
                    getIndex( ngrid, node_i, node_j - 1  < 0 ? node_j - 1 + ngrid : node_j - 1 , &index_W );
                    getIndex( ngrid, node_i, ( node_j + 1 ) % ngrid, &index_E );
                    getIndex( ngrid, ( node_i + 1 ) % ngrid, node_j, &index_S );
                    getIndex( ngrid, ( node_i - 1  < 0 ? node_i - 1  + ngrid : node_i - 1 ), node_j, &index_N );
                    /// harmonic mean of the permeabilitis and
                    double K_E = ( permeability_x[ index ] * porosity_x[ index ]) ;
                    double K_N = ( permeability_y[ index_N ] * porosity_y[ index_N ]) ;
                    double K_S = ( permeability_y[ index ] * porosity_y[ index ]) ;
                    double K_W = ( permeability_x[ index_W ] * porosity_x[ index_W ]) ;
                    double D;
                    if ( (K_E > K_N) && (K_E > K_S) && (K_E > K_W) ) { //the fluid will go to the upper site
                        D = 0.01 * (K_E) ; //1e-4 coming from = (Characteristic Permeability) / (viscosity * effective Compressibility)
                    }
                    
                    if ( (K_N > K_E) && (K_N > K_S) && (K_N > K_W) ) { //the fluid will go to the upper site
                        D = 0.01 * (K_N) ; //1e-4 coming from = (Characteristic Permeability) / (viscosity * effective Compressibility)
                    }
                    
                    if ( (K_S > K_N) && (K_S > K_N) && (K_S > K_W) ) { //the fluid will go to the upper site
                        D = 0.01 * (K_S) ; //1e-4 coming from = (Characteristic Permeability) / (viscosity * effective Compressibility)
                    }
                    if ( (K_W > K_N) && (K_W > K_S) && (K_W > K_E) ) { //the fluid will go to the upper site
                        D = 0.01 * (K_W) ; //1e-4 coming from = (Characteristic Permeability) / (viscosity * effective Compressibility)
                    }
                    rho[ index ] = exp((-1 * R2)/(4 * D * itime));
                    //---
                }
                index++;
            }
        }
        //--- yield function
        unsigned int time_internal = 0, avalanche_area = 0, avalanche_size = 0;
        double  energy_rel = 0.0;
        initialize( &failure, nsize );
        unsigned int signtaure = 0;//signtaure_of_saving_the_first_failure
        int X_hypocenter,Y_hypocenter; // the coordinates of rhe first initiation site
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
                        
                        if ( signtaure == 0 ){
                            X_hypocenter = ((ngrid - 1 ) / 2 ) - (Node_i) ;
                            Y_hypocenter = ((ngrid - 1 ) / 2 ) - (Node_j) ;
                            signtaure = 1;
                        }

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
                                permeability_y[ index_n ] = 100 * permeability_y[ index_n ] ;
                                K_y_change[ index_n ] = 1;
                                
                            }
                        }
                        if (tau[ index_e ] + rho[ index_e ] > S[ index_e ]) { //--for the eastern neighbour, the index of bond is the same as the index of the i site
                            if (K_x_change[ i ] == 0 ){
                                permeability_x[ i ] = 100 * permeability_x[ i ] ;
                                K_x_change[ i ] = 1;
                            }
                        }
                        if ( tau[ index_s ] + rho[ index_s ] > S[ index_s ]) { //-- Also for the southern neighbour, the index of bond is the same as the index of the i site
                            if ( K_y_change[ i ] == 0 ){
                                permeability_y[ i ] = 100 * permeability_y[ i ] ;
                                K_y_change[ i ] = 1;
                            }
                        }
                        if ( tau[ index_w ] + rho[ index_w ] > S[ index_w ]) {
                            if (K_x_change[ index_w ] == 0){
                                permeability_x[ index_w ] = 100 * permeability_x[ index_w ] ;
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
        //--- output
        for( unsigned int i = 0; i < nsize; i++ )
            avalanche_area += ( failure[ i ] == 1 );
        
        //--- /Distance to event centroid, MSD(mean squared distance), RMS(root mean squared) of events
        //--- /Distance to event centroid, MSD(mean squared distance), RMS(root mean squared) of events
        if( avalanche_area > 0 )
        {
            double X_com = 0.0, Y_com = 0.0,Number = 0.0 , r3 = 0.0 ,denominatorR = 0.0 , r3p = 0.0 , denominator = 0.0, R3 = 0.0 , deno = 0.0  ;
            unsigned int IND = 0;
            double dis_to_centroid , RMS_of_events , MSD_pressure, MSD_rad;
            for( int xx = 0; xx < ngrid; xx++ )
            {
                for( int yy = 0; yy < ngrid; yy++ )
                {
                    if (failure[ IND ] == 1 ){
                        int XXX = ((ngrid - 1 ) / 2 ) - (xx);
                        int YYY = ((ngrid - 1 ) / 2 ) - (yy);
                        X_com += XXX;
                        Y_com += YYY;
                        double R = sqrt(XXX * XXX + YYY * YYY);
                        r3 += R * R * R ;
                        denominatorR += R ;
                        Number += 1.0;
                        dis_to_centroid = sqrt((X_com * X_com)/(Number * Number)+(Y_com * Y_com)/(Number * Number));
                        RMS_of_events =  sqrt(r3/denominatorR); //root mean squared of events
                    }
                    if( rho[ IND ]  > 0.0 ){
                        int XX = ((ngrid - 1 ) / 2 ) - (xx);
                        int YY = ((ngrid - 1 ) / 2 ) - (yy);
                        double r = sqrt(XX * XX + YY * YY);
                        r3p += r * r * r * rho[ IND];
                        denominator += r * rho[ IND];
                        MSD_pressure = r3p / denominator;
                        R3 += r * r * r;
                        deno += r;
                        MSD_rad = R3 / deno;
                    }
                    IND++;
                }
            }
            fprintf( f_fracture, "%d\t%e\t%e\t%e\t%e\t%d\t%d\t%f\t%d\t%d\n", itime, dis_to_centroid, RMS_of_events , MSD_pressure, MSD_rad , avalanche_area, avalanche_size, energy_rel , X_hypocenter , Y_hypocenter );
            fflush( f_fracture );
        }
        
        //--- plot fluid pressure
        if( Number_ofevents == Total_events )
        {
            fptr = fopen( "vmd.xyz", "a" );
            fprintf( fptr, "%d\n", nsize );
            fprintf( fptr, "ITIME=%d\n", itime );
            for( unsigned int index = 0; index < nsize; index++ ){
                fprintf( fptr, "%d\t%d\t%d\t%e\n", index, index / ngrid , index % ngrid, fabs( rho[ index ] ) > 0.0 ? rho[ index ] : 0.0 );
                fprintf(f_pressure,"%e\n",rho[index]);
                fprintf(f_Tau,"%e\n",tau[index]);
                }
            fclose( fptr );
            fflush(f_pressure);
            fflush(f_Tau);
        }
        /*
        //--- plot fractured sites
        if( avalanche_area > 100 )
        {
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
        */
        itime += 1 ;
        }
}
