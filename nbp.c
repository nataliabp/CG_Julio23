


#include "Header.h"


//Stuff to compute division thresholds depending on the spatial location
double Compute_division_threshold_test(
            struct sim_parameters *params,
                 //array encoding that info vs spatial location
            double oxygen_concentration //rightmost slot entering the computation (aka nslots here)
){
    
    long x_slot;
    struct sim_parameters *p;
    double p6_over_p3;
    double c_crit;
    double aminus_hat;

    p= (struct sim_parameters *) params;
    p6_over_p3=(p->p6_over_p3);
    c_crit=(p->critical_oxy_hat);
    aminus_hat=(p->aminus_hat);
    double division_threshold; 
    
   

        
    if(oxygen_concentration<=c_crit){division_threshold=INFINITE_TIME;}
    else{
        division_threshold=aminus_hat*pow(oxygen_concentration/c_crit-1,-beta);
                //division_threshold[x_slot] = 100;
                //printf("fraccion %lf", oxygen_concentration[x_slot]/c_crit-1); 
        
        
    };
    
   return division_threshold; 
  
};// end Compute_division_threshold


//Computing the coarse-grained proliferation rate lambda
//Including a survival rate (see eq. 48 de la Cruez et al 2016)
//We removed tau_p since we perform with adimensional units (see implementation notes)
double Get_Eigenvalue_test(
                    double threshold, 
                  
                    double death_rate,
                    double tau_p,//deprecated
                    double survival_rate
                    ){
    
    double aux=0.0;
    double eigenvalue=0.0;
    int i=0;
    

    //Newton's procedure
    //printf("division threshold %lf\n", threshold[xslot]);

    if(threshold>=INFINITE_TIME){eigenvalue=-death_rate;}
    else{
        do{
            eigenvalue=aux;

            aux=eigenvalue
               -(eigenvalue + death_rate+1)
                 *(threshold*(eigenvalue + death_rate)
                    +log(1+eigenvalue + death_rate)-log(2*survival_rate)
                 )
                /(1 + threshold*(eigenvalue + death_rate+1));

            i++;
        }while((fabs(aux-eigenvalue)>TOL)&&(i<7));
        eigenvalue=aux;
        //eigenvalue = 10; 
        
    };//end if-else

    return(eigenvalue);
    
}; //end Get_Eigenvalue
