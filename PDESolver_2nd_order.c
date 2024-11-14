#include "Header.h"



void EvalPopulation(struct sim_parameters *params,
                    double *kpop, //rhs readout
                    double *density_previous,
                    double *division_threshold, double *lambda
                    ){
 
    long x_slot;
    double eigenvalue;
    //Note that the eigenvalue depends on space actually
    long n_xslots;
    double death_rate;
    double tau_p;
    double survival_rate;
    double diff_coef_pop;
    struct sim_parameters *p;

    p= (struct sim_parameters *) params;

    n_xslots=(p->n_xslots);
    death_rate=(p->death_rate_hat);
    tau_p=(p->tau_p);
    survival_rate=(p->survival_rate);
    diff_coef_pop=(p->diff_coef_eff);
    //printf("dif coef ef %e \n", diff_coef_pop); 
    //delta_x=(p->Delta_x_hat);

    
    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
            
            eigenvalue=Get_Eigenvalue(division_threshold,
                            x_slot,
                            death_rate,
                            tau_p,
                            survival_rate);
            
            //printf("eigenvalue %lf \n", eigenvalue);
            lambda[x_slot] = eigenvalue; 
            //printf("lambda %lf, division threshold %lf \n", eigenvalue, *division_threshold); 
            kpop[x_slot]=density_previous[x_slot]*eigenvalue
                +(density_previous[x_slot-1]+density_previous[x_slot+1]-2*density_previous[x_slot])
                    *diff_coef_pop;
            
           
    };
    
    eigenvalue=Get_Eigenvalue(division_threshold,
                            0,
                            death_rate,
                            tau_p,
                            survival_rate);
            //printf("eigenvalue %lf \n", eigenvalue);
    
    lambda[0] = eigenvalue; 
    kpop[0] = density_previous[0]*eigenvalue + 2*diff_coef_pop*(density_previous[1]-density_previous[0]); 
    //kpop[0] = density_previous[0]*eigenvalue + diff_coef_pop*((4*density_previous[1]-density_previous[2])/3); 

    eigenvalue=Get_Eigenvalue(division_threshold,
                            n_xslots-1,
                            death_rate,
                            tau_p,
                            survival_rate);
            //printf("eigenvalue %lf \n", eigenvalue);
    lambda[n_xslots-1] = eigenvalue; 
    kpop[n_xslots-1] = density_previous[n_xslots-1]*eigenvalue + 2*diff_coef_pop*(density_previous[n_xslots-2]-density_previous[n_xslots-1]); 
    //kpop[n_xslots-1] = density_previous[n_xslots-1]*eigenvalue + diff_coef_pop*((4*density_previous[n_xslots-2]-density_previous[n_xslots-3])/3); 

}; //end EvalPopulation(_targeted)




void EvalOxygen(struct sim_parameters *params,
                    double *koxy, //rhs readout
                    double *oxygen_concentration, //old status oxygen
                    double *density_of_individuals
                    ){
    
    double diffusion, reaction;
    long x_slot;

    struct sim_parameters *p;

    double k_consumption;
    double k_decay;
    double diff_coef;
    double source_oxygen;
    long n_xslots;
    double delta_x; //NABP
    double reaction0; 
    p= (struct sim_parameters *) params;
  
    k_consumption=(p->k_consumption_hat);
    k_decay=(p->k_decay_hat);
    diff_coef=(p->diff_coefs_ratio);

    source_oxygen=(p->source_oxygen_hat);
    n_xslots=(p->n_xslots);
    delta_x=(p->Delta_x_hat); //NABP
    double diffusion0; 
    double diffusionn; 
    double reactionn; 

    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
            
            reaction= source_oxygen 
              -k_consumption*oxygen_concentration[x_slot]*density_of_individuals[x_slot]
                -k_decay*oxygen_concentration[x_slot];
            
            diffusion= oxygen_concentration[x_slot-1]+oxygen_concentration[x_slot+1]-2*oxygen_concentration[x_slot];
            
            koxy[x_slot]=reaction+diff_coef*diffusion;//(delta_x*delta_x);
            //koxy[x_slot]=reaction+diff_coef*diffusion/(delta_x*delta_x); //NABP


    }; //end for
    reaction0 = source_oxygen 
              -k_consumption*oxygen_concentration[0]*density_of_individuals[0]
                -k_decay*oxygen_concentration[0];
    diffusion0 = oxygen_concentration[1]-oxygen_concentration[0];
    //diffusion0 = (4*oxygen_concentration[1]-oxygen_concentration[2])/3;
    koxy[0] = reaction0+diff_coef*2*diffusion0; 

    reactionn = source_oxygen 
              -k_consumption*oxygen_concentration[n_xslots-1]*density_of_individuals[n_xslots-1]
                -k_decay*oxygen_concentration[n_xslots-1];
    diffusionn = oxygen_concentration[n_xslots-2]-oxygen_concentration[n_xslots-1];
    //diffusionn = (4*oxygen_concentration[n_xslots-2]-oxygen_concentration[n_xslots-3])/3;

    koxy[n_xslots-1] = reactionn+2*diff_coef*diffusionn; 


}; //End EvalOxygen



void  Eulerpop_oxy( struct sim_parameters *params,
                double delta_t,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration, double *lambda
                ){

    
    double *stat1_pop;
    double *stat1_oxy;
    long i=0;

    struct sim_parameters *p;
    long n_xslots;

    p= (struct sim_parameters *) params;
    n_xslots=(p->n_xslots);

    //book memory  
    if((stat1_pop= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned (PDE evolution) \n");
        exit(1);
    };
    if((stat1_oxy= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned (PDE evolution) \n");
        exit(1);
    };
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)

    EvalPopulation(p,
                    stat1_pop, 
                    density_of_individuals,
                    division_threshold, lambda
                    );

    EvalOxygen(p,
                stat1_oxy,
                oxygen_concentration,
                density_of_individuals
                );
    //printf("delta_t %lf \n", delta_t); 
    //Step2: Euler's formula
    for(i=0;i<n_xslots;i++){
        
        //printf("stat1_pop %lf\n", stat1_pop[i]);
        density_of_individuals[i]=density_of_individuals[i]+delta_t*stat1_pop[i];
        
        oxygen_concentration[i]=oxygen_concentration[i]+delta_t*stat1_oxy[i];
    };
    
    //Boundary conditions treated after bulk advancement is made
    //Zero Neumann at both sides for both species
    //density_of_individuals[0]=(4*density_of_individuals[1]-density_of_individuals[2])/3.0;
    //density_of_individuals[n_xslots-1]=(4*density_of_individuals[n_xslots-2]-density_of_individuals[n_xslots-3])/3.0;
    //oxygen_concentration[0]=(4*oxygen_concentration[1]-oxygen_concentration[2])/3.0;
    //oxygen_concentration[n_xslots-1]=(4*oxygen_concentration[n_xslots-2]-oxygen_concentration[n_xslots-3])/3.0;
    //In case we want inflow conditions, use something like this: (BEWARE OF THE SPATIAL STEP,to be included)
    //oxygen_level[0]=(4*oxygen_level[1]-oxygen_level[2]+2*incoming_oxygen_flux)/3.0;
    //oxygen_level[n_xslots-1]=(4*oxygen_level[n_xslots-2]-oxygen_level[n_xslots-3]+2*incoming_oxygen_flux)/3.0;
    //The incoming flux needs to be rescaled before that (...)
    
   
    //Update division threshold and total number of cells 
    
    Compute_division_threshold(p,
                                division_threshold,
                                oxygen_concentration,
                                n_xslots-1);
    
    //free stuff
    free(stat1_pop);
    free(stat1_oxy);
    
}; //end Eulerpop_oxy




void GlobalPDE_Handler(	
                struct sim_parameters *params,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration,
                double *lambda){    
    double DeltaT; //refined time step
    double peak_density_of_cells;
    double stiff_consumption;
    long j=1;
    double TIME_INCREMENT=1.0;

    struct sim_parameters *p;

    double k_consumption;
    long n_xslots;
    double CFL;

    p= (struct sim_parameters *) params;
  
    TIME_INCREMENT=(p->Delta_t_hat); //Time interval to advance in this call (dimensionless)
     //Time interval to advance in this call (dimensionless)
    //TIME_INCREMENT= 0.0000001; 

    k_consumption=(p->k_consumption_hat);
    n_xslots=(p->n_xslots);
    CFL=(p->CFL_number);
    

    /////////////
    //Refine time step according to stability considerations
    peak_density_of_cells=GetMax(density_of_individuals,n_xslots);
    stiff_consumption=0.5/(1.0+k_consumption*peak_density_of_cells); 
    //we hope that during our time window "peak_number_of_cells" will not double, despite how fast it changes through
    // printf("stiff consumption %e \n", stiff_consumption); 
    // printf("CFL %e \n", CFL); 
   
    DeltaT=min(TIME_INCREMENT,min(CFL,stiff_consumption));
    //printf("Delta T %e \n", DeltaT); 

    //printf("Time increment %e \n", TIME_INCREMENT);
    //DeltaT =0.0000001;
    //printf("Delta t %e \n",DeltaT);
     //DeltaT =0.0000001;
    //////////////////////////////////
    //Perform temporal iteration
    //Create the sub-iteration structure (on each such step we perform Euler marching)
   
    while(j*DeltaT<TIME_INCREMENT){ 
        
        Eulerpop_oxy(p,
                DeltaT,
                density_of_individuals,
                division_threshold,
                oxygen_concentration, lambda
                );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet TIME_INCREMENT =1 (internal uds) exactly:
    DeltaT=TIME_INCREMENT-(j-1)*DeltaT;

    
        Eulerpop_oxy(p,
                DeltaT,
                density_of_individuals,
                division_threshold,
                oxygen_concentration, lambda
                );

    
}; //end GlobalPDE_Handler
