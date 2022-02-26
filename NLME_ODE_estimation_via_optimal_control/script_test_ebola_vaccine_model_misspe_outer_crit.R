###########################################################
######## POPULATION PARAMETERs ESTIMATION FOR SECTION 4.3.1
###########################################################

pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')


# test data retrieval
nb_subject =20
load(paste('data_gen_ebola_vaccine_model_misspe_pop',nb_subject,'.Rdata',sep=""))

####### Specify if theta_Si has to be estimated (=1 yes, no otherwise)
delta_S_estimation = 0

# specification of scaling parameter 
coeff_mult = list_Obs_subjects[[1]][[1]][[3]][[3]]


#########  Specification of the population parameter value (theta and delta) as first guess for the outer criteria optimization algorithm

#theta, sigma and delta specification
x0_mean <- c( 479.0007)
sd_x0 = log(258)

log_delta_S = log(log(2)/1.2)
log_delta_L = log(log(2)/(364*6))
log_phi_S = log(2755)
log_phi_L = log(16.6)
log_delta_Ab =log(log(2)/24)

theta_gen <- c(log_delta_S ,log_delta_L,log_phi_S,log_phi_L,log_delta_Ab )

if (delta_S_estimation ==1){
  theta_pop <- c(log_delta_S,log_phi_S,log_phi_L,log_delta_Ab )
}else{
  theta_pop <- c(log_phi_S,log_phi_L,log_delta_Ab )
}

std_log_mu_S = 0.92
std_log_mu_L = 0.85
std_log_delta_Ab = 0.30
std_param = c(std_log_mu_S,std_log_mu_L,std_log_delta_Ab)

std = 100*coeff_mult

Psi = diag(std_param^2)
Chol_R= chol(Psi)
Delta_mat= chol((std^2)*solve(Psi))


#number of tested trial
nb_trial = 100


# Specification of the required model element and dimension for parameter estimation
dim_sub = nrow(Psi)
dim_control = 1
dim_syst = 1
dim_obs = 1


if (delta_S_estimation ==1){
  mat_A_pop  = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_matA(t,param_sub,param_pop,exo_par)}
  vect_r_pop = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_vectR(t,param_sub,param_pop,exo_par)}
}else{
  mat_A_pop  = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_matA_deltaS_known(t,param_sub,param_pop,exo_par)}
  vect_r_pop = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_vectR_deltaS_known(t,param_sub,param_pop,exo_par)}
}

mat_B =diag(c(1),dim_syst,dim_control)
mat_C =  matrix(c(1),dim_obs,dim_syst)


log_prior_function = function(param_pop,Delta_mat){return(0)}

# precize the if the subject specific parameters have to be estimated or are fixed
est_pop_parameter_only=1

# precize the if the subject specific parameters variance has to be estimated 
delta_known =0

# precize if Asymptotic variance-covariance has to be estimated
est_var = 0

# Specification of the mesh size
mesh_iter=40

#Weighing parameter specification
lambda_seq = c(10^3,10^4,10^5)
nb_tested_U= length(lambda_seq)


#For each tested trial and subject in it: 
#1) Compute the estimation of the population paramters where the optimization algorithm starts from
# 1.1) the true parameter value (registered in list_res_trial)
# 1.2) a wrongly chosen parameter value to test practical identifiability issues (registered in list_res_trial_dc)

list_res_trial =  list()
list_res_trial_dc =  list()

# specify if the inner criteria optimization is parallelized (using snow package) among the subjects
cl_cur <- makeSOCKcluster(rep("localhost",nb_subject))
#cl_cur <- list()

for (nb_t in 1:nb_trial){
  
  list_res_est= list()
  list_res_est_dc= list()
  
  
  Obs_subjects_nb_t = list_Obs_subjects[[nb_t]]
  bi_list_est= list_seq_true_bi[[nb_t]]
  bi_list_est_dc= list_seq_true_bi_dc[[nb_t]]
  x0i_list_est = list_seq_true_x0i[[nb_t]]
  known_x0i_list = list_seq_known_x0i[[nb_t]] 
  
  theta_pop_ini = theta_pop 
  delta_ini = log(diag(Delta_mat))
  
  theta_pop_ini_dc = 0.8*theta_pop_ini
  delta_ini_dc =  0.8*delta_ini 
  
  
  for (lambda  in lambda_seq){
    T1<-Sys.time() 
    
    mat_U = lambda*diag(dim_control)
    out_outer_crit_dc = est_param_oca_outer_criteria_lincase_prof_par_version(Obs_subjects_nb_t,param_pop_ini=theta_pop_ini_dc,log_Delta_vect_ini=delta_ini_dc,
                                                                              mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,
                                                                              mesh_iter,bi_list_est,known_x0i_list,
                                                                              est_pop_parameter_only,delta_known,log_prior=log_prior_function,
                                                                              type_optim =0,est_var =est_var,type_optim_inner=0,nb_iter_max=c(1500),cl_cur)
    
    
    out_outer_crit = est_param_oca_outer_criteria_lincase_prof_par_version(Obs_subjects_nb_t,param_pop_ini=theta_pop_ini,log_Delta_vect_ini=delta_ini,
                                                                           mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,
                                                                           mesh_iter,bi_list_est,known_x0i_list,
                                                                           est_pop_parameter_only,delta_known,log_prior=log_prior_function,
                                                                           type_optim =0,est_var =est_var,type_optim_inner=0,nb_iter_max=c(1500),cl_cur)
    T2<-Sys.time() 
    time_parlapply = T2 -T1
     print("True population parameter values")
    print(c(theta_pop ,log(diag(Delta_mat))))
    
    print("parameter  estimation from true initial guess ")
    print(c(out_outer_crit$theta_est,out_outer_crit$delta_est))
    
    print("parameter  estimation from wrong initial guess ")
    print(c(out_outer_crit_dc$theta_est,out_outer_crit_dc$delta_est))
    
    if (est_var ==1){
      print("Variance estimation from true initial guess")
      print(diag(out_outer_crit$Variance_component$est_var_cov_Matrix))
      
      print("Variance estimation from wrong initial guess")
      print(diag(out_outer_crit_dc$Variance_component$est_var_cov_Matrix))
    }
    
    theta_pop_ini = out_outer_crit$theta_est
    delta_ini = out_outer_crit$delta_est
    
    theta_pop_ini_dc = out_outer_crit_dc$theta_est
    delta_ini_dc = out_outer_crit_dc$delta_est
    
    list_res_est = append(list_res_est,list(list(out_outer_crit,time_parlapply)))
    list_res_est_dc=  append(list_res_est_dc,list(list(out_outer_crit_dc,time_parlapply)))
    
  }
  list_res_trial = append(list_res_trial,list(list_res_est))
  list_res_trial_dc = append(list_res_trial_dc,list(list_res_est_dc))
}
stopCluster(cl_cur)
