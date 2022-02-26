##################################################################
######## SUBJECT SPECIFIC PARAMETERs ESTIMATION FOR SECTION 4.3.1
##################################################################

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
load(paste('data_gen_ebola_vaccine_model_pop',nb_subject,'.Rdata',sep=""))


# specification of scaling parameter 
coeff_mult = list_Obs_subjects[[1]][[1]][[3]][[3]]

#exogenous paramter specification
delta_S=  list_Obs_subjects[[1]][[1]][[3]][[1]]
delta_L = list_Obs_subjects[[1]][[1]][[3]][[2]]

#########  Specification of the population parameter value (theta and delta)


#theta, sigma and delta specification
x0_mean <- c( 479.0007)
sd_x0 = log(258)

log_delta_S = log(log(2)/1.2)
log_delta_L = log(log(2)/(364*6))
log_phi_S = log(2755)
log_phi_L = log(16.6)
log_delta_Ab =log(log(2)/24)

theta_pop <- c(log_delta_S,log_phi_S,log_phi_L,log_delta_Ab )

std_log_mu_S = 0.92
std_log_mu_L = 0.85
std_log_delta_Ab = 0.30
std_param = c(std_log_mu_S,std_log_mu_L,std_log_delta_Ab)

std = 100*coeff_mult

Psi = diag(std_param^2)
Chol_R= chol(Psi)
true_Delta_mat= chol((std^2)*solve(Psi))
true_theta = theta_pop


#number of tested trial
nb_trial = 1


# Specification of the required model element and dimension for parameter estimation
dim_sub = nrow(Psi)
dim_control = 1
dim_syst = 1
dim_obs = 1


mat_A_pop  = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_matA(t,param_sub,param_pop,exo_par)}
vect_r_pop = function(t,param_sub,param_pop,exo_par){Ebola_Vaccine_model_vectR(t,param_sub,param_pop,exo_par)}

mat_B =diag(c(1),dim_syst,dim_control)
mat_C =  matrix(c(1),dim_obs,dim_syst)

log_prior_function = function(param_pop,Delta_mat){return(0)}

# Specification of the mesh size
mesh_iter=40

#Weighing parameter specification
lambda_seq = c(10^3,10^4,10^5)
nb_tested_U= length(lambda_seq)


#For each tested trial and subject in it: 
#1) Compute the estimation of the b_i  where the optimization algorithm starts from
# 1.1) the true parameter value
# 1.2) a wrongly chosen parameter value to test practical identifiability issues

#2) Plot the solution of the related optimal control (in green), the true original ODE solution (in red) and the noisy data (circle)
for (nb_t  in 1:nb_trial){
  
  
  Obs_subjects_nb_t = list_Obs_subjects[[nb_t]]
  bi_list_est= list_seq_true_bi[[nb_t]]
  bi_list_est_dc= list_seq_true_bi_dc[[nb_t]]
  x0i_list_est = list_seq_true_x0i[[nb_t]]
  known_x0i_list = list_seq_known_x0i[[nb_t]] 
  
  for (nbs in 1:nb_subject){
    
    Times_obs_i =Obs_subjects_nb_t[[nbs]][[1]]
    Y_i =Obs_subjects_nb_t[[nbs]][[2]]
    true_bi =  bi_list_est[[nbs]]
    bi_dc = bi_list_est_dc[[nbs]]
    true_x0i = x0i_list_est[[nbs]]
    known_x0i =  known_x0i_list[[nbs]]
    
    nb_time_integ_i = mesh_iter*(length(Times_obs_i)-1)+1;
    
    Times_integ_i =  matrix(0, nrow = 1, ncol = nb_time_integ_i )
    weight_integ_i = matrix(0, nrow = 1, ncol = nb_time_integ_i )
    
    pseudo_Y_i =  matrix(0, nrow = dim_obs, ncol = nb_time_integ_i )
    
    for (obs_i in 1:(length(Times_obs_i)-1)){
      
      delta_i = Times_obs_i[obs_i+1] - Times_obs_i[obs_i] 
      
      refined_interval_i =  seq(Times_obs_i[obs_i],Times_obs_i[obs_i+1],by = (delta_i/mesh_iter))
      refined_interval_i[length( refined_interval_i)] = Times_obs_i[obs_i+1]
      Times_integ_i[seq(mesh_iter*(obs_i-1)+1,(mesh_iter*obs_i+1),by =1)] = refined_interval_i
      
      weight_integ_i[mesh_iter*(obs_i-1)+1] =  1/(refined_interval_i[2]-Times_obs_i[obs_i])
      
      
      pseudo_Y_i[,mesh_iter*(obs_i-1)+1]= Y_i[,obs_i]
      pseudo_Y_i[,mesh_iter*obs_i+1] =Y_i[,obs_i+1]
      
    }
    Times_integ_i = c( Times_integ_i)
    weight_integ_i = c(weight_integ_i)
    
    for (nb_l  in 1:nb_tested_U){
      
      mat_U = lambda_seq[nb_l]*diag(dim_control)
      print(mat_U)
      
      list_inner_arg = list(Times_integ_i, pseudo_Y_i, weight_integ_i,true_bi,true_Delta_mat,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C, known_x0i,2,0
                            ,true_theta,list(delta_S,delta_L,coeff_mult),type_algorithm_optim = 1,c(200,1))
      res_est_inner = est_param_oca_inner_criteria_lincase_par_version(list_inner_arg)
      
      
      list_inner_arg_dc =list(Times_integ_i, pseudo_Y_i, weight_integ_i,bi_dc,true_Delta_mat,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C, known_x0i,2,0
                              ,true_theta,list(delta_S,delta_L,coeff_mult),type_algorithm_optim = 1,c(200,1))
      res_est_inner_dc = est_param_oca_inner_criteria_lincase_par_version(list_inner_arg_dc)
      
      
      
      print(t(true_bi))
      print("estimation from true initial guess ")
      print(  res_est_inner$bi_est)
      print("estimation from wrong initial guess ")
      print(  res_est_inner_dc$bi_est)
      
      main_title = paste("Solution for subject ",nbs,", lambda = ", lambda_seq[nb_l],sep="")
      out <- ode( y =  true_x0i,times = Times_integ_i,func = Ebola_vaccine_model_sim, parms =c(theta_gen, true_bi))
      matplot(Times_integ_i,coeff_mult*out[,2],type="l" ,xlab = "Time", ylab = "X",main = main_title,col="red",ylim= c(0,1.2*coeff_mult*max(out[,2])))
      lines(Times_integ_i,res_est_inner$Xi_opt[1,],col="green")
      points(Times_obs_i, Y_i,pch=16,bg="black")
      legend("topright",c("True Trajectories","Est. Optimal Trajectories","Observations"),pch=c(NA,NA,16),lty=c(1,1,NA),bg=c(NA,NA,"black"),col=c("red", "green", "black"), lwd = c(1,1,1))
    }
    
    
  }
  
}
