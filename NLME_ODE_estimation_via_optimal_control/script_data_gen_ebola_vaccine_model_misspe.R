############################################
######## DATA GENERATION FOR SECTION 4.3.2
############################################


pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')

nb_trial = 100

# specification of scaled parameting  for Ab value, population size and observation numbers
nb_subject =20
coeff_mult = 0.01
Times_obs<- seq(0, 364, by = 364/10)
Times_ref <- seq(0, 364, by = 364/200)

nb_obs <- length(Times_obs)

# data file name
rep_file = paste('data_gen_ebola_vaccine_model_misspe_pop',nb_subject,'.Rdata',sep="")
rep_file_trajectory = paste('true_trajectory_ebola_vaccine_model_misspe_pop_',nb_subject,'.Rdata',sep="")

######### Parameter value specification

#theta, sigma and delta specification
x0_mean <- c( 479.0007)
sd_x0 = log(258)

log_delta_S = log(log(2)/1.2)
log_delta_L = log(log(2)/(364*6))
log_phi_S = log(2755)
log_phi_L = log(16.6)
log_delta_Ab =log(log(2)/24)

theta_gen <- c(log_delta_S ,log_delta_L,log_phi_S,log_phi_L,log_delta_Ab )

std_log_mu_S = 0.92
std_log_mu_L = 0.85
std_log_delta_Ab = 0.30
std_param = c(std_log_mu_S,std_log_mu_L,std_log_delta_Ab)

std = 100*coeff_mult

Psi = diag(std_param^2)
Chol_R= chol(Psi)
Delta_mat= chol((std^2)*solve(Psi))


list_Obs_subjects = list()
list_true_trajectory = list()
list_seq_true_bi = list()
list_seq_true_bi_dc = list()
list_seq_true_x0i = list()
list_seq_known_x0i = list()

# Data generation
for (nbt in 1:nb_trial){
  
  Obs_subjects = list()
  true_trajectory = list()
  seq_true_bi = list()
  seq_true_bi_dc = list()
  seq_true_x0i = list()
  seq_known_x0i = list()
  
  
  for (nbs in 1:nb_subject){
    true_bi =  Chol_R%*%rnorm(3 ,0);
    x0i = min(max(x0_mean+ sd_x0*rnorm(1,0),200),800)
    
    out_ref <- ode(func = Ebola_vaccine_model_sim , y = x0i, parms =c(theta_gen,true_bi), times = Times_ref)
    out_pert_ref = Ebola_vaccine_model_pert_sim(deb=0,h=Times_obs[length(Times_obs)]/1000,fin=Times_obs[length(Times_obs)],
                                        Pars = c(theta_gen,true_bi),x0i,sd_stoch_pert =10)
    out_pert = t(out_pert_ref[[2]][seq(1,length(out_pert_ref[[1]]),1001/10)])
     out_pert  =c( out_pert,  out_pert_ref[[2]][1001])
    Time_pert  =  t(out_pert_ref[[1]][seq(1,length(out_pert_ref[[1]]),1001/10)])
    Time_pert = c(cbind( Time_pert,out_pert_ref[[1]][1001]))
    
    Y_noise = t(coeff_mult*out_pert +std*rnorm(n = nb_obs))
    
    true_trajectory =  append(true_trajectory,list(out_pert_ref))
    
    Obs_subjects = append(Obs_subjects,list(list(Time_pert,Y_noise,list(log_delta_S,log_delta_L,coeff_mult))))
    seq_true_bi  = append(seq_true_bi, list(true_bi))
    seq_true_bi_dc = append(seq_true_bi_dc, list(1.3*true_bi))
    seq_true_x0i  = append(seq_true_x0i,list(x0i))
    seq_known_x0i  = append(seq_known_x0i,list(list()))
    
    if (nbt ==1){
      if (nbs==1){
        matplot(Times_ref ,coeff_mult*out_ref[,2], type = "l", xlab = "Time", ylab = "",main = "Perturbed Ebola vaccine model",ylim = c(0,1.4*coeff_mult*max(out_ref[,2])))
        lines(out_pert_ref[[1]],coeff_mult*out_pert_ref[[2]][1,],lty=2) 
        points(Times_obs,Y_noise[1,],pch=16,bg="black") 
        legend("topright",c("True stochastic curves","Deterministic approximation (ODE solution)","Observations"),pch=c(NA,NA,16),lty=c(2,1,NA),bg=c(NA,NA,"black"),col=c("black", "black", "black"), lwd = c(1,1,1))
      }
    }
  }
  
  list_true_trajectory = append(list_true_trajectory,list(true_trajectory))
  
  list_Obs_subjects  = append( list_Obs_subjects,list( Obs_subjects ))
  list_seq_true_bi = append(list_seq_true_bi,list( seq_true_bi))
  list_seq_true_bi_dc  = append(list_seq_true_bi_dc ,list(seq_true_bi_dc))
  list_seq_true_x0i = append(list_seq_true_x0i, list(seq_true_x0i) )
  list_seq_known_x0i = append(list_seq_known_x0i, list(seq_known_x0i) )
}

save(file=rep_file,list_Obs_subjects,list_seq_true_bi , list_seq_true_bi_dc, list_seq_true_x0i, list_seq_known_x0i )
save(file=rep_file_trajectory,list_true_trajectory)
