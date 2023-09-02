functions {
   vector expand_history_param_to_individuals(int option, real fixed_value, vector mu, vector sigma, vector rnd, int randomN, int link_function){
        // Uses option and parameters to compute value for each individual.
        //
        // Parameters:
        // -----------
        // int option : 1 (constant), 2 (single value), 3 (random independent), 4 (random pooled)
        // real fixed_value : value for option == 1 (oConstant)
        // vector mu : population mean, applicable for option == 2 or option == 4
        // vector sigma : population variance, applicable for option == 4
        // vector rnd : either individual means (option == 3) or individual z-scores (option == 4)
        // int randomN : number of random clusters
        // int link_function : 1 (identityt), 2 (log), 3 (logit)
        //
        // Returns:
        // -----------
        // vector[randomN] : value for each individual.

       vector[randomN] ind;
       if (option == 1) { // constant
           ind = rep_vector(fixed_value, randomN);
           return ind;
       }
       else if (option == 2) { // single
           ind = rep_vector(mu[1], randomN);
       }
       else if (option == 3) { // random independent
           ind = rnd;
       }
       else {
           // pooled: mean + variance * z-score
           ind = mu[1] + sigma[1] * rnd;
       }

       // using link function
       if (link_function == 2) { //log
           return exp(ind);
       }
       else if (link_function == 3) { // logit
           return inv_logit(ind);
       }

       return(ind);
   }
}
data{
    // --- Family choice ---
    // 1 - Gamma, linear model for both shape and rate
    // 2 - Log normal with linear model for the mean
    // 3 - Normal with linear model for the mean
    int<lower=1, upper=3> family;

    // --- Complete time-series ---
    int<lower=1> rowsN;   // Number of rows in the COMPLETE multi-timeseries table including mixed phase.
    array[rowsN] real duration; // Duration of a dominance/transition phase
    array[rowsN] int istate;     // Index of a dominance istate, 1 and 2 code for two competing clear states, 3 - transition/mixed.
    array[rowsN] int is_used;   // Whether history value must used to predict duration or ignored (mixed phases, warm-up period, last, etc.)
    array[rowsN] int run_start; // 1 marks a beginning of the new time-series (run/block/etc.)
    array[rowsN] real session_tmean;    // Mean dominance phase duration for both CLEAR percepts. Used to scale time-constant.

    // --- A shorter clear-states only time-series ---
    int clearN;                  // Number of rows in the clear-states only time-series
    array[clearN] real clear_duration; // Durations of clear percepts only.

    // --- Random effects ---
    int<lower=1> randomN;                               // Number of levels for random effects
    array[rowsN] int<lower=1, upper=randomN> irandom;         // Index of a random effect cluster (participant, display, etc.)

    // --- Fixed effects ---
    int<lower=0> fixedN;
    matrix[rowsN, fixedN > 0 ? fixedN : 0] fixed;

    // --- Cumulative history parameters
    array[2] real<lower=0, upper=1> history_starting_values; // Starting values for cumulative history at the beginning of the run

    // time constant
    int<lower=1, upper=4> tau_option;  // 1 - constant provided by user, 2 - fit single tau for all, 3 - independent taus, 4 - pooled (multilevel) taus
    real<lower=0> fixed_tau;           // a fixed option (tau_option == 1)
    int tau_mu_size;                   // dimensionality, 1 - sampled, 0 - unused
    int tau_sigma_size;                // dimensionality, 1 - sampled, 0 - unused
    int tau_rnd_size;                  // dimensionality, randomN - sampled, 0 - unused
    array[2] real tau_prior;                 // prior

    // Mixed state
    int<lower=1, upper=4> mixed_state_option; // 1 - constant provided by user, 2 - fit single tau for all, 3 - independent taus, 4 - pooled (multilevel) taus
    real<lower=0, upper=1> fixed_mixed_state; // a fixed option (mixed_state_option == 1)
    int mixed_state_mu_size;                  // dimensionality, 1 - sampled, 0 - unused
    int mixed_state_sigma_size;               // dimensionality, 1 - sampled, 0 - unused
    int mixed_state_rnd_size;                 // dimensionality, randomN - sampled, 0 - unused
    array[2] real mixed_state_prior;                // prior

    // --- Linear model ---
    int<lower=1> lmN; // number of linear models, i.e., distribution parameters that are modelled
    int<lower=0, upper=1> varianceN; // number of variance parameters. Effectively, whether variance is sampled at all(1) or not (0)

    // intercept: independent for each parameter and random cluster
    array[lmN, 2] real a_prior;

    // priors for effect of history
    array[lmN, 2] real bH_prior;

    // priors for fixed effects
    matrix[fixedN > 0 ? fixedN : 0, 2] fixed_priors;
}
transformed data {
    // Constants for likelihood index
    int lGamma = 1;
    int lLogNormal = 2;
    int lNormal = 3;

    // Constants for sampling options for cumulative history parameters
    int oConstant = 1;
    int oSingle = 2;
    int oIndependent = 3;
    int oPooled = 4;

    // Options for fitting fixed effects
    int fSingle = 1;
    int fPooled = 2;

    // Link function codes
    int lfIdentity = 1;
    int lfLog = 2;
    int lfLogit = 3;
}
parameters {
    // --- History ---
    // tau
    vector[tau_mu_size] tau_mu;       // population-level mean
    vector<lower=0>[tau_sigma_size] tau_sigma; // population-level variance
    vector[tau_rnd_size] tau_rnd;     // individuals

    // mixed state
    vector[mixed_state_mu_size] mixed_state_mu;       // population-level mean
    vector<lower=0>[mixed_state_sigma_size] mixed_state_sigma; // population-level variance
    vector[mixed_state_rnd_size] mixed_state_rnd;     // individuals

    // --- Linear model ---
    // Independent intercept for each random cluster
    array[lmN] vector[randomN] a;

    // history term, always multilevel but this makes sense only for randomN > 1
    array[lmN] real bH_mu;
    array[randomN > 1 ? lmN : 0] real<lower=0> bH_sigma;
    array[randomN > 1 ? lmN : 0] vector[randomN > 1 ? randomN : 0] bH_rnd;

    // fixed-effects
    array[fixedN > 0 ? lmN : 0] row_vector[fixedN] bF;

    // variance
    vector<lower=0>[varianceN] sigma;
}
transformed parameters {
    array[lmN] vector[clearN] lm_param;
    {
        // Service variables for computing cumulative history
        matrix[2, 3] level;
        array[2] real current_history;
        real tau;
        real dH;

        // Index of clear percepts
        int iC = 1;

        // Cumulative history parameters
        vector[randomN] tau_ind = expand_history_param_to_individuals(tau_option, fixed_tau, tau_mu, tau_sigma, tau_rnd, randomN, lfLog);
        vector[randomN] mixed_state_ind = expand_history_param_to_individuals(mixed_state_option, fixed_mixed_state, mixed_state_mu, mixed_state_sigma, mixed_state_rnd, randomN, lfLogit);
        array[lmN] vector[randomN] bH_ind;
        for(iLM in 1:lmN) {
            if (randomN == 1) {
                bH_ind[iLM] = rep_vector(bH_mu[iLM], randomN);
            }
            else {
                bH_ind[iLM] = bH_mu[iLM] + bH_sigma[iLM] * bH_rnd[iLM];
            }
        }

        for(iT in 1:rowsN){
            // new time-series, recompute absolute tau and reset history state
            if (run_start[iT]){
                current_history = history_starting_values;
                tau = session_tmean[iT] * tau_ind[irandom[iT]];

                // matrix with signal levels
                level = [[1, 0, mixed_state_ind[irandom[iT]]],
                         [0, 1, mixed_state_ind[irandom[iT]]]];
            }

            // for valid percepts, we use history for parameter computation
            if (is_used[iT] == 1){
                // history mixture
                dH = current_history[3-istate[iT]] - current_history[istate[iT]];

                // computing lm for parameters
                for(iLM in 1:lmN) {
                    lm_param[iLM][iC] = a[iLM][irandom[iT]] + bH_ind[iLM][irandom[iT]] * dH;
                    if (fixedN > 0) {
                        lm_param[iLM][iC] += sum(bF[iLM] .* fixed[iT]);
                    }
                }

                iC += 1;
            }

            // computing history for the NEXT episode
            for(iState in 1:2){
                current_history[iState] = level[iState, istate[iT]] + (current_history[iState] - level[iState, istate[iT]]) * exp(-duration[iT] / tau);
            }
        }
    }
}
model {
    {// tau
        if (tau_option == oSingle) {
            tau_mu ~ normal(tau_prior[1], tau_prior[2]);
        }
        else if (tau_option == oIndependent) {
            tau_rnd ~ normal(tau_prior[1], tau_prior[2]);
        }
        else if (tau_option == oPooled) {
            tau_mu ~ normal(tau_prior[1], tau_prior[2]);
            tau_sigma ~ exponential(10);
            tau_rnd ~ normal(0, 1);
        }
    }
    { // Mixed state
        if (mixed_state_option == oSingle) {
            mixed_state_mu ~ normal(mixed_state_prior[1], mixed_state_prior[2]);
        }
        else if (mixed_state_option == oIndependent) {
            mixed_state_rnd ~ normal(mixed_state_prior[1], mixed_state_prior[2]);
        }
        else if (mixed_state_option == oPooled) {
            mixed_state_mu ~ normal(mixed_state_prior[1], mixed_state_prior[2]);
            mixed_state_sigma ~ exponential(10);
            mixed_state_rnd ~ normal(0, 1);
        }
    }

    // linear model parameters
    for(iLM in 1:lmN) {
        // intercepts
        a[iLM] ~ normal(a_prior[iLM, 1], a_prior[iLM, 2]);

        // effect of history
        bH_mu[iLM] ~ normal(bH_prior[iLM, 1], bH_prior[iLM, 2]);
        if (randomN > 1) {
            bH_sigma[iLM] ~ exponential(1);
            bH_rnd[iLM] ~ normal(0, 1);
        }

        if (fixedN > 0) {
          for(iF in 1:fixedN) bF[iLM][iF] ~ normal(fixed_priors[iF, 1], fixed_priors[iF, 2]);
        }
    }

    // variance
    if (varianceN > 0) {
        sigma ~ exponential(1);
    }

    // predicting clear durations
    if (family == lGamma) {
        clear_duration ~ gamma(exp(lm_param[1]), 1 ./ exp(lm_param[2]));
    } else if (family == lLogNormal) {
        clear_duration ~ lognormal(exp(lm_param[1]), sigma[1]);
    } else if (family == lNormal) {
        clear_duration ~ normal(lm_param[1], sigma[1]);
    }
}
generated quantities{
    vector[clearN] log_lik;

    if (family == lGamma) {
        for(iC in 1:clearN) log_lik[iC] = gamma_lpdf(clear_duration[iC] | exp(lm_param[1][iC]), 1 / exp(lm_param[2][iC]));
    } else if (family == lLogNormal) {
        for(iC in 1:clearN) log_lik[iC] = lognormal_lpdf(clear_duration[iC] | exp(lm_param[1][iC]), sigma[1]);
    } else if (family == lNormal) {
        for(iC in 1:clearN) log_lik[iC] = normal_lpdf(clear_duration[iC] | lm_param[1][iC], sigma[1]);
    }
}
