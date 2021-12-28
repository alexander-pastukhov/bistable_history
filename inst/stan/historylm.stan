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
    int<lower=1, upper=3> family;
    // 1 - Gamma, linear model for both shape and rate
    // 2 - Log normal with linear model for the mean
    // 3 - Normal with linear model for the mean

    // --- Complete time-series ---
    int<lower=1> rowsN;   // Number of rows in the COMPLETE multi-timeseries table including mixed phase.
    real duration[rowsN]; // Duration of a dominance/transition phase
    int istate[rowsN];     // Index of a dominance istate, 1 and 2 code for two competing clear states, 3 - transition/mixed.
    int is_used[rowsN];   // Whether history value must used to predict duration or ignored (mixed phases, warm-up period, last, etc.)
    int run_start[rowsN]; // 1 marks a beginning of the new time-series (run/block/etc.)
    real session_tmean[rowsN];    // Mean dominance phase duration for both CLEAR percepts. Used to scale time-constant.

    // --- A shorter clear-states only time-series ---
    int clearN;                  // Number of rows in the clear-states only time-series
    real clear_duration[clearN]; // Durations of clear percepts only.

    // --- Random effects ---
    int<lower=1> randomN;                               // Number of levels for random effects
    int<lower=1, upper=randomN> irandom[rowsN];         // Index of a random effect cluster (participant, display, etc.)

    // --- Fixed effects ---
    int<lower=0> fixedN;
    matrix[clearN, fixedN > 0 ? fixedN : 1] fixed_clear;

    // --- Cumulative history parameters
    real<lower=0, upper=1> history_starting_values[2]; // Starting values for cumulative history at the beginning of the run

    // time constant
    int<lower=1, upper=4> tau_option;  // 1 - constant provided by user, 2 - fit single tau for all, 3 - independent taus, 4 - pooled (multilevel) taus
    real<lower=0> fixed_tau;           // a fixed option (tau_option == 1)
    int tau_mu_size;                   // dimensionality, 1 - sampled, 0 - unused
    int tau_sigma_size;                // dimensionality, 1 - sampled, 0 - unused
    int tau_rnd_size;                  // dimensionality, randomN - sampled, 0 - unused
    real tau_prior[2];                 // prior

    // Mixed state
    int<lower=1, upper=4> mixed_state_option; // 1 - constant provided by user, 2 - fit single tau for all, 3 - independent taus, 4 - pooled (multilevel) taus
    real<lower=0, upper=1> fixed_mixed_state; // a fixed option (mixed_state_option == 1)
    int mixed_state_mu_size;                  // dimensionality, 1 - sampled, 0 - unused
    int mixed_state_sigma_size;               // dimensionality, 1 - sampled, 0 - unused
    int mixed_state_rnd_size;                 // dimensionality, randomN - sampled, 0 - unused
    real mixed_state_prior[2];                // prior

    // History-mixing proportion, used as history_mix * history_same - (1-history_mix) * history_different
    int<lower=1, upper=4> history_mix_option; // 1 - constant provided by user, 2 - fit single tau for all, 3 - independent taus, 4 - pooled (multilevel) taus
    real<lower=0, upper=1> fixed_history_mix; // fixed proportion of history mixing (history_mix_option == 1)
    int history_mix_mu_size;                  // dimensionality, 1 - sampled, 0 - unused
    int history_mix_sigma_size;               // dimensionality, 1 - sampled, 0 - unused
    int history_mix_rnd_size;                 // dimensionality, randomN - sampled, 0 - unused
    real history_mix_prior[2];                // prior

    // --- Linear model ---
    int<lower=1> lmN; // number of linear models, i.e., distribution parameters that are modelled
    int<lower=0, upper=1> varianceN; // number of variance parameters. Effectively, whether variance is sampled at all(1) or not (0)

    // intercept: independent for each parameter and random cluster
    real a_prior[lmN, 2];

    // priors for effect of history
    real bH_prior[lmN, 2];

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

    // History mixture
    vector[history_mix_mu_size] history_mix_mu; // population-level mean
    vector<lower=0>[history_mix_sigma_size] history_mix_sigma; // population-level variance
    vector[history_mix_rnd_size] history_mix_rnd; // individuals

    // --- Linear model ---
    // Independent intercept for each random cluster
    vector[randomN] a[lmN];

    // history term, always multilevel but this makes sense only for randomN > 1
    real bH_mu[lmN];
    real<lower=0> bH_sigma[randomN > 1 ? lmN : 0];
    vector[randomN > 1 ? randomN : 0] bH_rnd[randomN > 1 ? lmN : 0];

    // fixed-effects
    row_vector[fixedN] bF[fixedN > 0 ? lmN : 0];

    // variance
    vector<lower=0>[varianceN] sigma;
}
transformed parameters {
    vector[clearN] lm_param[lmN];
    {
        // Service variables for computing cumulative history
        matrix[2, 3] level;
        real current_history[2];
        real tau;
        real hmix;

        // Index of clear percepts
        int iC = 1;

        // Cumulative history parameters
        vector[randomN] tau_ind = expand_history_param_to_individuals(tau_option, fixed_tau, tau_mu, tau_sigma, tau_rnd, randomN, lfLog);
        vector[randomN] mixed_state_ind = expand_history_param_to_individuals(mixed_state_option, fixed_mixed_state, mixed_state_mu, mixed_state_sigma, mixed_state_rnd, randomN, lfLogit);
        vector[randomN] history_mix_ind = expand_history_param_to_individuals(history_mix_option, fixed_history_mix, history_mix_mu, history_mix_sigma, history_mix_rnd, randomN, lfLogit);
        vector[randomN] bH_ind[lmN];
        for(iP in 1:lmN) {
            if (randomN == 1) {
                bH_ind[iP] = rep_vector(bH_mu[iP], randomN);
            }
            else {
                bH_ind[iP] = bH_mu[iP] + bH_sigma[iP] * bH_rnd[iP];
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
                hmix = history_mix_ind[irandom[iT]] * current_history[istate[iT]] +
                      (1 - history_mix_ind[irandom[iT]]) * current_history[3-istate[iT]];

                // computing lm for parameters
                for(iLM in 1:lmN) {
                    lm_param[iLM][iC] = a[iLM][irandom[iT]] + bH_ind[iLM][irandom[iT]] * hmix;
                    if (fixedN > 0) {
                        lm_param[iLM][iC] += sum(bF[iLM] .* fixed_clear[iC]);
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
    { // History mixture
        if (history_mix_option == oSingle) {
            history_mix_mu ~ normal(history_mix_prior[1], history_mix_prior[2]);
        }
        else if (history_mix_option == oIndependent) {
            history_mix_rnd ~ normal(history_mix_prior[1], history_mix_prior[2]);
        }
        else if (history_mix_option == oPooled) {
            history_mix_mu ~ normal(history_mix_prior[1], history_mix_prior[2]);
            history_mix_sigma ~ exponential(10);
            history_mix_rnd ~ normal(0, 1);
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
        clear_duration ~ gamma(exp(lm_param[1]), exp(lm_param[2]));
    } else if (family == lLogNormal) {
        clear_duration ~ lognormal(exp(lm_param[1]), sigma[1]);
    } else if (family == lNormal) {
        clear_duration ~ normal(lm_param[1], sigma[1]);
    }
}
generated quantities{
    vector[clearN] log_lik;

    if (family == lGamma) {
        for(iC in 1:clearN) log_lik[iC] = gamma_lpdf(clear_duration[iC] | exp(lm_param[1][iC]), exp(lm_param[2][iC]));
    } else if (family == lLogNormal) {
        for(iC in 1:clearN) log_lik[iC] = lognormal_lpdf(clear_duration[iC] | exp(lm_param[1][iC]), sigma[1]);
    } else if (family == lNormal) {
        for(iC in 1:clearN) log_lik[iC] = normal_lpdf(clear_duration[iC] | lm_param[1][iC], sigma[1]);
    }
}
