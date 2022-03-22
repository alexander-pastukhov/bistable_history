#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// distributions
const int L_GAMMA = 1;
const int L_LOGNORMAL = 2;
const int L_NORMAL = 3;

// history options
const int H_DIFFERENCE = -1; // H_OTHER - H_SAME
const int H_DOMINANT = -2;
const int H_SUPPRESSED = -3;
const int H_1 = -4; // history for the first percept (index 1)
const int H_2 = -5; // history for the second percept (index 2)


NumericMatrix state_to_signal_levels(double mixed_level){
  NumericMatrix signal_level(2, 3);

  signal_level(0, 0) = 1;
  signal_level(0, 1) = 0;
  signal_level(0, 2) = mixed_level;
  signal_level(1, 0) = 0;
  signal_level(1, 1) = 1;
  signal_level(1, 2) = mixed_level;

  return signal_level;
}

//' Computes cumulative history
//'
//' Computes cumulative history based on common \code{history} values and
//' \code{normalized_tau} and \code{mixed_state} that are defined for each
//' random cluster / individual.
//'
//' @param df DataFrame with \code{"state"} (integer, 1 and 2 clear state, 3 - mixed state), \code{"duration"} (double),
//'   \code{"irandom"} (integer, 1-based index of a random cluster), \code{"run_start"} (integer, 1 for the first entry of
//'   the run, 0 otherwise), \code{"session_tmean"} (double)
//' @param normalized_tau DoubleVector A normalized tau value for each random cluster / individual. Thus, its length must be
//'   equal to the number of unique indexes in \code{df["irandom"]}.
//' @param mixed_state DoubleVector A values used for the mixed state for each random cluster / individual.
//'   Thus, its length must be equal to the number of unique indexes in \code{df["irandom"]}.
//' @param history_init DoubleVector, size 2. Initial values of history for a run.
//' @return NumericMatrix, size \code{df.nrows()} × 2. Computed history values for each state.
//' @export
//'
//' @examples
//' df <- preprocess_data(br_singleblock, state="State", duration="Duration")
//' fast_history_compute(df, 1, 0.5, c(0, 0))
// [[Rcpp::export]]
NumericMatrix fast_history_compute(DataFrame df, DoubleVector normalized_tau, DoubleVector mixed_state, DoubleVector history_init){
  // output matrix
  NumericMatrix history(df.nrows(), 2);

  // table columns
  IntegerVector state = df["istate"];
  DoubleVector duration = df["duration"];
  IntegerVector irandom = df["irandom"];
  IntegerVector run_start = df["run_start"];
  DoubleVector session_tmean = df["session_tmean"];

  // temporary variables
  DoubleVector current_history = clone(history_init);
  double tau;
  NumericMatrix level(2, 3);

  for(int iR = 0; iR < df.nrows();iR++){
    // new time-series, recompute absolute tau and reset history state
    if (run_start[iR]){
      current_history = clone(history_init);
      tau = session_tmean[iR] * normalized_tau[irandom[iR] - 1];
      level = state_to_signal_levels(mixed_state[irandom[iR] - 1]);
    }

    // saving current history state
    for(int iState = 0; iState < 2; iState ++){
      history(iR, iState) = current_history[iState];
    }

    // computing history for the NEXT episode
    for(int iState = 0; iState < 2; iState ++){
      current_history[iState] = level(iState, state[iR]-1) + (current_history[iState] - level(iState, state[iR]-1)) * exp(-duration[iR] / tau);
    }
  }

  return history;
}


//' Computes prediction for a each sample.
//'
//' Computing prediction for each sample,
//' recomputing cumulative history and uses
//' fitted parameter values.
//'
//' @param family int, distribution family: gamma (1), lognormal(2), or
//' normal (3).
//' @param fixedN int, number of fixed parameters (>= 0).
//' @param randomN int, number of random factors (>= 1).
//' @param lmN int, number of linear models (>= 1).
//' @param istate IntegerVector, zero-based perceptual state 0 or 1,
//' 2 is mixed state.
//' @param duration DoubleVector, duration of a dominance phase.
//' @param is_used IntegerVector, whether dominance phase is used for
//' prediction (1) or not (0).
//' @param run_start IntegerVector, 1 whenever a new run starts.
//' @param session_tmean DoubleVector, average dominance phase duration.
//' @param irandom IntegerVector, zero-based index of a random effect.
//' @param fixed NumericMatrix, matrix with fixed effect values.
//' @param tau_ind NumericMatrix, matrix with samples of tau for each
//' random level.
//' @param mixed_state_ind NumericMatrix, matrix with samples of
//' mixed_state for each random level.
//' @param history_init DoubleVector, Initial values of history for a run
//' @param a NumericMatrix, matrix with samples of
//' a (intercept) for each random level.
//' @param bH NumericMatrix, matrix with sample of
//' bH for each linear model and random level.
//' @param bF NumericMatrix, matrix with sample of
//' bF for each linear model and fixed factor.
//' @param sigma DoubleVector, samples of sigma.
//'
//' @return NumericMatrix with predicted durations for each sample.
//' @export
// [[Rcpp::export]]
NumericMatrix predict_samples(int family,
                              int fixedN,
                              int randomN,
                              int lmN,
                              IntegerVector istate, // data
                              DoubleVector duration,
                              IntegerVector is_used,
                              IntegerVector run_start,
                              DoubleVector session_tmean,
                              IntegerVector irandom,
                              NumericMatrix fixed,
                              NumericMatrix tau_ind, // history parameters
                              NumericMatrix mixed_state_ind,
                              DoubleVector history_init,
                              NumericMatrix a, // linear model parameters
                              NumericMatrix bH,
                              NumericMatrix bF,
                              DoubleVector sigma)
{
  // initialize matrix to NA
  NumericMatrix predicted(sigma.length(), istate.length());
  std::fill(predicted.begin(), predicted.end(), NumericVector::get_na());

  // running history parameters
  DoubleVector current_history = clone(history_init);
  double tau;
  NumericMatrix level(2, 3);
  DoubleVector lm_param(lmN);

  for(int iSample = 0; iSample < sigma.length(); iSample++) {
    for(int iT = 0; iT < istate.length(); iT++) {
      // initialize new run
      if (run_start[iT]) {
        current_history = clone(history_init);
        tau = session_tmean[iT] * tau_ind(iSample, irandom[iT]);
        level = state_to_signal_levels(mixed_state_ind(iSample, irandom[iT]));
      }

      // we need to make a prediction
      if (is_used[iT]) {
        double dH = current_history[1 - istate[iT]] - current_history[istate[iT]];

        // computing parameter(s) from linear model(s)
        for(int iLM = 0; iLM < lmN; iLM++) {
          lm_param[iLM] = a(iSample, irandom[iT] + iLM * randomN) +
            bH(iSample, irandom[iT] + iLM * randomN) * dH;
          for(int iF = 0; iF < fixedN; iF++) {
            lm_param[iLM] += fixed(iT, iF) * bF(iSample, iF + iLM * fixedN);
          }
        }

        // computing mean of the distribution (from its parameters)
        switch (family) {
        // predicting duration (default mode)
        case L_GAMMA:
          predicted(iSample, iT) = exp(lm_param[0]) * exp(lm_param[1]);
          break;
        case L_LOGNORMAL:
          predicted(iSample, iT) = exp(exp(lm_param[0]) + sigma(iSample) / 2);
          break;
        case L_NORMAL:
          predicted(iSample, iT) = lm_param[0];
          break;

        // predicting history
        case H_DIFFERENCE:
          predicted(iSample, iT) = dH;
          break;
        case H_DOMINANT:
          predicted(iSample, iT) = current_history[istate[iT]];
          break;
        case H_SUPPRESSED:
          predicted(iSample, iT) = current_history[1 - istate[iT]];
          break;
        case H_1:
          predicted(iSample, iT) = current_history[0];
          break;
        case H_2:
          predicted(iSample, iT) = current_history[1];
          break;
        }
      }

      // computing history for the NEXT episode
      for(int iS = 0; iS < 2; iS ++){
        current_history[iS] = level(iS, istate[iT]) + (current_history[iS] - level(iS, istate[iT])) * exp(-duration[iT] / tau);
      }
    }
  }
  return predicted;
}
