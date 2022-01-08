#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

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
//'   the run, 0 othwerwise), \code{"session_tmean"} (double)
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
