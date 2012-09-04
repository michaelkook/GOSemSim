#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP infoContentMethod_cpp(SEXP ic1_, SEXP ic2_, SEXP mica_, SEXP mic_, SEXP method_) {
  double ic1 = as<double>(ic1_);
  double ic2 = as<double>(ic2_);
  double mica = as<double>(mica_);
  double mic = as<double>(mic_);
  std::string method = as<std::string>(method_);
  double sim;

  // Resnik does not consider how distant the terms are from their common ancestor.
  //  Lin and Jiang take that distance into account.
  if (method == "Resnik") {
    sim = mica;
  }
  if (method == "Lin") {
    sim = 2 * mica/(ic1+ic2);
  }
  if (method == "Jiang") {
    double d = -2 * mica + ic1 + ic2;
    if ( d >= 1)
      d = 1;
    sim = 1- d;
  }
  if (method == "Rel") {
    sim = 2 * mica/(ic1+ic2) * (1-exp(-mica*mic));
    // mica*mic equals to the original IC value.
    // and exp(-mica*mic) equals to the probability of the term's occurence.
  }
  return wrap(sim);
}
