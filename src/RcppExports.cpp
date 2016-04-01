// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// addBacCpp
List addBacCpp(arma::sp_mat occmat, DataFrame orgdat, int amount, double growth, int type, int ptype);
RcppExport SEXP BacArena_addBacCpp(SEXP occmatSEXP, SEXP orgdatSEXP, SEXP amountSEXP, SEXP growthSEXP, SEXP typeSEXP, SEXP ptypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::sp_mat >::type occmat(occmatSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type orgdat(orgdatSEXP);
    Rcpp::traits::input_parameter< int >::type amount(amountSEXP);
    Rcpp::traits::input_parameter< double >::type growth(growthSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type ptype(ptypeSEXP);
    __result = Rcpp::wrap(addBacCpp(occmat, orgdat, amount, growth, type, ptype));
    return __result;
END_RCPP
}
// diffuseGrajdeanuCpp
void diffuseGrajdeanuCpp(Rcpp::NumericMatrix y, double mu, bool donut);
RcppExport SEXP BacArena_diffuseGrajdeanuCpp(SEXP ySEXP, SEXP muSEXP, SEXP donutSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type donut(donutSEXP);
    diffuseGrajdeanuCpp(y, mu, donut);
    return R_NilValue;
END_RCPP
}
// diffuseNaiveCpp
void diffuseNaiveCpp(Rcpp::NumericMatrix y, bool donut);
RcppExport SEXP BacArena_diffuseNaiveCpp(SEXP ySEXP, SEXP donutSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type donut(donutSEXP);
    diffuseNaiveCpp(y, donut);
    return R_NilValue;
END_RCPP
}
// diffuseSteveCpp
void diffuseSteveCpp(Rcpp::NumericMatrix y, double D, double h, double tstep);
RcppExport SEXP BacArena_diffuseSteveCpp(SEXP ySEXP, SEXP DSEXP, SEXP hSEXP, SEXP tstepSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type tstep(tstepSEXP);
    diffuseSteveCpp(y, D, h, tstep);
    return R_NilValue;
END_RCPP
}
// updateSubmat
NumericMatrix updateSubmat(NumericMatrix submat, NumericMatrix sublb_red);
RcppExport SEXP BacArena_updateSubmat(SEXP submatSEXP, SEXP sublb_redSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type submat(submatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sublb_red(sublb_redSEXP);
    __result = Rcpp::wrap(updateSubmat(submat, sublb_red));
    return __result;
END_RCPP
}
// duplicateCpp
DataFrame duplicateCpp(DataFrame orgdat, int n, int m, List cellweight, IntegerMatrix occupyM);
RcppExport SEXP BacArena_duplicateCpp(SEXP orgdatSEXP, SEXP nSEXP, SEXP mSEXP, SEXP cellweightSEXP, SEXP occupyMSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< DataFrame >::type orgdat(orgdatSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< List >::type cellweight(cellweightSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type occupyM(occupyMSEXP);
    __result = Rcpp::wrap(duplicateCpp(orgdat, n, m, cellweight, occupyM));
    return __result;
END_RCPP
}
// movementCpp
void movementCpp(DataFrame orgdat, int n, int m, IntegerMatrix occupyM);
RcppExport SEXP BacArena_movementCpp(SEXP orgdatSEXP, SEXP nSEXP, SEXP mSEXP, SEXP occupyMSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< DataFrame >::type orgdat(orgdatSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type occupyM(occupyMSEXP);
    movementCpp(orgdat, n, m, occupyM);
    return R_NilValue;
END_RCPP
}