// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_bsvarTVPs_RCPPEXPORTS_H_GEN_
#define RCPP_bsvarTVPs_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace bsvarTVPs {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("bsvarTVPs", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("bsvarTVPs", "_bsvarTVPs_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in bsvarTVPs");
            }
        }
    }

    inline arma::field<arma::cube> bsvarTVPs_ir_ms(arma::field<arma::cube>& posterior_B, arma::cube& posterior_A, const int horizon, const int p) {
        typedef SEXP(*Ptr_bsvarTVPs_ir_ms)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarTVPs_ir_ms p_bsvarTVPs_ir_ms = NULL;
        if (p_bsvarTVPs_ir_ms == NULL) {
            validateSignature("arma::field<arma::cube>(*bsvarTVPs_ir_ms)(arma::field<arma::cube>&,arma::cube&,const int,const int)");
            p_bsvarTVPs_ir_ms = (Ptr_bsvarTVPs_ir_ms)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvarTVPs_ir_ms");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarTVPs_ir_ms(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(horizon)), Shield<SEXP>(Rcpp::wrap(p)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::cube> >(rcpp_result_gen);
    }

    inline arma::field<arma::cube> bsvarTVPs_ir(arma::cube& posterior_B, arma::cube& posterior_A, const int horizon, const int p) {
        typedef SEXP(*Ptr_bsvarTVPs_ir)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarTVPs_ir p_bsvarTVPs_ir = NULL;
        if (p_bsvarTVPs_ir == NULL) {
            validateSignature("arma::field<arma::cube>(*bsvarTVPs_ir)(arma::cube&,arma::cube&,const int,const int)");
            p_bsvarTVPs_ir = (Ptr_bsvarTVPs_ir)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvarTVPs_ir");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarTVPs_ir(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(horizon)), Shield<SEXP>(Rcpp::wrap(p)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::cube> >(rcpp_result_gen);
    }

    inline arma::cube bsvarTVPs_filter_forecast_smooth(Rcpp::List& posterior, const arma::mat& Y, const arma::mat& X, const bool forecasted, const bool smoothed) {
        typedef SEXP(*Ptr_bsvarTVPs_filter_forecast_smooth)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarTVPs_filter_forecast_smooth p_bsvarTVPs_filter_forecast_smooth = NULL;
        if (p_bsvarTVPs_filter_forecast_smooth == NULL) {
            validateSignature("arma::cube(*bsvarTVPs_filter_forecast_smooth)(Rcpp::List&,const arma::mat&,const arma::mat&,const bool,const bool)");
            p_bsvarTVPs_filter_forecast_smooth = (Ptr_bsvarTVPs_filter_forecast_smooth)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvarTVPs_filter_forecast_smooth");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarTVPs_filter_forecast_smooth(Shield<SEXP>(Rcpp::wrap(posterior)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(forecasted)), Shield<SEXP>(Rcpp::wrap(smoothed)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube bsvarTVPs_fitted_values(arma::cube& posterior_A, arma::mat& X) {
        typedef SEXP(*Ptr_bsvarTVPs_fitted_values)(SEXP,SEXP);
        static Ptr_bsvarTVPs_fitted_values p_bsvarTVPs_fitted_values = NULL;
        if (p_bsvarTVPs_fitted_values == NULL) {
            validateSignature("arma::cube(*bsvarTVPs_fitted_values)(arma::cube&,arma::mat&)");
            p_bsvarTVPs_fitted_values = (Ptr_bsvarTVPs_fitted_values)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvarTVPs_fitted_values");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarTVPs_fitted_values(Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube bsvarTVPs_structural_shocks(const arma::field<arma::cube>& posterior_B, const arma::cube& posterior_A, const arma::cube& posterior_xi, const arma::mat& Y, const arma::mat& X) {
        typedef SEXP(*Ptr_bsvarTVPs_structural_shocks)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarTVPs_structural_shocks p_bsvarTVPs_structural_shocks = NULL;
        if (p_bsvarTVPs_structural_shocks == NULL) {
            validateSignature("arma::cube(*bsvarTVPs_structural_shocks)(const arma::field<arma::cube>&,const arma::cube&,const arma::cube&,const arma::mat&,const arma::mat&)");
            p_bsvarTVPs_structural_shocks = (Ptr_bsvarTVPs_structural_shocks)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvarTVPs_structural_shocks");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarTVPs_structural_shocks(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(posterior_xi)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube bsvars_structural_shocks(const arma::cube& posterior_B, const arma::cube& posterior_A, const arma::mat& Y, const arma::mat& X) {
        typedef SEXP(*Ptr_bsvars_structural_shocks)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvars_structural_shocks p_bsvars_structural_shocks = NULL;
        if (p_bsvars_structural_shocks == NULL) {
            validateSignature("arma::cube(*bsvars_structural_shocks)(const arma::cube&,const arma::cube&,const arma::mat&,const arma::mat&)");
            p_bsvars_structural_shocks = (Ptr_bsvars_structural_shocks)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvars_structural_shocks");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvars_structural_shocks(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline Rcpp::List bsvar_mss_s4_sv_cpp(const int& SS, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB, const Rcpp::List& starting_values, const int thin = 100) {
        typedef SEXP(*Ptr_bsvar_mss_s4_sv_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvar_mss_s4_sv_cpp p_bsvar_mss_s4_sv_cpp = NULL;
        if (p_bsvar_mss_s4_sv_cpp == NULL) {
            validateSignature("Rcpp::List(*bsvar_mss_s4_sv_cpp)(const int&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&,const Rcpp::List&,const int)");
            p_bsvar_mss_s4_sv_cpp = (Ptr_bsvar_mss_s4_sv_cpp)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvar_mss_s4_sv_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvar_mss_s4_sv_cpp(Shield<SEXP>(Rcpp::wrap(SS)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(starting_values)), Shield<SEXP>(Rcpp::wrap(thin)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List bsvar_mss_sv_cpp(const int& SS, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB, const Rcpp::List& starting_values, const int thin = 100) {
        typedef SEXP(*Ptr_bsvar_mss_sv_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvar_mss_sv_cpp p_bsvar_mss_sv_cpp = NULL;
        if (p_bsvar_mss_sv_cpp == NULL) {
            validateSignature("Rcpp::List(*bsvar_mss_sv_cpp)(const int&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&,const Rcpp::List&,const int)");
            p_bsvar_mss_sv_cpp = (Ptr_bsvar_mss_sv_cpp)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvar_mss_sv_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvar_mss_sv_cpp(Shield<SEXP>(Rcpp::wrap(SS)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(starting_values)), Shield<SEXP>(Rcpp::wrap(thin)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List bsvar_s4_sv_cpp(const int& SS, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB, const Rcpp::List& starting_values, const int thin = 100) {
        typedef SEXP(*Ptr_bsvar_s4_sv_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvar_s4_sv_cpp p_bsvar_s4_sv_cpp = NULL;
        if (p_bsvar_s4_sv_cpp == NULL) {
            validateSignature("Rcpp::List(*bsvar_s4_sv_cpp)(const int&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&,const Rcpp::List&,const int)");
            p_bsvar_s4_sv_cpp = (Ptr_bsvar_s4_sv_cpp)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_bsvar_s4_sv_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvar_s4_sv_cpp(Shield<SEXP>(Rcpp::wrap(SS)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(starting_values)), Shield<SEXP>(Rcpp::wrap(thin)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::mat sample_B_heterosk1(arma::mat aux_B, const arma::mat& aux_A, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB) {
        typedef SEXP(*Ptr_sample_B_heterosk1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_B_heterosk1 p_sample_B_heterosk1 = NULL;
        if (p_sample_B_heterosk1 == NULL) {
            validateSignature("arma::mat(*sample_B_heterosk1)(arma::mat,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&)");
            p_sample_B_heterosk1 = (Ptr_sample_B_heterosk1)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_B_heterosk1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_B_heterosk1(Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::cube sample_B_mss(arma::cube aux_B, const arma::mat& aux_A, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& aux_xi, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB) {
        typedef SEXP(*Ptr_sample_B_mss)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_B_mss p_sample_B_mss = NULL;
        if (p_sample_B_mss == NULL) {
            validateSignature("arma::cube(*sample_B_mss)(arma::cube,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&)");
            p_sample_B_mss = (Ptr_sample_B_mss)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_B_mss");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_B_mss(Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(aux_xi)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline Rcpp::List sample_B_heterosk1_s4(arma::mat aux_B, arma::ivec aux_SL, const arma::mat& aux_A, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VBL) {
        typedef SEXP(*Ptr_sample_B_heterosk1_s4)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_B_heterosk1_s4 p_sample_B_heterosk1_s4 = NULL;
        if (p_sample_B_heterosk1_s4 == NULL) {
            validateSignature("Rcpp::List(*sample_B_heterosk1_s4)(arma::mat,arma::ivec,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&)");
            p_sample_B_heterosk1_s4 = (Ptr_sample_B_heterosk1_s4)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_B_heterosk1_s4");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_B_heterosk1_s4(Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_SL)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VBL)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List sample_B_mss_s4(arma::cube aux_B, arma::imat aux_SL, const arma::mat& aux_A, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& aux_xi, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior, const arma::field<arma::mat>& VB) {
        typedef SEXP(*Ptr_sample_B_mss_s4)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_B_mss_s4 p_sample_B_mss_s4 = NULL;
        if (p_sample_B_mss_s4 == NULL) {
            validateSignature("Rcpp::List(*sample_B_mss_s4)(arma::cube,arma::imat,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&)");
            p_sample_B_mss_s4 = (Ptr_sample_B_mss_s4)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_B_mss_s4");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_B_mss_s4(Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_SL)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(aux_xi)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(VB)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::mat sample_A_heterosk1(arma::mat aux_A, const arma::mat& aux_B, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior) {
        typedef SEXP(*Ptr_sample_A_heterosk1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_A_heterosk1 p_sample_A_heterosk1 = NULL;
        if (p_sample_A_heterosk1 == NULL) {
            validateSignature("arma::mat(*sample_A_heterosk1)(arma::mat,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&)");
            p_sample_A_heterosk1 = (Ptr_sample_A_heterosk1)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_A_heterosk1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_A_heterosk1(Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat sample_A_heterosk1_mss(arma::mat aux_A, const arma::cube& aux_B, const arma::mat& aux_xi, const arma::vec& aux_hyper, const arma::mat& aux_sigma, const arma::mat& Y, const arma::mat& X, const Rcpp::List& prior) {
        typedef SEXP(*Ptr_sample_A_heterosk1_mss)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_A_heterosk1_mss p_sample_A_heterosk1_mss = NULL;
        if (p_sample_A_heterosk1_mss == NULL) {
            validateSignature("arma::mat(*sample_A_heterosk1_mss)(arma::mat,const arma::cube&,const arma::mat&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const Rcpp::List&)");
            p_sample_A_heterosk1_mss = (Ptr_sample_A_heterosk1_mss)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_A_heterosk1_mss");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_A_heterosk1_mss(Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_xi)), Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(prior)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::vec sample_hyperparameters_s4(arma::vec aux_hyper, const arma::mat& aux_B, const arma::mat& aux_A, const arma::field<arma::mat>& VB, const arma::ivec& aux_SL, const Rcpp::List& prior) {
        typedef SEXP(*Ptr_sample_hyperparameters_s4)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_hyperparameters_s4 p_sample_hyperparameters_s4 = NULL;
        if (p_sample_hyperparameters_s4 == NULL) {
            validateSignature("arma::vec(*sample_hyperparameters_s4)(arma::vec,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const arma::ivec&,const Rcpp::List&)");
            p_sample_hyperparameters_s4 = (Ptr_sample_hyperparameters_s4)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_hyperparameters_s4");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_hyperparameters_s4(Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(aux_SL)), Shield<SEXP>(Rcpp::wrap(prior)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec sample_hyperparameters_mss(arma::vec aux_hyper, const arma::cube& aux_B, const arma::mat& aux_A, const arma::field<arma::mat>& VB, const Rcpp::List& prior) {
        typedef SEXP(*Ptr_sample_hyperparameters_mss)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_hyperparameters_mss p_sample_hyperparameters_mss = NULL;
        if (p_sample_hyperparameters_mss == NULL) {
            validateSignature("arma::vec(*sample_hyperparameters_mss)(arma::vec,const arma::cube&,const arma::mat&,const arma::field<arma::mat>&,const Rcpp::List&)");
            p_sample_hyperparameters_mss = (Ptr_sample_hyperparameters_mss)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_hyperparameters_mss");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_hyperparameters_mss(Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(prior)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::mat sample_hyperparameters_mss_s4(arma::vec aux_hyper, const arma::cube& aux_B, const arma::mat& aux_A, const arma::field<arma::mat>& VB, const arma::imat& aux_SL, const Rcpp::List& prior) {
        typedef SEXP(*Ptr_sample_hyperparameters_mss_s4)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_hyperparameters_mss_s4 p_sample_hyperparameters_mss_s4 = NULL;
        if (p_sample_hyperparameters_mss_s4 == NULL) {
            validateSignature("arma::mat(*sample_hyperparameters_mss_s4)(arma::vec,const arma::cube&,const arma::mat&,const arma::field<arma::mat>&,const arma::imat&,const Rcpp::List&)");
            p_sample_hyperparameters_mss_s4 = (Ptr_sample_hyperparameters_mss_s4)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_hyperparameters_mss_s4");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_hyperparameters_mss_s4(Shield<SEXP>(Rcpp::wrap(aux_hyper)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_A)), Shield<SEXP>(Rcpp::wrap(VB)), Shield<SEXP>(Rcpp::wrap(aux_SL)), Shield<SEXP>(Rcpp::wrap(prior)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline Rcpp::List svar_nc1(arma::rowvec aux_h_n, double aux_rho_n, double aux_omega_n, double aux_sigma2v_n, double aux_sigma2_omega_n, double aux_s_n, arma::urowvec aux_S_n, const arma::rowvec& u, const Rcpp::List& prior, bool sample_s_ = true) {
        typedef SEXP(*Ptr_svar_nc1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_svar_nc1 p_svar_nc1 = NULL;
        if (p_svar_nc1 == NULL) {
            validateSignature("Rcpp::List(*svar_nc1)(arma::rowvec,double,double,double,double,double,arma::urowvec,const arma::rowvec&,const Rcpp::List&,bool)");
            p_svar_nc1 = (Ptr_svar_nc1)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_svar_nc1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_svar_nc1(Shield<SEXP>(Rcpp::wrap(aux_h_n)), Shield<SEXP>(Rcpp::wrap(aux_rho_n)), Shield<SEXP>(Rcpp::wrap(aux_omega_n)), Shield<SEXP>(Rcpp::wrap(aux_sigma2v_n)), Shield<SEXP>(Rcpp::wrap(aux_sigma2_omega_n)), Shield<SEXP>(Rcpp::wrap(aux_s_n)), Shield<SEXP>(Rcpp::wrap(aux_S_n)), Shield<SEXP>(Rcpp::wrap(u)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(sample_s_)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List svar_nc1_mss(arma::rowvec& aux_h_n, double& aux_rho_n, arma::rowvec& aux_omega_n, double& aux_sigma2_omega_n, double& aux_s_n, arma::urowvec& aux_S_n, const arma::mat& aux_xi, const arma::rowvec& u, const Rcpp::List& prior, bool sample_s_ = true) {
        typedef SEXP(*Ptr_svar_nc1_mss)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_svar_nc1_mss p_svar_nc1_mss = NULL;
        if (p_svar_nc1_mss == NULL) {
            validateSignature("Rcpp::List(*svar_nc1_mss)(arma::rowvec&,double&,arma::rowvec&,double&,double&,arma::urowvec&,const arma::mat&,const arma::rowvec&,const Rcpp::List&,bool)");
            p_svar_nc1_mss = (Ptr_svar_nc1_mss)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_svar_nc1_mss");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_svar_nc1_mss(Shield<SEXP>(Rcpp::wrap(aux_h_n)), Shield<SEXP>(Rcpp::wrap(aux_rho_n)), Shield<SEXP>(Rcpp::wrap(aux_omega_n)), Shield<SEXP>(Rcpp::wrap(aux_sigma2_omega_n)), Shield<SEXP>(Rcpp::wrap(aux_s_n)), Shield<SEXP>(Rcpp::wrap(aux_S_n)), Shield<SEXP>(Rcpp::wrap(aux_xi)), Shield<SEXP>(Rcpp::wrap(u)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(sample_s_)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::mat filtering(const arma::cube& Z, const arma::mat& aux_PR_TR, const arma::vec& pi_0) {
        typedef SEXP(*Ptr_filtering)(SEXP,SEXP,SEXP);
        static Ptr_filtering p_filtering = NULL;
        if (p_filtering == NULL) {
            validateSignature("arma::mat(*filtering)(const arma::cube&,const arma::mat&,const arma::vec&)");
            p_filtering = (Ptr_filtering)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_filtering");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_filtering(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(aux_PR_TR)), Shield<SEXP>(Rcpp::wrap(pi_0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat smoothing(const arma::mat& filtered, const arma::mat& aux_PR_TR) {
        typedef SEXP(*Ptr_smoothing)(SEXP,SEXP);
        static Ptr_smoothing p_smoothing = NULL;
        if (p_smoothing == NULL) {
            validateSignature("arma::mat(*smoothing)(const arma::mat&,const arma::mat&)");
            p_smoothing = (Ptr_smoothing)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_smoothing");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_smoothing(Shield<SEXP>(Rcpp::wrap(filtered)), Shield<SEXP>(Rcpp::wrap(aux_PR_TR)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat sample_Markov_process_mss(arma::mat aux_xi, const arma::mat& E, const arma::cube& aux_B, const arma::mat& aux_sigma, const arma::mat& aux_PR_TR, const arma::vec& aux_pi_0, const bool finiteM = true) {
        typedef SEXP(*Ptr_sample_Markov_process_mss)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_Markov_process_mss p_sample_Markov_process_mss = NULL;
        if (p_sample_Markov_process_mss == NULL) {
            validateSignature("arma::mat(*sample_Markov_process_mss)(arma::mat,const arma::mat&,const arma::cube&,const arma::mat&,const arma::mat&,const arma::vec&,const bool)");
            p_sample_Markov_process_mss = (Ptr_sample_Markov_process_mss)R_GetCCallable("bsvarTVPs", "_bsvarTVPs_sample_Markov_process_mss");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_Markov_process_mss(Shield<SEXP>(Rcpp::wrap(aux_xi)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(aux_B)), Shield<SEXP>(Rcpp::wrap(aux_sigma)), Shield<SEXP>(Rcpp::wrap(aux_PR_TR)), Shield<SEXP>(Rcpp::wrap(aux_pi_0)), Shield<SEXP>(Rcpp::wrap(finiteM)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_bsvarTVPs_RCPPEXPORTS_H_GEN_
