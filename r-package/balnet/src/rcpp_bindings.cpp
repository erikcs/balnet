#include "bindings.hpp"

#include <memory>

#include <adelie_core/configs.hpp>
#include <adelie_core/util/omp.hpp>
#include <adelie_ext/solver_cbps_naive.hpp>


// [[Rcpp::export]]
Rcpp::List rcpp_solver(
    const Rcpp::List args
)
{
    // Setup state
    const Rcpp::NumericMatrix& X = args["X"];
    Eigen::Map<const dense_64F_t> mat_view(X.begin(), X.nrow(), X.ncol());
    matrix_naive_dense_64F_t X_ad(mat_view, args["n_threads"]);
    auto state = std::unique_ptr<state_glm_naive_64_t>(make_state(&X_ad, args));

    // Setup GLM
    const Eigen::Map<vec_value_t> y = args["y"];
    const Eigen::Map<vec_value_t> weights = args["weights"];
    value_t target_scale = args["target_scale"];
    auto glm = std::unique_ptr<glm_cbps_64_t>(new glm_cbps_64_t(y, weights, target_scale));

    // Invoke solver
    ad::Configs::pb_symbol = "\U0001F428";
    bool progress_bar = args["progress_bar"];
    auto pb = ad::util::tq::trange(0);
    pb.set_display(progress_bar);
    pb.set_ostream(Rcpp::Rcerr);
    const auto check_user_interrupt = [&]() {
        Rcpp::checkUserInterrupt();
    };

    std::string error;
    try {
        ad::solver::cbps::naive::solve(
            *state,
            *glm,
            pb,
            [](){ return false; },
            check_user_interrupt
        );
    } catch(const std::exception& e) {
        error = e.what();
    }

    // Stich together results and copy back to R
    const auto& lmda_path = state->lmda_path;
    const auto& lmdas = state->lmdas;
    const auto& devs = state->devs;
    const auto& intercepts = state->intercepts;
    auto betas = convert_betas(X.ncol(), state->betas);

    return Rcpp::List::create(
        Rcpp::Named("error") = error,
        Rcpp::Named("lmda_path") = Rcpp::wrap(lmda_path),
        Rcpp::Named("lmdas") = Rcpp::wrap(lmdas),
        Rcpp::Named("devs") = Rcpp::wrap(devs),
        Rcpp::Named("intercepts") = Rcpp::wrap(intercepts),
        Rcpp::Named("betas") = Rcpp::wrap(betas) // TODO-balnet: could avoid betas allocated twice (not critical unless huge)
    );
}

/* minor utilities */

// [[Rcpp::export]]
Rcpp::List rcpp_col_stats(
    const Rcpp::NumericMatrix& X, // n * p
    const Rcpp::NumericMatrix& weights, // n * L
    bool compute_sd,
    size_t n_threads
)
{
    auto n_threads_def = Eigen::nbThreads();
    Eigen::setNbThreads(n_threads);

    size_t n = X.nrow();
    size_t p = X.ncol();
    size_t L = weights.ncol();
    Rcpp::NumericMatrix center = Rcpp::no_init(L, p);
    Rcpp::NumericMatrix scale;

    Eigen::Map<const dense_64F_t> X_map(X.begin(), n, p);
    Eigen::Map<const dense_64F_t> weights_map(weights.begin(), n, L);
    Eigen::Map<dense_64F_t> center_map(center.begin(), L, p);

    auto weight_sum = weights_map.colwise().sum(); // 1 * L
    center_map.noalias() = weights_map.transpose() * X_map;
    center_map.array().colwise() /= weight_sum.transpose().array();

    if (compute_sd) {
        scale = Rcpp::no_init(L, p);
        Eigen::Map<dense_64F_t> scale_map(scale.begin(), L, p);
        scale_map.noalias() = weights_map.transpose() * X_map.array().square().matrix();
        scale_map.array().colwise() /= weight_sum.transpose().array();
        scale_map.noalias() -= center_map.array().square().matrix();
        // >= 0
        scale_map.array() = (scale_map.array() <= value_t(0)).select(value_t(0), scale_map.array());
        scale_map.array() = scale_map.array().sqrt();
    }
    Eigen::setNbThreads(n_threads_def);

    return Rcpp::List::create(
        Rcpp::Named("center") = center,
        Rcpp::Named("scale") = scale
    );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_standardize(
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& center,
    const Rcpp::NumericVector& scale,
    size_t n_threads
)
{
    size_t n = X.nrow();
    size_t p = X.ncol();
    Rcpp::NumericMatrix out = Rcpp::no_init(n, p);

    Eigen::Map<dense_64F_t> out_map(out.begin(), n, p);
    Eigen::Map<const dense_64F_t> X_map(X.begin(), n, p);

    ad::util::omp_parallel_for(
        [&](auto j) {
            const auto col_in = X_map.col(j).array();
            auto col_out = out_map.col(j).array();

            double mean_j = center[j];
            double sd_j = scale[j];
            if (sd_j <= std::numeric_limits<value_t>::epsilon()) {
                sd_j = 1.0;
            }

            col_out = (col_in - mean_j) / sd_j;
        },
        0, p, n_threads
    );

    return out;
}

// [[Rcpp::export]]
void rcpp_standardize_inplace(
    Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& center,
    const Rcpp::NumericVector& scale,
    size_t n_threads
)
{
    size_t n = X.nrow();
    size_t p = X.ncol();
    Eigen::Map<dense_64F_t> X_map(X.begin(), n, p);

    ad::util::omp_parallel_for(
        [&](auto j) {
            auto col = X_map.col(j).array();

            double mean_j = center[j];
            double sd_j = scale[j];
            if (sd_j <= std::numeric_limits<value_t>::epsilon()) {
                sd_j = 1.0;
            }

            col = (col - mean_j) / sd_j;
        },
        0, p, n_threads
    );
}
